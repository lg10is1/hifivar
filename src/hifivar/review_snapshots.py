"""Explicit IGV execution of existing external source plans, with technical QA."""
from __future__ import annotations

from dataclasses import replace
import html
import json
from pathlib import Path

import yaml

from hifivar.command import CommandRunner
from hifivar.context import AnalysisContext
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.igv import IgvWrapper
from hifivar.launch_inputs import IndexedBamPreflight, unique_bam_index
from hifivar.downstream_plan import _PlanLoader as _Loader
from hifivar.review import ReviewManifest, ReviewResult
from hifivar.review_cli import _output_path, _reserve, _write_text, _SEALED_FILES
from hifivar.review_sources import read_review_source_plan, resolve_review_targets
from hifivar.serialization import write_json_atomic, write_yaml_atomic, redact_sensitive_data, utc_now_iso8601
from hifivar.validation import compute_sha256

_DEFAULTS = {'view_mode': 'scientific', 'small_flank_bp': 40, 'breakpoint_flank_bp': 250,
             'long_variant_bp': 1000, 'max_panel_height': 1000, 'sleep_interval_ms': 1000,
             'timeout_seconds': 14400}


def _settings(path: Path | None) -> dict:
    settings = dict(_DEFAULTS)
    if path is not None:
        if not path.is_file() or not 0 < path.stat().st_size <= 65536:
            raise InputValidationError('Snapshot settings must be a nonempty YAML file <=64 KiB.')
        try:
            data = yaml.load(path.read_text(encoding='utf-8-sig'), Loader=_Loader)
        except (yaml.YAMLError, UnicodeError, RecursionError) as error:
            raise InputValidationError('Invalid snapshot settings YAML.') from error
        if not isinstance(data, dict) or set(data) - set(settings):
            raise InputValidationError('Unknown snapshot setting.')
        settings.update(data)
    if not isinstance(settings['view_mode'], str) or settings['view_mode'] not in {'scientific', 'legacy'}:
        raise InputValidationError('Snapshot view_mode must be scientific or legacy.')
    limits = {'small_flank_bp': (1, 100000), 'breakpoint_flank_bp': (1, 100000),
              'long_variant_bp': (1, 100000), 'max_panel_height': (64, 10000),
              'sleep_interval_ms': (0, 60000), 'timeout_seconds': (1, 86400)}
    for key, (low, high) in limits.items():
        if type(settings[key]) is not int or not low <= settings[key] <= high:
            raise InputValidationError(f'Snapshot {key} must be an integer in {low}..{high}.')
    return settings


def _metadata(paths: set[Path]) -> dict:
    return {str(path.resolve()): {'bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns}
            for path in paths}


def _program(value: str, base: Path) -> str:
    return str((base / value).resolve()) if ('/' in value or '\\' in value) and not Path(value).is_absolute() else value


def run_review_snapshots(config, output_directory: str | Path, *, settings_file: str | Path | None = None,
                         dry_run: bool = False, runner: CommandRunner | None = None) -> dict:
    """Run IGV only for an explicit external source plan; leave decisions pending.

    Generated sources remain on their existing owned-launcher execution path.
    No old stored batch command or unverified receipt is executed by this entry.
    """
    output = _output_path(output_directory)
    runner = runner or CommandRunner()
    user = config.sources.get('user')
    base = user.resolve().parent if user is not None else Path.cwd()
    values = config.to_dict()
    source = values['review'].get('source_plan')
    if not source or values['review'].get('selection_file'):
        raise InputValidationError('review snapshots requires review.source_plan and no selection_file.')
    source = Path(source).expanduser()
    source = (base / source if not source.is_absolute() else source).resolve()
    settings_path = Path(settings_file).expanduser().resolve() if settings_file is not None else None
    settings = _settings(settings_path)
    context = AnalysisContext.from_config(config)
    if not context.reference.build:
        raise InputValidationError('review snapshots requires explicit reference.build.')
    plan = read_review_source_plan(source, sample_ids=context.sample_ids,
                                  reference_build=context.reference.build, config=values)
    if any(item.source.origin != 'EXTERNAL' for item in plan.rows):
        raise InputValidationError('Direct snapshots accepts EXTERNAL sources; GENERATED sources use the owned launcher DAG.')
    resolved = resolve_review_targets(plan, config=values, context=context, output_directory=output)
    targets = tuple(replace(target, **{k: settings[k] for k in
                    ('view_mode', 'small_flank_bp', 'breakpoint_flank_bp', 'long_variant_bp')}) for target in resolved.targets)
    if not targets:
        raise InputValidationError('Snapshot plan has no selected targets.')
    wrapper = IgvWrapper(_program(values['review']['igv_executable'], base), runner=runner)
    batch = wrapper.plan(targets, batch_path=output / 'review.igv.batch',
                         snapshot_directory=output / 'screenshots',
                         max_panel_height=settings['max_panel_height'], sleep_interval_ms=settings['sleep_interval_ms'],
                         verify_session=True)
    if not batch.screenshots:
        raise InputValidationError('Snapshot selection has no renderable endpoints; no screenshot run was started.')
    inputs = {source, context.reference.fasta, context.reference.fai}
    if user is not None:
        inputs.add(user)
    if context.source_sample_sheet:
        inputs.add(context.source_sample_sheet)
    if settings_path is not None:
        inputs.add(settings_path)
    for target in targets:
        inputs.update((target.alignment_path, unique_bam_index(target.alignment_path), target.source_vcf,
                       Path(str(target.source_vcf) + '.tbi')))
        if target.bnd_bundle_manifest:
            inputs.add(Path(target.bnd_bundle_manifest))
        if target.trgt_visualization_path:
            inputs.add(target.trgt_visualization_path)
        if target._junction_record is not None and target._junction_record.bundle_binding is not None:
            binding = target._junction_record.bundle_binding
            inputs.add(Path(binding['manifest']['path']))
            inputs.update(Path(item['path']) for item in binding['artifacts'].values())
    for resource in inputs:
        if output.resolve() == resource.resolve() or output.resolve() in resource.resolve().parents:
            raise OutputValidationError('Snapshot output overlaps protected input.')
    before = _metadata(inputs)
    hashable = {p for p in inputs if p.suffix.lower() not in {'.bam', '.fa', '.fasta', '.fna'}}
    hashes = {str(p.resolve()): compute_sha256(p) for p in hashable}
    outcome = {'status': 'PLANNED' if dry_run else 'SNAPSHOTS_GENERATED_TECHNICAL_QA',
               'target_count': len(targets), 'expected_screenshots': len(batch.screenshots),
               'manual_review_status': 'NOT_REVIEWED', 'variant_accuracy': 'NOT_ASSESSED',
               'output_directory': str(output), 'display_settings': settings,
               'igv_preferences_status': 'WINDOWS_PANEL_DELAY_FIXED_OTHER_USER_PREFERENCES_NOT_FULLY_FROZEN',
               'alignment_hash_policy': 'STAT_PLUS_NATIVE_HEADER_EOF_INDEX_PROBE_NOT_FULL_BAM_HASH',
               'source_bindings': list(resolved.source_bindings)}
    if dry_run:
        outcome['screenshots_generated'] = False
        return outcome
    from hifivar.snapshot_validation import require_png_decoder
    require_png_decoder()
    _reserve(output)
    try:
        probe = IndexedBamPreflight(runner, _program(values['workflow']['samtools_executable'], base))
        by_sample = {target.sample_id: target.alignment_path for target in targets}
        native = [probe.run(path, unique_bam_index(path), sample, context.reference, output / 'logs')
                  for sample, path in by_sample.items()]
        igv = wrapper.run(batch, stdout_path=output / 'logs/igv.stdout.log',
                          stderr_path=output / 'logs/igv.stderr.log', timeout_seconds=settings['timeout_seconds'])
        manifest = ReviewManifest(tuple(ReviewResult(item.target, item) for item in igv.evidence))
        manifest.write_json(output / 'review_manifest.json')
        manifest.write_yaml(output / 'review_manifest.yaml')
        manifest.write_tsv(output / 'review_manifest.tsv')
        write_json_atomic(resolved.to_dict(), output / 'source_bindings.json')
        write_yaml_atomic(redact_sensitive_data(values), output / 'effective_config.yaml')
        write_yaml_atomic(settings, output / 'capture_settings.yaml')
        cards = []
        for target in targets:
            pictures = []
            for locus, shot in zip(target.loci, target.screenshot_paths, strict=True):
                relative = shot.relative_to(output).as_posix()
                pictures.append(f'<figure><a href="{html.escape(relative, quote=True)}"><img src="{html.escape(relative, quote=True)}" alt="{html.escape(locus.igv_locus, quote=True)}"></a><figcaption>{html.escape(locus.label+": "+locus.igv_locus)}</figcaption></figure>')
            cards.append('<section><h2>'+html.escape(target.review_id+' / '+target.sample_id+' / '+target.variant_id)+'</h2><p>NOT_REVIEWED</p>'+''.join(pictures)+'</section>')
        _write_text(output / 'index.html', '<!doctype html><meta charset="utf-8"><title>HiFiVar candidate review</title><style>body{font:15px Arial,sans-serif;color:#222;background:#fff;max-width:1400px;margin:24px auto}img{max-width:100%}section{border-top:1px solid #bbb;margin:24px 0}figure{margin:12px 0}</style><h1>HiFiVar candidate review</h1><p>Technical images. Human decisions are pending; images are not independent truth.</p>'+''.join(cards))
        decisions = 'review_id\tsample\tvariant_id\tstatus\treviewer\treviewed_at\tcomment\n'
        decisions += ''.join('\t'.join((target.review_id,target.sample_id,target.variant_id,'NOT_REVIEWED','','',''))+'\n' for target in targets)
        _write_text(output / 'manual_decisions.tsv', decisions)
        if before != _metadata(inputs) or hashes != {str(p.resolve()): compute_sha256(p) for p in hashable}:
            raise InputValidationError('Snapshot source inputs changed during execution; no completion receipt.')
        fresh = resolve_review_targets(plan, config=values, context=context, output_directory=output)
        if fresh.source_bindings != resolved.source_bindings:
            raise InputValidationError('Snapshot selected record identity changed during execution.')
        for item in igv.snapshot_validation:
            if compute_sha256(Path(item['path'])) != item['sha256']:
                raise OutputValidationError('Validated snapshot image changed before completion.')
        for item in igv.session_validation:
            if compute_sha256(Path(item['session_path'])) != item['sha256']:
                raise OutputValidationError('Validated IGV session changed before completion.')
        outcome.update(screenshots_generated=True, native_bam_preflight=native, igv=igv.to_dict(),
                       source_metadata=before, source_sha256=hashes)
        write_json_atomic(outcome, output / 'snapshot_result.json')
        # A snapshot bundle remains compatible with the independent decisions importer.
        write_json_atomic({'schema_version':1,'kind':'hifivar_review_snapshots',
            'created_at':utc_now_iso8601(), 'result':{'status':outcome['status'],'screenshots_generated':True},
            'sealed_artifacts':{name:compute_sha256(output/name) for name in _SEALED_FILES},
            'snapshot_result_sha256':compute_sha256(output/'snapshot_result.json'),
            'snapshot_sha256':{shot.relative_to(output).as_posix():compute_sha256(shot) for shot in batch.screenshots},
            'session_sha256':{path.relative_to(output).as_posix():compute_sha256(path) for path in batch.sessions},
            'session_validation_policy':'REFERENCE_LOCUS_AND_SOURCE_RESOURCES',
            'manual_review_status':'NOT_REVIEWED'}, output / 'review_bundle.json')
        return outcome
    except Exception as error:
        # Preserve partial native evidence; a failure receipt never claims completion.
        write_json_atomic({'status':'FAILED','error_type':type(error).__name__,
                           'manual_review_status':'NOT_REVIEWED'}, output / 'snapshot_failure.json')
        raise
