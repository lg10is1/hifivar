"""Owned review source controls, consumer stage and native receipt binding.

This module orchestrates existing review/IGV infrastructure; it never selects
candidates automatically, modifies calls or assigns human decisions.
"""
from __future__ import annotations

import copy
import csv
from dataclasses import asdict
import io
import json
from pathlib import Path
import re
from typing import Mapping

import yaml

from hifivar.attempt_receipts import artifact_identity, _path, _read
from hifivar.context import AnalysisContext
from hifivar.exceptions import ConfigurationError, InputValidationError, OutputValidationError
from hifivar.review import read_review_selection
from hifivar.review_sources import read_review_source_plan, resolve_review_targets
from hifivar.serialization import standardize_data, write_json_atomic, write_yaml_atomic

_LIMIT = 1024 * 1024
_CONTROL_NAMES = {'effective.yaml', 'selection.tsv', 'source_bindings.json'}
_FIELDS = ('review_id','sample','variant_id','variant_type','contig','start','end',
           'source_vcf','source_caller','evidence_class','flank_bp','bnd_bundle_manifest')


def review_enabled(config: Mapping) -> bool:
    return config.get('review', {}).get('enabled') is True


def _document(plan) -> dict:
    rows = []
    for request in plan.rows:
        row = standardize_data(asdict(request.source))
        row.pop('state'); row.pop('variant_ids')
        candidates = []
        for item in request.candidates:
            candidate = {'review_id':item.review_id,'variant_type':item.variant_type.value}
            if item.site is not None:
                candidate['site'] = asdict(item.site)
            else:
                candidate['variant_id'] = item.variant_id
            candidates.append(candidate)
        row['candidates'] = candidates
        rows.append(row)
    return {'schema_version':1,'review':rows}


def prepare_review_intent(config: dict, context, base: Path, root: Path) -> tuple[dict, list[Path], dict]:
    """Freeze explicit source intent before mkdir; do not run IGV or callers."""
    if not review_enabled(config):
        return {}, [], {}
    section = config['review']
    if not section.get('source_plan') or section.get('selection_file'):
        raise ConfigurationError('hifivar run review requires explicit review.source_plan, not selection_file mode.')
    from hifivar.launcher import _absolute, _ready_vcf
    source = Path(_absolute(section['source_plan'], base))
    if config.get('tr', {}).get('catalog'):
        config['tr']['catalog'] = _absolute(config['tr']['catalog'], base)
    plan = read_review_source_plan(source, sample_ids=context.sample_ids,
                                  reference_build=context.reference.build, config=config)
    inputs = [source]
    if config.get('tr', {}).get('catalog'):
        inputs.append(Path(config['tr']['catalog']))
    audits = []
    for request in plan.rows:
        row = request.source
        if row.origin != 'EXTERNAL':
            continue
        if not str(row.path).endswith('.vcf.gz') or row.index_path != Path(str(row.path)+'.tbi'):
            raise InputValidationError('Review sources require explicit BGZF VCF and conventional TBI.')
        header = _ready_vcf(row.path, row.sample_id, context)
        inputs.extend((row.path, row.index_path))
        if row.bnd_bundle_manifest is not None:
            from hifivar.junction_consumers import validate_bnd_consumer_bundle
            bound = validate_bnd_consumer_bundle(row.bnd_bundle_manifest,row.path,context.reference,row.sample_id)
            inputs.append(Path(bound['manifest']['path']))
            inputs.extend(Path(info['path']) for info in bound['artifacts'].values())
        audits.append({'sample_id':row.sample_id,'track':row.track,'header_audit':header,
                       'version_source':'USER_DECLARED_NOT_REPROBED'})
    config['review']['source_plan'] = str(root/'review_source_plan.yaml')
    return {'review_source_plan.yaml':yaml.safe_dump(_document(plan),sort_keys=False,allow_unicode=True)}, inputs, {
        'kind':'EXPLICIT_REVIEW_SOURCE_INTENT_NOT_EXECUTION','original_plan_sha256':plan.sha256,
        'external_headers':audits,'generated_sources':sum(r.source.origin=='GENERATED' for r in plan.rows),
        'human_review':'NOT_REVIEWED','accuracy':'NOT_ASSESSED'}


def without_review(config: dict) -> dict:
    result = copy.deepcopy(config)
    result['review'].update(enabled=False,source_plan=None,selection_file=None)
    return result


def _json(path: Path):
    if not path.is_file() or not 0 < path.stat().st_size <= _LIMIT:
        raise InputValidationError('Review metadata is missing or exceeds 1 MiB.')
    try:
        return json.loads(path.read_text(encoding='utf-8'))
    except (OSError,UnicodeError,ValueError) as error:
        raise InputValidationError('Invalid review metadata.') from error


def _alignment_witness(root: Path, attempt: Path, context, execution_context: dict, identity: dict) -> dict:
    """Bind prepared BAMs or this attempt's exact native header/EOF/index probes."""
    from hifivar.launch_inputs import unique_bam_index, validate_bam_header
    selected = {r.sample.sample_id:{'bam':artifact_identity(r.sample.input.files[0]),
                 'index':artifact_identity(unique_bam_index(r.sample.input.files[0]))} for r in context.samples}
    if identity['config']['alignment'].get('enabled'):
        from hifivar.alignment_preparation import verify_handoff
        handoff = execution_context.get('alignment_handoff')
        verified = verify_handoff(root,attempt,handoff,verify_bams=True)
        artifacts = {item['path']:item for item in verified['selected_artifacts']}
        if any(artifacts.get(info['path'])!=info for roles in selected.values() for info in roles.values()):
            raise InputValidationError('Review BAM/index differs from the actual prepared alignment handoff.')
        return {'mode':'THIS_ATTEMPT_ALIGNMENT_HANDOFF','handoff':handoff,'selected_bams':selected}
    frozen = {item['path']:item for item in identity['inputs']}
    if any(frozen.get(info['path'])!=info for roles in selected.values() for info in roles.values()):
        raise InputValidationError('Review BAM/index is not an exact frozen external alignment input.')
    probe = _path(attempt/'input_preflight.json',root=root)
    checks = _json(probe)
    if not isinstance(checks,list) or len(checks)!=len(context.sample_ids):
        raise InputValidationError('Review requires this attempt native BAM preflight for every declared sample.')
    if {x.get('sample_id') for x in checks if isinstance(x,dict)} != set(context.sample_ids):
        raise InputValidationError('Review BAM preflight sample set mismatch.')
    logs = {}
    for item in checks:
        sample = item['sample_id']; bam = selected[sample]['bam']['path']; index = selected[sample]['index']['path']
        if any(item.get(k)!='PASS' for k in ('header_identity','quickcheck','indexed_access')):
            raise InputValidationError('Review native BAM preflight is not successful.')
        steps = item.get('steps',[])
        if not isinstance(steps,list) or len(steps)!=4:
            raise InputValidationError('Review BAM preflight must retain its four explicit native commands.')
        for label,step in zip(('version','quickcheck','header','index_probe'),steps,strict=True):
            command = step.get('command')
            if step.get('returncode')!=0 or not isinstance(command,list):
                raise InputValidationError('Review BAM preflight command failed or lacks argv.')
            valid = (command[1:]==['--version'] if label=='version' else
                     command[1:]==['quickcheck','-v',bam] if label=='quickcheck' else
                     command[1:]==['view','--no-PG','-H',bam] if label=='header' else
                     command[1:]==['view','-c','-X',bam,index,f'{context.reference.contigs[0].name}:1-1'])
            if not valid:
                raise InputValidationError('Review BAM probe argv does not select the actual alignment/index.')
            for stream in ('stdout','stderr'):
                path = attempt/f'{sample}.{label}.{stream}.log'
                if step.get(stream+'_path')!=str(path):
                    raise InputValidationError('Review BAM preflight log path differs from its owned attempt.')
                logs[f'{sample}.{label}.{stream}'] = artifact_identity(path)
        header = (attempt/f'{sample}.header.stdout.log').read_text(encoding='utf-8')
        validate_bam_header(header,sample,context.reference)
        version = (attempt/f'{sample}.version.stdout.log').read_text(encoding='utf-8')
        match = re.search(r'(?m)^samtools\s+(\S+)',version)
        if match is None or match.group(1)!=item.get('samtools_version'):
            raise InputValidationError('Review BAM preflight version does not match its native log.')
        if not (attempt/f'{sample}.index_probe.stdout.log').read_text().strip().isdigit():
            raise InputValidationError('Review BAM index probe has no numeric result.')
    return {'mode':'THIS_ATTEMPT_NATIVE_HEADER_EOF_INDEX_PROBES_NOT_FULL_RECORD_QC',
            'record':artifact_identity(probe),'selected_bams':selected,'probe_logs':logs}


def review_inputs(targets) -> dict:
    from hifivar.launch_inputs import unique_bam_index
    result = {}
    for item in targets:
        roles = {'vcf':item.source_vcf,'index':Path(str(item.source_vcf)+'.tbi'),
                 'alignment':item.alignment_path,'alignment_index':unique_bam_index(item.alignment_path)}
        if item.bnd_bundle_manifest is not None:
            roles['bnd_bundle_manifest'] = item.bnd_bundle_manifest
        result[item.review_id] = roles
    return result


def review_outputs(targets, output: Path) -> dict:
    result = {name:output/filename for name,filename in {
        'manifest_json':'review_manifest.json','manifest_yaml':'review_manifest.yaml',
        'manifest_tsv':'review_manifest.tsv','batch':'review.igv.batch',
        'provenance':'review_provenance.json'}.items()}
    for number,path in enumerate((p for item in targets for p in item.screenshot_paths),1):
        result[f'screenshot_{number:04d}'] = path
    return result


def _selection(targets) -> str:
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream,fieldnames=_FIELDS,delimiter='\t',lineterminator='\n')
    writer.writeheader()
    for target in targets:
        data = target.to_dict()
        data['bnd_bundle_manifest'] = str(target.bnd_bundle_manifest) if target.bnd_bundle_manifest else ''
        writer.writerow({key:data[key] for key in _FIELDS})
    return stream.getvalue()


def materialize_review(config: dict, context, root: Path, attempt: Path, *, calling_config: dict) -> tuple[Path,dict]:
    """Write a fresh consumer-only configuration after exact producer resolution."""
    from hifivar.attempt_receipts import _context
    from hifivar.downstream_handoff import consumer_config
    execution,identity,_ = _context(attempt/'execution_context.json',live=False)
    if execution['dry_run']:
        raise InputValidationError('Dry-run cannot materialize review execution sources.')
    context = AnalysisContext.from_config(calling_config)
    alignment = _alignment_witness(root,attempt,context,execution,identity)
    planfile = root/'review_source_plan.yaml'
    plan = read_review_source_plan(planfile,sample_ids=context.sample_ids,
                                  reference_build=context.reference.build,config=config)
    output = root/'results/review'
    resolved = resolve_review_targets(plan,config=config,context=context,output_directory=output,
                                     run_root=root,attempt_dir=attempt)
    destination = _path(attempt/'review',root=root)
    if destination.exists():
        raise OutputValidationError('Review materialization exists; use a new explicit attempt.')
    effective = consumer_config(calling_config)
    for name in ('annotation','cohort'):
        effective[name]['enabled'] = False
    effective['review'] = copy.deepcopy(config['review'])
    effective['review'].update(source_plan=None,selection_file=str(destination/'selection.tsv'),overwrite=False)
    from hifivar.config import validate_config
    validate_config(effective,require_complete=True)
    # Reject oversized bindings before any review materialization writes.
    if len(json.dumps(resolved.to_dict(),ensure_ascii=False).encode('utf-8'))>_LIMIT:
        raise InputValidationError('Review source binding exceeds the 1 MiB metadata contract.')
    destination.mkdir(parents=True,exist_ok=False)
    with (destination/'selection.tsv').open('x',encoding='utf-8',newline='') as stream:
        stream.write(_selection(resolved.targets))
    write_json_atomic(resolved.to_dict(),destination/'source_bindings.json')
    write_yaml_atomic(effective,destination/'effective.yaml')
    record = {'schema_version':1,'kind':'hifivar_review_materialization','attempt_id':attempt.name,
              'context_sha256':artifact_identity(attempt/'execution_context.json')['sha256'],
              'plan_sha256':artifact_identity(planfile)['sha256'],
              'controls':{name:artifact_identity(destination/name) for name in sorted(_CONTROL_NAMES)},
              'alignment':alignment,'source_receipts':[b['source_receipt'] for b in resolved.source_bindings if b['source_receipt']],
              'selected_sample_order':[s for s in context.sample_ids if any(t.sample_id==s for t in resolved.targets)],
              'output_directory':str(output),'human_review':'NOT_REVIEWED','accuracy':'NOT_ASSESSED'}
    if len(json.dumps(record).encode('utf-8'))>_LIMIT:
        raise InputValidationError('Review materialization exceeds the bounded metadata contract.')
    write_json_atomic(record,destination/'materialization.json')
    return destination/'effective.yaml',record


def validate_review_materialization(root: Path, attempt: Path, *, live: bool=False):
    """Independently re-resolve current sources/alignments against frozen intent."""
    from hifivar.attempt_receipts import _context
    execution,identity,_ = _context(attempt/'execution_context.json',live=live)
    if execution['dry_run'] or not identity['config'].get('review',{}).get('source_plan'):
        raise InputValidationError('Review materialization requires a native explicitly configured source plan.')
    destination = _path(attempt/'review',root=root)
    record = _read(destination/'materialization.json')
    status = _read(attempt/'status.json')
    if (record.get('schema_version')!=1 or record.get('kind')!='hifivar_review_materialization'
        or record.get('attempt_id')!=attempt.name or record.get('context_sha256')!=artifact_identity(attempt/'execution_context.json')['sha256']
        or record.get('plan_sha256')!=artifact_identity(root/'review_source_plan.yaml')['sha256']
        or artifact_identity(destination/'materialization.json')['sha256']!=status.get('review_materialization_sha256')):
        raise InputValidationError('Review materialization is not bound to this exact attempt/stage/plan.')
    controls = record.get('controls')
    if not isinstance(controls,dict) or set(controls)!=_CONTROL_NAMES:
        raise InputValidationError('Review materialization controls are incomplete.')
    for name,info in controls.items():
        if info!=artifact_identity(destination/name):
            raise InputValidationError('Review control changed during execution.')
    effective = yaml.safe_load((destination/'effective.yaml').read_text(encoding='utf-8'))
    from hifivar.config import validate_config
    validate_config(effective,require_complete=True)
    if (effective['review'].get('source_plan') is not None
        or effective['review'].get('selection_file')!=str(destination/'selection.tsv')
        or not effective['review'].get('enabled') or effective['review'].get('overwrite')
        or any(effective[g].get('enabled') for g in ('alignment','small','sv','tr','phasing','assembly','assembly_sv','annotation','cohort'))
        or effective['sv']['harmonization'].get('enabled')):
        raise InputValidationError('Review consumer configuration reactivates producers or changes its source mode.')
    context = AnalysisContext.from_config(effective)
    if _alignment_witness(root,attempt,context,execution,identity)!=record.get('alignment'):
        raise InputValidationError('Review actual alignment/probe binding changed.')
    original = identity['config']
    plan = read_review_source_plan(root/'review_source_plan.yaml',sample_ids=context.sample_ids,
                                  reference_build=context.reference.build,config=original)
    output = root/'results/review'
    if record.get('output_directory')!=str(output):
        raise InputValidationError('Review output escaped its explicit owned result directory.')
    fresh = resolve_review_targets(plan,config=original,context=context,output_directory=output,run_root=root,attempt_dir=attempt)
    if _json(destination/'source_bindings.json')!=fresh.to_dict() or (destination/'selection.tsv').read_text(encoding='utf-8')!=_selection(fresh.targets):
        raise InputValidationError('Review selected records or normalized selection differ from frozen source intent.')
    expected_samples = tuple(s for s in context.sample_ids if any(t.sample_id==s for t in fresh.targets))
    if record.get('selected_sample_order')!=list(expected_samples):
        raise InputValidationError('Review selected sample order differs from the declared cohort context.')
    if record.get('source_receipts')!=[b['source_receipt'] for b in fresh.source_bindings if b['source_receipt']]:
        raise InputValidationError('Review generated source receipt links changed.')
    frozen_inputs = {item['path']:item for item in identity['inputs']}
    for binding in fresh.source_bindings:
        if binding['origin']=='EXTERNAL':
            for info in binding['source_identity'].values():
                if frozen_inputs.get(info['path'])!=info:
                    raise InputValidationError('External review VCF/index differs from its frozen launcher input.')
            for selected in binding['selected_records']:
                bundle = selected.get('bundle_binding')
                if bundle and any(frozen_inputs.get(info['path'])!={k:info[k] for k in ('path','bytes','sha256')} for info in [bundle['manifest'],*bundle['artifacts'].values()]):
                    raise InputValidationError('External review BND bundle differs from frozen source roles.')
    return record,effective,fresh,expected_samples


def review_receipt_binding(context_path,context,identity,root,samples,inputs,input_manifest,*,live):
    record,effective,resolved,expected = validate_review_materialization(root,context_path.parent,live=live)
    selected = context_path.parent/'review/selection.tsv'
    if tuple(samples)!=expected or input_manifest is None or Path(input_manifest).absolute()!=selected:
        raise InputValidationError('Actual review unit/manifest differs from its selected attempt controls.')
    actual = {key:{role:Path(value).absolute() if value is not None else None for role,value in roles.items()} for key,roles in (inputs or {}).items()}
    expected_inputs = review_inputs(resolved.targets)
    if actual!=expected_inputs:
        raise InputValidationError('Actual review VCF/index/alignment inputs differ from selected sources.')
    return {'schema_version':1,'mode':'EXPLICIT_REVIEW_SOURCE_PLAN','execution_stage':'REVIEW_DAG_RUNNING',
            'input_manifest':artifact_identity(selected),'source_inputs':{
                key:{role:artifact_identity(path) for role,path in roles.items()} for key,roles in actual.items()},
            'materialization':artifact_identity(context_path.parent/'review/materialization.json'),
            'controls':record['controls'],'alignment':record['alignment'],'producer_receipts':record['source_receipts'],
            'selection':{'samples':list(expected),'review_ids':[t.review_id for t in resolved.targets]},
            'expected_outputs':{role:str(path) for role,path in review_outputs(resolved.targets,root/'results/review').items()},
            'human_review':'NOT_REVIEWED','accuracy':'NOT_ASSESSED'}


def _png_header(path: Path) -> dict:
    """Compatibility helper now performs the shared complete PNG check."""
    from hifivar.snapshot_validation import validate_png_snapshot
    return validate_png_snapshot(path)


def run_review_step(config, context, targets, *, output: Path, selection: Path, igv_wrapper=None):
    """Execute existing IGV through CommandRunner and seal actual command evidence."""
    from hifivar.attempt_receipts import attempt_step
    from hifivar.igv import IgvWrapper, _protect_write_paths
    from hifivar.phase10 import run_phase10, Phase10Status
    wrapper = igv_wrapper or IgvWrapper(executable=config['review']['igv_executable'])
    targets = tuple(targets)
    samples = tuple(s for s in context.sample_ids if any(t.sample_id==s for t in targets))
    outputs = review_outputs(targets,output)
    # Match the collector's single named context hook. Its presence is not
    # execution evidence until attempt_step validates the complete context.
    import hifivar.attempt_receipts as receipts
    owned = bool(receipts.os.environ.get(receipts.ATTEMPT_CONTEXT_ENV))
    if owned and not targets:
        raise InputValidationError('Owned source plans cannot produce an empty review unit.')
    plan = wrapper.plan(targets,batch_path=outputs['batch'],snapshot_directory=output/'screenshots')
    _protect_write_paths(plan,executable=wrapper.executable,
        stdout_path=output/'logs/igv.stdout.log',stderr_path=output/'logs/igv.stderr.log',
        additional_outputs=(('review_provenance',outputs['provenance']),))
    if not owned:
        # Preserve the independent workflow's plain-VCF/CRAM and empty-selection
        # contracts. It emits no launcher receipt or generated-source claim.
        report = run_phase10(targets,output_directory=output,igv_wrapper=wrapper,overwrite=False)
        _protect_write_paths(plan,executable=wrapper.executable,
            stdout_path=output/'logs/igv.stdout.log',stderr_path=output/'logs/igv.stderr.log',
            additional_outputs=(('review_provenance',outputs['provenance']),))
        write_json_atomic(report.to_dict(),outputs['provenance'])
        return report
    with attempt_step('review',samples,outputs=outputs,inputs=review_inputs(targets),input_manifest=selection) as receipt:
        if receipt.context_path is not None:
            _,_,fresh,_ = validate_review_materialization(receipt.context_path.parent.parent.parent,receipt.context_path.parent,live=True)
            if targets != fresh.targets:
                raise InputValidationError("Actual IGV targets differ from the frozen source-record coordinates or provenance.")
        report = run_phase10(targets,output_directory=output,igv_wrapper=wrapper,overwrite=False)
        if report.status is not Phase10Status.COMPLETED or not report.igv.tool_version:
            raise OutputValidationError('Review workflow did not finish native IGV execution.')
        dimensions = {str(path):_png_header(path) for key,path in outputs.items() if key.startswith('screenshot_')}
        _protect_write_paths(plan,executable=wrapper.executable,
            stdout_path=output/'logs/igv.stdout.log',stderr_path=output/'logs/igv.stderr.log',
            additional_outputs=(('review_provenance',outputs['provenance']),))
        write_json_atomic(report.to_dict(),outputs['provenance'])
        receipt.complete(outputs=outputs,tool_versions={'igv':report.igv.tool_version},
            validation={'validated':True,'screenshots':dimensions,'scope':'PNG_CRC_AND_FULL_DECODE_NOT_TRACK_OR_HUMAN_QA'},
            metadata={'reference_build':context.reference.build,'human_review':'NOT_REVIEWED',
                      'accuracy':'NOT_ASSESSED','review_ids':[t.review_id for t in targets],
                      'tool_version_sources':{'igv':'COMMAND'},'caller_support_is_confidence':False})
    return report
