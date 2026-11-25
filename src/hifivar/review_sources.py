"""Explicit source/record binding for future unified IGV review orchestration.

This library does not start a workflow or IGV, select candidates automatically,
write a selection/manifest, or certify a human review. Separate launcher
orchestration lives in review_handoff; these bindings alone prove no execution.
"""
from __future__ import annotations

from dataclasses import dataclass
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re
from typing import Mapping, Sequence

import yaml

from hifivar.attempt_receipts import artifact_identity
from hifivar.downstream_plan import (
    DownstreamPlanRow, _PlanLoader, _ROW_KEYS, _build, _configured_catalog,
    _mapping, _plan_bytes, _row, _sample_ids, _text, _unknown_keys,
)
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.junction_consumers import (
    SelectedJunctionRecord, _info, _lines, read_selected_junction_records,
)
from hifivar.repeat_selection import SelectedRepeatRecord, read_selected_repeat_records
from hifivar.review import EvidenceClass, ReviewTarget, VariantClass, _SAFE_ID
from hifivar.sample import InputType


_LIMIT = 4096
_SEQUENCE = re.compile(r'[ACGTNacgtn]+\Z')


@dataclass(frozen=True, slots=True)
class ReviewSite:
    """Exact source-record alleles for unnamed ordinary records, never an ID rewrite."""
    contig: str
    position: int
    ref: str
    alt: str

    def __post_init__(self) -> None:
        for name in ('contig','ref','alt'):
            value = _text(getattr(self,name),name)
            if any(c.isspace() for c in value) or ';' in value:
                raise InputValidationError('Review site values must be explicit single-line fields.')
        if type(self.position) is not int or self.position < 1:
            raise InputValidationError('Review site position must be a positive integer.')

    @property
    def key(self) -> str:
        return 'SITE:' + json.dumps([self.contig,self.position,self.ref,self.alt],separators=(',',':'))


@dataclass(frozen=True, slots=True)
class SelectedReviewSite:
    source_vcf: Path
    source_sha256: str
    ordinal: int
    record_sha256: str
    selector: str
    raw_vcf_id: str
    contig: str
    position: int
    end: int
    ref: str
    alt: str
    junction: None = None

    def to_dict(self) -> dict[str, object]:
        return {'source_vcf':str(self.source_vcf),'source_sha256':self.source_sha256,
                'record_ordinal':self.ordinal,'selected_record_sha256':self.record_sha256,
                'source_record_key':f'sha256:{self.source_sha256}:record:{self.ordinal}',
                'selection_namespace':'EXACT_SITE_RECORD','selector':self.selector,
                'raw_vcf_id':self.raw_vcf_id,'contig':self.contig,'position':self.position,
                'end':self.end,'ref':self.ref,'alt':self.alt,'raw_vcf_modified':False}


@dataclass(frozen=True, slots=True)
class ReviewCandidate:
    """Explicit source-record selector and desired display class, not confidence."""
    review_id: str
    variant_id: str
    variant_type: VariantClass
    site: ReviewSite | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.review_id, str) or _SAFE_ID.fullmatch(self.review_id) is None:
            raise InputValidationError('Review source candidate review_id must be ASCII-safe.')
        if _text(self.variant_id, 'variant_id') == '.' or ';' in self.variant_id:
            raise InputValidationError('Review source candidates need explicit unambiguous variant IDs.')
        if not isinstance(self.variant_type, VariantClass):
            raise InputValidationError('Review source candidate variant_type must be VariantClass.')
        if self.site is not None and (not isinstance(self.site,ReviewSite) or self.variant_id!=self.site.key):
            raise InputValidationError('Review site and its canonical selector label differ.')


@dataclass(frozen=True, slots=True)
class ReviewSourceRequest:
    source: DownstreamPlanRow
    candidates: tuple[ReviewCandidate, ...]


@dataclass(frozen=True, slots=True)
class ReviewSourcePlan:
    path: Path
    sha256: str
    sample_ids: tuple[str, ...]
    reference_build: str
    rows: tuple[ReviewSourceRequest, ...]
    schema_version: int = 1


@dataclass(frozen=True, slots=True)
class ResolvedReviewTargets:
    """Read-only record/producer bindings; no snapshots or execution receipt."""
    plan: ReviewSourcePlan
    targets: tuple[ReviewTarget, ...]
    source_bindings: tuple[dict, ...]

    def to_dict(self) -> dict[str, object]:
        return {'schema_version': 1, 'source_plan': str(self.plan.path),
                'source_plan_sha256': self.plan.sha256,
                'targets': [target.to_dict() for target in self.targets],
                'source_bindings': deepcopy(list(self.source_bindings)),
                'review_status': 'NOT_REVIEWED', 'IGV_executed': False,
                'main_launcher_integration': 'NOT_IMPLEMENTED_BY_THIS_LIBRARY',
                'selection_unit': 'EXPLICIT_SOURCE_RECORD_ALL_ALLELES_RETAINED',
                'accuracy': 'NOT_ASSESSED', 'caller_support_is_confidence': False}


def read_review_source_plan(
    path: Path, *, sample_ids: Sequence[str], reference_build: str,
    config: Mapping,
) -> ReviewSourcePlan:
    """Read bounded strict YAML; generated selectors never invent future paths."""
    samples = _sample_ids(sample_ids)
    build = _build(reference_build, 'reference_build')
    if not isinstance(config, Mapping):
        raise InputValidationError('Review source-plan config must be a mapping.')
    source, raw = _plan_bytes(path)
    try:
        text = raw.decode('utf-8-sig')
        if len(text.splitlines()) > _LIMIT:
            raise InputValidationError('Review source plan exceeds the 4096-line limit.')
        doc = _mapping(yaml.load(text, Loader=_PlanLoader), 'Review source plan')
    except (UnicodeError, yaml.YAMLError, RecursionError) as error:
        raise InputValidationError('Review source plan requires bounded valid UTF-8 YAML.') from error
    _unknown_keys(doc, frozenset({'schema_version', 'review'}), 'Review source plan')
    if type(doc.get('schema_version')) is not int or doc['schema_version'] != 1:
        raise InputValidationError('Review source plan schema_version must be integer 1.')
    rows = doc.get('review')
    if not isinstance(rows, list) or not rows or len(rows) > _LIMIT:
        raise InputValidationError('Review source plan requires 1-4096 explicit source rows.')
    has_tr = any(isinstance(row, Mapping) and row.get('track') == 'tr' for row in rows)
    catalog = _configured_catalog(config, source.parent) if has_tr else None
    if has_tr and config.get('tr', {}).get('catalog_reference_build') is not None:
        if _build(config['tr']['catalog_reference_build'], 'TR catalog build') != build:
            raise InputValidationError('Review TR catalog reference build mismatch.')
    parsed = []
    seen_ids, seen_rows = set(), set()
    for value in rows:
        row = _mapping(value, 'Review source row')
        _unknown_keys(row, (_ROW_KEYS - {'variant_ids', 'state'}) | {'candidates'}, 'Review source row')
        entries = row.get('candidates')
        if not isinstance(entries, list) or not entries:
            raise InputValidationError('Review sources require explicit nonempty candidates.')
        candidates = []
        for value in entries:
            candidate = _mapping(value, 'Review candidate')
            _unknown_keys(candidate, frozenset({'review_id', 'variant_id', 'variant_type', 'site'}), 'Review candidate')
            try:
                kind = VariantClass(_text(candidate.get('variant_type'), 'variant_type'))
            except ValueError as error:
                raise InputValidationError('Unsupported review candidate variant_type.') from error
            if ('site' in candidate)==('variant_id' in candidate):
                raise InputValidationError('Review candidate must choose variant_id OR exact site, not both/neither.')
            site=None
            if 'site' in candidate:
                values=_mapping(candidate['site'],'Review site')
                if set(values)!={'contig','position','ref','alt'}:
                    raise InputValidationError('Review site needs exact contig/position/ref/alt fields.')
                site=ReviewSite(**values)
                if kind in {VariantClass.BND,VariantClass.TR}:
                    raise InputValidationError('Review BND/TR use unique record IDs or explicit TRID, not site coercion.')
            item = ReviewCandidate(candidate.get('review_id'),site.key if site else candidate.get('variant_id'),kind,site)
            if item.review_id in seen_ids or len(seen_ids) >= _LIMIT:
                raise InputValidationError('Review IDs must be globally unique; at most4096 candidates.')
            seen_ids.add(item.review_id)
            candidates.append(item)
        fields = {key: value for key, value in row.items() if key != 'candidates'}
        fields['variant_ids'] = [item.variant_id for item in candidates]
        selected = _row(fields, 'review', samples, build, config, source.parent, catalog)
        allowed = {'small': {VariantClass.SNV, VariantClass.INDEL},
                   'sv': {VariantClass.INS, VariantClass.DEL, VariantClass.DUP, VariantClass.INV, VariantClass.BND},
                   'tr': {VariantClass.TR}}[selected.track]
        if any(item.variant_type not in allowed for item in candidates):
            raise InputValidationError('Review variant_type does not match the source category.')
        if selected.track == 'tr' and any(not item.variant_id.startswith('TRID:') for item in candidates):
            raise InputValidationError('Review TR sources require explicit TRID selectors, no ordinary-ID fallback.')
        key = (selected.sample_id, selected.track, selected.origin, selected.producer,
               selected.output_role, selected.path)
        if key in seen_rows:
            raise InputValidationError('Repeated review source rows must combine their candidate lists.')
        seen_rows.add(key)
        parsed.append(ReviewSourceRequest(selected, tuple(candidates)))
    return ReviewSourcePlan(source, hashlib.sha256(raw).hexdigest(), samples, build, tuple(parsed))


def _selected_sites(source: Path, candidates: tuple[ReviewCandidate, ...], reference, sample_id: str) -> dict[str, SelectedReviewSite]:
    from hifivar.annotation import audit_annotation_vcf_identity
    audit_annotation_vcf_identity(source,sample_id,reference)
    before=artifact_identity(source);wanted={item.site.key:item.site for item in candidates if item.site is not None};found={}
    lengths={c.name:c.length for c in reference.contigs};ordinal=0
    for raw in _lines(source):
        if raw.startswith(b'#') or not raw.strip():continue
        ordinal+=1;fields=raw.decode('utf-8').rstrip('\r\n').split('\t')
        if len(fields)<8:raise InputValidationError('Malformed source record during exact-site review selection.')
        try:position=int(fields[1])
        except ValueError as error:raise InputValidationError('Invalid source POS during review selection.') from error
        key='SITE:'+json.dumps([fields[0],position,fields[3],fields[4]],separators=(',',':'))
        if key not in wanted:continue
        if key in found:raise InputValidationError('Exact review site is ambiguous: duplicate source records.')
        info=_info(fields[7]);alt=fields[4]
        if ('[' in alt or ']' in alt or info.get('SVTYPE')=='BND' or alt.startswith('.') or alt.endswith('.')):
            raise InputValidationError('Exact-site BND selection is not exposed; use a unique actual BND record ID.')
        try:end=int(info.get('END',str(position)))
        except (ValueError,TypeError) as error:raise InputValidationError('Review site has invalid source END.') from error
        if fields[0] not in lengths or not 1<=position<=end<=lengths[fields[0]]:
            raise InputValidationError('Selected review site contradicts reference bounds.')
        found[key]=SelectedReviewSite(source,before['sha256'],ordinal,hashlib.sha256(raw).hexdigest(),key,
                                     fields[2],fields[0],position,end,fields[3],alt)
    if set(found)!=set(wanted):raise InputValidationError('Exact review site is absent; no positional/allele fallback is allowed.')
    if artifact_identity(source)!=before:raise InputValidationError('Review site source changed during selection.')
    return found


def _fields_for_records(source: Path, records: tuple) -> dict[int, list[str]]:
    selected = {item.ordinal: item.record_sha256 for item in records}
    found = {}
    ordinal = 0
    for raw in _lines(source):
        if raw.startswith(b'#') or not raw.strip():
            continue
        ordinal += 1
        if ordinal in selected:
            if hashlib.sha256(raw).hexdigest() != selected[ordinal]:
                raise InputValidationError('Selected review record changed during resolution.')
            found[ordinal] = raw.decode('utf-8').rstrip('\r\n').split('\t')
    if set(found) != set(selected):
        raise InputValidationError('Selected review records disappeared during resolution.')
    return found


def _coordinates(candidate: ReviewCandidate, record, fields: list[str], lengths: Mapping) -> tuple[str, int, int]:
    kind = candidate.variant_type
    if isinstance(record, SelectedRepeatRecord):
        if kind is not VariantClass.TR:
            raise InputValidationError('TRID records must retain their independent TR category.')
        return record.contig, record.start, record.end
    if record.junction is not None:
        if kind is not VariantClass.BND:
            raise InputValidationError('BND ALT cannot become an ordinary interval review target.')
        return record.contig, record.position, record.position
    if kind is VariantClass.BND:
        raise InputValidationError('Declared BND review target has no BND ALT.')
    info = _info(fields[7])
    if kind in {VariantClass.SNV, VariantClass.INDEL}:
        ref, alts = fields[3], fields[4].split(',')
        if not _SEQUENCE.fullmatch(ref) or any(not _SEQUENCE.fullmatch(alt) or alt == ref for alt in alts):
            raise InputValidationError('Small review source needs explicit sequence alleles, not symbolic/gVCF/star alleles.')
        valid = (len(ref) == 1 and all(len(alt) == 1 for alt in alts)) if kind is VariantClass.SNV else any(len(alt) != len(ref) for alt in alts)
        if not valid:
            raise InputValidationError('Small review class contradicts source allele lengths; do not coerce MNV.')
        end = record.position + len(ref) - 1
        if end > lengths[record.contig] or ('END' in info and record.end != end):
            raise InputValidationError('Small review source REF span/END contradicts reference coordinates.')
        return record.contig, record.position, end
    if info.get('SVTYPE') != kind.value:
        raise InputValidationError('Review SV class differs from the selected source SVTYPE.')
    if kind is VariantClass.INS:
        return record.contig, record.position, record.position
    if 'END' not in info or record.end <= record.position:
        raise InputValidationError('DEL/DUP/INV review requires an explicit nonempty source END span.')
    return record.contig, record.position, record.end


def resolve_review_targets(
    plan: ReviewSourcePlan, *, config: Mapping, context,
    output_directory: Path, run_root: Path | None = None, attempt_dir: Path | None = None,
) -> ResolvedReviewTargets:
    """Bind explicit candidates to immutable sources without invoking IGV.

    Generated rows use only the selected attempt's validated native producer
    receipt. A matching path, older attempt, metadata declaration or dry-run is
    not execution evidence. Alignments must be actual BAM context; CRAM/raw
    input routing stays deferred. This does not perform native BAM-header QC.
    """
    if not isinstance(plan, ReviewSourcePlan):
        raise InputValidationError('Review resolution requires ReviewSourcePlan.')
    fresh = read_review_source_plan(plan.path, sample_ids=context.sample_ids,
                                   reference_build=context.reference.build, config=config)
    if fresh != plan:
        raise InputValidationError('Review source plan or its resolved intent changed.')
    generated = any(item.source.origin == 'GENERATED' for item in plan.rows)
    root = attempt = None
    if generated:
        if run_root is None or attempt_dir is None:
            raise InputValidationError('Generated review sources require an exact run_root and attempt_dir.')
        root, attempt = Path(run_root), Path(attempt_dir)
        if (not root.is_absolute() or not attempt.is_absolute() or '..' in root.parts
                or '..' in attempt.parts or attempt.parent != root/'attempts' or not attempt.is_dir()):
            raise InputValidationError('Generated review attempt must belong to this explicit owned run root.')
    alignments = {item.sample.sample_id: item.sample for item in context.samples}
    lengths = {item.name: item.length for item in context.reference.contigs}
    targets, bindings = [], []
    seen = set()
    for request in plan.rows:
        row = request.source
        sample = alignments.get(row.sample_id)
        if sample is None or sample.input.input_type is not InputType.BAM or len(sample.input.files) != 1:
            raise InputValidationError('Review source resolution needs this sample actual indexed BAM context; raw/CRAM deferred.')
        path, index, tool, version, bundle = row.path, row.index_path, row.source_tool, row.source_version, row.bnd_bundle_manifest
        linkage = None
        if row.origin == 'GENERATED':
            from hifivar.downstream_handoff import _generated_source
            path, index, tool, version, linkage, bundle = _generated_source(root, attempt, row, context.reference)
        assert path is not None and index is not None
        from hifivar.launcher import _ready_vcf
        if not str(path).lower().endswith('.vcf.gz') or index != Path(str(path)+'.tbi'):
            raise InputValidationError('Review sources require a BGZF VCF and its explicit conventional TBI.')
        _ready_vcf(path, row.sample_id, context)
        before = {'vcf': artifact_identity(path), 'index': artifact_identity(index)}
        stat = path.stat();key = (row.sample_id, row.track, stat.st_dev, stat.st_ino)
        if key in seen:
            raise InputValidationError('Review source rows alias one physical source; combine candidate lists.')
        seen.add(key)
        if row.track == 'tr':
            if tool.casefold() != 'trgt':
                raise InputValidationError('Explicit TRID review is reserved for declared TRGT sources.')
            records = read_selected_repeat_records(path, row.variant_ids, reference=context.reference, sample_id=row.sample_id)
        else:
            explicit_ids=tuple(item.variant_id for item in request.candidates if item.site is None)
            named=read_selected_junction_records(path,explicit_ids,reference=context.reference,
                                                sample_id=row.sample_id,bundle_manifest=bundle) if explicit_ids else ()
            by_key={item.record_id:item for item in named}
            if any(item.site is not None for item in request.candidates):
                if bundle is not None:
                    raise InputValidationError('Exact-site review does not detach a derived BND bundle; select its BND IDs.')
                by_key.update(_selected_sites(path,request.candidates,context.reference,row.sample_id))
            if bundle is not None and any(item.variant_type is not VariantClass.BND for item in request.candidates):
                raise InputValidationError('Derived BND review targets must retain the complete BND target binding.')
            records=tuple(by_key[item.variant_id] for item in request.candidates)
        if len({record.ordinal for record in records}) != len(records):
            raise InputValidationError('Review candidate modes select one source record more than once.')
        fields = _fields_for_records(path, records)
        evidence = EvidenceClass.EXPLICIT if row.origin == 'EXTERNAL' else (
            EvidenceClass.SMALL_VARIANT if row.track == 'small' else EvidenceClass.TANDEM_REPEAT if row.track == 'tr' else EvidenceClass.READ_SV)
        if row.origin == 'GENERATED' and row.producer.startswith('assembly_sv.'):
            evidence = EvidenceClass.ASSEMBLY_SV
        elif row.origin == 'GENERATED' and row.producer == 'sv.harmonization':
            evidence = EvidenceClass.HARMONIZED_SV
        elif row.origin == 'GENERATED' and row.producer == 'phasing.hiphase':
            evidence = EvidenceClass.PHASED_VARIANT
        for candidate, record in zip(request.candidates, records, strict=True):
            contig, start, end = _coordinates(candidate, record, fields[record.ordinal], lengths)
            targets.append(ReviewTarget(candidate.review_id, row.sample_id, candidate.variant_id,
                candidate.variant_type, contig, start, end, path, tool, evidence,
                sample.input.files[0], context.reference.fasta, Path(output_directory),
                flank_bp=config.get('review', {}).get('flank_bp',500),
                bnd_bundle_manifest=bundle if candidate.variant_type is VariantClass.BND else None,
                reference_build=context.reference.build))
        after = {'vcf': artifact_identity(path), 'index': artifact_identity(index)}
        if before != after or any(record.source_sha256 != before['vcf']['sha256'] for record in records):
            raise OutputValidationError('Review source/index changed while its records were resolved.')
        bindings.append({'sample_id':row.sample_id,'track':row.track,'origin':row.origin,
            'alignment_path':str(sample.input.files[0]),
            'alignment_validation_scope':'CONTEXT_PATH_AND_INDEX_PRESENCE_NOT_NATIVE_HEADER_OR_ATTEMPT_ALIGNMENT_QC',
            'source_tool':tool,'source_version':version,
            'version_source':'SELECTED_ATTEMPT_PRODUCER_RECEIPT' if linkage else 'USER_DECLARED_NOT_REPROBED',
            'source_identity':before,'source_receipt':linkage,
            'selected_records':[record.to_dict() for record in records],
            'candidate_review_ids':[item.review_id for item in request.candidates],
            'catalog_id':row.catalog_id,'biological_validity':'NOT_ASSESSED'})
    if artifact_identity(plan.path)['sha256'] != plan.sha256:
        raise InputValidationError('Review source plan changed during resolution.')
    return ResolvedReviewTargets(plan, tuple(targets), tuple(bindings))


__all__ = ['ReviewSite','SelectedReviewSite','ReviewCandidate','ReviewSourceRequest','ReviewSourcePlan','ResolvedReviewTargets',
           'read_review_source_plan','resolve_review_targets']
