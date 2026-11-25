"""Tool-neutral Phase 11 annotation and region-overlap contracts."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import os
import re
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from types import MappingProxyType
from typing import Mapping, Sequence

from hifivar.benchmark_inputs import _contig_attributes, _explicit_builds
from hifivar.exceptions import InputValidationError, OutputValidationError, ReferenceError
from hifivar.junction_consumers import SelectedJunctionRecord, read_selected_junction_records, revalidate_junction_record
from hifivar.reference import ReferenceGenome
from hifivar.repeat_selection import SelectedRepeatRecord, read_selected_repeat_records, revalidate_repeat_records
from hifivar.validation import validate_file, validate_output_file, validate_vcf


_SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
_INPUT_FIELDS = (
    "sample", "variant_category", "source_vcf", "source_tool",
    "source_variant_ids",
)


class VariantCategory(str, Enum):
    """Independent input classes retained throughout annotation."""

    SMALL = "small"
    SV = "sv"
    TR = "tr"


class AnnotationSource(str, Enum):
    """Independent annotation sources; values are never merged by overwrite."""

    ANNOVAR = "annovar"
    VEP = "vep"
    REGION_OVERLAP = "region_overlap"


class AnnotationRunStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


class RegionCategory(str, Enum):
    GENE = "gene"
    EXON = "exon"
    REGULATORY = "regulatory"
    REPEAT = "repeat"
    SEGDUP = "segdup"


@dataclass(frozen=True, slots=True)
class AnnotationInput:
    """One immutable caller VCF handed to an annotation source."""

    sample_id: str
    variant_category: VariantCategory
    source_vcf: Path
    source_tool: str
    reference: ReferenceGenome
    source_variant_ids: tuple[str, ...] = ()
    bnd_bundle_manifest: Path | None = None
    source_identity: Mapping[str, object] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not isinstance(self.sample_id, str) or _SAFE_ID.fullmatch(self.sample_id) is None:
            raise InputValidationError("Annotation sample_id must be ASCII-safe.")
        if not isinstance(self.variant_category, VariantCategory):
            raise InputValidationError("Annotation variant_category must be VariantCategory.")
        if not isinstance(self.source_tool, str) or not self.source_tool.strip():
            raise InputValidationError("Annotation source_tool must be non-empty.")
        if not isinstance(self.reference, ReferenceGenome) or self.reference.build is None:
            raise InputValidationError("Annotation requires a reference with an explicit build.")
        identity = audit_annotation_vcf_identity(
            self.source_vcf, self.sample_id, self.reference,
        )
        source = Path(str(identity["vcf"]))
        if self.bnd_bundle_manifest is not None and self.variant_category is not VariantCategory.SV:
            raise InputValidationError('BND bundle annotation requires the independent SV category.')
        ids = tuple(self.source_variant_ids)
        if any(not isinstance(item, str) or not item.strip() for item in ids):
            raise InputValidationError("Annotation source variant IDs must be non-empty strings.")
        if len(ids) != len(set(ids)):
            raise InputValidationError("Annotation source variant IDs must be unique.")
        if bool(identity.get('explicit_bnd_derived_header')) != (self.bnd_bundle_manifest is not None):
            raise InputValidationError('BND derived annotation requires its explicit bundle manifest; native inputs must not attach a derived role.')
        if self.bnd_bundle_manifest is not None:
            from hifivar.junction_consumers import validate_bnd_consumer_bundle
            manifest_path = validate_file(self.bnd_bundle_manifest).absolute()
            identity['bnd_bundle_binding'] = validate_bnd_consumer_bundle(manifest_path, source, self.reference, self.sample_id)
            object.__setattr__(self, 'bnd_bundle_manifest', manifest_path)
        object.__setattr__(self, "source_vcf", source)
        object.__setattr__(self, "source_variant_ids", ids)
        object.__setattr__(self, "source_identity", MappingProxyType(identity))

    def to_dict(self) -> dict[str, object]:
        return {
            "sample": self.sample_id,
            "variant_category": self.variant_category.value,
            "source_vcf": str(self.source_vcf),
            "source_tool": self.source_tool,
            "source_variant_ids": list(self.source_variant_ids),
            "annovar_vep_input_scope": "whole_source_vcf",
            "source_variant_ids_role": "provenance_and_region_overlap_selection_not_annovar_vep_filter",
            "source_identity": dict(self.source_identity),
            "bnd_bundle_manifest": str(self.bnd_bundle_manifest) if self.bnd_bundle_manifest else None,
            "reference": self.reference.to_dict(),
            "raw_source_modified": False,
        }


def audit_annotation_vcf_identity(
    vcf: str | Path,
    sample_id: str,
    reference: ReferenceGenome,
) -> dict[str, object]:
    """Validate one small/SV/TR source header without requiring INFO/SVTYPE.

    Every declared contig must have its exact reference name and positive,
    matching length. Subsets are allowed, but absent lengths fail closed.
    Explicit GRCh37/hg19 versus GRCh38/hg38 contradictions fail; absent or
    unrecognized build headers are recorded as unverified, never inferred.
    This guard does not scan records or establish REF-allele sequence identity.
    It never rewrites, normalizes, indexes, or filters the source VCF.
    """
    if (
        not isinstance(sample_id, str)
        or _SAFE_ID.fullmatch(sample_id) is None
        or sample_id.casefold() in {"unknown", "unknownsamplename"}
    ):
        raise InputValidationError("Annotation requires an exact non-placeholder sample_id.")
    if not isinstance(reference, ReferenceGenome) or reference.build is None:
        raise InputValidationError("Annotation requires a reference with an explicit build.")
    path = validate_vcf(vcf, require_index=str(vcf).lower().endswith(".vcf.gz"))
    lengths: dict[str, int] = {}
    reference_headers: list[str] = []
    contig_assemblies: dict[str, str] = {}
    samples: tuple[str, ...] | None = None
    header_digest = hashlib.sha256()
    explicit_bnd_derived_header = False
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    try:
        with opener(path, "rt", encoding="utf-8", newline="") as handle:
            for raw in handle:
                header_digest.update(raw.encode("utf-8"))
                line = raw.rstrip("\r\n")
                if line.startswith("##contig="):
                    attributes = _contig_attributes(line, path)
                    name = attributes.get("ID", "")
                    if not name or any(character.isspace() for character in name):
                        raise InputValidationError("Annotation VCF contains an invalid contig ID.")
                    if name in lengths:
                        raise InputValidationError(f"Duplicate VCF contig declaration '{name}'.")
                    length = attributes.get("length")
                    if length is None:
                        raise InputValidationError(
                            f"CONTIG_LENGTH_UNVERIFIED: annotation VCF '{name}' lacks a length."
                        )
                    if not length.isascii() or not length.isdecimal() or int(length) <= 0:
                        raise InputValidationError(f"Invalid annotation VCF contig length for '{name}'.")
                    lengths[name] = int(length)
                    if attributes.get("assembly"):
                        contig_assemblies[name] = attributes["assembly"]
                elif line.startswith('##HIFIVAR_BND_'):
                    explicit_bnd_derived_header = True
                elif line.startswith("##reference="):
                    reference_headers.append(line.partition("=")[2])
                elif line == "#CHROM" or line.startswith("#CHROM\t"):
                    columns = line.split("\t")
                    expected = ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
                    if tuple(columns[:9]) != expected:
                        raise InputValidationError(
                            "SAMPLE_IDENTITY_UNVERIFIED: annotation requires a genotyped VCF header."
                        )
                    samples = tuple(columns[9:])
                    break
    except (OSError, EOFError, UnicodeError) as error:
        raise InputValidationError(f"Unable to read annotation VCF header '{path}': {error}") from error
    if not samples or any(not name or any(character.isspace() for character in name) for name in samples):
        raise InputValidationError("Annotation VCF has missing or empty sample names.")
    if len(samples) != len(set(samples)):
        raise InputValidationError("Annotation VCF has duplicate sample names.")
    if samples != (sample_id,):
        raise InputValidationError(
            f"ANNOTATION_SAMPLE_MISMATCH: expected exactly {(sample_id,)!r}, observed {samples!r}."
        )
    if not lengths:
        raise InputValidationError("CONTIG_LENGTH_UNVERIFIED: annotation VCF lacks contig declarations.")
    reference_lengths = {item.name: item.length for item in reference.contigs}
    for name, length in lengths.items():
        if name not in reference_lengths:
            raise ReferenceError(f"REFERENCE_CONTIG_MISMATCH: annotation contig '{name}' is absent from reference.")
        if length != reference_lengths[name]:
            raise ReferenceError(
                f"CONTIG_LENGTH_MISMATCH: annotation VCF '{name}' length {length} != {reference_lengths[name]}."
            )
    vcf_builds = _explicit_builds(reference_headers + list(contig_assemblies.values()))
    declared_builds = _explicit_builds([reference.build])
    if len(vcf_builds) > 1 or (vcf_builds and declared_builds and vcf_builds != declared_builds):
        raise ReferenceError(
            f"REFERENCE_BUILD_MISMATCH: annotation VCF identifies {sorted(vcf_builds)!r}, "
            f"declared reference is {reference.build!r}."
        )
    return {
        "scope": "HEADER_IDENTITY_AND_EXISTING_INDEX_PRESENCE_ONLY",
        "vcf": str(path), "sample_check": "PASS_EXACT_SINGLE_SAMPLE",
        "sample_ids": samples, "contig_name_check": "PASS_EXACT_SUBSET",
        "contig_length_check": "PASS", "contig_lengths": lengths,
        "reference_headers": reference_headers, "contig_assemblies": contig_assemblies,
        "explicit_vcf_builds": sorted(vcf_builds),
        "declared_reference_build": reference.build,
        "reference_build_check": "PASS" if vcf_builds and declared_builds else "UNVERIFIED_MISSING_OR_UNRECOGNIZED_METADATA",
        "header_sha256": header_digest.hexdigest(),
        "record_validation": "NOT_RUN_HEADER_ONLY",
        "explicit_bnd_derived_header": explicit_bnd_derived_header,
        "reference_allele_sequence_check": "NOT_RUN",
        "input_modified": False,
    }


@dataclass(frozen=True, slots=True)
class AnnotationDatabase:
    """Versioned local annotation-data provenance."""

    name: str
    version: str
    path: Path
    reference_build: str
    checksum: str | None = None

    def __post_init__(self) -> None:
        for name in ("name", "version", "reference_build"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise InputValidationError(f"Annotation database {name} must be non-empty.")
        path = Path(self.path).expanduser()
        if not path.exists():
            raise InputValidationError(f"Annotation database path is missing: '{path}'.")
        object.__setattr__(self, "path", path)

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "version": self.version,
            "path": str(self.path),
            "reference_build": self.reference_build,
            "checksum": self.checksum,
        }


@dataclass(frozen=True, slots=True)
class AnnotationArtifact:
    """Validated derived annotation with immutable-source provenance."""

    input: AnnotationInput
    source: AnnotationSource
    output_path: Path
    output_format: str
    tool_version: str
    databases: tuple[AnnotationDatabase, ...]
    command: tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.source, AnnotationSource):
            raise InputValidationError("Annotation artifact source is invalid.")
        if self.output_format not in {"tsv", "native_vcf"}:
            raise InputValidationError("Annotation output format must be tsv or native_vcf.")
        if not isinstance(self.tool_version, str) or not self.tool_version.strip():
            raise InputValidationError("Annotation tool version must be non-empty.")
        output = validate_output_file(self.output_path)
        object.__setattr__(self, "output_path", output)
        object.__setattr__(self, "databases", tuple(self.databases))

    def to_dict(self) -> dict[str, object]:
        return {
            "input": self.input.to_dict(),
            "annotation_source": self.source.value,
            "output": str(self.output_path),
            "output_format": self.output_format,
            "tool_version": self.tool_version,
            "databases": [item.to_dict() for item in self.databases],
            "command": list(self.command) if self.command else None,
            "scientific_policy": {
                "functional_impact_is_call_confidence": False,
                "raw_source_modified": False,
            },
        }


@dataclass(frozen=True, slots=True)
class AnnotationResult:
    input: AnnotationInput
    source: AnnotationSource
    status: AnnotationRunStatus
    command: tuple[str, ...]
    tool_version: str | None = None
    duration_seconds: float = 0.0
    artifact: AnnotationArtifact | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "input": self.input.to_dict(),
            "annotation_source": self.source.value,
            "status": self.status.value,
            "command": list(self.command),
            "tool_version": self.tool_version,
            "duration_seconds": self.duration_seconds,
            "artifact": self.artifact.to_dict() if self.artifact else None,
        }


@dataclass(frozen=True, slots=True)
class VariantLocus:
    """One source variant for explicit SV/TR region-overlap annotation."""

    sample_id: str
    variant_id: str
    variant_category: VariantCategory
    contig: str
    start: int
    end: int
    source_vcf: Path
    source_tool: str

    junction_record: SelectedJunctionRecord | None = None
    repeat_record: SelectedRepeatRecord | None = None

    def __post_init__(self) -> None:
        if self.variant_category not in {VariantCategory.SV, VariantCategory.TR}:
            raise InputValidationError("Region overlap supports SV or TR loci only.")
        for name in ("sample_id", "variant_id", "contig", "source_tool"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise InputValidationError(f"Variant locus {name} must be non-empty.")
        record = self.junction_record
        repeat = self.repeat_record
        if repeat is not None and (not isinstance(repeat, SelectedRepeatRecord) or record is not None
                or self.variant_category is not VariantCategory.TR
                or (self.variant_id, self.contig, self.start, self.end) != (repeat.selector, repeat.contig, repeat.start, repeat.end)
                or Path(self.source_vcf).absolute() != repeat.source_vcf):
            raise InputValidationError('TRGT annotation locus differs from its explicit selected repeat record.')
        if record is not None and (not isinstance(record, SelectedJunctionRecord) or record.junction is None
                or self.variant_category is not VariantCategory.SV
                or (self.variant_id, self.contig, self.start, self.end) != (record.record_id, record.contig, record.position, record.position)
                or Path(self.source_vcf).absolute() != record.source_vcf):
            raise InputValidationError('BND annotation locus differs from its selected source record.')
        if type(self.start) is not int or type(self.end) is not int or self.start < (0 if record else 1) or self.end < self.start:
            raise InputValidationError("Variant locus requires positive ordered coordinates or an explicit BND anchor.")
        object.__setattr__(self, "source_vcf", validate_vcf(self.source_vcf))


@dataclass(frozen=True, slots=True)
class RegionDatabase:
    category: RegionCategory
    path: Path
    version: str
    reference_build: str

    def __post_init__(self) -> None:
        if not isinstance(self.category, RegionCategory):
            raise InputValidationError("Region database category is invalid.")
        if not isinstance(self.version, str) or not self.version.strip():
            raise InputValidationError("Region database version must be non-empty.")
        if not isinstance(self.reference_build, str) or not self.reference_build.strip():
            raise InputValidationError("Region database reference build must be non-empty.")
        object.__setattr__(self, "path", validate_file(self.path))


@dataclass(frozen=True, slots=True)
class RegionOverlap:
    variant: VariantLocus
    category: RegionCategory
    feature_id: str
    feature_start: int
    feature_end: int
    database_version: str
    endpoint: str = 'SPAN'
    query_contig: str | None = None
    query_position: int | None = None
    endpoint_side: str | None = None

    def to_dict(self) -> dict[str, object]:
        result = {
            "sample": self.variant.sample_id,
            "source_variant_id": self.variant.variant_id,
            "variant_category": self.variant.variant_category.value,
            "source_vcf": str(self.variant.source_vcf),
            "source_tool": self.variant.source_tool,
            "contig": self.query_contig or self.variant.contig,
            "start": self.query_position if self.query_position is not None else self.variant.start,
            "end": self.query_position if self.query_position is not None else self.variant.end,
            "region_category": self.category.value,
            "feature_id": self.feature_id,
            "feature_start": self.feature_start,
            "feature_end": self.feature_end,
            "database_version": self.database_version,
            "breakpoint_modified": False,
            "functional_impact_is_call_confidence": False,
            "endpoint": self.endpoint, "endpoint_side": self.endpoint_side or '',
            "coordinate_semantics": ('TRGT_REPEAT_1_BASED_INCLUSIVE_EXCLUDING_LEFT_ANCHOR' if self.variant.repeat_record else
                                     'VCF_ANCHOR_NOT_EVENT_SPAN' if self.variant.junction_record else 'VCF_DECLARED_SPAN'),
            "source_record_key": (self.variant.repeat_record.to_dict()['source_record_key'] if self.variant.repeat_record else
                                  self.variant.junction_record.to_dict()['source_record_key'] if self.variant.junction_record else ''),
            "junction_json": json.dumps(self.variant.junction_record.to_dict(), sort_keys=True) if self.variant.junction_record else '',
        }
        if self.variant.repeat_record:
            result['repeat_json'] = json.dumps(self.variant.repeat_record.to_dict(), sort_keys=True)
        return result


@dataclass(frozen=True, slots=True)
class RegionOverlapResult:
    overlaps: tuple[RegionOverlap, ...]
    output_tsv: Path
    databases: tuple[RegionDatabase, ...]
    junction_annotations: tuple[dict, ...] = ()
    repeat_annotations: tuple[dict, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "output_tsv": str(self.output_tsv),
            "databases": [
                {"category": item.category.value, "path": str(item.path),
                 "version": item.version, "reference_build": item.reference_build}
                for item in self.databases
            ],
            "overlaps": [item.to_dict() for item in self.overlaps],
            "junction_annotations": list(self.junction_annotations),
            "repeat_annotations": list(self.repeat_annotations),
            "junction_annotation_is_fusion_or_confidence": False,
        }


def annotate_region_overlaps(
    variants: Sequence[VariantLocus],
    databases: Sequence[RegionDatabase],
    *,
    reference: ReferenceGenome,
    output_tsv: str | Path,
    overwrite: bool = False,
) -> RegionOverlapResult:
    """Stream BED rows for explicit SV/TR loci without changing breakpoints."""
    ordered_variants = tuple(variants)
    ordered_databases = tuple(databases)
    if reference.build is None:
        raise InputValidationError("Region overlap requires an explicit reference build.")
    for database in ordered_databases:
        if database.reference_build != reference.build:
            raise InputValidationError(
                f"Region database build '{database.reference_build}' does not match '{reference.build}'."
            )
    contigs = set(reference.contig_names)
    if any(item.contig not in contigs for item in ordered_variants):
        raise InputValidationError("REFERENCE_CONTIG_MISMATCH in region-overlap variants.")

    destination = Path(output_tsv).expanduser().absolute()
    protected = [reference.fasta, reference.fai, *(v.source_vcf for v in ordered_variants), *(d.path for d in ordered_databases)]
    repeat_records = tuple(v.repeat_record for v in ordered_variants if v.repeat_record)
    revalidate_repeat_records(repeat_records)
    for variant in ordered_variants:
        if variant.junction_record:
            revalidate_junction_record(variant.junction_record)
            binding = variant.junction_record.bundle_binding
            if binding:
                protected.extend(Path(info['path']) for info in (binding['manifest'], *binding['artifacts'].values()))
    for path in protected:
        if destination.resolve() == path.resolve() or (destination.exists() and destination.samefile(path)):
            raise OutputValidationError('Region overlap output cannot replace an immutable source/reference/database/bundle input.')
    overlaps: list[RegionOverlap] = []
    ordinary = tuple(item for item in ordered_variants if item.junction_record is None)
    endpoints = []
    for variant in ordered_variants:
        if variant.junction_record is None:
            continue
        junction = variant.junction_record.junction
        for label, contig, position, side, virtual in (
            ('LOCAL', junction.contig, junction.position, junction.local_side, junction.local_virtual),
            ('ALT_MATE', junction.mate_contig, junction.mate_position, junction.mate_side, junction.mate_virtual),
        ):
            endpoints.append({'variant': variant, 'endpoint': label, 'contig': contig,
                'position': position, 'side': side, 'virtual': virtual, 'overlap_count': 0,
                'status': 'UNKNOWN_SINGLE_BREAKEND' if contig is None else 'NO_REFERENCE_BASE_VIRTUAL_ANCHOR' if virtual else 'NO_OVERLAP'})
    by_contig: dict[str, tuple[VariantLocus, ...]] = {}
    for contig in reference.contig_names:
        by_contig[contig] = tuple(item for item in ordinary if item.contig == contig)
    for database in ordered_databases:
        try:
            with database.path.open("rt", encoding="utf-8", newline="") as handle:
                for line_number, raw in enumerate(handle, 1):
                    if not raw.strip() or raw.startswith("#"):
                        continue
                    fields = raw.rstrip("\r\n").split("\t")
                    if len(fields) < 3:
                        raise InputValidationError(
                            f"Malformed BED row at {database.path}:{line_number}."
                        )
                    try:
                        start0, end0 = int(fields[1]), int(fields[2])
                    except ValueError as error:
                        raise InputValidationError(
                            f"Non-integer BED coordinate at {database.path}:{line_number}."
                        ) from error
                    if start0 < 0 or end0 <= start0:
                        raise InputValidationError(
                            f"Invalid BED interval at {database.path}:{line_number}."
                        )
                    feature = fields[3].strip() if len(fields) > 3 and fields[3].strip() else f"line-{line_number}"
                    for endpoint in endpoints:
                        if endpoint['contig'] != fields[0] or endpoint['virtual'] or endpoint['position'] is None:
                            continue
                        if start0 < endpoint['position'] <= end0:
                            overlaps.append(RegionOverlap(endpoint['variant'], database.category, feature,
                                start0 + 1, end0, database.version, endpoint['endpoint'], fields[0], endpoint['position'], endpoint['side']))
                            endpoint['overlap_count'] += 1
                            endpoint['status'] = 'OVERLAP'
                    for variant in by_contig.get(fields[0], ()):
                        if variant.start <= end0 and variant.end >= start0 + 1:
                            overlaps.append(
                                RegionOverlap(
                                    variant, database.category, feature,
                                    start0 + 1, end0, database.version,
                                )
                            )
        except (OSError, UnicodeError) as error:
            raise InputValidationError(f"Unable to read region database '{database.path}': {error}") from error

    destination = Path(output_tsv).expanduser()
    fields = (
        "sample", "source_variant_id", "variant_category", "source_vcf",
        "source_tool", "contig", "start", "end", "region_category",
        "feature_id", "feature_start", "feature_end", "database_version",
        "breakpoint_modified", "functional_impact_is_call_confidence",
        "endpoint", "endpoint_side", "coordinate_semantics", "source_record_key", "junction_json",
    )
    has_repeats = any(v.repeat_record is not None for v in ordered_variants)
    if has_repeats:
        fields = (*fields, 'repeat_json')
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for overlap in overlaps:
        row = overlap.to_dict()
        if has_repeats:
            row.setdefault('repeat_json', '')
        writer.writerow(row)
    revalidate_repeat_records(repeat_records)
    for variant in ordered_variants:
        if variant.junction_record:
            revalidate_junction_record(variant.junction_record)
    _write_text_atomic(destination, stream.getvalue(), overwrite=overwrite)
    junction_annotations = tuple({**variant.junction_record.to_dict(), 'endpoints': [
        {key: value for key, value in endpoint.items() if key != 'variant'}
        for endpoint in endpoints if endpoint['variant'] is variant]} for variant in ordered_variants if variant.junction_record)
    repeat_annotations = tuple({**v.repeat_record.to_dict(),
        'overlap_count': sum(item.variant is v for item in overlaps)} for v in ordered_variants if v.repeat_record)
    return RegionOverlapResult(tuple(overlaps), destination, ordered_databases, junction_annotations, repeat_annotations)


def read_annotation_input_rows(path: str | Path) -> tuple[dict[str, str], ...]:
    """Read metadata for DAG construction without validating future VCF outputs."""
    source = validate_file(path)
    try:
        with source.open("rt", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = tuple(reader.fieldnames or ())
            if fields not in (_INPUT_FIELDS, (*_INPUT_FIELDS, 'bnd_bundle_manifest')):
                raise InputValidationError(
                    f"Annotation input columns must be exactly {_INPUT_FIELDS!r}."
                )
            rows = tuple(dict(row) for row in reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise InputValidationError(f"Unable to read annotation inputs '{source}': {error}") from error
    keys: list[tuple[str, str]] = []
    for line_number, row in enumerate(rows, 2):
        if set(row) not in (set(_INPUT_FIELDS), set((*_INPUT_FIELDS, 'bnd_bundle_manifest'))) or any(value is None for value in row.values()):
            raise InputValidationError(f"Annotation input row {line_number} has an invalid column count.")
        try:
            sample = row["sample"].strip()
            category = VariantCategory(row["variant_category"].strip().lower())
            source_vcf = row["source_vcf"].strip()
            source_tool = row["source_tool"].strip()
        except (AttributeError, KeyError, ValueError) as error:
            raise InputValidationError(f"Invalid annotation input row {line_number}: {error}.") from error
        if not sample or not source_vcf or not source_tool:
            raise InputValidationError(f"Annotation input row {line_number} has blank required values.")
        keys.append((sample, category.value))
    if len(keys) != len(set(keys)):
        raise InputValidationError("Annotation input manifest duplicates a sample/category pair.")
    if not rows:
        raise InputValidationError("Annotation input manifest must contain at least one row.")
    return rows


def read_annotation_inputs(
    path: str | Path,
    *,
    reference: ReferenceGenome,
    sample_id: str | None = None,
    variant_category: VariantCategory | str | None = None,
) -> tuple[AnnotationInput, ...]:
    """Validate only selected metadata rows, or all rows when filters are absent.

    Filtering occurs before constructing any ``AnnotationInput``: an unrelated
    future VCF may be absent while the current row is ready. The complete plan
    must still have valid, unique sample/category metadata. Filters must match
    at least one row; they never turn an empty or misspelled plan into success.
    """
    source = Path(path).expanduser()
    rows = read_annotation_input_rows(source)
    if sample_id is not None and (
        not isinstance(sample_id, str) or _SAFE_ID.fullmatch(sample_id) is None
    ):
        raise InputValidationError("Annotation sample_id filter must be ASCII-safe.")
    if variant_category is not None:
        try:
            variant_category = VariantCategory(variant_category)
        except (TypeError, ValueError) as error:
            raise InputValidationError("Invalid annotation variant_category filter.") from error
    rows = tuple(
        row for row in rows
        if (sample_id is None or row["sample"].strip() == sample_id)
        and (variant_category is None or row["variant_category"].strip().lower() == variant_category.value)
    )
    if not rows:
        raise InputValidationError("Annotation selection matches no sample/category row.")
    inputs: list[AnnotationInput] = []
    for row in rows:
        vcf = Path(row["source_vcf"].strip()).expanduser()
        if not vcf.is_absolute():
            vcf = source.parent / vcf
        ids = tuple(
            value.strip() for value in row["source_variant_ids"].split(";")
            if value.strip()
        )
        inputs.append(
            AnnotationInput(
                row["sample"].strip(),
                VariantCategory(row["variant_category"].strip().lower()),
                vcf,
                row["source_tool"].strip(),
                reference,
                ids,
                ((source.parent / Path(row['bnd_bundle_manifest'].strip())).absolute()
                 if row.get('bnd_bundle_manifest', '').strip() else None),
            )
        )
    return tuple(inputs)


def read_selected_variant_loci(annotation_input: AnnotationInput) -> tuple[VariantLocus, ...]:
    """Stream a VCF and materialize only explicitly named SV/TR records."""
    if annotation_input.variant_category not in {VariantCategory.SV, VariantCategory.TR}:
        raise InputValidationError("Region-overlap selection requires an SV or TR input.")
    selected = set(annotation_input.source_variant_ids)
    if not selected:
        raise InputValidationError(
            "Region-overlap annotation requires explicit source_variant_ids; full-VCF loading is disabled."
        )
    if annotation_input.variant_category is VariantCategory.TR and any(x.startswith('TRID:') for x in annotation_input.source_variant_ids):
        if annotation_input.source_tool.casefold() != 'trgt':
            raise InputValidationError('Explicit TRID selection is reserved for declared TRGT inputs in the independent TR category.')
        records = read_selected_repeat_records(annotation_input.source_vcf, annotation_input.source_variant_ids,
            reference=annotation_input.reference, sample_id=annotation_input.sample_id)
        return tuple(VariantLocus(annotation_input.sample_id, record.selector, VariantCategory.TR,
            record.contig, record.start, record.end, annotation_input.source_vcf, annotation_input.source_tool,
            repeat_record=record) for record in records)
    records = read_selected_junction_records(annotation_input.source_vcf, annotation_input.source_variant_ids,
        reference=annotation_input.reference, sample_id=annotation_input.sample_id,
        bundle_manifest=annotation_input.bnd_bundle_manifest)
    return tuple(VariantLocus(annotation_input.sample_id, record.record_id,
        annotation_input.variant_category, record.contig, record.position, record.end,
        annotation_input.source_vcf, annotation_input.source_tool,
        record if record.junction else None) for record in records)


def _write_text_atomic(path: Path, text: str, *, overwrite: bool) -> Path:
    if path.exists() and not overwrite:
        raise OutputValidationError(f"Annotation output already exists: '{path}'.")
    temporary = path.with_name(f".{path.name}.hifivar.tmp")
    if temporary.exists() and not overwrite:
        raise OutputValidationError(f"Annotation temporary output exists: '{temporary}'.")
    temporary.unlink(missing_ok=True)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary.write_text(text, encoding="utf-8", newline="\n")
        os.replace(temporary, path)
    except OSError as error:
        temporary.unlink(missing_ok=True)
        raise OutputValidationError(f"Unable to write annotation output '{path}': {error}") from error
    return path


__all__ = [
    "AnnotationArtifact", "AnnotationDatabase", "AnnotationInput",
    "AnnotationResult", "AnnotationRunStatus", "AnnotationSource",
    "RegionCategory", "RegionDatabase", "RegionOverlap", "RegionOverlapResult",
    "VariantCategory", "VariantLocus", "annotate_region_overlaps",
    "audit_annotation_vcf_identity",
    "read_annotation_input_rows", "read_annotation_inputs",
    "read_selected_variant_loci",
]
