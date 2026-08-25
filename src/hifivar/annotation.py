"""Tool-neutral Phase 11 annotation and region-overlap contracts."""

from __future__ import annotations

import csv
import gzip
import io
import os
import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Mapping, Sequence

from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.reference import ReferenceGenome
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

    def __post_init__(self) -> None:
        if not isinstance(self.sample_id, str) or _SAFE_ID.fullmatch(self.sample_id) is None:
            raise InputValidationError("Annotation sample_id must be ASCII-safe.")
        if not isinstance(self.variant_category, VariantCategory):
            raise InputValidationError("Annotation variant_category must be VariantCategory.")
        if not isinstance(self.source_tool, str) or not self.source_tool.strip():
            raise InputValidationError("Annotation source_tool must be non-empty.")
        if not isinstance(self.reference, ReferenceGenome) or self.reference.build is None:
            raise InputValidationError("Annotation requires a reference with an explicit build.")
        source = validate_vcf(
            self.source_vcf,
            require_index=str(self.source_vcf).lower().endswith(".vcf.gz"),
        )
        ids = tuple(self.source_variant_ids)
        if any(not isinstance(item, str) or not item.strip() for item in ids):
            raise InputValidationError("Annotation source variant IDs must be non-empty strings.")
        if len(ids) != len(set(ids)):
            raise InputValidationError("Annotation source variant IDs must be unique.")
        object.__setattr__(self, "source_vcf", source)
        object.__setattr__(self, "source_variant_ids", ids)

    def to_dict(self) -> dict[str, object]:
        return {
            "sample": self.sample_id,
            "variant_category": self.variant_category.value,
            "source_vcf": str(self.source_vcf),
            "source_tool": self.source_tool,
            "source_variant_ids": list(self.source_variant_ids),
            "reference": self.reference.to_dict(),
            "raw_source_modified": False,
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

    def __post_init__(self) -> None:
        if self.variant_category not in {VariantCategory.SV, VariantCategory.TR}:
            raise InputValidationError("Region overlap supports SV or TR loci only.")
        for name in ("sample_id", "variant_id", "contig", "source_tool"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise InputValidationError(f"Variant locus {name} must be non-empty.")
        if not isinstance(self.start, int) or not isinstance(self.end, int) or self.start < 1 or self.end < self.start:
            raise InputValidationError("Variant locus requires positive ordered coordinates.")
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

    def to_dict(self) -> dict[str, object]:
        return {
            "sample": self.variant.sample_id,
            "source_variant_id": self.variant.variant_id,
            "variant_category": self.variant.variant_category.value,
            "source_vcf": str(self.variant.source_vcf),
            "source_tool": self.variant.source_tool,
            "contig": self.variant.contig,
            "start": self.variant.start,
            "end": self.variant.end,
            "region_category": self.category.value,
            "feature_id": self.feature_id,
            "feature_start": self.feature_start,
            "feature_end": self.feature_end,
            "database_version": self.database_version,
            "breakpoint_modified": False,
            "functional_impact_is_call_confidence": False,
        }


@dataclass(frozen=True, slots=True)
class RegionOverlapResult:
    overlaps: tuple[RegionOverlap, ...]
    output_tsv: Path
    databases: tuple[RegionDatabase, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "output_tsv": str(self.output_tsv),
            "databases": [
                {"category": item.category.value, "path": str(item.path),
                 "version": item.version, "reference_build": item.reference_build}
                for item in self.databases
            ],
            "overlaps": [item.to_dict() for item in self.overlaps],
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

    overlaps: list[RegionOverlap] = []
    by_contig: dict[str, tuple[VariantLocus, ...]] = {}
    for contig in reference.contig_names:
        by_contig[contig] = tuple(item for item in ordered_variants if item.contig == contig)
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
    )
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for overlap in overlaps:
        writer.writerow(overlap.to_dict())
    _write_text_atomic(destination, stream.getvalue(), overwrite=overwrite)
    return RegionOverlapResult(tuple(overlaps), destination, ordered_databases)


def read_annotation_input_rows(path: str | Path) -> tuple[dict[str, str], ...]:
    """Read metadata for DAG construction without validating future VCF outputs."""
    source = validate_file(path)
    try:
        with source.open("rt", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = tuple(reader.fieldnames or ())
            if fields != _INPUT_FIELDS:
                raise InputValidationError(
                    f"Annotation input columns must be exactly {_INPUT_FIELDS!r}."
                )
            rows = tuple(dict(row) for row in reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise InputValidationError(f"Unable to read annotation inputs '{source}': {error}") from error
    keys: list[tuple[str, str]] = []
    for line_number, row in enumerate(rows, 2):
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
    return rows


def read_annotation_inputs(
    path: str | Path,
    *,
    reference: ReferenceGenome,
) -> tuple[AnnotationInput, ...]:
    """Construct validated independent small/SV/TR annotation inputs."""
    source = Path(path).expanduser()
    rows = read_annotation_input_rows(source)
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
    opener = gzip.open if str(annotation_input.source_vcf).lower().endswith(".gz") else open
    loci: list[VariantLocus] = []
    try:
        with opener(annotation_input.source_vcf, "rt", encoding="utf-8", newline="") as handle:
            for raw in handle:
                if raw.startswith("#") or not raw.strip():
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) < 8 or fields[2] not in selected:
                    continue
                try:
                    position = int(fields[1])
                    info = {
                        item.split("=", 1)[0]: item.split("=", 1)[1]
                        for item in fields[7].split(";") if "=" in item
                    }
                    end = int(info.get("END", position))
                except ValueError as error:
                    raise InputValidationError(
                        f"Invalid selected VCF coordinates for variant '{fields[2]}'."
                    ) from error
                loci.append(
                    VariantLocus(
                        annotation_input.sample_id, fields[2],
                        annotation_input.variant_category, fields[0], position, end,
                        annotation_input.source_vcf, annotation_input.source_tool,
                    )
                )
    except (OSError, UnicodeError) as error:
        raise InputValidationError(
            f"Unable to read annotation VCF '{annotation_input.source_vcf}': {error}"
        ) from error
    found = {item.variant_id for item in loci}
    missing = sorted(selected - found)
    if missing:
        raise InputValidationError(f"Selected source variant IDs are absent from VCF: {missing!r}.")
    return tuple(loci)


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
    "read_annotation_input_rows", "read_annotation_inputs",
    "read_selected_variant_loci",
]
