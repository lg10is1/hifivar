"""Variant-centered manual-review contracts and manifest serialization."""

from __future__ import annotations

import csv
import io
import json
import os
import re
from dataclasses import dataclass, field, replace
from enum import Enum
from pathlib import Path
from typing import Mapping, Sequence

from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.junction_consumers import SelectedJunctionRecord, read_selected_junction_records
from hifivar.reference import ReferenceGenome
from hifivar.serialization import (
    utc_now_iso8601,
    write_json_atomic,
    write_yaml_atomic,
)
from hifivar.validation import (
    read_fai_contigs,
    validate_alignment_file,
    validate_fasta,
    validate_file,
    validate_vcf,
)


_SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
_SELECTION_REQUIRED = (
    "review_id",
    "sample",
    "variant_id",
    "variant_type",
    "contig",
    "start",
    "end",
    "source_vcf",
    "source_caller",
    "evidence_class",
)
_SELECTION_OPTIONAL = (
    "mate_contig",
    "mate_position",
    "flank_bp",
    "trgt_visualization",
    "bnd_bundle_manifest",
)


class VariantClass(str, Enum):
    """Variant classes with explicit Phase 10 window semantics."""

    SNV = "SNV"
    INDEL = "INDEL"
    DEL = "DEL"
    DUP = "DUP"
    INV = "INV"
    INS = "INS"
    BND = "BND"
    TR = "TR"


class ReviewStatus(str, Enum):
    """Human observation only; never a truth or pathogenicity label."""

    NOT_REVIEWED = "NOT_REVIEWED"
    SUPPORT = "SUPPORT"
    NOT_SUPPORT = "NOT_SUPPORT"
    UNCERTAIN = "UNCERTAIN"


class EvidenceClass(str, Enum):
    """Origin of the immutable variant evidence under review."""

    SMALL_VARIANT = "small_variant"
    READ_SV = "read_sv"
    ASSEMBLY_SV = "assembly_sv"
    HARMONIZED_SV = "harmonized_sv"
    PHASED_VARIANT = "phased_variant"
    TANDEM_REPEAT = "tandem_repeat"
    EXPLICIT = "explicit"


@dataclass(frozen=True, slots=True)
class ReviewLocus:
    """One 1-based closed IGV review interval."""

    label: str
    contig: str
    variant_start: int
    variant_end: int
    window_start: int
    window_end: int

    def __post_init__(self) -> None:
        if not self.label or not self.contig:
            raise InputValidationError("Review locus label and contig must be non-empty.")
        for name in ("variant_start", "variant_end", "window_start", "window_end"):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value < 1:
                raise InputValidationError(f"Review locus {name} must be a positive integer.")
        if self.variant_end < self.variant_start:
            raise InputValidationError("Review locus variant_end precedes variant_start.")
        if self.window_start > self.variant_start or self.window_end < self.variant_end:
            raise InputValidationError("Review window must contain the variant locus.")

    @property
    def igv_locus(self) -> str:
        return f"{self.contig}:{self.window_start}-{self.window_end}"

    def to_dict(self) -> dict[str, object]:
        return {
            "label": self.label,
            "contig": self.contig,
            "variant_start": self.variant_start,
            "variant_end": self.variant_end,
            "window_start": self.window_start,
            "window_end": self.window_end,
            "igv_locus": self.igv_locus,
        }


@dataclass(frozen=True, slots=True)
class ReviewTarget:
    """One explicitly selected variant and its immutable source evidence."""

    review_id: str
    sample_id: str
    variant_id: str
    variant_class: VariantClass
    contig: str
    start: int
    end: int
    source_vcf: Path
    source_caller: str
    evidence_class: EvidenceClass
    alignment_path: Path
    reference_fasta: Path
    output_directory: Path
    flank_bp: int = 500
    mate_contig: str | None = None
    mate_position: int | None = None
    trgt_visualization_path: Path | None = None
    bnd_bundle_manifest: Path | None = None
    reference_build: str | None = None
    view_mode: str = "legacy"
    small_flank_bp: int = 40
    breakpoint_flank_bp: int = 250
    long_variant_bp: int = 1000
    _junction_record: SelectedJunctionRecord | None = field(default=None, init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not isinstance(self.view_mode, str) or self.view_mode not in {"legacy", "scientific"}:
            raise InputValidationError("Review view_mode must be legacy or scientific.")
        for name in ("small_flank_bp", "breakpoint_flank_bp", "long_variant_bp"):
            value = getattr(self, name)
            if type(value) is not int or not 1 <= value <= 100000:
                raise InputValidationError(f"Review {name} must be an integer in 1..100000.")
        for name in ("review_id", "sample_id"):
            value = getattr(self, name)
            if not isinstance(value, str) or _SAFE_ID.fullmatch(value) is None:
                raise InputValidationError(
                    f"Review {name} must be ASCII-safe: letters, digits, '.', '_' or '-'."
                )
        for name in ("variant_id", "contig", "source_caller"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip() or any(c in value for c in "\r\n"):
                raise InputValidationError(f"Review {name} must be non-empty and single-line.")
        if not isinstance(self.variant_class, VariantClass):
            raise InputValidationError("Review variant_class must be VariantClass.")
        if not isinstance(self.evidence_class, EvidenceClass):
            raise InputValidationError("Review evidence_class must be EvidenceClass.")
        for name in ("start", "end"):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value < (0 if self.variant_class is VariantClass.BND else 1):
                raise InputValidationError(f"Review {name} must be a positive integer or explicit virtual BND anchor.")
        if self.end < self.start:
            raise InputValidationError("Review end must not precede start.")
        if not isinstance(self.flank_bp, int) or isinstance(self.flank_bp, bool) or self.flank_bp < 0:
            raise InputValidationError("Review flank_bp must be a non-negative integer.")
        if self.variant_class is VariantClass.BND:
            if (self.mate_contig is None) != (self.mate_position is None):
                raise InputValidationError('BND mate declarations must supply both fields or neither.')
            if self.mate_position is not None and (type(self.mate_position) is not int or self.mate_position < 0):
                raise InputValidationError('BND mate position must be an explicit nonnegative integer.')
        elif self.mate_contig is not None or self.mate_position is not None:
            raise InputValidationError("Breakend mate coordinates are valid only for BND review.")
        if self.variant_class is not VariantClass.TR and self.trgt_visualization_path is not None:
            raise InputValidationError("TRGT visualization metadata is valid only for TR targets.")

        source = validate_vcf(self.source_vcf, require_index=str(self.source_vcf).lower().endswith(".vcf.gz"))
        alignment = validate_alignment_file(self.alignment_path, require_index=True)
        reference = validate_fasta(self.reference_fasta, require_fai=True)
        trgt = validate_file(self.trgt_visualization_path) if self.trgt_visualization_path else None
        contigs = read_fai_contigs(Path(f"{reference}.fai"))
        requested = {self.contig}
        if self.mate_contig:
            requested.add(self.mate_contig)
        missing = sorted(requested.difference(contigs))
        if missing:
            raise InputValidationError(
                f"REFERENCE_CONTIG_MISMATCH for review '{self.review_id}': {missing!r}."
            )
        if self.variant_class is VariantClass.BND:
            genome = ReferenceGenome.from_fasta(reference, build=self.reference_build or 'UNVERIFIED_REVIEW_BUILD')
            record = read_selected_junction_records(source, (self.variant_id,), reference=genome,
                sample_id=self.sample_id, bundle_manifest=self.bnd_bundle_manifest)[0]
            junction = record.junction
            if junction is None or (self.contig, self.start, self.end) != (record.contig, record.position, record.position):
                raise InputValidationError('BND review coordinates must match the exact selected VCF anchor, not END.')
            if self.mate_contig is not None and (self.mate_contig, self.mate_position) != (junction.mate_contig, junction.mate_position):
                raise InputValidationError('Declared BND mate differs from the exact selected ALT.')
            object.__setattr__(self, 'mate_contig', junction.mate_contig)
            object.__setattr__(self, 'mate_position', junction.mate_position)
            object.__setattr__(self, '_junction_record', record)
        elif self.bnd_bundle_manifest is not None:
            raise InputValidationError('Explicit BND bundle review requires a BND target.')
        object.__setattr__(self, "source_vcf", source)
        object.__setattr__(self, "alignment_path", alignment)
        object.__setattr__(self, "reference_fasta", reference)
        object.__setattr__(self, "output_directory", Path(self.output_directory).expanduser())
        object.__setattr__(self, "trgt_visualization_path", trgt)

    @property
    def loci(self) -> tuple[ReviewLocus, ...]:
        """Generate transparent windows; INS uses its anchor, BND uses two loci."""
        if self.variant_class is VariantClass.BND:
            record = self._junction_record
            assert record is not None and record.junction is not None
            junction = record.junction
            lengths = read_fai_contigs(Path(str(self.reference_fasta) + '.fai'))
            result = []
            for label, contig, position, virtual in (
                ('primary', junction.contig, junction.position, junction.local_virtual),
                ('mate', junction.mate_contig, junction.mate_position, junction.mate_virtual)):
                if contig is None or virtual:
                    continue
                flank = self.breakpoint_flank_bp if self.view_mode == 'scientific' else self.flank_bp
                result.append(ReviewLocus(label, contig, position, position,
                    max(1, position - flank), min(lengths[contig], position + flank)))
            return tuple(result)
        primary_end = self.start if self.variant_class in {VariantClass.SNV, VariantClass.INS} else self.end
        if self.view_mode == "scientific":
            length = read_fai_contigs(Path(str(self.reference_fasta) + '.fai'))[self.contig]
            if self.end > length:
                raise InputValidationError('Review target extends beyond its reference contig.')
            def window(label, start, end, flank):
                return ReviewLocus(label, self.contig, start, end,
                    max(1, start-flank), min(length, end+flank))
            if self.variant_class in {VariantClass.SNV, VariantClass.INDEL, VariantClass.INS}:
                return (window('primary', self.start, primary_end, self.small_flank_bp),)
            overview = window('overview', self.start, primary_end, self.flank_bp)
            if primary_end-self.start+1 >= self.long_variant_bp:
                return (overview,
                        window('left_breakpoint', self.start, self.start, self.breakpoint_flank_bp),
                        window('right_breakpoint', primary_end, primary_end, self.breakpoint_flank_bp))
            return (overview,)
        return (_make_locus('primary', self.contig, self.start, primary_end, self.flank_bp),)

    @property
    def unrendered_bnd_endpoints(self) -> tuple[dict, ...]:
        if self._junction_record is None:
            return ()
        junction = self._junction_record.junction
        return tuple({'endpoint': label, 'contig': contig, 'position': position,
            'status': 'UNKNOWN_SINGLE_BREAKEND' if contig is None else 'NO_REFERENCE_BASE_VIRTUAL_ANCHOR'}
            for label, contig, position, virtual in (
                ('LOCAL', junction.contig, junction.position, junction.local_virtual),
                ('ALT_MATE', junction.mate_contig, junction.mate_position, junction.mate_virtual))
            if contig is None or virtual)

    @property
    def screenshot_paths(self) -> tuple[Path, ...]:
        root = self.output_directory / "screenshots"
        return tuple(root / f"{self.review_id}.locus{index:02d}.png" for index, _ in enumerate(self.loci, 1))

    def to_dict(self) -> dict[str, object]:
        return {
            "review_id": self.review_id,
            "sample": self.sample_id,
            "variant_id": self.variant_id,
            "variant_type": self.variant_class.value,
            "contig": self.contig,
            "start": self.start,
            "end": self.end,
            "flank_bp": self.flank_bp,
            "mate_contig": self.mate_contig,
            "mate_position": self.mate_position,
            "junction_record": self._junction_record.to_dict() if self._junction_record else None,
            "unrendered_bnd_endpoints": list(self.unrendered_bnd_endpoints),
            "view_mode": self.view_mode,
            "loci": [locus.to_dict() for locus in self.loci],
            "source_vcf": str(self.source_vcf),
            "source_caller": self.source_caller,
            "evidence_class": self.evidence_class.value,
            "bam": str(self.alignment_path),
            "reference": str(self.reference_fasta),
            "screenshots": [str(path) for path in self.screenshot_paths],
            "trgt_visualization": str(self.trgt_visualization_path) if self.trgt_visualization_path else None,
        }


@dataclass(frozen=True, slots=True)
class ReviewEvidence:
    """Traceable IGV/TR evidence bundle without a scientific interpretation."""

    target: ReviewTarget
    batch_path: Path
    batch_command: tuple[str, ...]
    screenshots: tuple[Path, ...]
    generated: bool

    def to_dict(self) -> dict[str, object]:
        return {
            "target": self.target.to_dict(),
            "batch_path": str(self.batch_path),
            "batch_command": list(self.batch_command),
            "screenshots": [str(path) for path in self.screenshots],
            "generated": self.generated,
        }


@dataclass(frozen=True, slots=True)
class ReviewResult:
    """Manual-review row; status is explicitly not truth or pathogenicity."""

    target: ReviewTarget
    evidence: ReviewEvidence
    status: ReviewStatus = ReviewStatus.NOT_REVIEWED
    notes: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.status, ReviewStatus):
            raise InputValidationError("Review result status must be ReviewStatus.")
        if not isinstance(self.notes, str) or any(c in self.notes for c in "\r\n"):
            raise InputValidationError("Review notes must be a single-line string.")
        if self.evidence.target != self.target:
            raise InputValidationError("Review evidence does not match its target.")

    @property
    def is_truth(self) -> bool:
        return False

    def with_manual_review(self, status: ReviewStatus, notes: str = "") -> ReviewResult:
        return replace(self, status=status, notes=notes)

    def to_dict(self) -> dict[str, object]:
        return {
            "review_id": self.target.review_id,
            "sample": self.target.sample_id,
            "variant_id": self.target.variant_id,
            "variant_type": self.target.variant_class.value,
            "locus": [item.igv_locus for item in self.target.loci],
            "evidence_class": self.target.evidence_class.value,
            "source_caller": self.target.source_caller,
            "source_vcf": str(self.target.source_vcf),
            "source_bam": str(self.target.alignment_path),
            "source_reference": str(self.target.reference_fasta),
            "screenshots": [str(path) for path in self.evidence.screenshots],
            "trgt_visualization": str(self.target.trgt_visualization_path) if self.target.trgt_visualization_path else None,
            "status": self.status.value,
            "notes": self.notes,
            "status_is_truth": False,
            "junction_record": self.target._junction_record.to_dict() if self.target._junction_record else None,
            "unrendered_bnd_endpoints": list(self.target.unrendered_bnd_endpoints),
        }


@dataclass(frozen=True, slots=True)
class ReviewManifest:
    """Ordered JSON/YAML/TSV-friendly manual-review manifest."""

    results: tuple[ReviewResult, ...]
    created_at: str = ""
    schema_version: int = 1

    def __post_init__(self) -> None:
        ids = [result.target.review_id for result in self.results]
        if len(ids) != len(set(ids)):
            raise InputValidationError("Review manifest contains duplicate review IDs.")
        if not self.created_at:
            object.__setattr__(self, "created_at", utc_now_iso8601())

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": self.schema_version,
            "created_at": self.created_at,
            "scientific_policy": {
                "manual_status_is_truth": False,
                "manual_status_is_pathogenicity": False,
                "raw_variant_artifacts_modified": False,
            },
            "results": [result.to_dict() for result in self.results],
        }

    def write_json(self, path: str | Path, *, overwrite: bool = False) -> Path:
        return write_json_atomic(self.to_dict(), path, overwrite=overwrite, artifact_name="Review manifest")

    def write_yaml(self, path: str | Path, *, overwrite: bool = False) -> Path:
        return write_yaml_atomic(self.to_dict(), path, overwrite=overwrite, artifact_name="Review manifest")

    def write_tsv(self, path: str | Path, *, overwrite: bool = False) -> Path:
        fields = (
            "review_id", "sample", "variant_id", "variant_type", "locus",
            "evidence_class", "source_caller", "source_vcf", "source_bam",
            "source_reference", "screenshot", "trgt_visualization", "status", "notes",
            "junction_json", "unrendered_bnd_endpoints_json",
        )
        stream = io.StringIO(newline="")
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for result in self.results:
            row = result.to_dict()
            writer.writerow({
                "review_id": row["review_id"], "sample": row["sample"],
                "variant_id": row["variant_id"], "variant_type": row["variant_type"],
                "locus": ";".join(row["locus"]), "evidence_class": row["evidence_class"],
                "source_caller": row["source_caller"], "source_vcf": row["source_vcf"],
                "source_bam": row["source_bam"], "source_reference": row["source_reference"],
                "screenshot": ";".join(row["screenshots"]),
                "trgt_visualization": row["trgt_visualization"] or "",
                "status": row["status"], "notes": row["notes"],
                "junction_json": json.dumps(row["junction_record"], sort_keys=True),
                "unrendered_bnd_endpoints_json": json.dumps(row["unrendered_bnd_endpoints"], sort_keys=True),
            })
        return _write_text_atomic(Path(path), stream.getvalue(), overwrite=overwrite, label="Review TSV manifest")


def read_review_selection(
    path: str | Path,
    *,
    alignments: Mapping[str, Path],
    reference_fasta: Path,
    output_directory: Path,
    default_flank_bp: int = 500,
    reference_build: str | None = None,
) -> tuple[ReviewTarget, ...]:
    """Read an explicit TSV selection without inferring confidence or truth."""
    selection = Path(path).expanduser()
    rows = read_review_selection_rows(selection)
    base = selection.parent
    targets: list[ReviewTarget] = []
    for line_number, row in enumerate(rows, start=2):
        sample = (row.get("sample") or "").strip()
        if sample not in alignments:
            raise InputValidationError(
                f"Review selection line {line_number} has no alignment for sample '{sample}'."
            )
        try:
            start = int(row["start"])
            end = int(row["end"])
            flank_text = (row.get("flank_bp") or "").strip()
            flank = int(flank_text) if flank_text else default_flank_bp
            mate_text = (row.get("mate_position") or "").strip()
            mate_position = int(mate_text) if mate_text else None
            variant_class = VariantClass(row["variant_type"].strip().upper())
            evidence_class = EvidenceClass(row["evidence_class"].strip().lower())
        except (ValueError, TypeError, KeyError) as error:
            raise InputValidationError(
                f"Review selection line {line_number} has invalid enum/integer fields: {error}."
            ) from error
        target = ReviewTarget(
            review_id=row["review_id"].strip(), sample_id=sample,
            variant_id=row["variant_id"].strip(), variant_class=variant_class,
            contig=row["contig"].strip(), start=start, end=end,
            source_vcf=_relative_path(row["source_vcf"], base),
            source_caller=row["source_caller"].strip(), evidence_class=evidence_class,
            alignment_path=alignments[sample], reference_fasta=reference_fasta,
            output_directory=output_directory, flank_bp=flank,
            mate_contig=(row.get("mate_contig") or "").strip() or None,
            mate_position=mate_position,
            trgt_visualization_path=_optional_relative_path(row.get("trgt_visualization"), base),
            bnd_bundle_manifest=_optional_relative_path(row.get("bnd_bundle_manifest"), base),
            reference_build=reference_build,
        )
        targets.append(target)
    ids = [target.review_id for target in targets]
    if len(ids) != len(set(ids)):
        raise InputValidationError("Review selection contains duplicate review IDs.")
    return tuple(targets)


def read_review_selection_rows(path: str | Path) -> tuple[dict[str, str], ...]:
    """Read selection metadata without requiring not-yet-built DAG artifacts."""
    selection = validate_file(path, require_nonempty=False)
    try:
        with selection.open("rt", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = tuple(reader.fieldnames or ())
            missing = [name for name in _SELECTION_REQUIRED if name not in fields]
            unknown = sorted(set(fields).difference((*_SELECTION_REQUIRED, *_SELECTION_OPTIONAL)))
            if missing or unknown:
                raise InputValidationError(
                    f"Review selection columns invalid; missing={missing!r}, unknown={unknown!r}."
                )
            rows = tuple(dict(row) for row in reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise InputValidationError(f"Unable to read review selection '{selection}': {error}") from error

    return rows


def _make_locus(label: str, contig: str, start: int, end: int, flank: int) -> ReviewLocus:
    return ReviewLocus(label, contig, start, end, max(1, start - flank), end + flank)


def _relative_path(value: str, base: Path) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise InputValidationError("Review selection path must be non-empty.")
    path = Path(value).expanduser()
    return path if path.is_absolute() else base / path


def _optional_relative_path(value: str | None, base: Path) -> Path | None:
    return _relative_path(value, base) if isinstance(value, str) and value.strip() else None


def _write_text_atomic(path: Path, text: str, *, overwrite: bool, label: str) -> Path:
    if path.exists() and not overwrite:
        raise OutputValidationError(f"{label} already exists: '{path}'.")
    temporary = path.with_name(f".{path.name}.hifivar.tmp")
    if temporary.exists():
        if not overwrite:
            raise OutputValidationError(f"{label} temporary file already exists: '{temporary}'.")
        temporary.unlink()
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary.write_text(text, encoding="utf-8", newline="\n")
        os.replace(temporary, path)
    except OSError as error:
        temporary.unlink(missing_ok=True)
        raise OutputValidationError(f"Unable to write {label} '{path}': {error}") from error
    return path


__all__ = [
    "EvidenceClass", "ReviewEvidence", "ReviewLocus", "ReviewManifest",
    "ReviewResult", "ReviewStatus", "ReviewTarget", "VariantClass",
    "read_review_selection", "read_review_selection_rows",
]
