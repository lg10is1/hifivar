"""Tandem-repeat catalog and result models for TRGT."""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path

from hifivar import validation
from hifivar.exceptions import InputValidationError, OutputValidationError, ReferenceError
from hifivar.reference import ReferenceGenome


@dataclass(frozen=True, slots=True)
class TandemRepeatCatalog:
    """One reference-specific TRGT BED catalog."""

    path: Path
    reference_build: str | None = None

    def __post_init__(self) -> None:
        path = Path(self.path).expanduser()
        if path.suffix.lower() != ".bed":
            raise InputValidationError("TRGT catalog must use the .bed suffix.")
        object.__setattr__(self, "path", path)
        if self.reference_build is not None and not self.reference_build.strip():
            raise InputValidationError("TRGT catalog reference_build must be non-empty or None.")

    def validate(self, reference: ReferenceGenome) -> "TandemRepeatCatalog":
        """Stream catalog rows and validate the required TRGT fields."""
        validation.validate_file(self.path, require_nonempty=True)
        reference_names = {contig.name for contig in reference.contigs}
        records = 0
        identifiers: set[str] = set()
        try:
            with self.path.open("r", encoding="utf-8") as handle:
                for line_number, raw_line in enumerate(handle, start=1):
                    line = raw_line.strip()
                    if not line or line.startswith("#"):
                        continue
                    fields = line.split("\t")
                    if len(fields) < 4:
                        raise InputValidationError(
                            f"TRGT catalog '{self.path}' line {line_number} requires at least four tab-separated fields."
                        )
                    contig, start_text, end_text, annotations = fields[:4]
                    if contig not in reference_names:
                        raise ReferenceError(
                            f"REFERENCE_CONTIG_MISMATCH in TRGT catalog '{self.path}' line {line_number}: '{contig}'."
                        )
                    try:
                        start, end = int(start_text), int(end_text)
                    except ValueError as error:
                        raise InputValidationError(
                            f"TRGT catalog '{self.path}' line {line_number} has non-integer coordinates."
                        ) from error
                    if start < 0 or end <= start:
                        raise InputValidationError(
                            f"TRGT catalog '{self.path}' line {line_number} has invalid BED coordinates {start}-{end}."
                        )
                    tags = _parse_catalog_tags(annotations, self.path, line_number)
                    identifier = tags["ID"]
                    if identifier in identifiers:
                        raise InputValidationError(
                            f"TRGT catalog '{self.path}' contains duplicate ID '{identifier}'."
                        )
                    identifiers.add(identifier)
                    records += 1
        except UnicodeDecodeError as error:
            raise InputValidationError(f"TRGT catalog '{self.path}' is not valid UTF-8.") from error
        if records == 0:
            raise InputValidationError(f"TRGT catalog '{self.path}' contains no repeat records.")
        if self.reference_build is not None and self.reference_build != reference.build:
            raise ReferenceError(
                f"TRGT catalog build '{self.reference_build}' conflicts with reference build '{reference.build}'."
            )
        return self

    def to_dict(self) -> dict[str, object]:
        return {"path": str(self.path), "reference_build": self.reference_build}


@dataclass(frozen=True, slots=True)
class TandemRepeatArtifact:
    """Validated, sorted/indexed single-sample TRGT outputs."""

    sample_id: str
    catalog: TandemRepeatCatalog
    reference_build: str | None
    vcf_path: Path
    vcf_index_path: Path
    spanning_bam_path: Path
    spanning_bam_index_path: Path
    trgt_version: str
    commands: tuple[tuple[str, ...], ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "catalog": self.catalog.to_dict(),
            "reference_build": self.reference_build,
            "vcf_path": str(self.vcf_path),
            "vcf_index_path": str(self.vcf_index_path),
            "spanning_bam_path": str(self.spanning_bam_path),
            "spanning_bam_index_path": str(self.spanning_bam_index_path),
            "trgt_version": self.trgt_version,
            "commands": [list(command) for command in self.commands],
        }


def validate_tandem_repeat_outputs(
    *,
    sample_id: str,
    reference: ReferenceGenome,
    catalog: TandemRepeatCatalog,
    vcf_path: Path,
    vcf_index_path: Path,
    spanning_bam_path: Path,
    spanning_bam_index_path: Path,
    trgt_version: str,
    commands: tuple[tuple[str, ...], ...],
) -> TandemRepeatArtifact:
    """Validate lightweight TRGT headers and post-processed output files."""
    for path in (vcf_path, vcf_index_path, spanning_bam_path, spanning_bam_index_path):
        validation.validate_output_file(path)
    _validate_bgzf(vcf_path, "TRGT VCF")
    _validate_tabix_index(vcf_index_path)
    _validate_trgt_vcf_header(vcf_path, sample_id, {item.name for item in reference.contigs})
    return TandemRepeatArtifact(
        sample_id,
        catalog,
        reference.build,
        vcf_path,
        vcf_index_path,
        spanning_bam_path,
        spanning_bam_index_path,
        trgt_version,
        commands,
    )


def _parse_catalog_tags(annotations: str, path: Path, line_number: int) -> dict[str, str]:
    tags: dict[str, str] = {}
    for item in annotations.split(";"):
        key, separator, value = item.partition("=")
        if separator and key and value:
            tags[key] = value
    missing = sorted({"ID", "MOTIFS", "STRUC"}.difference(tags))
    if missing:
        raise InputValidationError(
            f"TRGT catalog '{path}' line {line_number} is missing annotation tags {missing!r}."
        )
    return tags


def _validate_trgt_vcf_header(path: Path, sample_id: str, reference_contigs: set[str]) -> None:
    fileformat = False
    trgt_fields: set[str] = set()
    contigs: set[str] = set()
    columns: str | None = None
    try:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                line = line.rstrip("\r\n")
                if line.startswith("##fileformat=VCFv4"):
                    fileformat = True
                for field in ("TRID", "MOTIFS", "STRUC"):
                    if line.startswith(f"##INFO=<ID={field},"):
                        trgt_fields.add(field)
                if line.startswith("##contig=<"):
                    match = re.search(r"(?:^|,)ID=([^,>]+)", line[10:])
                    if match:
                        contigs.add(match.group(1))
                if line.startswith("#CHROM\t"):
                    columns = line
                    break
    except (OSError, EOFError, UnicodeDecodeError) as error:
        raise OutputValidationError(f"Unable to read TRGT VCF header '{path}': {error}") from error
    if not fileformat or columns is None:
        raise OutputValidationError(f"TRGT VCF header is incomplete: '{path}'.")
    samples = columns.split("\t")[9:]
    if samples != [sample_id]:
        raise OutputValidationError(
            f"TRGT VCF sample mismatch in '{path}': expected '{sample_id}', observed {samples!r}."
        )
    missing = sorted({"TRID", "MOTIFS", "STRUC"}.difference(trgt_fields))
    if missing:
        raise OutputValidationError(f"TRGT VCF is missing INFO declarations {missing!r}: '{path}'.")
    unexpected = sorted(contigs.difference(reference_contigs))
    if unexpected:
        raise ReferenceError(f"REFERENCE_CONTIG_MISMATCH in TRGT VCF '{path}': {unexpected!r}.")


def _validate_bgzf(path: Path, label: str) -> None:
    try:
        with path.open("rb") as handle:
            fixed = handle.read(12)
            if len(fixed) != 12 or fixed[:4] != b"\x1f\x8b\x08\x04":
                raise OutputValidationError(f"{label} is not BGZF-compressed: '{path}'.")
            extra_length = int.from_bytes(fixed[10:12], "little")
            extra = handle.read(extra_length)
    except OSError as error:
        raise OutputValidationError(f"Unable to inspect {label} '{path}': {error}") from error
    offset = 0
    while offset + 4 <= len(extra):
        subfield_length = int.from_bytes(extra[offset + 2 : offset + 4], "little")
        if extra[offset : offset + 2] == b"BC" and subfield_length == 2:
            return
        offset += 4 + subfield_length
    raise OutputValidationError(f"{label} is not BGZF-compressed: '{path}'.")


def _validate_tabix_index(path: Path) -> None:
    _validate_bgzf(path, "tabix index")
    try:
        with gzip.open(path, "rb") as handle:
            magic = handle.read(4)
    except (OSError, EOFError) as error:
        raise OutputValidationError(f"Unable to read tabix index '{path}': {error}") from error
    if magic != b"TBI\x01":
        raise OutputValidationError(f"Invalid tabix index magic in '{path}'.")


__all__ = [
    "TandemRepeatArtifact",
    "TandemRepeatCatalog",
    "validate_tandem_repeat_outputs",
]
