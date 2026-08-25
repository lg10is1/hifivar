"""Tool-neutral read-based phasing contracts and lightweight validation."""

from __future__ import annotations

import gzip
import os
import re
from dataclasses import dataclass
from pathlib import Path

from hifivar import validation
from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import AlignmentArtifact
from hifivar.exceptions import InputValidationError, OutputValidationError, ReferenceError
from hifivar.small import SmallVariantArtifact


@dataclass(frozen=True, slots=True)
class PhasingResources:
    """Scheduler-neutral resources for one HiPhase invocation."""

    threads: int = 16
    memory_mb: int = 32000
    runtime_minutes: int = 1440

    def __post_init__(self) -> None:
        for name in ("threads", "memory_mb", "runtime_minutes"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise InputValidationError(f"Phasing {name} must be a positive integer.")

    def to_dict(self) -> dict[str, int]:
        return {name: getattr(self, name) for name in ("threads", "memory_mb", "runtime_minutes")}


@dataclass(frozen=True, slots=True)
class PhasingRequest:
    """One indexed alignment and small-variant VCF to phase without mutation."""

    alignment: AlignmentArtifact
    small_variants: SmallVariantArtifact
    output_vcf: Path
    resources: PhasingResources = PhasingResources()
    overwrite: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.alignment, AlignmentArtifact):
            raise InputValidationError("Phasing requires an AlignmentArtifact, not raw FASTQ.")
        if self.alignment.output_format is not AlignmentOutputFormat.BAM:
            raise InputValidationError("Phase 6 HiPhase currently requires indexed BAM input.")
        if self.alignment.index_path is None:
            raise InputValidationError("Phase 6 HiPhase requires an indexed BAM.")
        if not isinstance(self.small_variants, SmallVariantArtifact):
            raise InputValidationError("Phasing requires a validated SmallVariantArtifact.")
        if self.alignment.sample_id != self.small_variants.sample_id:
            raise InputValidationError("Phasing BAM and small VCF sample IDs differ.")
        if self.alignment.reference.build != self.small_variants.reference_build:
            raise ReferenceError("Phasing BAM and small VCF reference builds differ.")
        if not isinstance(self.resources, PhasingResources):
            raise InputValidationError("Phasing resources must be PhasingResources.")
        if not isinstance(self.overwrite, bool):
            raise InputValidationError("Phasing overwrite must be boolean.")
        output = Path(self.output_vcf).expanduser()
        if output.name != f"{self.sample_id}.phased.vcf.gz":
            raise InputValidationError("Phased VCF must follow {sample}.phased.vcf.gz.")
        protected = {
            _identity(self.alignment.path),
            _identity(self.alignment.index_path),
            _identity(self.small_variants.vcf_path),
            _identity(self.small_variants.vcf_index_path),
            _identity(self.alignment.reference.fasta),
        }
        if _identity(output) in protected:
            raise InputValidationError("Phased output must not replace a source artifact.")
        for path in (output, Path(f"{output}.tbi")):
            if path.exists() and not self.overwrite:
                raise OutputValidationError(f"Phasing output already exists: '{path}'.")
        object.__setattr__(self, "output_vcf", output)

    @property
    def sample_id(self) -> str:
        return self.alignment.sample_id

    @property
    def output_vcf_index(self) -> Path:
        return Path(f"{self.output_vcf}.tbi")

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "alignment": self.alignment.to_dict(),
            "small_variants": self.small_variants.to_dict(),
            "output_vcf": str(self.output_vcf),
            "output_vcf_index": str(self.output_vcf_index),
            "resources": self.resources.to_dict(),
            "overwrite": self.overwrite,
        }


@dataclass(frozen=True, slots=True)
class PhasedVariantArtifact:
    """Validated phased VCF plus immutable source/provenance references."""

    sample_id: str
    reference_build: str | None
    vcf_path: Path
    vcf_index_path: Path
    source_bam: Path
    source_small_vcf: Path
    hiphase_version: str
    command: tuple[str, ...]
    execution_backend: str = "native"

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "reference_build": self.reference_build,
            "vcf_path": str(self.vcf_path),
            "vcf_index_path": str(self.vcf_index_path),
            "source_bam": str(self.source_bam),
            "source_small_vcf": str(self.source_small_vcf),
            "hiphase_version": self.hiphase_version,
            "command": list(self.command),
            "execution_backend": self.execution_backend,
        }


def validate_phased_variant_output(
    request: PhasingRequest, *, hiphase_version: str, command: tuple[str, ...]
) -> PhasedVariantArtifact:
    """Stream only the phased VCF header and validate its tabix index."""
    validation.validate_output_file(request.output_vcf)
    validation.validate_output_file(request.output_vcf_index)
    _validate_bgzf(request.output_vcf, "phased VCF")
    _validate_bgzf(request.output_vcf_index, "tabix index")
    try:
        with gzip.open(request.output_vcf_index, "rb") as handle:
            if handle.read(4) != b"TBI\x01":
                raise OutputValidationError("Phased VCF tabix index has invalid magic.")
    except (OSError, EOFError) as error:
        raise OutputValidationError(f"Unable to read phased VCF index: {error}") from error
    fileformat = False
    phase_set = False
    contigs: set[str] = set()
    columns: list[str] | None = None
    try:
        with gzip.open(request.output_vcf, "rt", encoding="utf-8") as handle:
            for number, line in enumerate(handle, start=1):
                if number > 100000:
                    raise OutputValidationError("Phased VCF header exceeds 100000 lines.")
                fileformat |= line.startswith("##fileformat=VCFv4")
                phase_set |= line.startswith("##FORMAT=<ID=PS,")
                if line.startswith("##contig=<"):
                    match = re.search(r"(?:^|,)ID=([^,>]+)", line[10:])
                    if match:
                        contigs.add(match.group(1))
                if line.startswith("#CHROM\t"):
                    columns = line.rstrip("\r\n").split("\t")
                    break
    except (OSError, EOFError, UnicodeError) as error:
        raise OutputValidationError(f"Unable to read phased VCF header: {error}") from error
    if not fileformat or not phase_set or columns is None:
        raise OutputValidationError("Phased VCF header is incomplete or lacks FORMAT/PS.")
    if columns[9:] != [request.sample_id]:
        raise OutputValidationError(f"Phased VCF sample mismatch: observed {columns[9:]!r}.")
    reference_contigs = {item.name for item in request.alignment.reference.contigs}
    unexpected = sorted(contigs.difference(reference_contigs))
    if unexpected:
        raise ReferenceError(f"REFERENCE_CONTIG_MISMATCH in phased VCF: {unexpected!r}.")
    return PhasedVariantArtifact(
        request.sample_id,
        request.alignment.reference.build,
        request.output_vcf,
        request.output_vcf_index,
        request.alignment.path,
        request.small_variants.vcf_path,
        hiphase_version,
        command,
    )


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
        length = int.from_bytes(extra[offset + 2 : offset + 4], "little")
        if extra[offset : offset + 2] == b"BC" and length == 2:
            return
        offset += 4 + length
    raise OutputValidationError(f"{label} is not BGZF-compressed: '{path}'.")


def _identity(path: Path) -> str:
    return os.path.normcase(os.path.normpath(str(Path(path).absolute())))


__all__ = [
    "PhasedVariantArtifact",
    "PhasingRequest",
    "PhasingResources",
    "validate_phased_variant_output",
]
