"""Tool-neutral single-sample small-variant calling models."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import gzip
import os
from pathlib import Path
import re

from hifivar.alignment_postprocess import AlignmentArtifact
from hifivar.exceptions import InputValidationError, OutputValidationError


class SmallVariantModelType(str, Enum):
    """DeepVariant model types supported by the current HiFi-only phase."""

    PACBIO = "PACBIO"


class SmallVariantResultStatus(str, Enum):
    """Execution states for a single-sample small-variant call."""

    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class SmallVariantResources:
    """Scheduler-neutral resources for one DeepVariant invocation."""

    threads: int = 8
    memory_mb: int = 32_000
    runtime_minutes: int = 1_440

    def __post_init__(self) -> None:
        for name in ("threads", "memory_mb", "runtime_minutes"):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
                raise InputValidationError(
                    f"Small-variant {name} must be a positive integer."
                )

    def to_dict(self) -> dict[str, int]:
        return {
            "threads": self.threads,
            "memory_mb": self.memory_mb,
            "runtime_minutes": self.runtime_minutes,
        }


@dataclass(frozen=True, slots=True)
class DeepVariantRequest:
    """Immutable request for one aligned sample and separate VCF/gVCF outputs."""

    artifact: AlignmentArtifact
    output_vcf: Path
    output_gvcf: Path
    resources: SmallVariantResources = SmallVariantResources()
    model_type: SmallVariantModelType = SmallVariantModelType.PACBIO
    overwrite: bool = False
    intermediate_directory: Path | None = None
    logging_directory: Path | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.artifact, AlignmentArtifact):
            raise InputValidationError(
                "DeepVariant request artifact must be AlignmentArtifact."
            )
        if not isinstance(self.resources, SmallVariantResources):
            raise InputValidationError(
                "DeepVariant request resources must be SmallVariantResources."
            )
        if self.model_type is not SmallVariantModelType.PACBIO:
            raise InputValidationError(
                "HiFiVar Phase 3 requires the DeepVariant PACBIO model."
            )
        if not isinstance(self.overwrite, bool):
            raise InputValidationError("DeepVariant overwrite must be boolean.")

        output_vcf = _coerce_path(self.output_vcf, "small VCF")
        output_gvcf = _coerce_path(self.output_gvcf, "gVCF")
        _require_suffix(output_vcf, ".small.vcf.gz", "small VCF")
        _require_suffix(output_gvcf, ".g.vcf.gz", "gVCF")
        if _path_identity(output_vcf) == _path_identity(output_gvcf):
            raise InputValidationError("DeepVariant VCF and gVCF must be distinct.")
        protected = {
            _path_identity(self.artifact.path),
            _path_identity(self.artifact.reference.fasta),
        }
        if self.artifact.index_path is not None:
            protected.add(_path_identity(self.artifact.index_path))
        if _path_identity(output_vcf) in protected or _path_identity(output_gvcf) in protected:
            raise InputValidationError(
                "DeepVariant outputs must not replace an input alignment, index, or reference."
            )
        for output in (output_vcf, output_gvcf):
            if output.exists() and output.is_dir():
                raise OutputValidationError(
                    f"DeepVariant output path is a directory: '{output}'."
                )
            if output.exists() and not self.overwrite:
                raise OutputValidationError(
                    f"DeepVariant output already exists: '{output}'."
                )
        object.__setattr__(self, "output_vcf", output_vcf)
        object.__setattr__(self, "output_gvcf", output_gvcf)

        intermediate = self.intermediate_directory or (
            output_vcf.parent / f"{self.sample_id}.deepvariant-intermediate"
        )
        logging = self.logging_directory or (
            output_vcf.parent / "logs" / self.sample_id
        )
        object.__setattr__(
            self,
            "intermediate_directory",
            _coerce_path(intermediate, "DeepVariant intermediate directory"),
        )
        object.__setattr__(
            self,
            "logging_directory",
            _coerce_path(logging, "DeepVariant logging directory"),
        )

    @property
    def sample_id(self) -> str:
        return self.artifact.sample_id

    @property
    def alignment_path(self) -> Path:
        return self.artifact.path

    @property
    def reference_fasta(self) -> Path:
        return self.artifact.reference.fasta

    @property
    def output_vcf_index(self) -> Path:
        return Path(f"{self.output_vcf}.tbi")

    @property
    def output_gvcf_index(self) -> Path:
        return Path(f"{self.output_gvcf}.tbi")

    @classmethod
    def create(
        cls,
        artifact: AlignmentArtifact,
        output_directory: str | Path,
        *,
        resources: SmallVariantResources | None = None,
        model_type: SmallVariantModelType = SmallVariantModelType.PACBIO,
        overwrite: bool = False,
    ) -> DeepVariantRequest:
        output_root = _coerce_path(output_directory, "small-variant output directory")
        sample_id = artifact.sample_id if isinstance(artifact, AlignmentArtifact) else "sample"
        return cls(
            artifact=artifact,
            output_vcf=output_root / f"{sample_id}.small.vcf.gz",
            output_gvcf=output_root / f"{sample_id}.g.vcf.gz",
            resources=resources or SmallVariantResources(),
            model_type=model_type,
            overwrite=overwrite,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "alignment": self.artifact.to_dict(),
            "reference_fasta": str(self.reference_fasta),
            "model_type": self.model_type.value,
            "output_vcf": str(self.output_vcf),
            "output_gvcf": str(self.output_gvcf),
            "output_vcf_index": str(self.output_vcf_index),
            "output_gvcf_index": str(self.output_gvcf_index),
            "intermediate_directory": str(self.intermediate_directory),
            "logging_directory": str(self.logging_directory),
            "resources": self.resources.to_dict(),
            "overwrite": self.overwrite,
        }


@dataclass(frozen=True, slots=True)
class SmallVariantCommandPlan:
    """Serializable shell-free representation of a DeepVariant command."""

    args: tuple[str, ...]
    display_command: str
    shell: bool = False

    def to_dict(self) -> dict[str, object]:
        return {
            "args": list(self.args),
            "display_command": self.display_command,
            "shell": self.shell,
        }


@dataclass(frozen=True, slots=True)
class SmallVariantResult:
    """One planned or completed DeepVariant result with provenance."""

    request: DeepVariantRequest
    status: SmallVariantResultStatus
    command: SmallVariantCommandPlan
    tool_version: str | None = None
    duration_seconds: float = 0.0
    artifact: SmallVariantArtifact | None = None

    @property
    def executed(self) -> bool:
        return self.status is SmallVariantResultStatus.COMPLETED

    def to_dict(self) -> dict[str, object]:
        return {
            "request": self.request.to_dict(),
            "status": self.status.value,
            "executed": self.executed,
            "command": self.command.to_dict(),
            "tool_version": self.tool_version,
            "duration_seconds": self.duration_seconds,
            "artifact": self.artifact.to_dict() if self.artifact else None,
        }


@dataclass(frozen=True, slots=True)
class SmallVariantArtifact:
    """Validated DeepVariant VCF/gVCF pair and their tabix indexes."""

    sample_id: str
    reference_build: str | None
    vcf_path: Path
    gvcf_path: Path
    vcf_index_path: Path
    gvcf_index_path: Path
    tool: str = "deepvariant"
    tool_version: str | None = None
    reference_compatibility: str = "declared_not_header_verified"

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "reference_build": self.reference_build,
            "vcf_path": str(self.vcf_path),
            "gvcf_path": str(self.gvcf_path),
            "vcf_index_path": str(self.vcf_index_path),
            "gvcf_index_path": str(self.gvcf_index_path),
            "tool": self.tool,
            "tool_version": self.tool_version,
            "reference_compatibility": self.reference_compatibility,
        }


def validate_small_variant_outputs(
    request: DeepVariantRequest,
    *,
    tool_version: str | None = None,
) -> SmallVariantArtifact:
    """Validate BGZF VCF/gVCF headers, sample, contigs, and tabix indexes.

    This streams only VCF headers. It does not load variant records or claim
    that the input BAM/CRAM header matches the declared reference.
    """
    if not isinstance(request, DeepVariantRequest):
        raise InputValidationError(
            "Small-variant output validation requires DeepVariantRequest."
        )
    reference_contigs = {contig.name for contig in request.artifact.reference.contigs}
    _validate_vcf_header(
        request.output_vcf,
        sample_id=request.sample_id,
        reference_contigs=reference_contigs,
        require_deepvariant_gvcf=False,
    )
    _validate_vcf_header(
        request.output_gvcf,
        sample_id=request.sample_id,
        reference_contigs=reference_contigs,
        require_deepvariant_gvcf=True,
    )
    _validate_tabix_index(request.output_vcf_index)
    _validate_tabix_index(request.output_gvcf_index)
    return SmallVariantArtifact(
        sample_id=request.sample_id,
        reference_build=request.artifact.reference.build,
        vcf_path=request.output_vcf,
        gvcf_path=request.output_gvcf,
        vcf_index_path=request.output_vcf_index,
        gvcf_index_path=request.output_gvcf_index,
        tool_version=tool_version,
    )


def _validate_vcf_header(
    path: Path,
    *,
    sample_id: str,
    reference_contigs: set[str],
    require_deepvariant_gvcf: bool,
) -> None:
    _validate_bgzf(path, "VCF")
    fileformat = False
    deepvariant_version = False
    ref_call_filter = False
    reference_block_depth = False
    contigs: set[str] = set()
    column_header: str | None = None
    try:
        with gzip.open(path, mode="rt", encoding="utf-8", errors="strict") as handle:
            for line_number, line in enumerate(handle, start=1):
                if line_number > 100_000:
                    raise OutputValidationError(
                        f"VCF header exceeds 100000 lines: '{path}'."
                    )
                if line.startswith("##fileformat=VCFv4"):
                    fileformat = True
                if line.startswith("##DeepVariant_version="):
                    deepvariant_version = True
                if line.startswith("##FILTER=<ID=RefCall"):
                    ref_call_filter = True
                if line.startswith((
                    "##FORMAT=<ID=MIN_DP",
                    "##FORMAT=<ID=MED_DP",
                )):
                    reference_block_depth = True
                if line.startswith("##contig=<"):
                    match = re.search(r"(?:^|,)ID=([^,>]+)", line[10:])
                    if match is not None:
                        contigs.add(match.group(1))
                if line.startswith("#CHROM\t"):
                    column_header = line.rstrip("\r\n")
                    break
    except (OSError, UnicodeError, EOFError) as error:
        raise OutputValidationError(f"Unable to read VCF header '{path}': {error}") from error
    if not fileformat or column_header is None:
        raise OutputValidationError(
            f"VCF header is incomplete or invalid: '{path}'."
        )
    columns = column_header.split("\t")
    if len(columns) != 10 or columns[9] != sample_id:
        observed = columns[9:] if len(columns) >= 10 else []
        raise OutputValidationError(
            f"VCF sample mismatch for '{path}': expected '{sample_id}', "
            f"observed {observed!r}."
        )
    if not contigs:
        raise OutputValidationError(f"VCF has no contig header records: '{path}'.")
    unexpected = sorted(contigs.difference(reference_contigs))
    if unexpected:
        raise OutputValidationError(
            f"REFERENCE_CONTIG_MISMATCH in '{path}': unexpected contigs "
            f"{unexpected!r}."
        )
    if require_deepvariant_gvcf and not (
        deepvariant_version and ref_call_filter and reference_block_depth
    ):
        missing: list[str] = []
        if not deepvariant_version:
            missing.append("##DeepVariant_version")
        if not ref_call_filter:
            missing.append("##FILTER=<ID=RefCall")
        if not reference_block_depth:
            missing.append("##FORMAT=<ID=MIN_DP or MED_DP")
        raise OutputValidationError(
            f"gVCF is missing DeepVariant reference-block header markers "
            f"{missing!r}: '{path}'."
        )


def _validate_bgzf(path: Path, label: str) -> None:
    try:
        with path.open("rb") as handle:
            fixed = handle.read(12)
            if len(fixed) != 12 or fixed[:4] != b"\x1f\x8b\x08\x04":
                raise OutputValidationError(f"{label} is not BGZF-compressed: '{path}'.")
            extra_length = int.from_bytes(fixed[10:12], "little")
            extra = handle.read(extra_length)
    except OutputValidationError:
        raise
    except OSError as error:
        raise OutputValidationError(f"Unable to read {label} '{path}': {error}") from error
    offset = 0
    while offset + 4 <= len(extra):
        subfield_id = extra[offset : offset + 2]
        subfield_length = int.from_bytes(extra[offset + 2 : offset + 4], "little")
        offset += 4
        if subfield_id == b"BC" and subfield_length == 2:
            return
        offset += subfield_length
    raise OutputValidationError(f"{label} is not BGZF-compressed: '{path}'.")


def _validate_tabix_index(path: Path) -> None:
    _validate_bgzf(path, "tabix index")
    try:
        with gzip.open(path, "rb") as handle:
            magic = handle.read(4)
    except (OSError, EOFError) as error:
        raise OutputValidationError(f"Unable to read tabix index '{path}': {error}") from error
    if magic != b"TBI\x01":
        raise OutputValidationError(f"Invalid tabix index header: '{path}'.")


def _coerce_path(value: str | Path, label: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise InputValidationError(f"{label} path must be string or Path.")
    if isinstance(value, str) and not value.strip():
        raise InputValidationError(f"{label} path must not be empty.")
    return Path(value).expanduser()


def _require_suffix(path: Path, suffix: str, label: str) -> None:
    if not str(path).lower().endswith(suffix):
        raise InputValidationError(f"DeepVariant {label} must end with '{suffix}'.")


def _path_identity(path: Path) -> str:
    return os.path.normcase(os.path.normpath(str(path.absolute())))


__all__ = [
    "DeepVariantRequest",
    "SmallVariantCommandPlan",
    "SmallVariantArtifact",
    "SmallVariantModelType",
    "SmallVariantResources",
    "SmallVariantResult",
    "SmallVariantResultStatus",
    "validate_small_variant_outputs",
]
