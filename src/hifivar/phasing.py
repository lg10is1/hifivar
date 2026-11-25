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
from hifivar.joint_phasing import PhasedTrackArtifact, PhasingTrackInput
from hifivar.reference import ReferenceGenome
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
    """One indexed alignment and separate supported VCF tracks, without source mutation."""

    alignment: AlignmentArtifact
    small_variants: SmallVariantArtifact
    output_vcf: Path
    resources: PhasingResources = PhasingResources()
    overwrite: bool = False
    additional_inputs: tuple[PhasingTrackInput, ...] = ()
    global_realignment: str = "auto"

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
        if self.global_realignment not in ("auto", "enabled", "disabled"):
            raise InputValidationError("Phasing global_realignment must be auto, enabled, or disabled.")
        if not isinstance(self.additional_inputs, tuple) or any(
            not isinstance(item, PhasingTrackInput) for item in self.additional_inputs
        ):
            raise InputValidationError("Phasing additional_inputs must be a tuple of PhasingTrackInput.")
        categories = [item.category for item in self.additional_inputs]
        if len(categories) != len(set(categories)):
            raise InputValidationError("Joint phasing permits one explicitly selected source per category.")
        for item in self.additional_inputs:
            if item.sample_id != self.sample_id:
                raise InputValidationError("Joint phasing VCF sample IDs differ.")
            if item.reference_build != self.alignment.reference.build:
                raise ReferenceError("Joint phasing source reference builds differ.")
        object.__setattr__(self, "additional_inputs", tuple(sorted(
            self.additional_inputs, key=lambda item: {"sv": 0, "tr": 1}[item.category],
        )))
        output = Path(self.output_vcf).expanduser()
        if output.name != f"{self.sample_id}.phased.vcf.gz":
            raise InputValidationError("Phased VCF must follow {sample}.phased.vcf.gz.")
        object.__setattr__(self, "output_vcf", output)
        self.validate_output_safety()

    @property
    def sample_id(self) -> str:
        return self.alignment.sample_id

    @property
    def output_vcf_index(self) -> Path:
        return Path(f"{self.output_vcf}.tbi")

    @property
    def use_global_realignment(self) -> bool:
        return self.global_realignment == "enabled" or (
            self.global_realignment == "auto" and bool(self.additional_inputs)
        )

    @property
    def vcf_pairs(self) -> tuple[tuple[str, Path, Path, Path, Path, str], ...]:
        """Canonical small, SV, TR order for matching input/output CLI flags."""
        pairs = [("small", self.small_variants.vcf_path, self.small_variants.vcf_index_path,
                  self.output_vcf, self.output_vcf_index, "deepvariant")]
        for item in self.additional_inputs:
            output = self.output_vcf.parent / f"{self.sample_id}.{item.category}.phased.vcf.gz"
            pairs.append((item.category, item.vcf_path, item.index_path, output,
                          Path(f"{output}.tbi"), item.source_tool))
        return tuple(pairs)

    @property
    def protected_input_paths(self) -> tuple[Path, ...]:
        """Explicit sources/indexes/reference used for write-path guards."""
        protected = [self.alignment.path, self.alignment.index_path,
                     self.alignment.reference.fasta, self.alignment.reference.fai,
                     self.small_variants.gvcf_path, self.small_variants.gvcf_index_path]
        protected.extend(p for _, source, index, _, _, _ in self.vcf_pairs for p in (source, index))
        return tuple(Path(p) for p in protected if p is not None)

    def validate_output_safety(self) -> None:
        """Recheck existing/aliased outputs before execution, including every index."""
        protected = self.protected_input_paths
        source_ids = {_identity(p) for p in protected}
        outputs = [p for _, _, _, output, index, _ in self.vcf_pairs for p in (output, index)]
        if len({_identity(p) for p in outputs}) != len(outputs):
            raise InputValidationError("Joint phased output paths must be distinct.")
        sources = [source for _, source, _, _, _, _ in self.vcf_pairs]
        if len({_identity(p) for p in sources}) != len(sources):
            raise InputValidationError("Joint phasing source VCFs must be distinct.")
        if any(a.exists() and b.exists() and a.samefile(b)
               for i, a in enumerate(sources) for b in sources[i + 1:]):
            raise InputValidationError("Joint phasing source VCFs must not be hardlink aliases.")
        for path in outputs:
            if _identity(path) in source_ids or (
                path.exists() and any(p.exists() and path.samefile(p) for p in protected)
            ):
                raise InputValidationError("Phased output must not replace a source artifact.")
            if (path.exists() or path.is_symlink()) and not self.overwrite:
                raise OutputValidationError(f"Phasing output already exists: '{path}'.")

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "alignment": self.alignment.to_dict(),
            "small_variants": self.small_variants.to_dict(),
            "output_vcf": str(self.output_vcf),
            "output_vcf_index": str(self.output_vcf_index),
            "resources": self.resources.to_dict(),
            "overwrite": self.overwrite,
            "additional_inputs": [item.to_dict() for item in self.additional_inputs],
            "global_realignment": self.global_realignment,
            "effective_global_realignment": self.use_global_realignment,
            "phased_outputs": {category: str(out) for category, _, _, out, _, _ in self.vcf_pairs},
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
    additional_artifacts: tuple[PhasedTrackArtifact, ...] = ()
    records_identity_and_genotype_checked: int | None = None

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
            "additional_artifacts": [item.to_dict() for item in self.additional_artifacts],
            "records_identity_and_genotype_checked": self.records_identity_and_genotype_checked,
        }


def validate_phased_variant_output(
    request: PhasingRequest, *, hiphase_version: str, command: tuple[str, ...]
) -> PhasedVariantArtifact:
    """Validate all separate outputs; joint mode additionally checks ordered records."""
    from hifivar.joint_phasing import validate_record_preservation
    from hifivar.annotation import audit_annotation_vcf_identity

    _validate_phased_path(request.output_vcf, request.output_vcf_index,
                         request.sample_id, request.alignment.reference)
    if request.additional_inputs:
        audit_annotation_vcf_identity(request.output_vcf, request.sample_id, request.alignment.reference)
    small_count = (validate_record_preservation(request.small_variants.vcf_path, request.output_vcf)
                   if request.additional_inputs else None)
    additional = []
    for category, source, _, output, index, tool in request.vcf_pairs[1:]:
        _validate_phased_path(output, index, request.sample_id, request.alignment.reference)
        audit_annotation_vcf_identity(output, request.sample_id, request.alignment.reference)
        count = validate_record_preservation(source, output)
        additional.append(PhasedTrackArtifact(category, tool, source, output, index, count))
    return PhasedVariantArtifact(
        request.sample_id, request.alignment.reference.build,
        request.output_vcf, request.output_vcf_index, request.alignment.path,
        request.small_variants.vcf_path, hiphase_version, command,
        additional_artifacts=tuple(additional), records_identity_and_genotype_checked=small_count,
    )


def _validate_phased_path(vcf_path: Path, index_path: Path, sample_id: str, reference: ReferenceGenome) -> None:
    """Keep the legacy BGZF/TBI/PS header contract for each output category."""
    validation.validate_output_file(vcf_path)
    validation.validate_output_file(index_path)
    _validate_bgzf(vcf_path, "phased VCF")
    _validate_bgzf(index_path, "tabix index")
    try:
        with gzip.open(index_path, "rb") as handle:
            if handle.read(4) != b"TBI\x01":
                raise OutputValidationError("Phased VCF tabix index has invalid magic.")
    except (OSError, EOFError) as error:
        raise OutputValidationError(f"Unable to read phased VCF index: {error}") from error
    fileformat = False
    phase_set = False
    contigs: set[str] = set()
    columns: list[str] | None = None
    try:
        with gzip.open(vcf_path, "rt", encoding="utf-8") as handle:
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
    if columns[9:] != [sample_id]:
        raise OutputValidationError(f"Phased VCF sample mismatch: observed {columns[9:]!r}.")
    reference_contigs = {item.name for item in reference.contigs}
    unexpected = sorted(contigs.difference(reference_contigs))
    if unexpected:
        raise ReferenceError(f"REFERENCE_CONTIG_MISMATCH in phased VCF: {unexpected!r}.")


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
    return os.path.normcase(str(Path(path).expanduser().resolve()))


__all__ = [
    "PhasedVariantArtifact",
    "PhasingRequest",
    "PhasingResources",
    "validate_phased_variant_output",
]
