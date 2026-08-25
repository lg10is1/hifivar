"""Assembly-derived structural-variant contracts and validation."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hifivar.assembly import AssemblyRole, HaplotypeAssemblyArtifact
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.reference import ReferenceGenome
from hifivar.sv import validate_sv_vcf
from hifivar.validation import validate_output_file


class SVEvidenceSource(str, Enum):
    READ = "read"
    ASSEMBLY = "assembly"


class AssemblySvCaller(str, Enum):
    PAV = "pav"
    SVIM_ASM = "svim_asm"


@dataclass(frozen=True, slots=True)
class AssemblySvResources:
    threads: int = 16
    memory_mb: int = 64000
    runtime_minutes: int = 2880

    def __post_init__(self) -> None:
        for name in ("threads", "memory_mb", "runtime_minutes"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise InputValidationError(f"Assembly-SV {name} must be a positive integer.")

    def to_dict(self) -> dict[str, int]:
        return {name: getattr(self, name) for name in ("threads", "memory_mb", "runtime_minutes")}


@dataclass(frozen=True, slots=True)
class AssemblySvRequest:
    """One caller-specific assembly/reference comparison without implicit merging."""

    sample_id: str
    caller: AssemblySvCaller
    reference: ReferenceGenome
    assemblies: tuple[HaplotypeAssemblyArtifact, ...]
    work_directory: Path
    output_vcf: Path
    resources: AssemblySvResources = AssemblySvResources()
    overwrite: bool = False

    def __post_init__(self) -> None:
        if not self.sample_id:
            raise InputValidationError("Assembly-SV sample_id must be non-empty.")
        if not isinstance(self.caller, AssemblySvCaller):
            raise InputValidationError("Assembly-SV caller must be AssemblySvCaller.")
        if not isinstance(self.reference, ReferenceGenome):
            raise InputValidationError("Assembly-SV request requires ReferenceGenome.")
        if not self.assemblies or len(self.assemblies) > 2:
            raise InputValidationError("Assembly-SV requires one or two haplotype assemblies.")
        roles = tuple(item.role for item in self.assemblies)
        if len(set(roles)) != len(roles):
            raise InputValidationError("Assembly-SV haplotype roles must be unique.")
        if any(item.sample_id != self.sample_id for item in self.assemblies):
            raise InputValidationError("Assembly-SV sample and assembly sample IDs differ.")
        allowed = {AssemblyRole.HAPLOTYPE1, AssemblyRole.HAPLOTYPE2}
        if any(role not in allowed for role in roles):
            raise InputValidationError("Assembly-SV accepts explicit haplotype1/haplotype2 FASTAs only.")
        if len(roles) == 2 and set(roles) != allowed:
            raise InputValidationError("Diploid assembly-SV requires haplotype1 and haplotype2.")
        expected = f"{self.sample_id}.{self.caller.value}.assembly.sv.vcf.gz"
        output = Path(self.output_vcf)
        if output.name != expected:
            raise InputValidationError(f"Assembly-SV output must be '{expected}'.")
        if not isinstance(self.resources, AssemblySvResources):
            raise InputValidationError("Assembly-SV resources must be AssemblySvResources.")
        if not isinstance(self.overwrite, bool):
            raise InputValidationError("Assembly-SV overwrite must be boolean.")
        object.__setattr__(self, "work_directory", Path(self.work_directory))
        object.__setattr__(self, "output_vcf", output)
        for path in (output, Path(f"{output}.tbi")):
            if path.exists() and not self.overwrite:
                raise OutputValidationError(f"Assembly-SV output already exists: '{path}'.")

    @property
    def output_index(self) -> Path:
        return Path(f"{self.output_vcf}.tbi")

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "caller": self.caller.value,
            "evidence_source": SVEvidenceSource.ASSEMBLY.value,
            "reference": self.reference.to_dict(),
            "assemblies": [item.to_dict() for item in self.assemblies],
            "work_directory": str(self.work_directory),
            "output_vcf": str(self.output_vcf),
            "output_index": str(self.output_index),
            "resources": self.resources.to_dict(),
            "overwrite": self.overwrite,
        }

@dataclass(frozen=True, slots=True)
class AssemblySvArtifact:
    caller: AssemblySvCaller
    sample_id: str
    reference_fasta: Path
    reference_build: str | None
    assemblies: tuple[HaplotypeAssemblyArtifact, ...]
    raw_vcf: Path
    vcf_path: Path
    index_path: Path
    intermediate_files: tuple[Path, ...]
    caller_version: str
    backend: str
    commands: tuple[tuple[str, ...], ...]
    evidence_source: SVEvidenceSource = SVEvidenceSource.ASSEMBLY
    harmonized: bool = False

    def __post_init__(self) -> None:
        if self.evidence_source is not SVEvidenceSource.ASSEMBLY or self.harmonized:
            raise InputValidationError("Raw assembly caller artifacts must remain assembly evidence.")
        if self.index_path != Path(f"{self.vcf_path}.tbi"):
            raise InputValidationError("Assembly-SV index must be '<vcf>.tbi'.")

    def to_dict(self) -> dict[str, object]:
        return {
            "caller": self.caller.value,
            "sample_id": self.sample_id,
            "evidence_source": self.evidence_source.value,
            "reference_fasta": str(self.reference_fasta),
            "reference_build": self.reference_build,
            "assemblies": [item.to_dict() for item in self.assemblies],
            "raw_vcf": str(self.raw_vcf),
            "vcf_path": str(self.vcf_path),
            "index_path": str(self.index_path),
            "intermediate_files": [str(path) for path in self.intermediate_files],
            "caller_version": self.caller_version,
            "backend": self.backend,
            "commands": [list(item) for item in self.commands],
            "harmonized": self.harmonized,
        }


@dataclass(frozen=True, slots=True)
class AssemblySvCollection:
    sample_id: str
    artifacts: tuple[AssemblySvArtifact, ...]

    def __post_init__(self) -> None:
        if any(item.sample_id != self.sample_id for item in self.artifacts):
            raise InputValidationError("Assembly-SV collection contains another sample.")
        callers = tuple(item.caller for item in self.artifacts)
        if len(callers) != len(set(callers)):
            raise InputValidationError("Assembly-SV collection contains duplicate callers.")

    def to_dict(self) -> dict[str, object]:
        return {"sample_id": self.sample_id, "artifacts": [item.to_dict() for item in self.artifacts]}


def create_assembly_sv_artifact(
    request: AssemblySvRequest,
    *,
    raw_vcf: Path,
    intermediate_files: tuple[Path, ...],
    caller_version: str,
    backend: str,
    commands: tuple[tuple[str, ...], ...],
) -> AssemblySvArtifact:
    """Validate and retain one caller's raw and finalized evidence."""
    validate_output_file(raw_vcf)
    validate_sv_vcf(
        request.output_vcf,
        request.output_index,
        sample_id=request.sample_id,
        reference=request.reference,
    )
    return AssemblySvArtifact(
        request.caller,
        request.sample_id,
        request.reference.fasta,
        request.reference.build,
        request.assemblies,
        Path(raw_vcf),
        request.output_vcf,
        request.output_index,
        tuple(Path(path) for path in intermediate_files),
        caller_version,
        backend,
        commands,
    )


__all__ = [
    "AssemblySvArtifact",
    "AssemblySvCaller",
    "AssemblySvCollection",
    "AssemblySvRequest",
    "AssemblySvResources",
    "SVEvidenceSource",
    "create_assembly_sv_artifact",
]
