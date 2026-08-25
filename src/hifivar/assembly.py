"""Reference-independent assembly contracts and explicit GFA conversion."""

from __future__ import annotations

import os
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hifivar import validation
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.sample import InputType, Sample


class AssemblyRole(str, Enum):
    PRIMARY = "primary"
    HAPLOTYPE1 = "haplotype1"
    HAPLOTYPE2 = "haplotype2"


@dataclass(frozen=True, slots=True)
class AssemblyResources:
    threads: int = 32
    memory_mb: int = 128000
    runtime_minutes: int = 4320

    def __post_init__(self) -> None:
        for name in ("threads", "memory_mb", "runtime_minutes"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise InputValidationError(f"Assembly {name} must be a positive integer.")

    def to_dict(self) -> dict[str, int]:
        return {
            name: getattr(self, name)
            for name in ("threads", "memory_mb", "runtime_minutes")
        }


@dataclass(frozen=True, slots=True)
class AssemblyRequest:
    """One ordered HiFi FASTQ dataset and deterministic hifiasm paths."""

    sample: Sample
    output_prefix: Path
    assembly_directory: Path
    resources: AssemblyResources = AssemblyResources()
    overwrite: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.sample, Sample):
            raise InputValidationError("Assembly request requires Sample.")
        if self.sample.input.input_type is not InputType.FASTQ:
            raise InputValidationError(
                "Phase 7 assembly requires primary HiFi FASTQ; BAM/CRAM are not converted."
            )
        if not isinstance(self.resources, AssemblyResources):
            raise InputValidationError("Assembly resources must be AssemblyResources.")
        if not isinstance(self.overwrite, bool):
            raise InputValidationError("Assembly overwrite must be boolean.")
        prefix = Path(self.output_prefix).expanduser()
        directory = Path(self.assembly_directory).expanduser()
        if prefix.name != f"{self.sample_id}.asm":
            raise InputValidationError("hifiasm output prefix must follow {sample}.asm.")
        if prefix.parent == directory:
            raise InputValidationError(
                "Raw hifiasm work files and final FASTA directory must be distinct."
            )
        object.__setattr__(self, "output_prefix", prefix)
        object.__setattr__(self, "assembly_directory", directory)
        for path in (*self.raw_gfa_paths.values(), *self.fasta_paths.values()):
            if path.exists() and not self.overwrite:
                raise OutputValidationError(f"Assembly output already exists: '{path}'.")

    @property
    def sample_id(self) -> str:
        return self.sample.sample_id

    @property
    def fastq_files(self) -> tuple[Path, ...]:
        return self.sample.input.files

    @property
    def raw_gfa_paths(self) -> dict[AssemblyRole, Path]:
        return {
            AssemblyRole.PRIMARY: Path(f"{self.output_prefix}.bp.p_ctg.gfa"),
            AssemblyRole.HAPLOTYPE1: Path(f"{self.output_prefix}.bp.hap1.p_ctg.gfa"),
            AssemblyRole.HAPLOTYPE2: Path(f"{self.output_prefix}.bp.hap2.p_ctg.gfa"),
        }

    @property
    def fasta_paths(self) -> dict[AssemblyRole, Path]:
        return {
            AssemblyRole.PRIMARY: self.assembly_directory
            / f"{self.sample_id}.primary.fa",
            AssemblyRole.HAPLOTYPE1: self.assembly_directory
            / f"{self.sample_id}.hap1.fa",
            AssemblyRole.HAPLOTYPE2: self.assembly_directory
            / f"{self.sample_id}.hap2.fa",
        }

    def to_dict(self) -> dict[str, object]:
        return {
            "sample": self.sample.to_dict(),
            "fastq_files": [str(path) for path in self.fastq_files],
            "output_prefix": str(self.output_prefix),
            "assembly_directory": str(self.assembly_directory),
            "raw_gfa_paths": {
                role.value: str(path) for role, path in self.raw_gfa_paths.items()
            },
            "fasta_paths": {
                role.value: str(path) for role, path in self.fasta_paths.items()
            },
            "resources": self.resources.to_dict(),
            "overwrite": self.overwrite,
        }


@dataclass(frozen=True, slots=True)
class AssemblyConversion:
    source_gfa: Path
    output_fasta: Path
    method: str = "hifivar_gfa_segment_stream_v1"

    def to_dict(self) -> dict[str, str]:
        return {
            "source_gfa": str(self.source_gfa),
            "output_fasta": str(self.output_fasta),
            "method": self.method,
        }


@dataclass(frozen=True, slots=True)
class HaplotypeAssemblyArtifact:
    sample_id: str
    role: AssemblyRole
    path: Path
    source_gfa: Path
    hifiasm_version: str
    file_size: int
    reference_independent: bool = True

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "role": self.role.value,
            "path": str(self.path),
            "source_gfa": str(self.source_gfa),
            "hifiasm_version": self.hifiasm_version,
            "file_size": self.file_size,
            "reference_independent": self.reference_independent,
        }


@dataclass(frozen=True, slots=True)
class AssemblyArtifact:
    """Raw hifiasm graph outputs and explicit derived FASTA artifacts."""

    sample_id: str
    raw_gfas: tuple[Path, ...]
    assemblies: tuple[HaplotypeAssemblyArtifact, ...]
    conversions: tuple[AssemblyConversion, ...]
    hifiasm_version: str
    command: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "raw_gfas": [str(path) for path in self.raw_gfas],
            "assemblies": [item.to_dict() for item in self.assemblies],
            "conversions": [item.to_dict() for item in self.conversions],
            "hifiasm_version": self.hifiasm_version,
            "command": list(self.command),
        }


def convert_gfa_to_fasta(
    source: Path,
    destination: Path,
    *,
    overwrite: bool = False,
) -> AssemblyConversion:
    """Stream GFA S records to FASTA with an atomic, auditable conversion."""
    validation.validate_output_file(source)
    destination = Path(destination)
    if destination.exists() and not overwrite:
        raise OutputValidationError(f"Assembly FASTA already exists: '{destination}'.")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.hifivar.tmp")
    if temporary.exists():
        raise OutputValidationError(f"Owned conversion temporary already exists: '{temporary}'.")
    records = 0
    try:
        with Path(source).open("r", encoding="utf-8") as reader, temporary.open(
            "x", encoding="utf-8", newline="\n"
        ) as writer:
            for number, line in enumerate(reader, start=1):
                if not line.startswith("S\t"):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 3 or not fields[1] or fields[2] in {"", "*"}:
                    raise OutputValidationError(
                        f"GFA segment at line {number} lacks an embedded sequence: '{source}'."
                    )
                writer.write(f">{fields[1]}\n{fields[2]}\n")
                records += 1
        if records == 0:
            raise OutputValidationError(f"GFA contains no sequence-bearing segments: '{source}'.")
        os.replace(temporary, destination)
    except Exception:
        if temporary.exists():
            temporary.unlink()
        raise
    validation.validate_output_file(destination)
    return AssemblyConversion(Path(source), destination)


def build_assembly_artifact(
    request: AssemblyRequest,
    *,
    hifiasm_version: str,
    command: tuple[str, ...],
) -> AssemblyArtifact:
    """Validate expected HiFi-only GFAs and derive three explicit FASTAs."""
    raw_paths = request.raw_gfa_paths
    fasta_paths = request.fasta_paths
    conversions: list[AssemblyConversion] = []
    assemblies: list[HaplotypeAssemblyArtifact] = []
    for role in (
        AssemblyRole.PRIMARY,
        AssemblyRole.HAPLOTYPE1,
        AssemblyRole.HAPLOTYPE2,
    ):
        source = raw_paths[role]
        conversion = convert_gfa_to_fasta(
            source,
            fasta_paths[role],
            overwrite=request.overwrite,
        )
        conversions.append(conversion)
        assemblies.append(
            HaplotypeAssemblyArtifact(
                request.sample_id,
                role,
                conversion.output_fasta,
                source,
                hifiasm_version,
                conversion.output_fasta.stat().st_size,
            )
        )
    return AssemblyArtifact(
        request.sample_id,
        tuple(raw_paths.values()),
        tuple(assemblies),
        tuple(conversions),
        hifiasm_version,
        command,
    )


__all__ = [
    "AssemblyArtifact",
    "AssemblyConversion",
    "AssemblyRequest",
    "AssemblyResources",
    "AssemblyRole",
    "HaplotypeAssemblyArtifact",
    "build_assembly_artifact",
    "convert_gfa_to_fasta",
]
