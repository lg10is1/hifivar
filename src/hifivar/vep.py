"""Independent Ensembl VEP offline/cache wrapper."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from hifivar.annotation import (
    AnnotationArtifact,
    AnnotationDatabase,
    AnnotationInput,
    AnnotationResult,
    AnnotationRunStatus,
    AnnotationSource,
)
from hifivar.command import CommandRunner
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar.validation import validate_output_file


_VERSION = re.compile(
    r"(?:ensembl-vep|\bvep\b|version|release)\s*(?:version\s*)?[:= ]\s*v?(\d+(?:\.\d+)*)",
    re.IGNORECASE,
)


@dataclass(frozen=True, slots=True)
class VepRequest:
    input: AnnotationInput
    cache_directory: Path
    cache_version: str
    species: str
    assembly: str
    output_tsv: Path
    threads: int = 1
    overwrite: bool = False

    def __post_init__(self) -> None:
        cache = Path(self.cache_directory).expanduser()
        if not cache.is_dir():
            raise InputValidationError(f"VEP cache directory is missing: '{cache}'.")
        for name in ("cache_version", "species", "assembly"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise InputValidationError(f"VEP {name} must be non-empty.")
        if self.input.reference.build != self.assembly:
            raise InputValidationError(
                f"VEP assembly '{self.assembly}' does not match reference build "
                f"'{self.input.reference.build}'."
            )
        if not isinstance(self.threads, int) or isinstance(self.threads, bool) or self.threads <= 0:
            raise InputValidationError("VEP threads must be a positive integer.")
        if not isinstance(self.overwrite, bool):
            raise InputValidationError("VEP overwrite must be boolean.")
        object.__setattr__(self, "cache_directory", cache)
        object.__setattr__(self, "output_tsv", Path(self.output_tsv).expanduser())

    def to_dict(self) -> dict[str, object]:
        return {
            "input": self.input.to_dict(),
            "cache_directory": str(self.cache_directory),
            "cache_version": self.cache_version,
            "species": self.species,
            "assembly": self.assembly,
            "reference_fasta": str(self.input.reference.fasta),
            "output_tsv": str(self.output_tsv),
            "threads": self.threads,
            "overwrite": self.overwrite,
            "offline": True,
            "cache_download_performed": False,
        }


class VepWrapper:
    """Run VEP in explicit offline cache mode through CommandRunner."""

    def __init__(self, executable: str = "vep", *, runner: CommandRunner | None = None) -> None:
        if not isinstance(executable, str) or not executable.strip():
            raise InputValidationError("VEP executable must be non-empty.")
        self.executable = executable
        self.runner = runner or CommandRunner()

    def build_command(self, request: VepRequest) -> tuple[str, ...]:
        return (
            self.executable,
            "--input_file", str(request.input.source_vcf.absolute()),
            "--output_file", str(request.output_tsv.absolute()),
            "--format", "vcf",
            "--tab",
            "--offline",
            "--cache",
            "--dir_cache", str(request.cache_directory.absolute()),
            "--cache_version", request.cache_version,
            "--species", request.species,
            "--assembly", request.assembly,
            "--fasta", str(request.input.reference.fasta.absolute()),
            "--fork", str(request.threads),
        )

    def detect_version(self) -> str:
        self.runner.require_executable(self.executable)
        result = self.runner.run((self.executable, "--help"))
        output = "\n".join(value for value in (result.stdout, result.stderr) if value)
        match = _VERSION.search(output)
        if match is None:
            raise ToolVersionError("Unable to parse Ensembl VEP version from --help output.")
        return match.group(1)

    def run(
        self,
        request: VepRequest,
        *,
        dry_run: bool = False,
        stdout_path: Path | None = None,
        stderr_path: Path | None = None,
    ) -> AnnotationResult:
        command = self.build_command(request)
        if dry_run:
            self.runner.run(command, dry_run=True)
            return AnnotationResult(
                request.input, AnnotationSource.VEP,
                AnnotationRunStatus.PLANNED, command,
            )
        if request.output_tsv.exists() and not request.overwrite:
            raise OutputValidationError(f"VEP output already exists: '{request.output_tsv}'.")
        if request.output_tsv.exists() and request.output_tsv.is_dir():
            raise OutputValidationError(f"VEP output is a directory: '{request.output_tsv}'.")
        if request.output_tsv.exists():
            request.output_tsv.unlink()
        version = self.detect_version()
        request.output_tsv.parent.mkdir(parents=True, exist_ok=True)
        result = self.runner.run(
            command,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
        )
        validate_output_file(request.output_tsv)
        database = AnnotationDatabase(
            name=f"ensembl-vep-cache:{request.species}",
            version=request.cache_version,
            path=request.cache_directory,
            reference_build=request.assembly,
        )
        artifact = AnnotationArtifact(
            request.input,
            AnnotationSource.VEP,
            request.output_tsv,
            "tsv",
            version,
            (database,),
            command,
        )
        return AnnotationResult(
            request.input,
            AnnotationSource.VEP,
            AnnotationRunStatus.COMPLETED,
            command,
            version,
            result.duration_seconds,
            artifact,
        )


__all__ = ["VepRequest", "VepWrapper"]
