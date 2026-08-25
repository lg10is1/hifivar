"""ANNOVAR table_annovar.pl adapter executed only through CommandRunner."""

from __future__ import annotations

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
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.validation import validate_output_file, validate_vcf


_OPERATIONS = frozenset({"g", "gx", "r", "f"})


@dataclass(frozen=True, slots=True)
class AnnovarRequest:
    input: AnnotationInput
    database_root: Path
    database_version: str
    protocols: tuple[str, ...]
    operations: tuple[str, ...]
    output_prefix: Path
    tool_version: str
    overwrite: bool = False

    def __post_init__(self) -> None:
        root = Path(self.database_root).expanduser()
        if not root.is_dir():
            raise InputValidationError(f"ANNOVAR database root is missing: '{root}'.")
        if not isinstance(self.database_version, str) or not self.database_version.strip():
            raise InputValidationError("ANNOVAR database version must be explicit.")
        if not isinstance(self.tool_version, str) or not self.tool_version.strip():
            raise InputValidationError("ANNOVAR release/version must be explicit.")
        protocols, operations = tuple(self.protocols), tuple(self.operations)
        if not protocols or len(protocols) != len(operations):
            raise InputValidationError("ANNOVAR protocols and operations must be non-empty and paired.")
        if any(not isinstance(item, str) or not item.strip() or "," in item for item in protocols):
            raise InputValidationError("ANNOVAR protocol names must be non-empty and comma-free.")
        if any(item not in _OPERATIONS for item in operations):
            raise InputValidationError("ANNOVAR operations must be one of g, gx, r, or f.")
        if not isinstance(self.overwrite, bool):
            raise InputValidationError("ANNOVAR overwrite must be boolean.")
        object.__setattr__(self, "database_root", root)
        object.__setattr__(self, "output_prefix", Path(self.output_prefix).expanduser())
        object.__setattr__(self, "protocols", protocols)
        object.__setattr__(self, "operations", operations)

    @property
    def annovar_build(self) -> str:
        build = self.input.reference.build
        assert build is not None
        aliases = {"GRCh37": "hg19", "GRCh38": "hg38"}
        return aliases.get(build, build)

    @property
    def output_tsv(self) -> Path:
        return Path(f"{self.output_prefix}.{self.annovar_build}_multianno.txt")

    @property
    def output_vcf(self) -> Path:
        return Path(f"{self.output_prefix}.{self.annovar_build}_multianno.vcf")

    def to_dict(self) -> dict[str, object]:
        return {
            "input": self.input.to_dict(),
            "database_root": str(self.database_root),
            "database_version": self.database_version,
            "protocols": list(self.protocols),
            "operations": list(self.operations),
            "output_prefix": str(self.output_prefix),
            "output_tsv": str(self.output_tsv),
            "native_output_vcf": str(self.output_vcf),
            "tool_version": self.tool_version,
            "overwrite": self.overwrite,
            "database_download_performed": False,
        }


class AnnovarWrapper:
    """Minimal table_annovar.pl wrapper with externally managed databases."""

    def __init__(self, executable: str = "table_annovar.pl", *, runner: CommandRunner | None = None) -> None:
        if not isinstance(executable, str) or not executable.strip():
            raise InputValidationError("ANNOVAR executable must be non-empty.")
        self.executable = executable
        self.runner = runner or CommandRunner()

    def build_command(self, request: AnnovarRequest) -> tuple[str, ...]:
        return (
            self.executable,
            str(request.input.source_vcf.absolute()),
            str(request.database_root.absolute()),
            "-buildver", request.annovar_build,
            "-out", str(request.output_prefix.absolute()),
            "-protocol", ",".join(request.protocols),
            "-operation", ",".join(request.operations),
            "-nastring", ".",
            "-vcfinput",
            "-polish",
        )

    def run(
        self,
        request: AnnovarRequest,
        *,
        dry_run: bool = False,
        stdout_path: Path | None = None,
        stderr_path: Path | None = None,
    ) -> AnnotationResult:
        command = self.build_command(request)
        if dry_run:
            self.runner.run(command, dry_run=True)
            return AnnotationResult(
                request.input, AnnotationSource.ANNOVAR,
                AnnotationRunStatus.PLANNED, command,
            )
        self.runner.require_executable(self.executable)
        _protect_outputs(request)
        request.output_prefix.parent.mkdir(parents=True, exist_ok=True)
        result = self.runner.run(
            command,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
        )
        validate_output_file(request.output_tsv)
        validate_vcf(request.output_vcf)
        database = AnnotationDatabase(
            name="annovar:" + ",".join(request.protocols),
            version=request.database_version,
            path=request.database_root,
            reference_build=request.input.reference.build or request.annovar_build,
        )
        artifact = AnnotationArtifact(
            request.input,
            AnnotationSource.ANNOVAR,
            request.output_tsv,
            "tsv",
            request.tool_version,
            (database,),
            command,
        )
        return AnnotationResult(
            request.input,
            AnnotationSource.ANNOVAR,
            AnnotationRunStatus.COMPLETED,
            command,
            request.tool_version,
            result.duration_seconds,
            artifact,
        )


def _protect_outputs(request: AnnovarRequest) -> None:
    for path in (request.output_tsv, request.output_vcf):
        if path.exists() and not request.overwrite:
            raise OutputValidationError(f"ANNOVAR output already exists: '{path}'.")
        if path.exists() and path.is_dir():
            raise OutputValidationError(f"ANNOVAR output is a directory: '{path}'.")
        if path.exists():
            path.unlink()


__all__ = ["AnnovarRequest", "AnnovarWrapper"]
