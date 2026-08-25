"""Sniffles2 single-sample structural-variant wrapper."""

from __future__ import annotations

import re
from collections.abc import Collection
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hifivar import validation
from hifivar.alignment_postprocess import AlignmentArtifact, validate_alignment_artifact
from hifivar.command import CommandRunner, format_command
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar.logging_utils import get_logger


_LOGGER = get_logger(__name__)
_VERSION_PATTERN = re.compile(r"(?:sniffles2?|version)\D*v?(\d+(?:\.\d+)+)", re.I)


@dataclass(frozen=True, slots=True)
class Sniffles2Resources:
    threads: int = 8
    memory_mb: int = 16_000
    runtime_minutes: int = 720

    def __post_init__(self) -> None:
        for name in ("threads", "memory_mb", "runtime_minutes"):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
                raise InputValidationError(f"Sniffles2 {name} must be a positive integer.")

    def to_dict(self) -> dict[str, int]:
        return {name: getattr(self, name) for name in ("threads", "memory_mb", "runtime_minutes")}


@dataclass(frozen=True, slots=True)
class Sniffles2Request:
    artifact: AlignmentArtifact
    output_vcf: Path
    resources: Sniffles2Resources = Sniffles2Resources()
    minimum_support: int | None = None
    minimum_sv_length: int = 50
    overwrite: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.artifact, AlignmentArtifact):
            raise InputValidationError("Sniffles2 artifact must be AlignmentArtifact.")
        if not isinstance(self.resources, Sniffles2Resources):
            raise InputValidationError("Sniffles2 resources must be Sniffles2Resources.")
        if self.minimum_support is not None and (
            not isinstance(self.minimum_support, int) or isinstance(self.minimum_support, bool) or self.minimum_support <= 0
        ):
            raise InputValidationError("Sniffles2 minimum_support must be null or a positive integer.")
        if not isinstance(self.minimum_sv_length, int) or isinstance(self.minimum_sv_length, bool) or self.minimum_sv_length <= 0:
            raise InputValidationError("Sniffles2 minimum_sv_length must be a positive integer.")
        if not isinstance(self.overwrite, bool):
            raise InputValidationError("Sniffles2 overwrite must be boolean.")
        output = _coerce_path(self.output_vcf, "Sniffles2 output VCF")
        if not str(output).lower().endswith(".sniffles2.sv.vcf.gz"):
            raise InputValidationError("Sniffles2 output must end with '.sniffles2.sv.vcf.gz'.")
        object.__setattr__(self, "output_vcf", output)

    @property
    def sample_id(self) -> str:
        return self.artifact.sample_id

    @property
    def output_index(self) -> Path:
        return Path(f"{self.output_vcf}.tbi")

    @classmethod
    def create(cls, artifact: AlignmentArtifact, output_directory: str | Path, **kwargs: object) -> "Sniffles2Request":
        root = _coerce_path(output_directory, "Sniffles2 output directory")
        sample = artifact.sample_id if isinstance(artifact, AlignmentArtifact) else "sample"
        return cls(artifact=artifact, output_vcf=root / f"{sample}.sniffles2.sv.vcf.gz", **kwargs)  # type: ignore[arg-type]

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "alignment": self.artifact.to_dict(),
            "output_vcf": str(self.output_vcf),
            "output_index": str(self.output_index),
            "resources": self.resources.to_dict(),
            "minimum_support": self.minimum_support,
            "minimum_sv_length": self.minimum_sv_length,
            "overwrite": self.overwrite,
        }


@dataclass(frozen=True, slots=True)
class Sniffles2CommandPlan:
    args: tuple[str, ...]
    display_command: str
    shell: bool = False

    def to_dict(self) -> dict[str, object]:
        return {"args": list(self.args), "display_command": self.display_command, "shell": self.shell}


class Sniffles2ResultStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class Sniffles2Result:
    request: Sniffles2Request
    status: Sniffles2ResultStatus
    command: Sniffles2CommandPlan
    tool_version: str | None = None
    duration_seconds: float = 0.0

    def to_dict(self) -> dict[str, object]:
        return {
            "request": self.request.to_dict(),
            "status": self.status.value,
            "executed": self.status is Sniffles2ResultStatus.COMPLETED,
            "command": self.command.to_dict(),
            "tool_version": self.tool_version,
            "duration_seconds": self.duration_seconds,
        }


class Sniffles2Wrapper:
    def __init__(self, *, runner: CommandRunner | None = None, executable: str = "sniffles") -> None:
        if not isinstance(executable, str) or not executable.strip():
            raise InputValidationError("Sniffles2 executable must be non-empty.")
        self.runner = runner or CommandRunner()
        self.executable = executable

    def build_command(self, request: Sniffles2Request) -> list[str]:
        self._validate_contract(request)
        command = [
            self.executable,
            "--input",
            str(request.artifact.path.absolute()),
            "--vcf",
            str(request.output_vcf.absolute()),
            "--reference",
            str(request.artifact.reference.fasta.absolute()),
            "--threads",
            str(request.resources.threads),
            "--sample-id",
            request.sample_id,
            "--minsvlen",
            str(request.minimum_sv_length),
        ]
        if request.minimum_support is not None:
            command.extend(("--minsupport", str(request.minimum_support)))
        return command

    def plan_command(self, request: Sniffles2Request, *, redact_values: Collection[str] | None = None) -> Sniffles2CommandPlan:
        command = self.build_command(request)
        return Sniffles2CommandPlan(tuple(command), format_command(command, redact_values=redact_values))

    def detect_version(self) -> str:
        self.runner.require_executable(self.executable)
        result = self.runner.run([self.executable, "--version"])
        output = "\n".join(part for part in (result.stdout, result.stderr) if part)
        match = _VERSION_PATTERN.search(output)
        if match is None:
            raise ToolVersionError(f"Unable to parse Sniffles2 version from {output.strip()!r}.")
        return match.group(1)

    def run(self, request: Sniffles2Request, *, dry_run: bool = False, timeout: float | None = None, redact_values: Collection[str] | None = None, stderr_path: str | Path | None = None) -> Sniffles2Result:
        self._validate_inputs(request)
        command = self.plan_command(request, redact_values=redact_values)
        if dry_run:
            self.runner.run(command.args, dry_run=True, timeout=timeout, redact_values=redact_values, stderr_path=stderr_path)
            return Sniffles2Result(request, Sniffles2ResultStatus.PLANNED, command)
        version = self.detect_version()
        self._prepare_outputs(request)
        result = self.runner.run(command.args, timeout=timeout, redact_values=redact_values, stderr_path=stderr_path)
        validation.validate_output_file(request.output_vcf)
        validation.validate_output_file(request.output_index)
        _LOGGER.info("Sniffles2 completed sample=%s version=%s runtime=%.3fs", request.sample_id, version, result.duration_seconds)
        return Sniffles2Result(request, Sniffles2ResultStatus.COMPLETED, command, version, result.duration_seconds)

    @staticmethod
    def _validate_contract(request: Sniffles2Request) -> None:
        if not isinstance(request, Sniffles2Request):
            raise InputValidationError("Sniffles2 wrapper requires Sniffles2Request.")

    def _validate_inputs(self, request: Sniffles2Request) -> None:
        self._validate_contract(request)
        validation.validate_fasta(request.artifact.reference.fasta, require_fai=True)
        validate_alignment_artifact(request.artifact, require_index=True)

    @staticmethod
    def _prepare_outputs(request: Sniffles2Request) -> None:
        for output in (request.output_vcf, request.output_index):
            if output.exists() and output.is_dir():
                raise OutputValidationError(f"Sniffles2 output is a directory: '{output}'.")
            if output.exists() and not request.overwrite:
                raise OutputValidationError(f"Sniffles2 output already exists: '{output}'.")
            if output.exists():
                output.unlink()
        request.output_vcf.parent.mkdir(parents=True, exist_ok=True)


def _coerce_path(value: str | Path, label: str) -> Path:
    if not isinstance(value, (str, Path)) or (isinstance(value, str) and not value.strip()):
        raise InputValidationError(f"{label} must be a non-empty string or Path.")
    return Path(value).expanduser()


__all__ = ["Sniffles2CommandPlan", "Sniffles2Request", "Sniffles2Resources", "Sniffles2Result", "Sniffles2ResultStatus", "Sniffles2Wrapper"]
