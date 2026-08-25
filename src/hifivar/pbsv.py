"""PacBio pbsv discover/call wrapper for one HiFi sample."""

from __future__ import annotations

import re
from collections.abc import Collection
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hifivar import validation
from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import AlignmentArtifact, validate_alignment_artifact
from hifivar.command import CommandRunner, format_command
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar.logging_utils import get_logger


_LOGGER = get_logger(__name__)
_VERSION_PATTERN = re.compile(r"pbsv(?:\s+version)?\s+v?(\d+(?:\.\d+)+)", re.I)


@dataclass(frozen=True, slots=True)
class PbsvResources:
    threads: int = 8
    memory_mb: int = 32_000
    runtime_minutes: int = 1_440

    def __post_init__(self) -> None:
        for name in ("threads", "memory_mb", "runtime_minutes"):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
                raise InputValidationError(f"pbsv {name} must be a positive integer.")

    def to_dict(self) -> dict[str, int]:
        return {name: getattr(self, name) for name in ("threads", "memory_mb", "runtime_minutes")}


@dataclass(frozen=True, slots=True)
class PbsvRequest:
    artifact: AlignmentArtifact
    signatures_path: Path
    raw_vcf: Path
    resources: PbsvResources = PbsvResources()
    overwrite: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.artifact, AlignmentArtifact):
            raise InputValidationError("pbsv artifact must be AlignmentArtifact.")
        if self.artifact.output_format is not AlignmentOutputFormat.BAM:
            raise InputValidationError("pbsv discover currently requires BAM input; CRAM must be converted explicitly upstream.")
        if not isinstance(self.resources, PbsvResources):
            raise InputValidationError("pbsv resources must be PbsvResources.")
        if not isinstance(self.overwrite, bool):
            raise InputValidationError("pbsv overwrite must be boolean.")
        signatures = _coerce_path(self.signatures_path, "pbsv signatures")
        raw_vcf = _coerce_path(self.raw_vcf, "pbsv raw VCF")
        if not str(signatures).lower().endswith(".pbsv.svsig.gz"):
            raise InputValidationError("pbsv signatures must end with '.pbsv.svsig.gz'.")
        if not str(raw_vcf).lower().endswith(".pbsv.raw.vcf"):
            raise InputValidationError("pbsv native output must end with '.pbsv.raw.vcf'.")
        object.__setattr__(self, "signatures_path", signatures)
        object.__setattr__(self, "raw_vcf", raw_vcf)

    @property
    def sample_id(self) -> str:
        return self.artifact.sample_id

    @classmethod
    def create(cls, artifact: AlignmentArtifact, output_directory: str | Path, work_directory: str | Path, **kwargs: object) -> "PbsvRequest":
        output_root = _coerce_path(output_directory, "pbsv output directory")
        work_root = _coerce_path(work_directory, "pbsv work directory")
        sample = artifact.sample_id if isinstance(artifact, AlignmentArtifact) else "sample"
        return cls(artifact, work_root / f"{sample}.pbsv.svsig.gz", output_root / f"{sample}.pbsv.raw.vcf", **kwargs)  # type: ignore[arg-type]

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "alignment": self.artifact.to_dict(),
            "signatures_path": str(self.signatures_path),
            "raw_vcf": str(self.raw_vcf),
            "resources": self.resources.to_dict(),
            "overwrite": self.overwrite,
        }


@dataclass(frozen=True, slots=True)
class PbsvCommandPlan:
    args: tuple[str, ...]
    display_command: str
    step: str
    shell: bool = False

    def to_dict(self) -> dict[str, object]:
        return {"step": self.step, "args": list(self.args), "display_command": self.display_command, "shell": self.shell}


class PbsvResultStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class PbsvResult:
    request: PbsvRequest
    status: PbsvResultStatus
    commands: tuple[PbsvCommandPlan, PbsvCommandPlan]
    tool_version: str | None = None
    duration_seconds: float = 0.0

    def to_dict(self) -> dict[str, object]:
        return {
            "request": self.request.to_dict(),
            "status": self.status.value,
            "executed": self.status is PbsvResultStatus.COMPLETED,
            "commands": [command.to_dict() for command in self.commands],
            "tool_version": self.tool_version,
            "duration_seconds": self.duration_seconds,
        }


class PbsvWrapper:
    """Run `pbsv discover` and `pbsv call` without a shell pipeline."""

    def __init__(self, *, runner: CommandRunner | None = None, executable: str = "pbsv") -> None:
        if not isinstance(executable, str) or not executable.strip():
            raise InputValidationError("pbsv executable must be non-empty.")
        self.runner = runner or CommandRunner()
        self.executable = executable

    def build_commands(self, request: PbsvRequest) -> tuple[list[str], list[str]]:
        self._validate_contract(request)
        discover = [self.executable, "discover", str(request.artifact.path.absolute()), str(request.signatures_path.absolute())]
        call = [
            self.executable,
            "call",
            "--ccs",
            "-j",
            str(request.resources.threads),
            str(request.artifact.reference.fasta.absolute()),
            str(request.signatures_path.absolute()),
            str(request.raw_vcf.absolute()),
        ]
        return discover, call

    def plan_commands(self, request: PbsvRequest, *, redact_values: Collection[str] | None = None) -> tuple[PbsvCommandPlan, PbsvCommandPlan]:
        discover, call = self.build_commands(request)
        return (
            PbsvCommandPlan(tuple(discover), format_command(discover, redact_values=redact_values), "discover"),
            PbsvCommandPlan(tuple(call), format_command(call, redact_values=redact_values), "call"),
        )

    def detect_version(self) -> str:
        self.runner.require_executable(self.executable)
        result = self.runner.run([self.executable, "--version"])
        output = "\n".join(part for part in (result.stdout, result.stderr) if part)
        match = _VERSION_PATTERN.search(output)
        if match is None:
            raise ToolVersionError(f"Unable to parse pbsv version from {output.strip()!r}.")
        return match.group(1)

    def run(self, request: PbsvRequest, *, dry_run: bool = False, timeout: float | None = None, redact_values: Collection[str] | None = None, stderr_path: str | Path | None = None) -> PbsvResult:
        self._validate_inputs(request)
        commands = self.plan_commands(request, redact_values=redact_values)
        if dry_run:
            for command in commands:
                self.runner.run(command.args, dry_run=True, timeout=timeout, redact_values=redact_values, stderr_path=_step_log(stderr_path, command.step))
            return PbsvResult(request, PbsvResultStatus.PLANNED, commands)
        version = self.detect_version()
        self._prepare_outputs(request)
        first = self.runner.run(commands[0].args, timeout=timeout, redact_values=redact_values, stderr_path=_step_log(stderr_path, "discover"))
        validation.validate_output_file(request.signatures_path)
        second = self.runner.run(commands[1].args, timeout=timeout, redact_values=redact_values, stderr_path=_step_log(stderr_path, "call"))
        validation.validate_output_file(request.raw_vcf)
        duration = first.duration_seconds + second.duration_seconds
        _LOGGER.info("pbsv completed sample=%s version=%s runtime=%.3fs", request.sample_id, version, duration)
        return PbsvResult(request, PbsvResultStatus.COMPLETED, commands, version, duration)

    @staticmethod
    def _validate_contract(request: PbsvRequest) -> None:
        if not isinstance(request, PbsvRequest):
            raise InputValidationError("pbsv wrapper requires PbsvRequest.")

    def _validate_inputs(self, request: PbsvRequest) -> None:
        self._validate_contract(request)
        validation.validate_fasta(request.artifact.reference.fasta, require_fai=True)
        validate_alignment_artifact(request.artifact, require_index=True)

    @staticmethod
    def _prepare_outputs(request: PbsvRequest) -> None:
        for output in (request.signatures_path, request.raw_vcf):
            if output.exists() and output.is_dir():
                raise OutputValidationError(f"pbsv output is a directory: '{output}'.")
            if output.exists() and not request.overwrite:
                raise OutputValidationError(f"pbsv output already exists: '{output}'.")
            if output.exists():
                output.unlink()
            output.parent.mkdir(parents=True, exist_ok=True)


def _coerce_path(value: str | Path, label: str) -> Path:
    if not isinstance(value, (str, Path)) or (isinstance(value, str) and not value.strip()):
        raise InputValidationError(f"{label} must be a non-empty string or Path.")
    return Path(value).expanduser()


def _step_log(path: str | Path | None, step: str) -> Path | None:
    if path is None:
        return None
    log = Path(path)
    return log.with_name(f"{log.stem}.{step}{log.suffix}")


__all__ = ["PbsvCommandPlan", "PbsvRequest", "PbsvResources", "PbsvResult", "PbsvResultStatus", "PbsvWrapper"]
