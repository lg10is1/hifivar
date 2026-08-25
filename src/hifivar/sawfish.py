"""Sawfish single-sample structural-variant wrapper."""

from __future__ import annotations

import os
import re
import shutil
import tempfile
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
_VERSION_PATTERN = re.compile(r"sawfish(?:\s+version)?\s+v?(\d+(?:\.\d+)+)", re.I)


@dataclass(frozen=True, slots=True)
class SawfishResources:
    """Scheduler-neutral resources used by both Sawfish steps."""

    threads: int = 16
    memory_mb: int = 32_000
    runtime_minutes: int = 1_440

    def __post_init__(self) -> None:
        for name in ("threads", "memory_mb", "runtime_minutes"):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
                raise InputValidationError(f"Sawfish {name} must be a positive integer.")

    def to_dict(self) -> dict[str, int]:
        return {name: getattr(self, name) for name in ("threads", "memory_mb", "runtime_minutes")}


@dataclass(frozen=True, slots=True)
class SawfishRequest:
    """Inputs and deterministic paths for one Sawfish discover/joint-call run."""

    artifact: AlignmentArtifact
    output_vcf: Path
    work_directory: Path
    resources: SawfishResources = SawfishResources()
    overwrite: bool = False
    disable_cnv: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.artifact, AlignmentArtifact):
            raise InputValidationError("Sawfish artifact must be AlignmentArtifact.")
        if not isinstance(self.resources, SawfishResources):
            raise InputValidationError("Sawfish resources must be SawfishResources.")
        if not isinstance(self.overwrite, bool) or not isinstance(self.disable_cnv, bool):
            raise InputValidationError("Sawfish overwrite and disable_cnv must be boolean.")
        output = _coerce_path(self.output_vcf, "Sawfish output VCF")
        work = _coerce_path(self.work_directory, "Sawfish work directory")
        if not str(output).lower().endswith(".sawfish.sv.vcf.gz"):
            raise InputValidationError("Sawfish output must end with '.sawfish.sv.vcf.gz'.")
        protected = {self.artifact.path.absolute(), self.artifact.reference.fasta.absolute()}
        if self.artifact.index_path is not None:
            protected.add(self.artifact.index_path.absolute())
        if output.absolute() in protected:
            raise InputValidationError("Sawfish output must not replace an input file.")
        object.__setattr__(self, "output_vcf", output)
        object.__setattr__(self, "work_directory", work)

    @property
    def sample_id(self) -> str:
        return self.artifact.sample_id

    @property
    def output_index(self) -> Path:
        return Path(f"{self.output_vcf}.tbi")

    @property
    def discover_directory(self) -> Path:
        return self.work_directory / "discover"

    @property
    def joint_call_directory(self) -> Path:
        return self.work_directory / "joint-call"

    @property
    def native_vcf(self) -> Path:
        return self.joint_call_directory / "genotyped.sv.vcf.gz"

    @property
    def native_index(self) -> Path:
        return Path(f"{self.native_vcf}.tbi")

    @classmethod
    def create(
        cls,
        artifact: AlignmentArtifact,
        output_directory: str | Path,
        work_directory: str | Path,
        *,
        resources: SawfishResources | None = None,
        overwrite: bool = False,
        disable_cnv: bool = False,
    ) -> "SawfishRequest":
        output_root = _coerce_path(output_directory, "Sawfish output directory")
        work_root = _coerce_path(work_directory, "Sawfish work root")
        sample = artifact.sample_id if isinstance(artifact, AlignmentArtifact) else "sample"
        return cls(
            artifact=artifact,
            output_vcf=output_root / f"{sample}.sawfish.sv.vcf.gz",
            work_directory=work_root / sample,
            resources=resources or SawfishResources(),
            overwrite=overwrite,
            disable_cnv=disable_cnv,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "alignment": self.artifact.to_dict(),
            "output_vcf": str(self.output_vcf),
            "output_index": str(self.output_index),
            "work_directory": str(self.work_directory),
            "resources": self.resources.to_dict(),
            "overwrite": self.overwrite,
            "disable_cnv": self.disable_cnv,
        }


@dataclass(frozen=True, slots=True)
class SawfishCommandPlan:
    args: tuple[str, ...]
    display_command: str
    step: str
    shell: bool = False

    def to_dict(self) -> dict[str, object]:
        return {"step": self.step, "args": list(self.args), "display_command": self.display_command, "shell": self.shell}


class SawfishResultStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class SawfishResult:
    request: SawfishRequest
    status: SawfishResultStatus
    commands: tuple[SawfishCommandPlan, SawfishCommandPlan]
    tool_version: str | None = None
    duration_seconds: float = 0.0

    @property
    def executed(self) -> bool:
        return self.status is SawfishResultStatus.COMPLETED

    def to_dict(self) -> dict[str, object]:
        return {
            "request": self.request.to_dict(),
            "status": self.status.value,
            "executed": self.executed,
            "commands": [command.to_dict() for command in self.commands],
            "tool_version": self.tool_version,
            "duration_seconds": self.duration_seconds,
        }


class SawfishWrapper:
    """Execute official Sawfish discover and single-sample joint-call steps."""

    def __init__(self, *, runner: CommandRunner | None = None, executable: str = "sawfish") -> None:
        if not isinstance(executable, str) or not executable.strip():
            raise InputValidationError("Sawfish executable must be non-empty.")
        self.runner = runner or CommandRunner()
        self.executable = executable

    def build_commands(self, request: SawfishRequest) -> tuple[list[str], list[str]]:
        self._validate_contract(request)
        discover = [
            self.executable,
            "discover",
            "--threads",
            str(request.resources.threads),
            "--ref",
            str(request.artifact.reference.fasta.absolute()),
            "--bam",
            str(request.artifact.path.absolute()),
            "--output-dir",
            str(request.discover_directory.absolute()),
        ]
        if request.disable_cnv:
            discover.append("--disable-cnv")
        joint = [
            self.executable,
            "joint-call",
            "--threads",
            str(request.resources.threads),
            "--ref",
            str(request.artifact.reference.fasta.absolute()),
            "--sample",
            str(request.discover_directory.absolute()),
            "--output-dir",
            str(request.joint_call_directory.absolute()),
        ]
        if request.overwrite:
            discover.append("--clobber")
            joint.append("--clobber")
        return discover, joint

    def plan_commands(self, request: SawfishRequest, *, redact_values: Collection[str] | None = None) -> tuple[SawfishCommandPlan, SawfishCommandPlan]:
        discover, joint = self.build_commands(request)
        return (
            SawfishCommandPlan(tuple(discover), format_command(discover, redact_values=redact_values), "discover"),
            SawfishCommandPlan(tuple(joint), format_command(joint, redact_values=redact_values), "joint-call"),
        )

    def detect_version(self) -> str:
        self.runner.require_executable(self.executable)
        result = self.runner.run([self.executable, "--version"])
        output = "\n".join(part for part in (result.stdout, result.stderr) if part)
        match = _VERSION_PATTERN.search(output)
        if match is None:
            raise ToolVersionError(f"Unable to parse Sawfish version from {output.strip()!r}.")
        return match.group(1)

    def run(
        self,
        request: SawfishRequest,
        *,
        dry_run: bool = False,
        timeout: float | None = None,
        redact_values: Collection[str] | None = None,
        stderr_path: str | Path | None = None,
    ) -> SawfishResult:
        self._validate_inputs(request)
        commands = self.plan_commands(request, redact_values=redact_values)
        if dry_run:
            for command in commands:
                self.runner.run(command.args, dry_run=True, timeout=timeout, redact_values=redact_values, stderr_path=_step_log(stderr_path, command.step))
            return SawfishResult(request, SawfishResultStatus.PLANNED, commands)
        version = self.detect_version()
        self._prepare_outputs(request)
        duration = 0.0
        for command in commands:
            result = self.runner.run(command.args, timeout=timeout, redact_values=redact_values, stderr_path=_step_log(stderr_path, command.step))
            duration += result.duration_seconds
        validation.validate_output_file(request.native_vcf)
        validation.validate_output_file(request.native_index)
        _copy_atomic(request.native_vcf, request.output_vcf)
        _copy_atomic(request.native_index, request.output_index)
        validation.validate_output_file(request.output_vcf)
        validation.validate_output_file(request.output_index)
        _LOGGER.info("Sawfish completed sample=%s version=%s runtime=%.3fs", request.sample_id, version, duration)
        return SawfishResult(request, SawfishResultStatus.COMPLETED, commands, version, duration)

    @staticmethod
    def _validate_contract(request: SawfishRequest) -> None:
        if not isinstance(request, SawfishRequest):
            raise InputValidationError("Sawfish wrapper requires SawfishRequest.")

    def _validate_inputs(self, request: SawfishRequest) -> None:
        self._validate_contract(request)
        validation.validate_fasta(request.artifact.reference.fasta, require_fai=True)
        validate_alignment_artifact(request.artifact, require_index=True)

    @staticmethod
    def _prepare_outputs(request: SawfishRequest) -> None:
        for path in (request.output_vcf, request.output_index):
            if path.exists() and not request.overwrite:
                raise OutputValidationError(f"Sawfish output already exists: '{path}'.")
            if path.exists() and path.is_dir():
                raise OutputValidationError(f"Sawfish output is a directory: '{path}'.")
        for directory in (request.discover_directory, request.joint_call_directory):
            if directory.exists() and not request.overwrite:
                raise OutputValidationError(f"Sawfish work directory already exists: '{directory}'.")
        request.output_vcf.parent.mkdir(parents=True, exist_ok=True)
        request.work_directory.parent.mkdir(parents=True, exist_ok=True)


def _copy_atomic(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent, delete=False) as handle:
            temporary = Path(handle.name)
            with source.open("rb") as source_handle:
                shutil.copyfileobj(source_handle, handle)
        os.replace(temporary, destination)
    except OSError as error:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
        raise OutputValidationError(f"Unable to materialize Sawfish output '{destination}': {error}") from error


def _coerce_path(value: str | Path, label: str) -> Path:
    if not isinstance(value, (str, Path)) or (isinstance(value, str) and not value.strip()):
        raise InputValidationError(f"{label} must be a non-empty string or Path.")
    return Path(value).expanduser()


def _step_log(path: str | Path | None, step: str) -> Path | None:
    if path is None:
        return None
    log = Path(path)
    return log.with_name(f"{log.stem}.{step}{log.suffix}")


__all__ = [
    "SawfishCommandPlan",
    "SawfishRequest",
    "SawfishResources",
    "SawfishResult",
    "SawfishResultStatus",
    "SawfishWrapper",
]
