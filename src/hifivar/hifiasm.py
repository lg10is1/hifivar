"""hifiasm HiFi-only single-sample assembly wrapper."""

from __future__ import annotations

import re
from collections.abc import Collection
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hifivar import validation
from hifivar.assembly import AssemblyArtifact, AssemblyRequest, build_assembly_artifact
from hifivar.command import CommandRunner, format_command
from hifivar.exceptions import OutputValidationError, ToolVersionError
from hifivar.logging_utils import get_logger


_LOGGER = get_logger(__name__)
_HIFIASM_VERSION = re.compile(
    r"(?:hifiasm(?:\s+version)?\s+)?v?(\d+(?:\.\d+)+(?:-r\d+)?)",
    re.I,
)


class AssemblyResultStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class AssemblyCommandPlan:
    args: tuple[str, ...]
    display_command: str

    def to_dict(self) -> dict[str, object]:
        return {
            "args": list(self.args),
            "display_command": self.display_command,
            "shell": False,
        }


@dataclass(frozen=True, slots=True)
class AssemblyResult:
    request: AssemblyRequest
    status: AssemblyResultStatus
    command: AssemblyCommandPlan
    hifiasm_version: str | None = None
    runtime_seconds: float = 0.0
    artifact: AssemblyArtifact | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "request": self.request.to_dict(),
            "status": self.status.value,
            "command": self.command.to_dict(),
            "hifiasm_version": self.hifiasm_version,
            "runtime_seconds": self.runtime_seconds,
            "artifact": self.artifact.to_dict() if self.artifact else None,
        }


class HifiasmWrapper:
    """Run official hifiasm HiFi-only mode through CommandRunner."""

    def __init__(
        self,
        *,
        executable: str = "hifiasm",
        runner: CommandRunner | None = None,
    ) -> None:
        self.executable = executable
        self.runner = runner or CommandRunner()

    def plan_command(
        self,
        request: AssemblyRequest,
        *,
        redact_values: Collection[str] | None = None,
    ) -> AssemblyCommandPlan:
        self._validate_inputs(request)
        args = [
            self.executable,
            "-o",
            str(request.output_prefix.absolute()),
            "-t",
            str(request.resources.threads),
        ]
        if request.overwrite:
            args.append("-i")
        args.extend(str(path.absolute()) for path in request.fastq_files)
        command = tuple(args)
        return AssemblyCommandPlan(
            command,
            format_command(command, redact_values=redact_values),
        )

    def detect_version(self) -> str:
        self.runner.require_executable(self.executable)
        result = self.runner.run([self.executable, "--version"])
        output = "\n".join(item for item in (result.stdout, result.stderr) if item)
        match = _HIFIASM_VERSION.search(output)
        if match is None:
            raise ToolVersionError(
                f"Unable to parse hifiasm version from {output!r}."
            )
        return match.group(1)

    def run(
        self,
        request: AssemblyRequest,
        *,
        dry_run: bool = False,
        stderr_path: Path | None = None,
    ) -> AssemblyResult:
        command = self.plan_command(request)
        if dry_run:
            self.runner.run(command.args, dry_run=True)
            return AssemblyResult(
                request,
                AssemblyResultStatus.PLANNED,
                command,
            )

        version = self.detect_version()
        self._validate_inputs(request)
        existing = tuple(
            request.output_prefix.parent.glob(f"{request.output_prefix.name}.*")
        )
        if existing and not request.overwrite:
            raise OutputValidationError(
                f"hifiasm prefix already owns outputs: {[str(path) for path in existing]!r}."
            )
        request.output_prefix.parent.mkdir(parents=True, exist_ok=True)
        request.assembly_directory.mkdir(parents=True, exist_ok=True)
        result = self.runner.run(command.args, stderr_path=stderr_path)
        for path in request.raw_gfa_paths.values():
            validation.validate_output_file(path)
        artifact = build_assembly_artifact(
            request,
            hifiasm_version=version,
            command=command.args,
        )
        _LOGGER.info("hifiasm completed for sample %s", request.sample_id)
        return AssemblyResult(
            request,
            AssemblyResultStatus.COMPLETED,
            command,
            version,
            result.duration_seconds,
            artifact,
        )

    @staticmethod
    def _validate_inputs(request: AssemblyRequest) -> None:
        for path in request.fastq_files:
            validation.validate_fastq(path)


__all__ = [
    "AssemblyCommandPlan",
    "AssemblyResult",
    "AssemblyResultStatus",
    "HifiasmWrapper",
]
