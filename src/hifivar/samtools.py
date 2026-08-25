"""Minimal samtools wrapper for explicit alignment indexing only."""

from __future__ import annotations

import re
from collections.abc import Collection
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hifivar import validation
from hifivar.alignment_postprocess import (
    AlignmentArtifact,
    AlignmentIndexFormat,
    AlignmentIndexRequest,
)
from hifivar.command import CommandRunner
from hifivar.exceptions import (
    InputValidationError,
    OutputValidationError,
    ToolVersionError,
)
from hifivar.logging_utils import get_logger


_LOGGER = get_logger(__name__)
_VERSION_PATTERN = re.compile(
    r"\bsamtools\s+([0-9]+(?:\.[0-9]+){1,2}(?:[-+._A-Za-z0-9]*)?)",
    re.IGNORECASE,
)


class IndexResultStatus(str, Enum):
    """Successful indexing states."""

    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class SamtoolsIndexResult:
    """Dry-run or completed samtools index result."""

    request: AlignmentIndexRequest
    status: IndexResultStatus
    command: tuple[str, ...]
    artifact: AlignmentArtifact
    tool_version: str | None = None
    duration_seconds: float | None = None

    def __post_init__(self) -> None:
        """Validate result state and immutable command metadata."""
        if not isinstance(self.request, AlignmentIndexRequest):
            raise InputValidationError(
                "Samtools index result requires AlignmentIndexRequest."
            )
        if not isinstance(self.status, IndexResultStatus):
            raise InputValidationError(
                "Samtools index result status must be IndexResultStatus."
            )
        command = tuple(self.command)
        if not command or any(not isinstance(arg, str) for arg in command):
            raise InputValidationError(
                "Samtools index result command must contain strings."
            )
        object.__setattr__(self, "command", command)
        if not isinstance(self.artifact, AlignmentArtifact):
            raise InputValidationError(
                "Samtools index result artifact must be AlignmentArtifact."
            )
        if self.status is IndexResultStatus.PLANNED:
            if self.artifact.index_path is not None:
                raise InputValidationError(
                    "Planned indexing cannot claim an index artifact."
                )
        elif self.artifact.index_path != self.request.output_path:
            raise InputValidationError(
                "Completed indexing result must attach its requested index."
            )

    @property
    def executed(self) -> bool:
        """Return whether samtools index completed."""
        return self.status is IndexResultStatus.COMPLETED

    def to_dict(self) -> dict[str, object]:
        """Return JSON/YAML-friendly indexing provenance."""
        return {
            "request": self.request.to_dict(),
            "status": self.status.value,
            "command": list(self.command),
            "artifact": self.artifact.to_dict(),
            "tool_version": self.tool_version,
            "duration_seconds": self.duration_seconds,
            "executed": self.executed,
        }


class SamtoolsWrapper:
    """Execute only `samtools index` through CommandRunner."""

    def __init__(
        self,
        *,
        runner: CommandRunner | None = None,
        executable: str | Path = "samtools",
    ) -> None:
        """Create a lazy wrapper without requiring samtools on Windows."""
        if not isinstance(executable, (str, Path)) or not str(executable).strip():
            raise InputValidationError(
                "samtools executable must be a non-empty string or Path."
            )
        self.runner = runner or CommandRunner()
        self.executable = str(executable)

    def detect_version(self) -> str:
        """Require samtools and parse its version output."""
        executable_path = self.runner.require_executable(self.executable)
        _LOGGER.debug("Detected samtools executable: %s", executable_path)
        result = self.runner.run([self.executable, "--version"])
        output = "\n".join(
            text for text in (result.stdout, result.stderr) if text is not None
        )
        match = _VERSION_PATTERN.search(output)
        if match is None:
            raise ToolVersionError(
                "Unable to parse samtools version from `samtools --version`."
            )
        version = match.group(1)
        _LOGGER.info("Detected samtools version: %s", version)
        return version

    def build_index_command(self, request: AlignmentIndexRequest) -> list[str]:
        """Build one explicit-output indexing command."""
        if not isinstance(request, AlignmentIndexRequest):
            raise InputValidationError(
                "samtools index requires AlignmentIndexRequest."
            )
        command = [
            self.executable,
            "index",
            "-@",
            str(request.threads),
        ]
        if request.index_format is AlignmentIndexFormat.BAI:
            command.append("-b")
        elif request.index_format is AlignmentIndexFormat.CSI:
            command.append("-c")
        command.extend((str(request.artifact.path), str(request.output_path)))
        return command

    def plan_index_command(
        self,
        request: AlignmentIndexRequest,
    ) -> tuple[str, ...]:
        """Return the shell-free dry-run argument tuple."""
        return tuple(self.build_index_command(request))

    def run_index(
        self,
        request: AlignmentIndexRequest,
        *,
        dry_run: bool = False,
        timeout: float | None = None,
        redact_values: Collection[str] | None = None,
        stderr_path: str | Path | None = None,
    ) -> SamtoolsIndexResult:
        """Validate, preview or execute indexing, then validate the index."""
        if not isinstance(request, AlignmentIndexRequest):
            raise InputValidationError(
                "samtools index requires AlignmentIndexRequest."
            )
        command = self.plan_index_command(request)

        if dry_run:
            result = self.runner.run(
                command,
                dry_run=True,
                timeout=timeout,
                redact_values=redact_values,
                stderr_path=stderr_path,
            )
            if result.executed:
                raise OutputValidationError(
                    "samtools index dry-run unexpectedly reported execution."
                )
            return SamtoolsIndexResult(
                request=request,
                status=IndexResultStatus.PLANNED,
                command=command,
                artifact=request.artifact,
            )

        validation.validate_output_file(request.artifact.path)
        version = self.detect_version()
        self._prepare_output(request)
        result = self.runner.run(
            command,
            timeout=timeout,
            redact_values=redact_values,
            stderr_path=stderr_path,
        )
        if not result.executed or result.returncode != 0:
            raise OutputValidationError(
                "samtools index returned without successful execution."
            )
        validation.validate_output_file(request.output_path)
        artifact = request.artifact.with_index(request.output_path)
        _LOGGER.info(
            "samtools indexing completed: sample=%s version=%s runtime=%.3fs "
            "index=%s",
            artifact.sample_id,
            version,
            result.duration_seconds,
            request.output_path,
        )
        return SamtoolsIndexResult(
            request=request,
            status=IndexResultStatus.COMPLETED,
            command=command,
            artifact=artifact,
            tool_version=version,
            duration_seconds=result.duration_seconds,
        )

    def _prepare_output(self, request: AlignmentIndexRequest) -> None:
        """Recheck output races and create only the index parent directory."""
        output = request.output_path
        if output.exists() and output.is_dir():
            raise OutputValidationError(
                f"samtools index output is a directory: '{output}'."
            )
        if output.exists() and not request.overwrite:
            raise OutputValidationError(
                f"samtools index output already exists: '{output}'."
            )
        try:
            if output.exists():
                output.unlink()
            output.parent.mkdir(parents=True, exist_ok=True)
        except OSError as error:
            raise OutputValidationError(
                f"Unable to prepare samtools index output '{output}': {error}"
            ) from error


__all__ = [
    "IndexResultStatus",
    "SamtoolsIndexResult",
    "SamtoolsWrapper",
]
