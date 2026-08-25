"""HiPhase single-sample small-variant phasing wrapper."""

from __future__ import annotations

import re
from collections.abc import Collection
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hifivar import validation
from hifivar.command import CommandRunner, format_command
from hifivar.exceptions import OutputValidationError, ToolVersionError
from hifivar.logging_utils import get_logger
from hifivar.phasing import (
    PhasedVariantArtifact,
    PhasingRequest,
    validate_phased_variant_output,
)


_LOGGER = get_logger(__name__)
_HIPHASE_VERSION = re.compile(
    r"(?:HiPhase|hiphase)(?:\s+version)?\s+v?(\d+(?:\.\d+)+)",
    re.I,
)


class PhasingResultStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class PhasingCommandPlan:
    args: tuple[str, ...]
    display_command: str

    def to_dict(self) -> dict[str, object]:
        return {
            "args": list(self.args),
            "display_command": self.display_command,
            "shell": False,
        }


@dataclass(frozen=True, slots=True)
class PhasingResult:
    request: PhasingRequest
    status: PhasingResultStatus
    commands: tuple[PhasingCommandPlan, ...]
    hiphase_version: str | None = None
    runtime_seconds: float = 0.0
    artifact: PhasedVariantArtifact | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "request": self.request.to_dict(),
            "status": self.status.value,
            "commands": [item.to_dict() for item in self.commands],
            "hiphase_version": self.hiphase_version,
            "runtime_seconds": self.runtime_seconds,
            "artifact": self.artifact.to_dict() if self.artifact else None,
        }


class HiPhaseWrapper:
    """Plan and run the official HiPhase CLI through CommandRunner."""

    def __init__(
        self,
        *,
        executable: str = "hiphase",
        tabix_executable: str = "tabix",
        runner: CommandRunner | None = None,
    ) -> None:
        self.executable = executable
        self.tabix_executable = tabix_executable
        self.runner = runner or CommandRunner()

    def plan_commands(
        self,
        request: PhasingRequest,
        *,
        redact_values: Collection[str] | None = None,
    ) -> tuple[PhasingCommandPlan, ...]:
        self._validate_inputs(request)
        hiphase = (
            self.executable,
            "--bam",
            str(request.alignment.path.absolute()),
            "--vcf",
            str(request.small_variants.vcf_path.absolute()),
            "--output-vcf",
            str(request.output_vcf.absolute()),
            "--reference",
            str(request.alignment.reference.fasta.absolute()),
            "--threads",
            str(request.resources.threads),
            "--sample-name",
            request.sample_id,
            "--disable-global-realignment",
        )
        tabix = (
            self.tabix_executable,
            "-f",
            "-p",
            "vcf",
            str(request.output_vcf.absolute()),
        )
        return tuple(
            PhasingCommandPlan(args, format_command(args, redact_values=redact_values))
            for args in (hiphase, tabix)
        )

    def detect_version(self) -> str:
        self.runner.require_executable(self.executable)
        result = self.runner.run([self.executable, "--version"])
        output = "\n".join(item for item in (result.stdout, result.stderr) if item)
        match = _HIPHASE_VERSION.search(output)
        if match is None:
            raise ToolVersionError(
                f"Unable to parse HiPhase version from {output!r}."
            )
        return match.group(1)

    def run(
        self,
        request: PhasingRequest,
        *,
        dry_run: bool = False,
        stderr_path: Path | None = None,
    ) -> PhasingResult:
        commands = self.plan_commands(request)
        if dry_run:
            for command in commands:
                self.runner.run(command.args, dry_run=True)
            return PhasingResult(request, PhasingResultStatus.PLANNED, commands)

        version = self.detect_version()
        self.runner.require_executable(self.tabix_executable)
        self._validate_inputs(request)
        request.output_vcf.parent.mkdir(parents=True, exist_ok=True)
        total = 0.0
        for index, command in enumerate(commands):
            result = self.runner.run(
                command.args,
                stderr_path=_numbered_log(stderr_path, index, len(commands)),
            )
            total += result.duration_seconds
            if index == 0:
                validation.validate_output_file(request.output_vcf)

        artifact = validate_phased_variant_output(
            request,
            hiphase_version=version,
            command=commands[0].args,
        )
        _LOGGER.info("HiPhase completed for sample %s", request.sample_id)
        return PhasingResult(
            request,
            PhasingResultStatus.COMPLETED,
            commands,
            version,
            total,
            artifact,
        )

    @staticmethod
    def _validate_inputs(request: PhasingRequest) -> None:
        validation.validate_file(request.alignment.path, require_nonempty=True)
        if request.alignment.index_path is None:
            raise OutputValidationError("HiPhase BAM index is missing.")
        validation.validate_file(
            request.alignment.index_path,
            require_nonempty=True,
        )
        validation.validate_file(
            request.small_variants.vcf_path,
            require_nonempty=True,
        )
        validation.validate_file(
            request.small_variants.vcf_index_path,
            require_nonempty=True,
        )
        validation.validate_file(request.alignment.reference.fasta, require_nonempty=True)
        validation.validate_file(request.alignment.reference.fai, require_nonempty=True)


def _numbered_log(path: Path | None, index: int, count: int) -> Path | None:
    if path is None or count == 1:
        return path
    path = Path(path)
    return path.with_name(f"{path.stem}.{index + 1}{path.suffix}")


__all__ = [
    "HiPhaseWrapper",
    "PhasingCommandPlan",
    "PhasingResult",
    "PhasingResultStatus",
]
