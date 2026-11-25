"""HiPhase single-sample phasing with optional, separately retained SV/TR VCFs."""

from __future__ import annotations

import re
from collections.abc import Collection
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hifivar import validation
from hifivar.command import CommandRunner, _default_stdout_log_path, format_command
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
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
        hiphase = [
            self.executable,
            "--bam",
            str(request.alignment.path.absolute()),
        ]
        for _, source, _, output, _, _ in request.vcf_pairs:
            hiphase.extend(("--vcf", str(source.absolute()), "--output-vcf", str(output.absolute())))
        hiphase.extend((
            "--reference",
            str(request.alignment.reference.fasta.absolute()),
            "--threads",
            str(request.resources.threads),
            "--sample-name",
            request.sample_id,
        ))
        if not request.use_global_realignment:
            hiphase.append("--disable-global-realignment")
        commands = [tuple(hiphase)]
        commands.extend((self.tabix_executable, "-f", "-p", "vcf", str(output.absolute()))
                        for _, _, _, output, _, _ in request.vcf_pairs)
        return tuple(
            PhasingCommandPlan(args, format_command(args, redact_values=redact_values))
            for args in commands
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
        self._validate_log_paths(request, stderr_path, len(commands))
        if dry_run:
            for command in commands:
                self.runner.run(command.args, dry_run=True)
            return PhasingResult(request, PhasingResultStatus.PLANNED, commands)

        version = self.detect_version()
        self.runner.require_executable(self.tabix_executable)
        self._validate_inputs(request)
        self._validate_log_paths(request, stderr_path, len(commands))
        request.output_vcf.parent.mkdir(parents=True, exist_ok=True)
        total = 0.0
        for index, command in enumerate(commands):
            result = self.runner.run(
                command.args,
                stderr_path=_numbered_log(stderr_path, index, len(commands)),
            )
            total += result.duration_seconds
            if index == 0:
                for _, _, _, output, _, _ in request.vcf_pairs:
                    validation.validate_output_file(output)

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
        request.validate_output_safety()
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
        if request.additional_inputs:
            from hifivar.joint_phasing import validate_joint_source_headers
            validate_joint_source_headers(request)

    def _validate_log_paths(self, request: PhasingRequest, stderr_path: Path | None, count: int) -> None:
        """Protect source, output and executable paths, including implicit stdout logs."""
        if stderr_path is None:
            return
        protected = list(request.protected_input_paths)
        protected.extend(p for _, _, _, output, index, _ in request.vcf_pairs for p in (output, index))
        for executable in (self.executable, self.tabix_executable):
            path = Path(executable).expanduser()
            if path.is_file():
                protected.append(path)
        actual = []
        for index in range(count):
            log = _numbered_log(stderr_path, index, count)
            actual.extend((log, _default_stdout_log_path(log)))
        for path in [Path(stderr_path), *actual]:
            if any(path.resolve() == p.resolve() or (path.exists() and p.exists() and path.samefile(p))
                   for p in protected):
                raise InputValidationError("HiPhase log path aliases a protected source, output, or executable.")
        if len({p.resolve() for p in actual}) != len(actual) or any(
            a.exists() and b.exists() and a.samefile(b)
            for i, a in enumerate(actual) for b in actual[i + 1:]
        ):
            raise InputValidationError("HiPhase log paths must be distinct.")


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
