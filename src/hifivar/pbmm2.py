"""pbmm2 alignment wrapper executed exclusively through CommandRunner."""

from __future__ import annotations

import os
import re
import tempfile
from collections.abc import Collection, Mapping
from dataclasses import dataclass
from pathlib import Path

from hifivar import validation
from hifivar.alignment import (
    AlignmentAction,
    AlignmentCommandPlan,
    AlignmentOutputFormat,
    AlignmentPlan,
    AlignmentRequest,
    AlignmentResult,
    AlignmentResultStatus,
    AlignmentTool,
)
from hifivar.command import CommandRunner
from hifivar.exceptions import (
    ConfigurationError,
    InputValidationError,
    OutputValidationError,
    ToolVersionError,
)
from hifivar.logging_utils import get_logger


_LOGGER = get_logger(__name__)
_VERSION_PATTERN = re.compile(
    r"\bpbmm2\s+([0-9]+(?:\.[0-9]+){1,2}(?:[-+._A-Za-z0-9]*)?)",
    re.IGNORECASE,
)
_HIFI_PRESETS = frozenset({"CCS", "HIFI"})
_LOG_LEVELS = frozenset({"TRACE", "DEBUG", "INFO", "WARN", "FATAL"})


@dataclass(frozen=True, slots=True)
class Pbmm2Options:
    """Validated HiFi-specific pbmm2 command settings."""

    preset: str = "CCS"
    log_level: str = "INFO"

    def __post_init__(self) -> None:
        """Restrict the wrapper to explicit HiFi presets and log levels."""
        if not isinstance(self.preset, str):
            raise ConfigurationError("pbmm2 preset must be a string.")
        preset = self.preset.upper()
        if preset not in _HIFI_PRESETS:
            allowed = ", ".join(sorted(_HIFI_PRESETS))
            raise ConfigurationError(
                f"pbmm2 preset must be a HiFi preset: {allowed}."
            )
        if not isinstance(self.log_level, str):
            raise ConfigurationError("pbmm2 log_level must be a string.")
        log_level = self.log_level.upper()
        if log_level not in _LOG_LEVELS:
            allowed = ", ".join(sorted(_LOG_LEVELS))
            raise ConfigurationError(
                f"pbmm2 log_level must be one of: {allowed}."
            )
        object.__setattr__(self, "preset", preset)
        object.__setattr__(self, "log_level", log_level)

    @classmethod
    def from_config(cls, config: Mapping[str, object]) -> Pbmm2Options:
        """Load tool settings from the validated alignment config section."""
        if not isinstance(config, Mapping):
            raise ConfigurationError("pbmm2 options require a config mapping.")
        section = config.get("alignment")
        if not isinstance(section, Mapping):
            raise ConfigurationError(
                "pbmm2 options require configuration section alignment."
            )
        return cls(
            preset=section.get("pbmm2_preset", "CCS"),  # type: ignore[arg-type]
            log_level=section.get("pbmm2_log_level", "INFO"),  # type: ignore[arg-type]
        )

    def to_dict(self) -> dict[str, str]:
        """Return standard tool-option provenance."""
        return {"preset": self.preset, "log_level": self.log_level}


class Pbmm2Wrapper:
    """Build and execute sorted HiFi FASTQ-to-BAM pbmm2 alignments."""

    tool = AlignmentTool.PBMM2

    def __init__(
        self,
        *,
        runner: CommandRunner | None = None,
        executable: str | Path = "pbmm2",
        options: Pbmm2Options | None = None,
    ) -> None:
        """Create a wrapper without checking the executable until execution."""
        if not isinstance(executable, (str, Path)) or not str(executable).strip():
            raise InputValidationError(
                "pbmm2 executable must be a non-empty string or Path."
            )
        self.runner = runner or CommandRunner()
        self.executable = str(executable)
        self.options = options or Pbmm2Options()

    def detect_version(self) -> str:
        """Require pbmm2 and return its parsed version string."""
        executable_path = self.runner.require_executable(self.executable)
        _LOGGER.debug("Detected pbmm2 executable: %s", executable_path)
        result = self.runner.run([self.executable, "--version"])
        output = "\n".join(
            text for text in (result.stdout, result.stderr) if text is not None
        )
        match = _VERSION_PATTERN.search(output)
        if match is None:
            raise ToolVersionError(
                "Unable to parse pbmm2 version from `pbmm2 --version` output."
            )
        version = match.group(1)
        _LOGGER.info("Detected pbmm2 version: %s", version)
        return version

    def build_command(self, request: AlignmentRequest) -> list[str]:
        """Build a sorted BAM command as an argument list without execution."""
        self._validate_request_contract(request)
        read_group = (
            f"@RG\tID:{request.sample.sample_id}\t"
            f"SM:{request.sample.sample_id}\tPL:PACBIO"
        )
        if request.resources is None:  # pragma: no cover - model invariant
            raise InputValidationError("pbmm2 request resources are missing.")
        return [
            self.executable,
            "align",
            str(request.reference.fasta),
            str(self.input_argument(request)),
            str(request.output_path),
            "--preset",
            self.options.preset,
            "--sort",
            "--bam-index",
            "NONE",
            "--rg",
            read_group,
            "-j",
            str(request.resources.threads),
            "--log-level",
            self.options.log_level,
        ]

    def plan_command(self, request: AlignmentRequest) -> AlignmentCommandPlan:
        """Return a reproducible dry-run command plan."""
        return AlignmentCommandPlan(self.tool, tuple(self.build_command(request)))

    def input_argument(self, request: AlignmentRequest) -> Path:
        """Return a FASTQ path or deterministic multi-FASTQ FOFN path."""
        self._validate_request_contract(request)
        if len(request.input_paths) == 1:
            return request.input_paths[0]
        return request.output_path.with_name(
            f"{request.sample.sample_id}.fastq.fofn"
        )

    def run(
        self,
        request: AlignmentRequest,
        *,
        dry_run: bool = False,
        timeout: float | None = None,
        redact_values: Collection[str] | None = None,
        stderr_path: str | Path | None = None,
    ) -> AlignmentResult:
        """Validate, preview or execute pbmm2, then validate its BAM output."""
        self._validate_inputs(request)
        command_plan = self.plan_command(request)
        plan = _plan_for_request(request)

        if dry_run:
            command_result = self.runner.run(
                command_plan.args,
                dry_run=True,
                timeout=timeout,
                redact_values=redact_values,
                stderr_path=stderr_path,
            )
            if command_result.executed:
                raise OutputValidationError(
                    "pbmm2 dry-run unexpectedly reported command execution."
                )
            _LOGGER.info(
                "Planned pbmm2 alignment: sample=%s output=%s",
                request.sample.sample_id,
                request.output_path,
            )
            return AlignmentResult(
                plan=plan,
                status=AlignmentResultStatus.PLANNED,
                command=command_plan,
            )

        version = self.detect_version()
        self._prepare_output(request)
        if len(request.input_paths) > 1:
            _write_fofn(
                self.input_argument(request),
                request.input_paths,
                overwrite=request.overwrite,
            )

        command_result = self.runner.run(
            command_plan.args,
            timeout=timeout,
            redact_values=redact_values,
            stderr_path=stderr_path,
        )
        if not command_result.executed or command_result.returncode != 0:
            raise OutputValidationError(
                "pbmm2 returned without a successful executed command result."
            )
        validation.validate_output_file(request.output_path)
        _LOGGER.info(
            "pbmm2 alignment completed: sample=%s version=%s runtime=%.3fs "
            "output=%s",
            request.sample.sample_id,
            version,
            command_result.duration_seconds,
            request.output_path,
        )
        return AlignmentResult(
            plan=plan,
            status=AlignmentResultStatus.COMPLETED,
            command=command_plan,
            tool_version=version,
            duration_seconds=command_result.duration_seconds,
        )

    def _validate_request_contract(self, request: AlignmentRequest) -> None:
        """Require the generic request subset supported by pbmm2 Phase 2.3."""
        if not isinstance(request, AlignmentRequest):
            raise InputValidationError(
                "pbmm2 requires an AlignmentRequest instance."
            )
        if request.tool is not AlignmentTool.PBMM2:
            raise InputValidationError(
                f"pbmm2 cannot handle alignment tool {request.tool.value}."
            )
        if request.output_format is not AlignmentOutputFormat.BAM:
            raise InputValidationError("pbmm2 Phase 2.3 output must be BAM.")
        if request.resources is None:  # pragma: no cover - model invariant
            raise InputValidationError("pbmm2 request resources are missing.")

    def _validate_inputs(self, request: AlignmentRequest) -> None:
        """Recheck files immediately before planning or execution."""
        self._validate_request_contract(request)
        validation.validate_fasta(request.reference.fasta, require_fai=True)
        for fastq in request.input_paths:
            validation.validate_fastq(fastq)

    def _prepare_output(self, request: AlignmentRequest) -> None:
        """Recheck overwrite policy and create only the requested parent."""
        output = request.output_path
        if output.exists() and output.is_dir():
            raise OutputValidationError(
                f"pbmm2 output path is a directory: '{output}'."
            )
        if output.exists() and not request.overwrite:
            raise OutputValidationError(
                f"pbmm2 output already exists: '{output}'."
            )
        try:
            if output.exists():
                output.unlink()
            output.parent.mkdir(parents=True, exist_ok=True)
        except OSError as error:
            raise OutputValidationError(
                f"Unable to prepare pbmm2 output '{output}': {error}"
            ) from error


def _plan_for_request(request: AlignmentRequest) -> AlignmentPlan:
    """Build the generic plan represented by one pbmm2 request."""
    return AlignmentPlan(
        sample=request.sample,
        reference=request.reference,
        action=AlignmentAction.ALIGN,
        alignment_path=request.output_path,
        output_format=request.output_format,
        request=request,
    )


def _write_fofn(
    path: Path,
    fastq_paths: tuple[Path, ...],
    *,
    overwrite: bool,
) -> None:
    """Atomically write absolute ordered FASTQ paths for pbmm2."""
    if path.exists() and not overwrite:
        raise OutputValidationError(
            f"pbmm2 FASTQ FOFN already exists: '{path}'."
        )
    content = "".join(f"{item.absolute()}\n" for item in fastq_paths)
    temporary_path: Path | None = None
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            mode="wt",
            encoding="utf-8",
            newline="\n",
            prefix=f".{path.name}.",
            suffix=".tmp",
            dir=path.parent,
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    except OSError as error:
        if temporary_path is not None:
            try:
                temporary_path.unlink(missing_ok=True)
            except OSError:
                pass
        raise OutputValidationError(
            f"Unable to write pbmm2 FASTQ FOFN '{path}': {error}"
        ) from error


__all__ = ["Pbmm2Options", "Pbmm2Wrapper"]
