"""Minimap2 alignment wrapper executed exclusively through CommandRunner.

Implements the full minimap2 -> samtools sort -> samtools index chain with:

R2-MAJOR-003: Atomic output replacement using .partial paths.
R2-MAJOR-004: Comprehensive output validation (quickcheck, sort, RG/SM, contigs, index).
R2-MAJOR-006: Structured provenance for each step.
R2-MAJOR-007: Separate logs per step, unique task-owned intermediate paths.
"""

from __future__ import annotations

import os
import re
import time
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from uuid import uuid4

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
    CommandExecutionError,
    ConfigurationError,
    InputValidationError,
    OutputValidationError,
    ToolVersionError,
)
from hifivar.logging_utils import get_logger


_LOGGER = get_logger(__name__)
_VERSION_PATTERN = re.compile(
    r"minimap2\s+([0-9]+(?:\.[0-9]+){1,2}(?:[-+._A-Za-z0-9]*)?)",
    re.IGNORECASE,
)
_VERSION_BARE_PATTERN = re.compile(
    r"^([0-9]+(?:\.[0-9]+){1,2}(?:[-+._A-Za-z0-9]*)?)\s*$",
    re.MULTILINE,
)
_HIFI_PRESETS = frozenset({"map-hifi", "asm20", "asm10", "asm5"})
_OUTPUT_FORMATS = frozenset({"bam"})
_LOG_LEVELS = frozenset({"DEBUG", "INFO", "WARNING", "ERROR"})


@dataclass(frozen=True, slots=True)
class Minimap2Options:
    """Validated minimap2 command settings."""

    preset: str = "map-hifi"
    output_format: str = "bam"
    log_level: str = "INFO"
    extra_flags: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Restrict the wrapper to explicit presets and output formats."""
        if not isinstance(self.preset, str):
            raise ConfigurationError("minimap2 preset must be a string.")
        preset = self.preset.lower()
        if preset not in _HIFI_PRESETS:
            allowed = ", ".join(sorted(_HIFI_PRESETS))
            raise ConfigurationError(
                f"minimap2 preset must be one of: {allowed}."
            )
        if not isinstance(self.output_format, str):
            raise ConfigurationError("minimap2 output_format must be a string.")
        output_format = self.output_format.lower()
        if output_format not in _OUTPUT_FORMATS:
            allowed = ", ".join(sorted(_OUTPUT_FORMATS))
            raise ConfigurationError(
                f"minimap2 output_format must be one of: {allowed}."
            )
        if not isinstance(self.log_level, str):
            raise ConfigurationError("minimap2 log_level must be a string.")
        log_level = self.log_level.upper()
        if log_level not in _LOG_LEVELS:
            allowed = ", ".join(sorted(_LOG_LEVELS))
            raise ConfigurationError(
                f"minimap2 log_level must be one of: {allowed}."
            )
        if not isinstance(self.extra_flags, tuple):
            raise ConfigurationError("minimap2 extra_flags must be a tuple.")
        for flag in self.extra_flags:
            if not isinstance(flag, str) or not flag.strip():
                raise ConfigurationError(
                    "minimap2 extra_flags must contain non-empty strings."
                )
        object.__setattr__(self, "preset", preset)
        object.__setattr__(self, "output_format", output_format)
        object.__setattr__(self, "log_level", log_level)

    @classmethod
    def from_config(cls, config: Mapping[str, object]) -> Minimap2Options:
        """Load tool settings from the validated alignment config section."""
        if not isinstance(config, Mapping):
            raise ConfigurationError("minimap2 options require a config mapping.")
        section = config.get("alignment")
        if not isinstance(section, Mapping):
            raise ConfigurationError(
                "minimap2 options require configuration section alignment."
            )
        extra = section.get("minimap2_extra_flags", ())
        if isinstance(extra, str):
            extra = tuple(extra.split()) if extra.strip() else ()
        elif not isinstance(extra, (list, tuple)):
            extra = ()
        return cls(
            preset=section.get("minimap2_preset", "map-hifi"),  # type: ignore[arg-type]
            output_format=section.get("minimap2_output_format", "bam"),  # type: ignore[arg-type]
            log_level=section.get("minimap2_log_level", "INFO"),  # type: ignore[arg-type]
            extra_flags=tuple(str(f) for f in extra),
        )

    def to_dict(self) -> dict[str, object]:
        """Return standard tool-option provenance."""
        return {
            "preset": self.preset,
            "output_format": self.output_format,
            "log_level": self.log_level,
            "extra_flags": list(self.extra_flags),
        }


@dataclass(frozen=True, slots=True)
class StepProvenance:
    """Structured provenance for one external command step."""

    step_name: str
    command: tuple[str, ...]
    returncode: int
    tool: str
    tool_version: str | None
    duration_seconds: float
    stdout_path: Path | None
    stderr_path: Path | None

    def to_dict(self) -> dict[str, object]:
        """Return JSON-friendly provenance for one step."""
        return {
            "step_name": self.step_name,
            "command": list(self.command),
            "returncode": self.returncode,
            "tool": self.tool,
            "tool_version": self.tool_version,
            "duration_seconds": self.duration_seconds,
            "stdout_path": str(self.stdout_path) if self.stdout_path else None,
            "stderr_path": str(self.stderr_path) if self.stderr_path else None,
        }


@dataclass(frozen=True, slots=True)
class AlignmentChainProvenance:
    """Complete structured provenance for the minimap2 alignment chain."""

    sample_id: str
    minimap2: StepProvenance
    samtools_sort: StepProvenance
    samtools_index: StepProvenance

    def to_dict(self) -> dict[str, object]:
        """Return JSON-friendly chain provenance."""
        return {
            "sample_id": self.sample_id,
            "chain": "minimap2 -> samtools sort -> samtools index",
            "minimap2": self.minimap2.to_dict(),
            "samtools_sort": self.samtools_sort.to_dict(),
            "samtools_index": self.samtools_index.to_dict(),
        }


class Minimap2Wrapper:
    """Build and execute sorted HiFi FASTQ-to-BAM minimap2 alignments.

    This wrapper does not assume minimap2 can directly consume uBAM input.
    uBAM must be converted to FASTQ upstream before alignment.

    R2-MAJOR-003 fix: Uses atomic replacement with .partial paths.
    R2-MAJOR-004 fix: Validates quickcheck, sort, RG/SM, contigs, index.
    R2-MAJOR-006 fix: Returns structured provenance for each step.
    R2-MAJOR-007 fix: Separate logs per step, unique intermediate paths.
    """

    tool = AlignmentTool.MINIMAP2

    def __init__(
        self,
        *,
        runner: CommandRunner | None = None,
        executable: str | Path = "minimap2",
        samtools_executable: str | Path = "samtools",
        options: Minimap2Options | None = None,
    ) -> None:
        """Create a wrapper without checking the executable until execution."""
        if not isinstance(executable, (str, Path)) or not str(executable).strip():
            raise InputValidationError(
                "minimap2 executable must be a non-empty string or Path."
            )
        if not isinstance(samtools_executable, (str, Path)) or not str(samtools_executable).strip():
            raise InputValidationError(
                "samtools executable must be a non-empty string or Path."
            )
        self.runner = runner or CommandRunner()
        self.executable = str(executable)
        self.samtools_executable = str(samtools_executable)
        self.options = options or Minimap2Options()

    def detect_version(self) -> str:
        """Require minimap2 and return its parsed version string."""
        executable_path = self.runner.require_executable(self.executable)
        _LOGGER.debug("Detected minimap2 executable: %s", executable_path)
        result = self.runner.run([self.executable, "--version"])
        output = "\n".join(
            text for text in (result.stdout, result.stderr) if text is not None
        )
        match = _VERSION_PATTERN.search(output)
        if match is None:
            match = _VERSION_BARE_PATTERN.search(output.strip())
        if match is None:
            raise ToolVersionError(
                "Unable to parse minimap2 version from `minimap2 --version` output."
            )
        version = match.group(1)
        _LOGGER.info("Detected minimap2 version: %s", version)
        return version

    def detect_samtools_version(self) -> str:
        """Require samtools and return its parsed version string."""
        from hifivar.alignment_header import SamtoolsHeaderParser
        parser = SamtoolsHeaderParser(
            runner=self.runner,
            executable=self.samtools_executable,
        )
        return parser.detect_version()

    def build_command(self, request: AlignmentRequest) -> list[str]:
        """Build a sorted BAM command as an argument list without execution.

        The command uses samtools-style piping for sorting:
        minimap2 ... | samtools sort -o output.bam

        Since we use CommandRunner with shell=False, we construct
        the minimap2 command with -a flag and explicit FASTQ input paths.
        """
        self._validate_request_contract(request)
        # minimap2 parses literal backslash-t escapes in its -R argument.
        # Actual tab characters belong in the resulting SAM header only.
        read_group = (
            f"@RG\\tID:{request.sample.sample_id}\\t"
            f"SM:{request.sample.sample_id}\\tPL:PACBIO"
        )
        if request.resources is None:
            raise InputValidationError("minimap2 request resources are missing.")

        args = [
            self.executable,
            "-a",
            "-x",
            self.options.preset,
            "-t",
            str(request.resources.threads),
            "--eqx",
            "-Y",
            "--MD",
            "-R",
            read_group,
            str(request.reference.fasta),
        ]

        args.extend(self.options.extra_flags)

        for fastq_path in request.input_paths:
            args.append(str(fastq_path))

        return args

    def build_samtools_sort_command(
        self,
        request: AlignmentRequest,
        input_mode: str = "file",
        input_path: Path | None = None,
        output_path: Path | None = None,
    ) -> list[str]:
        """Build the samtools sort command for coordinate sorting."""
        if request.resources is None:
            raise InputValidationError("minimap2 request resources are missing.")

        actual_output = output_path or request.output_path
        args = [
            self.samtools_executable,
            "sort",
            "-@",
            str(max(1, request.resources.threads - 1)),
            "-o",
            str(actual_output),
        ]

        if input_mode == "pipe":
            args.append("-")
        elif input_mode == "file":
            if input_path is None:
                raise InputValidationError(
                    "samtools sort input_mode='file' requires input_path."
                )
            args.append(str(input_path))
        else:
            raise InputValidationError(
                f"Unsupported samtools sort input_mode: '{input_mode}'."
            )

        return args

    def plan_command(self, request: AlignmentRequest) -> AlignmentCommandPlan:
        """Return a reproducible dry-run command plan."""
        return AlignmentCommandPlan(self.tool, tuple(self.build_command(request)))

    def run(
        self,
        request: AlignmentRequest,
        *,
        dry_run: bool = False,
        timeout: float | None = None,
        redact_values: Collection[str] | None = None,
        stderr_path: str | Path | None = None,
    ) -> AlignmentResult:
        """Validate, preview or execute minimap2, then validate its BAM output.

        R2-MAJOR-003: Uses .partial paths for atomic replacement.
        R2-MAJOR-004: Validates all properties before returning COMPLETED.
        R2-MAJOR-006: Returns structured provenance for each step.
        R2-MAJOR-007: Separate logs per step, consistent timeout/redaction.
        """
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
                    "minimap2 dry-run unexpectedly reported command execution."
                )
            _LOGGER.info(
                "Planned minimap2 alignment: sample=%s output=%s",
                request.sample.sample_id,
                request.output_path,
            )
            return AlignmentResult(
                plan=plan,
                status=AlignmentResultStatus.PLANNED,
                command=command_plan,
            )

        version = self.detect_version()
        samtools_version = self.detect_samtools_version()
        self._prepare_output(request)

        # R2-MAJOR-007: Unique task-owned intermediate paths
        output_stem = request.output_path.stem
        task_id = (
            f"{request.sample.sample_id}_{output_stem}_{uuid4().hex}"
        )
        intermediate_sam = request.output_path.parent / f".{task_id}.intermediate.sam"

        # R2-MAJOR-007: Separate stderr logs per step. When the caller supplies
        # a log path, keep minimap2 at that exact path and derive the two
        # samtools logs beside it.
        if stderr_path is not None:
            minimap2_stderr = Path(stderr_path)
            stderr_base = minimap2_stderr.parent
            log_prefix = minimap2_stderr.name
            samtools_sort_stderr = stderr_base / f"{log_prefix}.samtools_sort"
            samtools_index_stderr = stderr_base / f"{log_prefix}.samtools_index"
        else:
            stderr_base = request.output_path.parent / ".logs"
            minimap2_stderr = stderr_base / f".{task_id}.minimap2.stderr"
            samtools_sort_stderr = stderr_base / f".{task_id}.samtools_sort.stderr"
            samtools_index_stderr = stderr_base / f".{task_id}.samtools_index.stderr"
        stderr_base.mkdir(parents=True, exist_ok=True)

        partial_bam = request.output_path.with_suffix(f".{task_id}.partial.bam")
        partial_bai = Path(f"{partial_bam}.bai")

        try:
            provenance = self._execute_pipeline(
                request,
                intermediate_sam=intermediate_sam,
                partial_bam=partial_bam,
                partial_bai=partial_bai,
                timeout=timeout,
                redact_values=redact_values,
                minimap2_stderr=minimap2_stderr,
                samtools_sort_stderr=samtools_sort_stderr,
                samtools_index_stderr=samtools_index_stderr,
                version=version,
                samtools_version=samtools_version,
            )

            # R2-MAJOR-004: Comprehensive output validation
            self._validate_completed_output(
                partial_bam, partial_bai, request
            )

            # R2-MAJOR-003: Atomic replacement
            self._atomic_replace(
                partial_bam, partial_bai, request
            )

        except Exception:
            # R2-MAJOR-007: Cleanup partials on failure, preserve existing final
            for partial in (partial_bam, partial_bai, intermediate_sam):
                if partial.exists():
                    try:
                        partial.unlink()
                    except OSError:
                        pass
            raise

        _LOGGER.info(
            "minimap2 alignment completed: sample=%s minimap2_version=%s "
            "samtools_version=%s output=%s",
            request.sample.sample_id,
            version,
            samtools_version,
            request.output_path,
        )
        return AlignmentResult(
            plan=plan,
            status=AlignmentResultStatus.COMPLETED,
            command=command_plan,
            tool_version=f"minimap2={version};samtools={samtools_version}",
            duration_seconds=sum(
                step.duration_seconds
                for step in (
                    provenance.minimap2,
                    provenance.samtools_sort,
                    provenance.samtools_index,
                )
            ),
            provenance=provenance.to_dict(),
        )

    def _execute_pipeline(
        self,
        request: AlignmentRequest,
        *,
        intermediate_sam: Path,
        partial_bam: Path,
        partial_bai: Path,
        timeout: float | None,
        redact_values: Collection[str] | None,
        minimap2_stderr: Path,
        samtools_sort_stderr: Path,
        samtools_index_stderr: Path,
        version: str,
        samtools_version: str,
    ) -> AlignmentChainProvenance:
        """Execute minimap2 -> samtools sort -> samtools index with structured provenance.

        R2-MAJOR-006: Each step returns structured provenance.
        R2-MAJOR-007: Separate logs per step, consistent timeout.
        """
        # Step 1: Run minimap2 to produce SAM
        minimap2_args = self.build_command(request)
        start = time.perf_counter()
        result_minimap2 = self.runner.run(
            minimap2_args,
            stdout_path=str(intermediate_sam),
            stderr_path=str(minimap2_stderr),
            timeout=timeout,
            redact_values=redact_values,
            check=True,
        )
        minimap2_duration = time.perf_counter() - start

        minimap2_provenance = StepProvenance(
            step_name="minimap2",
            command=tuple(minimap2_args),
            returncode=result_minimap2.returncode or 0,
            tool="minimap2",
            tool_version=version,
            duration_seconds=minimap2_duration,
            stdout_path=intermediate_sam,
            stderr_path=minimap2_stderr,
        )

        # Step 2: Run samtools sort to produce coordinate-sorted BAM
        samtools_args = self.build_samtools_sort_command(
            request, input_mode="file", input_path=intermediate_sam,
            output_path=partial_bam,
        )
        start = time.perf_counter()
        result_sort = self.runner.run(
            samtools_args,
            stderr_path=str(samtools_sort_stderr),
            timeout=timeout,
            redact_values=redact_values,
            check=True,
        )
        sort_duration = time.perf_counter() - start

        sort_provenance = StepProvenance(
            step_name="samtools_sort",
            command=tuple(samtools_args),
            returncode=result_sort.returncode or 0,
            tool="samtools",
            tool_version=samtools_version,
            duration_seconds=sort_duration,
            stdout_path=None,
            stderr_path=samtools_sort_stderr,
        )

        # Step 3: Index the output BAM
        index_args = [
            self.samtools_executable,
            "index",
            str(partial_bam),
            str(partial_bai),
        ]
        start = time.perf_counter()
        result_index = self.runner.run(
            index_args,
            stderr_path=str(samtools_index_stderr),
            timeout=timeout,
            redact_values=redact_values,
            check=True,
        )
        index_duration = time.perf_counter() - start

        index_provenance = StepProvenance(
            step_name="samtools_index",
            command=tuple(index_args),
            returncode=result_index.returncode or 0,
            tool="samtools",
            tool_version=samtools_version,
            duration_seconds=index_duration,
            stdout_path=None,
            stderr_path=samtools_index_stderr,
        )

        # Clean up intermediate SAM
        if intermediate_sam.exists():
            try:
                intermediate_sam.unlink()
            except OSError:
                _LOGGER.warning(
                    "Failed to remove intermediate SAM: %s", intermediate_sam
                )

        return AlignmentChainProvenance(
            sample_id=request.sample.sample_id,
            minimap2=minimap2_provenance,
            samtools_sort=sort_provenance,
            samtools_index=index_provenance,
        )

    def _validate_completed_output(
        self,
        bam_path: Path,
        bai_path: Path,
        request: AlignmentRequest,
    ) -> None:
        """R2-MAJOR-004: Comprehensive output validation.

        Validates: file existence, non-empty, quickcheck, sort order,
        RG/SM presence, contigs, correct index type/readability.
        """
        # File existence and non-empty
        validation.validate_output_file(bam_path)
        validation.validate_output_file(bai_path)

        # Quickcheck
        quickcheck_args = [
            self.samtools_executable,
            "quickcheck",
            "-v",
            str(bam_path),
        ]
        result = self.runner.run(quickcheck_args, check=False)
        if result.returncode != 0:
            raise OutputValidationError(
                f"samtools quickcheck failed for '{bam_path}': "
                f"{result.stderr or result.stdout or 'unknown error'}"
            )

        # Sort order validation via header
        from hifivar.alignment_header import SamtoolsHeaderParser
        from hifivar.alignment_input import SortOrder
        parser = SamtoolsHeaderParser(
            runner=self.runner,
            executable=self.samtools_executable,
        )
        header = parser.parse_header(bam_path)

        if header.sort_order != SortOrder.COORDINATE:
            raise OutputValidationError(
                f"Output BAM '{bam_path}' sort order is "
                f"'{header.sort_order.value}', expected 'coordinate'."
            )

        # RG/SM presence
        if not header.read_group_ids:
            raise OutputValidationError(
                f"Output BAM '{bam_path}' is missing read group (RG) header lines."
            )
        if not header.sample_names:
            raise OutputValidationError(
                f"Output BAM '{bam_path}' is missing sample name (SM) in header."
            )

        # Sample identity match
        expected_sample = request.sample.sample_id
        if expected_sample not in header.sample_names:
            raise OutputValidationError(
                f"Output BAM '{bam_path}' sample names {header.sample_names} "
                f"do not include expected sample '{expected_sample}'."
            )

        # Reference contig validation
        reference_contigs = set(request.reference.contig_names)
        bam_contigs = set(header.reference_contigs)
        missing = bam_contigs - reference_contigs
        if missing:
            raise OutputValidationError(
                f"Output BAM '{bam_path}' contains contigs not in reference: "
                f"{', '.join(sorted(missing)[:5])}."
            )

        if not header.reference_contigs:
            raise OutputValidationError(
                f"Output BAM '{bam_path}' has no reference contigs."
            )
        index_query = [
            self.samtools_executable,
            "view",
            "-c",
            str(bam_path),
            f"{header.reference_contigs[0]}:1-1",
        ]
        index_result = self.runner.run(index_query, check=False)
        if index_result.returncode != 0:
            raise OutputValidationError(
                f"Output BAM index is not readable for '{bam_path}': "
                f"{index_result.stderr or index_result.stdout or 'unknown error'}"
            )

        _LOGGER.debug(
            "Output validation passed: bam=%s quickcheck=OK sort=coordinate "
            "rg=%s sm=%s contigs=%d",
            bam_path,
            header.read_group_ids,
            header.sample_names,
            len(header.reference_contigs),
        )

    def _atomic_replace(
        self,
        partial_bam: Path,
        partial_bai: Path,
        request: AlignmentRequest,
    ) -> None:
        """R2-MAJOR-003: Atomically replace final outputs with validated partials.

        On failure, existing final artifacts are preserved.
        """
        final_bam = request.output_path
        final_bai = Path(f"{final_bam}.bai")

        backup_token = uuid4().hex
        backup_bam = final_bam.parent / f".{final_bam.name}.{backup_token}.backup"
        backup_bai = final_bai.parent / f".{final_bai.name}.{backup_token}.backup"
        had_bam = final_bam.exists()
        had_bai = final_bai.exists()
        backed_up_bam = False
        backed_up_bai = False
        published_bam = False
        published_bai = False

        try:
            if had_bam:
                self._replace_path(final_bam, backup_bam)
                backed_up_bam = True
            if had_bai:
                self._replace_path(final_bai, backup_bai)
                backed_up_bai = True
            self._replace_path(partial_bam, final_bam)
            published_bam = True
            self._replace_path(partial_bai, final_bai)
            published_bai = True
        except OSError as error:
            rollback_errors: list[str] = []
            for was_published, published in (
                (published_bam, final_bam),
                (published_bai, final_bai),
            ):
                if was_published and published.exists():
                    try:
                        published.unlink()
                    except OSError as rollback_error:
                        rollback_errors.append(str(rollback_error))
            for was_backed_up, backup, destination in (
                (backed_up_bam, backup_bam, final_bam),
                (backed_up_bai, backup_bai, final_bai),
            ):
                if was_backed_up and backup.exists():
                    try:
                        self._replace_path(backup, destination)
                    except OSError as rollback_error:
                        rollback_errors.append(str(rollback_error))
            details = (
                f" Rollback errors: {'; '.join(rollback_errors)}"
                if rollback_errors
                else ""
            )
            raise OutputValidationError(
                "Failed to publish validated BAM/BAI pair atomically: "
                f"{error}.{details}"
            ) from error
        else:
            for backup in (backup_bam, backup_bai):
                if backup.exists():
                    try:
                        backup.unlink()
                    except OSError:
                        _LOGGER.warning("Failed to remove output backup: %s", backup)

        _LOGGER.debug(
            "Atomic replacement completed: bam=%s bai=%s",
            final_bam,
            final_bai,
        )

    @staticmethod
    def _replace_path(source: Path, destination: Path) -> None:
        """Replace one path through a testable filesystem boundary."""
        source.replace(destination)

    def _validate_request_contract(self, request: AlignmentRequest) -> None:
        """Require the generic request subset supported by minimap2."""
        if not isinstance(request, AlignmentRequest):
            raise InputValidationError(
                "minimap2 requires an AlignmentRequest instance."
            )
        if request.tool is not AlignmentTool.MINIMAP2:
            raise InputValidationError(
                f"minimap2 cannot handle alignment tool {request.tool.value}."
            )
        if request.output_format is not AlignmentOutputFormat.BAM:
            raise InputValidationError(
                "minimap2 CRAM output is deferred; request BAM output."
            )
        if request.resources is None:
            raise InputValidationError("minimap2 request resources are missing.")

    def _validate_inputs(self, request: AlignmentRequest) -> None:
        """Recheck files immediately before planning or execution."""
        self._validate_request_contract(request)
        validation.validate_fasta(request.reference.fasta, require_fai=True)
        for fastq in request.input_paths:
            validation.validate_fastq(fastq)

    def _prepare_output(self, request: AlignmentRequest) -> None:
        """Recheck overwrite policy and create only the requested parent."""
        output = request.output_path
        index = Path(f"{output}.bai")
        if output.exists() and output.is_dir():
            raise OutputValidationError(
                f"minimap2 output path is a directory: '{output}'."
            )
        if output.exists() and not request.overwrite:
            raise OutputValidationError(
                f"minimap2 output already exists: '{output}'."
            )
        if index.exists() and not request.overwrite:
            raise OutputValidationError(
                f"minimap2 output index already exists: '{index}'."
            )
        try:
            output.parent.mkdir(parents=True, exist_ok=True)
        except OSError as error:
            raise OutputValidationError(
                f"Unable to prepare minimap2 output '{output}': {error}"
            ) from error


def _plan_for_request(request: AlignmentRequest) -> AlignmentPlan:
    """Build the generic plan represented by one minimap2 request."""
    return AlignmentPlan(
        sample=request.sample,
        reference=request.reference,
        action=AlignmentAction.ALIGN,
        alignment_path=request.output_path,
        output_format=request.output_format,
        request=request,
    )


__all__ = [
    "AlignmentChainProvenance",
    "Minimap2Options",
    "Minimap2Wrapper",
    "StepProvenance",
]
