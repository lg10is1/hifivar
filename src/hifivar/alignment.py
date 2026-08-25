"""Tool-neutral alignment planning models and backend contract."""

from __future__ import annotations

import os
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Protocol

from hifivar.command import format_command
from hifivar.context import AnalysisContext
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.reference import ReferenceGenome
from hifivar.sample import InputType, Sample


PathInput = str | Path


class AlignmentTool(str, Enum):
    """Alignment implementations supported or planned by HiFiVar."""

    PBMM2 = "pbmm2"
    MINIMAP2 = "minimap2"


class AlignmentOutputFormat(str, Enum):
    """Supported primary alignment containers."""

    BAM = "bam"
    CRAM = "cram"

    @property
    def suffix(self) -> str:
        """Return the conventional filename suffix including its dot."""
        return f".{self.value}"


class AlignmentAction(str, Enum):
    """Whether a sample needs alignment or reuses an existing alignment."""

    ALIGN = "align"
    REUSE = "reuse"


class AlignmentResultStatus(str, Enum):
    """Successful alignment lifecycle states represented by this phase."""

    PLANNED = "planned"
    COMPLETED = "completed"
    REUSED = "reused"


@dataclass(frozen=True, slots=True)
class AlignmentResources:
    """Portable resource request without scheduler-specific syntax."""

    threads: int = 1
    memory_mb: int | None = None
    runtime_minutes: int | None = None

    def __post_init__(self) -> None:
        """Require positive resource values when specified."""
        _require_positive_integer(self.threads, "Alignment resource threads")
        if self.memory_mb is not None:
            _require_positive_integer(
                self.memory_mb,
                "Alignment resource memory_mb",
            )
        if self.runtime_minutes is not None:
            _require_positive_integer(
                self.runtime_minutes,
                "Alignment resource runtime_minutes",
            )

    def to_dict(self) -> dict[str, int | None]:
        """Return scheduler-neutral resource metadata."""
        return {
            "threads": self.threads,
            "memory_mb": self.memory_mb,
            "runtime_minutes": self.runtime_minutes,
        }


@dataclass(frozen=True, slots=True)
class AlignmentRequest:
    """Immutable request for aligning one HiFi FASTQ sample.

    The canonical resource field is :attr:`resources`. ``threads`` remains an
    accepted construction shorthand for the original Phase 2.2 API.
    """

    sample: Sample
    reference: ReferenceGenome
    output_path: Path
    tool: AlignmentTool
    output_format: AlignmentOutputFormat = AlignmentOutputFormat.BAM
    threads: int = 1
    resources: AlignmentResources | None = None
    overwrite: bool = False

    def __post_init__(self) -> None:
        """Normalize the output path and enforce the alignment boundary."""
        if not isinstance(self.sample, Sample):
            raise InputValidationError(
                "Alignment request sample must be a Sample instance."
            )
        if not isinstance(self.reference, ReferenceGenome):
            raise InputValidationError(
                "Alignment request reference must be a ReferenceGenome instance."
            )
        _require_positive_integer(self.threads, "Alignment request threads")
        if self.resources is None:
            resources = AlignmentResources(threads=self.threads)
        elif not isinstance(self.resources, AlignmentResources):
            raise InputValidationError(
                "Alignment request resources must be AlignmentResources or None."
            )
        else:
            resources = self.resources
            if self.threads != 1 and resources.threads != self.threads:
                raise InputValidationError(
                    "Alignment request threads conflicts with resources.threads."
                )
        object.__setattr__(self, "threads", resources.threads)
        object.__setattr__(self, "resources", resources)
        _validate_options(self.tool, self.output_format, resources)
        if not isinstance(self.overwrite, bool):
            raise InputValidationError(
                "Alignment request overwrite must be a boolean."
            )
        if self.sample.input.input_type is not InputType.FASTQ:
            raise InputValidationError(
                f"Alignment request for sample '{self.sample.sample_id}' requires "
                f"FASTQ input; received {self.sample.input.input_type.value}."
            )

        output_path = _coerce_path(self.output_path, "Alignment output")
        if output_path.suffix.lower() != self.output_format.suffix:
            raise InputValidationError(
                f"Alignment output '{output_path}' does not match requested "
                f"format {self.output_format.value.upper()}; expected suffix "
                f"'{self.output_format.suffix}'."
            )
        protected_paths = (
            *self.sample.input.files,
            self.reference.fasta,
            self.reference.fai,
        )
        output_identity = _path_identity(output_path)
        if any(output_identity == _path_identity(path) for path in protected_paths):
            raise InputValidationError(
                f"Alignment output '{output_path}' conflicts with an input or "
                "reference path."
            )
        if output_path.exists() and output_path.is_dir():
            raise OutputValidationError(
                f"Alignment output path is a directory: '{output_path}'."
            )
        if output_path.exists() and not self.overwrite:
            raise OutputValidationError(
                f"Alignment output already exists: '{output_path}'. Explicitly "
                "set overwrite=True to replace it."
            )
        object.__setattr__(self, "output_path", output_path)

    @property
    def input_paths(self) -> tuple[Path, ...]:
        """Return ordered FASTQ paths without copying file contents."""
        return self.sample.input.files

    def to_dict(self) -> dict[str, object]:
        """Return JSON/YAML-friendly request metadata."""
        if self.resources is None:  # pragma: no cover - post-init invariant
            raise InputValidationError("Alignment request resources are missing.")
        return {
            "sample_id": self.sample.sample_id,
            "tool": self.tool.value,
            "input": self.sample.input.to_dict(),
            "reference": self.reference.to_dict(),
            "output": {
                "path": str(self.output_path),
                "format": self.output_format.value,
            },
            "threads": self.threads,
            "resources": self.resources.to_dict(),
            "overwrite": self.overwrite,
        }


@dataclass(frozen=True, slots=True)
class AlignmentCommandPlan:
    """A reproducible shell-free command preview produced by a backend."""

    tool: AlignmentTool
    args: tuple[str, ...]

    def __post_init__(self) -> None:
        """Copy and validate one executable argument sequence."""
        if not isinstance(self.tool, AlignmentTool):
            raise InputValidationError(
                "Alignment command tool must be an AlignmentTool value."
            )
        if isinstance(self.args, (str, Path)):
            raise InputValidationError(
                "Alignment command args must be a sequence of strings."
            )
        args = tuple(self.args)
        if not args or any(not isinstance(arg, str) for arg in args):
            raise InputValidationError(
                "Alignment command args must contain one or more strings."
            )
        if not args[0].strip():
            raise InputValidationError(
                "Alignment command executable must not be empty."
            )
        object.__setattr__(self, "args", args)

    @property
    def display(self) -> str:
        """Return a display-only command string that is never shell input."""
        return format_command(self.args)

    def to_dict(self) -> dict[str, object]:
        """Return standard reproducibility metadata."""
        return {
            "tool": self.tool.value,
            "args": list(self.args),
            "display": self.display,
            "shell": False,
        }


@dataclass(frozen=True, slots=True)
class AlignmentPlan:
    """Per-sample decision to align FASTQ or reuse an existing BAM/CRAM."""

    sample: Sample
    reference: ReferenceGenome
    action: AlignmentAction
    alignment_path: Path
    output_format: AlignmentOutputFormat
    request: AlignmentRequest | None = None

    def __post_init__(self) -> None:
        """Ensure action, input type, and optional request agree."""
        if not isinstance(self.sample, Sample):
            raise InputValidationError("Alignment plan sample must be a Sample.")
        if not isinstance(self.reference, ReferenceGenome):
            raise InputValidationError(
                "Alignment plan reference must be a ReferenceGenome."
            )
        if not isinstance(self.action, AlignmentAction):
            raise InputValidationError(
                "Alignment plan action must be an AlignmentAction value."
            )
        if not isinstance(self.output_format, AlignmentOutputFormat):
            raise InputValidationError(
                "Alignment plan output_format must be an AlignmentOutputFormat."
            )
        path = _coerce_path(self.alignment_path, "Alignment plan")
        if path.suffix.lower() != self.output_format.suffix:
            raise InputValidationError(
                f"Alignment plan path '{path}' does not match format "
                f"{self.output_format.value.upper()}."
            )

        if self.action is AlignmentAction.ALIGN:
            if not isinstance(self.request, AlignmentRequest):
                raise InputValidationError(
                    "ALIGN plans require an AlignmentRequest."
                )
            if self.sample.input.input_type is not InputType.FASTQ:
                raise InputValidationError("ALIGN plans require FASTQ input.")
            if (
                self.request.sample != self.sample
                or self.request.reference != self.reference
                or _path_identity(self.request.output_path) != _path_identity(path)
                or self.request.output_format is not self.output_format
            ):
                raise InputValidationError(
                    "Alignment plan fields conflict with its AlignmentRequest."
                )
        else:
            if self.request is not None:
                raise InputValidationError(
                    "REUSE plans must not contain an AlignmentRequest."
                )
            if self.sample.input.input_type not in {InputType.BAM, InputType.CRAM}:
                raise InputValidationError("REUSE plans require BAM or CRAM input.")
            if _path_identity(self.sample.input.files[0]) != _path_identity(path):
                raise InputValidationError(
                    "REUSE plan path must be the existing primary alignment."
                )
        object.__setattr__(self, "alignment_path", path)

    @property
    def sample_id(self) -> str:
        """Return the planned sample identifier."""
        return self.sample.sample_id

    @property
    def requires_alignment(self) -> bool:
        """Return whether an external aligner must run."""
        return self.action is AlignmentAction.ALIGN

    @property
    def resources(self) -> AlignmentResources | None:
        """Return resources only for newly generated alignments."""
        return self.request.resources if self.request is not None else None

    def to_dict(self) -> dict[str, object]:
        """Return standard planning metadata."""
        return {
            "sample_id": self.sample_id,
            "input_type": self.sample.input.input_type.value,
            "action": self.action.value,
            "alignment_path": str(self.alignment_path),
            "output_format": self.output_format.value,
            "requires_alignment": self.requires_alignment,
            "request": self.request.to_dict() if self.request is not None else None,
        }


@dataclass(frozen=True, slots=True)
class AlignmentResult:
    """Tool-neutral successful or dry-run alignment outcome."""

    plan: AlignmentPlan
    status: AlignmentResultStatus
    command: AlignmentCommandPlan | None = None
    tool_version: str | None = None
    duration_seconds: float | None = None

    def __post_init__(self) -> None:
        """Reject contradictory action and result states."""
        if not isinstance(self.plan, AlignmentPlan):
            raise InputValidationError(
                "Alignment result plan must be an AlignmentPlan."
            )
        if not isinstance(self.status, AlignmentResultStatus):
            raise InputValidationError(
                "Alignment result status must be an AlignmentResultStatus."
            )
        if self.command is not None and not isinstance(
            self.command,
            AlignmentCommandPlan,
        ):
            raise InputValidationError(
                "Alignment result command must be an AlignmentCommandPlan or None."
            )
        if self.tool_version is not None and (
            not isinstance(self.tool_version, str) or not self.tool_version.strip()
        ):
            raise InputValidationError(
                "Alignment result tool_version must be a non-empty string or None."
            )
        if self.duration_seconds is not None and (
            isinstance(self.duration_seconds, bool)
            or not isinstance(self.duration_seconds, (int, float))
            or self.duration_seconds < 0
        ):
            raise InputValidationError(
                "Alignment result duration_seconds must be non-negative or None."
            )

        if self.plan.action is AlignmentAction.REUSE:
            if self.status is not AlignmentResultStatus.REUSED or self.command:
                raise InputValidationError(
                    "REUSE plans require REUSED status and no command."
                )
        elif self.status is AlignmentResultStatus.REUSED:
            raise InputValidationError("ALIGN plans cannot have REUSED status.")
        elif self.command is None:
            raise InputValidationError(
                "Planned or completed ALIGN results require a command."
            )

    @property
    def alignment_path(self) -> Path:
        """Return the existing or planned primary alignment path."""
        return self.plan.alignment_path

    @property
    def executed(self) -> bool:
        """Return whether a new alignment command completed."""
        return self.status is AlignmentResultStatus.COMPLETED

    def to_dict(self) -> dict[str, object]:
        """Return JSON/YAML-friendly result/provenance metadata."""
        return {
            "plan": self.plan.to_dict(),
            "status": self.status.value,
            "alignment_path": str(self.alignment_path),
            "executed": self.executed,
            "command": self.command.to_dict() if self.command is not None else None,
            "tool_version": self.tool_version,
            "duration_seconds": self.duration_seconds,
        }


class AlignmentBackend(Protocol):
    """Minimal command-construction contract for alignment wrappers."""

    @property
    def tool(self) -> AlignmentTool:
        """Return the tool implemented by this backend."""
        ...

    def build_command(self, request: AlignmentRequest) -> list[str]:
        """Build one shell-free argument list without executing it."""
        ...


def plan_alignment_command(
    backend: AlignmentBackend,
    request: AlignmentRequest,
) -> AlignmentCommandPlan:
    """Ask a backend for a dry-run command without executing external code."""
    if backend.tool is not request.tool:
        raise InputValidationError(
            f"Alignment backend {backend.tool.value} cannot handle request tool "
            f"{request.tool.value}."
        )
    command = backend.build_command(request)
    return AlignmentCommandPlan(tool=backend.tool, args=tuple(command))


def build_alignment_plans(
    context: AnalysisContext,
    output_directory: PathInput,
    *,
    tool: AlignmentTool,
    output_format: AlignmentOutputFormat = AlignmentOutputFormat.BAM,
    resources: AlignmentResources | None = None,
    threads: int | None = None,
    overwrite: bool = False,
) -> tuple[AlignmentPlan, ...]:
    """Plan FASTQ alignment and explicit BAM/CRAM reuse in context order.

    Existing BAM/CRAM inputs retain their original paths and never invoke an
    aligner. FASTQ outputs use deterministic sample names under the supplied
    output directory. The function does not create, copy, or execute anything.
    """
    if not isinstance(context, AnalysisContext):
        raise InputValidationError(
            "Alignment planning requires an AnalysisContext instance."
        )
    selected_resources = _select_resources(resources, threads)
    _validate_options(tool, output_format, selected_resources)
    if not isinstance(overwrite, bool):
        raise InputValidationError("Alignment overwrite must be a boolean.")
    output_root = _coerce_path(output_directory, "Alignment output directory")

    plans: list[AlignmentPlan] = []
    for record in context.samples:
        sample = record.sample
        if sample.input.input_type is InputType.FASTQ:
            request = AlignmentRequest(
                sample=sample,
                reference=context.reference,
                output_path=(
                    output_root / f"{sample.sample_id}.aligned{output_format.suffix}"
                ),
                tool=tool,
                output_format=output_format,
                resources=selected_resources,
                overwrite=overwrite,
            )
            plans.append(
                AlignmentPlan(
                    sample=sample,
                    reference=context.reference,
                    action=AlignmentAction.ALIGN,
                    alignment_path=request.output_path,
                    output_format=output_format,
                    request=request,
                )
            )
            continue

        existing_format = (
            AlignmentOutputFormat.BAM
            if sample.input.input_type is InputType.BAM
            else AlignmentOutputFormat.CRAM
        )
        plans.append(
            AlignmentPlan(
                sample=sample,
                reference=context.reference,
                action=AlignmentAction.REUSE,
                alignment_path=sample.input.files[0],
                output_format=existing_format,
            )
        )
    return tuple(plans)


def build_alignment_requests(
    context: AnalysisContext,
    output_directory: PathInput,
    *,
    tool: AlignmentTool,
    output_format: AlignmentOutputFormat = AlignmentOutputFormat.BAM,
    resources: AlignmentResources | None = None,
    threads: int | None = None,
    overwrite: bool = False,
) -> tuple[AlignmentRequest, ...]:
    """Build FASTQ-only requests using the original Phase 2.2 public API.

    Use :func:`build_alignment_plans` when a context may contain existing
    BAM/CRAM inputs.
    """
    if not isinstance(context, AnalysisContext):
        raise InputValidationError(
            "Alignment request planning requires an AnalysisContext instance."
        )
    incompatible = [
        f"{record.sample.sample_id} ({record.sample.input.input_type.value})"
        for record in context.samples
        if record.sample.input.input_type is not InputType.FASTQ
    ]
    if incompatible:
        details = ", ".join(incompatible)
        raise InputValidationError(
            "Alignment request planning requires FASTQ input for every sample; "
            f"incompatible samples: {details}."
        )
    plans = build_alignment_plans(
        context,
        output_directory,
        tool=tool,
        output_format=output_format,
        resources=resources,
        threads=threads,
        overwrite=overwrite,
    )
    return tuple(plan.request for plan in plans if plan.request is not None)


def _coerce_path(value: PathInput, label: str) -> Path:
    """Normalize a portable path spelling without resolving or creating it."""
    if not isinstance(value, (str, Path)):
        raise InputValidationError(f"{label} path must be a string or Path.")
    if isinstance(value, str) and not value.strip():
        raise InputValidationError(f"{label} path must not be empty.")
    return Path(value).expanduser()


def _select_resources(
    resources: AlignmentResources | None,
    threads: int | None,
) -> AlignmentResources:
    """Resolve canonical resources and the compatibility thread shorthand."""
    if resources is not None and not isinstance(resources, AlignmentResources):
        raise InputValidationError(
            "Alignment resources must be AlignmentResources or None."
        )
    if threads is not None:
        _require_positive_integer(threads, "Alignment threads")
        if resources is not None and resources.threads != threads:
            raise InputValidationError(
                "Alignment threads conflicts with resources.threads."
            )
    return resources or AlignmentResources(threads=threads or 1)


def _validate_options(
    tool: AlignmentTool,
    output_format: AlignmentOutputFormat,
    resources: AlignmentResources,
) -> None:
    """Validate shared planning options."""
    if not isinstance(tool, AlignmentTool):
        raise InputValidationError(
            "Alignment request tool must be an AlignmentTool value."
        )
    if not isinstance(output_format, AlignmentOutputFormat):
        raise InputValidationError(
            "Alignment request output_format must be an AlignmentOutputFormat "
            "value."
        )
    if not isinstance(resources, AlignmentResources):
        raise InputValidationError(
            "Alignment request resources must be AlignmentResources."
        )


def _require_positive_integer(value: object, label: str) -> None:
    """Require one non-boolean positive integer."""
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise InputValidationError(f"{label} must be a positive integer.")


def _path_identity(path: Path) -> str:
    """Return a platform-aware absolute identity without resolving symlinks."""
    return os.path.normcase(os.path.normpath(str(path.absolute())))


__all__ = [
    "AlignmentAction",
    "AlignmentBackend",
    "AlignmentCommandPlan",
    "AlignmentOutputFormat",
    "AlignmentPlan",
    "AlignmentRequest",
    "AlignmentResources",
    "AlignmentResult",
    "AlignmentResultStatus",
    "AlignmentTool",
    "build_alignment_plans",
    "build_alignment_requests",
    "plan_alignment_command",
]
