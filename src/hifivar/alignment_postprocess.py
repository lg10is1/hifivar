"""Alignment artifacts, explicit indexing plans, and lightweight QC."""

from __future__ import annotations

import os
from dataclasses import dataclass, replace
from enum import Enum
from pathlib import Path

from hifivar import validation
from hifivar.alignment import (
    AlignmentAction,
    AlignmentOutputFormat,
    AlignmentResult,
    AlignmentResultStatus,
    AlignmentTool,
)
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.qc import QCIssue, QCIssueLevel, QCMetric, QCResult, QCStatus
from hifivar.reference import ReferenceGenome
from hifivar.sample import validate_sample_id


_BAI_MAX_CONTIG_LENGTH = 2**29


class AlignmentSource(str, Enum):
    """How an alignment entered the Phase 2 post-processing layer."""

    GENERATED = "generated"
    EXISTING = "existing"


class AlignmentSortOrder(str, Enum):
    """Only sort-order claims that Phase 2 can establish without parsing."""

    COORDINATE = "coordinate"
    UNKNOWN = "unknown"


class AlignmentIndexFormat(str, Enum):
    """Supported BAM/CRAM random-access index formats."""

    BAI = "bai"
    CSI = "csi"
    CRAI = "crai"

    @property
    def suffix(self) -> str:
        """Return the appended index suffix."""
        return f".{self.value}"


@dataclass(frozen=True, slots=True)
class AlignmentArtifact:
    """Validated metadata boundary for an existing or generated alignment."""

    sample_id: str
    path: Path
    output_format: AlignmentOutputFormat
    reference: ReferenceGenome
    source: AlignmentSource
    sort_order: AlignmentSortOrder
    index_path: Path | None = None
    tool: AlignmentTool | None = None
    tool_version: str | None = None

    def __post_init__(self) -> None:
        """Normalize paths and reject unsupported metadata combinations."""
        validate_sample_id(self.sample_id)
        if not isinstance(self.output_format, AlignmentOutputFormat):
            raise InputValidationError(
                "Alignment artifact output_format must be AlignmentOutputFormat."
            )
        if not isinstance(self.reference, ReferenceGenome):
            raise InputValidationError(
                "Alignment artifact reference must be a ReferenceGenome."
            )
        if not isinstance(self.source, AlignmentSource):
            raise InputValidationError(
                "Alignment artifact source must be an AlignmentSource."
            )
        if not isinstance(self.sort_order, AlignmentSortOrder):
            raise InputValidationError(
                "Alignment artifact sort_order must be AlignmentSortOrder."
            )
        path = _coerce_path(self.path, "Alignment artifact")
        if path.suffix.lower() != self.output_format.suffix:
            raise InputValidationError(
                f"Alignment artifact '{path}' does not match format "
                f"{self.output_format.value.upper()}."
            )
        object.__setattr__(self, "path", path)

        if self.index_path is not None:
            index_path = _coerce_path(self.index_path, "Alignment index")
            _validate_index_suffix(self.output_format, index_path)
            object.__setattr__(self, "index_path", index_path)
        if self.tool is not None and not isinstance(self.tool, AlignmentTool):
            raise InputValidationError(
                "Alignment artifact tool must be AlignmentTool or None."
            )
        if self.source is AlignmentSource.GENERATED and self.tool is None:
            raise InputValidationError(
                "Generated alignment artifacts require a tool."
            )
        if self.tool_version is not None and (
            not isinstance(self.tool_version, str) or not self.tool_version.strip()
        ):
            raise InputValidationError(
                "Alignment artifact tool_version must be non-empty or None."
            )

    @classmethod
    def from_result(cls, result: AlignmentResult) -> AlignmentArtifact:
        """Create artifact metadata from a completed or reused result."""
        if not isinstance(result, AlignmentResult):
            raise InputValidationError(
                "AlignmentArtifact.from_result requires AlignmentResult."
            )
        if result.status is AlignmentResultStatus.PLANNED:
            raise OutputValidationError(
                "A dry-run alignment result cannot become an artifact."
            )
        if result.plan.action is AlignmentAction.REUSE:
            source = AlignmentSource.EXISTING
            sort_order = AlignmentSortOrder.UNKNOWN
            tool = None
        else:
            source = AlignmentSource.GENERATED
            tool = result.plan.request.tool if result.plan.request is not None else None
            sort_order = (
                AlignmentSortOrder.COORDINATE
                if tool is AlignmentTool.PBMM2
                else AlignmentSortOrder.UNKNOWN
            )
        return cls(
            sample_id=result.plan.sample_id,
            path=result.alignment_path,
            output_format=result.plan.output_format,
            reference=result.plan.reference,
            source=source,
            sort_order=sort_order,
            index_path=find_alignment_index(result.alignment_path),
            tool=tool,
            tool_version=result.tool_version,
        )

    def with_index(self, index_path: str | Path) -> AlignmentArtifact:
        """Return a new artifact retaining all provenance plus an index path."""
        return replace(self, index_path=Path(index_path).expanduser())

    def to_dict(self) -> dict[str, object]:
        """Return standard alignment artifact provenance."""
        return {
            "sample_id": self.sample_id,
            "path": str(self.path),
            "format": self.output_format.value,
            "reference_build": self.reference.build,
            "source": self.source.value,
            "sort_order": self.sort_order.value,
            "index_path": str(self.index_path) if self.index_path else None,
            "tool": self.tool.value if self.tool else None,
            "tool_version": self.tool_version,
        }


@dataclass(frozen=True, slots=True)
class AlignmentIndexRequest:
    """Explicit plan to index one known coordinate-sorted alignment."""

    artifact: AlignmentArtifact
    index_format: AlignmentIndexFormat
    output_path: Path
    threads: int = 1
    overwrite: bool = False

    def __post_init__(self) -> None:
        """Require compatible format, sort order, resources, and overwrite."""
        if not isinstance(self.artifact, AlignmentArtifact):
            raise InputValidationError(
                "Alignment index request artifact must be AlignmentArtifact."
            )
        if self.artifact.sort_order is not AlignmentSortOrder.COORDINATE:
            raise InputValidationError(
                "Alignment indexing requires explicitly coordinate-sorted input; "
                "UNKNOWN sort order is never assumed safe."
            )
        if not isinstance(self.index_format, AlignmentIndexFormat):
            raise InputValidationError(
                "Alignment index format must be AlignmentIndexFormat."
            )
        _validate_index_format(self.artifact.output_format, self.index_format)
        if (
            not isinstance(self.threads, int)
            or isinstance(self.threads, bool)
            or self.threads <= 0
        ):
            raise InputValidationError(
                "Alignment index threads must be a positive integer."
            )
        if not isinstance(self.overwrite, bool):
            raise InputValidationError(
                "Alignment index overwrite must be a boolean."
            )
        output_path = _coerce_path(self.output_path, "Alignment index output")
        if output_path.suffix.lower() != self.index_format.suffix:
            raise InputValidationError(
                f"Index output '{output_path}' does not match "
                f"{self.index_format.value.upper()}."
            )
        if output_path.exists() and output_path.is_dir():
            raise OutputValidationError(
                f"Alignment index output is a directory: '{output_path}'."
            )
        if output_path.exists() and not self.overwrite:
            raise OutputValidationError(
                f"Alignment index output already exists: '{output_path}'."
            )
        object.__setattr__(self, "output_path", output_path)

    @classmethod
    def create(
        cls,
        artifact: AlignmentArtifact,
        *,
        index_format: AlignmentIndexFormat | None = None,
        threads: int = 1,
        overwrite: bool = False,
    ) -> AlignmentIndexRequest:
        """Select a safe default index format and conventional output path."""
        if not isinstance(artifact, AlignmentArtifact):
            raise InputValidationError(
                "Alignment index planning requires AlignmentArtifact."
            )
        selected = index_format or choose_index_format(artifact)
        return cls(
            artifact=artifact,
            index_format=selected,
            output_path=index_path_for(artifact.path, selected),
            threads=threads,
            overwrite=overwrite,
        )

    def to_dict(self) -> dict[str, object]:
        """Return standard indexing plan metadata."""
        return {
            "artifact": self.artifact.to_dict(),
            "index_format": self.index_format.value,
            "output_path": str(self.output_path),
            "threads": self.threads,
            "overwrite": self.overwrite,
        }


def choose_index_format(artifact: AlignmentArtifact) -> AlignmentIndexFormat:
    """Choose CRAI for CRAM and BAI/CSI for BAM using the BAI contig limit."""
    if artifact.output_format is AlignmentOutputFormat.CRAM:
        return AlignmentIndexFormat.CRAI
    if any(
        contig.length > _BAI_MAX_CONTIG_LENGTH
        for contig in artifact.reference.contigs
    ):
        return AlignmentIndexFormat.CSI
    return AlignmentIndexFormat.BAI


def index_path_for(path: str | Path, index_format: AlignmentIndexFormat) -> Path:
    """Return an explicit adjacent index path without checking the filesystem."""
    alignment_path = _coerce_path(path, "Alignment")
    if not isinstance(index_format, AlignmentIndexFormat):
        raise InputValidationError(
            "index_format must be an AlignmentIndexFormat value."
        )
    return Path(f"{alignment_path}{index_format.suffix}")


def find_alignment_index(path: str | Path) -> Path | None:
    """Return the first readable conventional BAI/CSI/CRAI path."""
    alignment_path = _coerce_path(path, "Alignment")
    if alignment_path.suffix.lower() == ".bam":
        suffixes = (".bai", ".csi")
    elif alignment_path.suffix.lower() == ".cram":
        suffixes = (".crai",)
    else:
        raise InputValidationError(
            f"Unsupported alignment suffix for '{alignment_path}'."
        )
    candidates: list[Path] = []
    for suffix in suffixes:
        candidates.extend(
            (Path(f"{alignment_path}{suffix}"), alignment_path.with_suffix(suffix))
        )
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.R_OK):
            return candidate
    return None


def validate_alignment_artifact(
    artifact: AlignmentArtifact,
    *,
    require_index: bool = False,
) -> AlignmentArtifact:
    """Validate lightweight output existence and optional index presence."""
    if not isinstance(artifact, AlignmentArtifact):
        raise InputValidationError(
            "Alignment artifact validation requires AlignmentArtifact."
        )
    validation.validate_output_file(artifact.path)
    index_path = artifact.index_path or find_alignment_index(artifact.path)
    if require_index and index_path is None:
        raise OutputValidationError(
            f"Alignment index missing for '{artifact.path}'."
        )
    if index_path is not None:
        validation.validate_output_file(index_path)
    return artifact if index_path == artifact.index_path else replace(
        artifact,
        index_path=index_path,
    )


def run_alignment_qc(artifact: AlignmentArtifact) -> QCResult:
    """Run path/index/sort metadata QC without reading BAM/CRAM records."""
    validated = validate_alignment_artifact(artifact, require_index=False)
    index_path = validated.index_path
    metrics = (
        QCMetric("alignment_format", validated.output_format.value),
        QCMetric("alignment_path", str(validated.path.absolute())),
        QCMetric("file_size_bytes", validated.path.stat().st_size, unit="bytes"),
        QCMetric("alignment_source", validated.source.value),
        QCMetric("sort_order", validated.sort_order.value),
        QCMetric("index_present", index_path is not None),
        QCMetric("index_path", str(index_path.absolute()) if index_path else None),
    )
    issues: list[QCIssue] = []
    if validated.sort_order is AlignmentSortOrder.UNKNOWN:
        issues.append(
            QCIssue(
                code="ALIGNMENT_SORT_ORDER_UNKNOWN",
                level=QCIssueLevel.WARNING,
                message=(
                    "Alignment sort order is unknown; HiFiVar does not infer it "
                    "from the filename or silently sort the file."
                ),
            )
        )
    if index_path is None:
        issues.append(
            QCIssue(
                code="ALIGNMENT_INDEX_MISSING",
                level=QCIssueLevel.WARNING,
                message="Alignment index is missing.",
            )
        )
    status = QCStatus.WARN if issues else QCStatus.PASS
    return QCResult(
        status=status,
        metrics=metrics,
        issues=tuple(issues),
        sample_id=validated.sample_id,
    )


def _validate_index_format(
    output_format: AlignmentOutputFormat,
    index_format: AlignmentIndexFormat,
) -> None:
    """Require BAI/CSI for BAM and CRAI for CRAM."""
    if output_format is AlignmentOutputFormat.BAM:
        if index_format not in {AlignmentIndexFormat.BAI, AlignmentIndexFormat.CSI}:
            raise InputValidationError("BAM indexing requires BAI or CSI.")
    elif index_format is not AlignmentIndexFormat.CRAI:
        raise InputValidationError("CRAM indexing requires CRAI.")


def _validate_index_suffix(
    output_format: AlignmentOutputFormat,
    index_path: Path,
) -> None:
    """Validate an attached index suffix for its alignment container."""
    allowed = (
        {".bai", ".csi"}
        if output_format is AlignmentOutputFormat.BAM
        else {".crai"}
    )
    if index_path.suffix.lower() not in allowed:
        expected = ", ".join(sorted(allowed))
        raise InputValidationError(
            f"Alignment index '{index_path}' has an incompatible suffix; "
            f"expected one of: {expected}."
        )


def _coerce_path(value: str | Path, label: str) -> Path:
    """Normalize one portable path without resolving it."""
    if not isinstance(value, (str, Path)):
        raise InputValidationError(f"{label} path must be string or Path.")
    if isinstance(value, str) and not value.strip():
        raise InputValidationError(f"{label} path must not be empty.")
    return Path(value).expanduser()


__all__ = [
    "AlignmentArtifact",
    "AlignmentIndexFormat",
    "AlignmentIndexRequest",
    "AlignmentSortOrder",
    "AlignmentSource",
    "choose_index_format",
    "find_alignment_index",
    "index_path_for",
    "run_alignment_qc",
    "validate_alignment_artifact",
]
