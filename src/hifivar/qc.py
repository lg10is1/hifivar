"""Unified QC result models and lightweight primary-input inspection."""

from __future__ import annotations

import math
import os
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Final

from hifivar import __version__, validation
from hifivar.context import AnalysisContext
from hifivar.exceptions import InputValidationError
from hifivar.logging_utils import get_logger
from hifivar.sample import InputDataset, InputType, validate_sample_id
from hifivar.serialization import (
    standardize_data,
    utc_now_iso8601,
    write_json_atomic,
    write_yaml_atomic,
)


MetricValue = str | int | float | bool | None

_METRIC_NAME_PATTERN: Final[re.Pattern[str]] = re.compile(r"[a-z][a-z0-9_]*\Z")
_ISSUE_CODE_PATTERN: Final[re.Pattern[str]] = re.compile(r"[A-Z][A-Z0-9_]*\Z")
_MODULE_PATTERN: Final[re.Pattern[str]] = re.compile(r"[a-z][a-z0-9_]*\Z")
_LOGGER = get_logger(__name__)


class QCStatus(str, Enum):
    """Deterministic overall state for one QC result or report."""

    PASS = "pass"
    WARN = "warn"
    FAIL = "fail"
    NOT_CHECKED = "not_checked"


class QCIssueLevel(str, Enum):
    """Minimal actionable level for an individual QC finding."""

    WARNING = "warning"
    ERROR = "error"


@dataclass(frozen=True, slots=True)
class QCMetric:
    """One machine-readable QC measurement with optional human context."""

    name: str
    value: MetricValue
    unit: str | None = None
    description: str | None = None

    def __post_init__(self) -> None:
        """Require a stable name and a portable scalar value."""
        if (
            not isinstance(self.name, str)
            or _METRIC_NAME_PATTERN.fullmatch(self.name) is None
        ):
            raise InputValidationError(
                "QC metric name must use lowercase snake_case."
            )
        if isinstance(self.value, Enum) or not (
            self.value is None
            or isinstance(self.value, (str, int, float, bool))
        ):
            raise InputValidationError(
                "QC metric value must be a string, integer, float, boolean, or null."
            )
        if isinstance(self.value, float) and not math.isfinite(self.value):
            raise InputValidationError("QC metric float value must be finite.")
        _validate_optional_text(self.unit, "QC metric unit")
        _validate_optional_text(self.description, "QC metric description")

    def to_dict(self) -> dict[str, MetricValue]:
        """Return a deterministic standard-type metric mapping."""
        return {
            "name": self.name,
            "value": self.value,
            "unit": self.unit,
            "description": self.description,
        }


@dataclass(frozen=True, slots=True)
class QCIssue:
    """One stable machine code plus a human-readable QC finding."""

    code: str
    level: QCIssueLevel
    message: str

    def __post_init__(self) -> None:
        """Validate the stable code, issue level, and explanatory message."""
        if (
            not isinstance(self.code, str)
            or _ISSUE_CODE_PATTERN.fullmatch(self.code) is None
        ):
            raise InputValidationError(
                "QC issue code must use uppercase snake case."
            )
        if not isinstance(self.level, QCIssueLevel):
            raise InputValidationError("QC issue level must be a QCIssueLevel value.")
        if not isinstance(self.message, str) or not self.message.strip():
            raise InputValidationError("QC issue message must be non-empty text.")

    def to_dict(self) -> dict[str, str]:
        """Return a deterministic standard-type issue mapping."""
        return {
            "code": self.code,
            "level": self.level.value,
            "message": self.message,
        }


@dataclass(frozen=True, slots=True)
class QCResult:
    """One module/sample QC outcome with ordered metrics and findings."""

    status: QCStatus
    metrics: tuple[QCMetric, ...] = ()
    issues: tuple[QCIssue, ...] = ()
    module: str = "input_qc"
    sample_id: str | None = None

    def __post_init__(self) -> None:
        """Detach ordered values and prevent status/issue contradictions."""
        if not isinstance(self.status, QCStatus):
            raise InputValidationError("QC result status must be a QCStatus value.")
        if (
            not isinstance(self.module, str)
            or _MODULE_PATTERN.fullmatch(self.module) is None
        ):
            raise InputValidationError(
                "QC result module must use lowercase snake_case."
            )
        if self.sample_id is not None:
            validate_sample_id(self.sample_id)

        metrics = tuple(self.metrics)
        issues = tuple(self.issues)
        if any(not isinstance(metric, QCMetric) for metric in metrics):
            raise InputValidationError(
                "QC result metrics must contain only QCMetric objects."
            )
        if any(not isinstance(issue, QCIssue) for issue in issues):
            raise InputValidationError(
                "QC result issues must contain only QCIssue objects."
            )
        _validate_unique_metric_names(metrics, "QC result")
        _validate_status_issue_consistency(self.status, issues)
        object.__setattr__(self, "metrics", metrics)
        object.__setattr__(self, "issues", issues)

    def get_metric(self, name: str) -> QCMetric:
        """Return one metric by exact machine name or raise ``KeyError``."""
        for metric in self.metrics:
            if metric.name == name:
                return metric
        raise KeyError(name)

    def to_dict(self) -> dict[str, object]:
        """Return only standard JSON/YAML containers and scalar values."""
        return {
            "module": self.module,
            "sample_id": self.sample_id,
            "status": self.status.value,
            "metrics": [metric.to_dict() for metric in self.metrics],
            "issues": [issue.to_dict() for issue in self.issues],
        }


@dataclass(frozen=True, slots=True)
class RunQCReport:
    """Immutable run-level aggregation of ordered sample QC results."""

    sample_results: tuple[QCResult, ...]
    metrics: tuple[QCMetric, ...] = ()
    created_at: str = field(default_factory=utc_now_iso8601)
    hifivar_version: str = __version__
    overall_status: QCStatus = field(init=False)

    def __post_init__(self) -> None:
        """Detach nested values and derive the only allowed overall status."""
        sample_results = tuple(self.sample_results)
        metrics = tuple(self.metrics)
        if any(not isinstance(result, QCResult) for result in sample_results):
            raise InputValidationError(
                "Run QC sample_results must contain only QCResult objects."
            )
        if any(not isinstance(metric, QCMetric) for metric in metrics):
            raise InputValidationError(
                "Run QC metrics must contain only QCMetric objects."
            )
        _validate_unique_metric_names(metrics, "Run QC report")
        if not isinstance(self.created_at, str) or not self.created_at.strip():
            raise InputValidationError("Run QC created_at must be non-empty text.")
        if not isinstance(self.hifivar_version, str) or not self.hifivar_version:
            raise InputValidationError(
                "Run QC hifivar_version must be non-empty text."
            )
        object.__setattr__(self, "sample_results", sample_results)
        object.__setattr__(self, "metrics", metrics)
        object.__setattr__(
            self,
            "overall_status",
            aggregate_qc_status(result.status for result in sample_results),
        )

    @property
    def status_counts(self) -> dict[str, int]:
        """Return deterministic counts for every status, including zeros."""
        return {
            status.value: sum(
                result.status is status for result in self.sample_results
            )
            for status in QCStatus
        }

    def get_metric(self, name: str) -> QCMetric:
        """Return one run metric by exact machine name or raise ``KeyError``."""
        for metric in self.metrics:
            if metric.name == name:
                return metric
        raise KeyError(name)

    def to_dict(self) -> dict[str, object]:
        """Return an independent report payload with standard data types."""
        payload = {
            "hifivar_version": self.hifivar_version,
            "created_at": self.created_at,
            "overall_status": self.overall_status.value,
            "status_counts": self.status_counts,
            "metrics": [metric.to_dict() for metric in self.metrics],
            "sample_results": [result.to_dict() for result in self.sample_results],
        }
        standardized = standardize_data(payload, context="QC report value")
        if not isinstance(standardized, dict):  # pragma: no cover - invariant guard
            raise InputValidationError("QC report serialization produced invalid data.")
        return standardized

    def write_json(self, path: str | Path, *, overwrite: bool = False) -> Path:
        """Atomically write UTF-8 JSON, refusing replacement by default."""
        return write_json_atomic(
            self.to_dict(),
            path,
            overwrite=overwrite,
            artifact_name="QC report",
        )

    def write_yaml(self, path: str | Path, *, overwrite: bool = False) -> Path:
        """Atomically write UTF-8 YAML, refusing replacement by default."""
        return write_yaml_atomic(
            self.to_dict(),
            path,
            overwrite=overwrite,
            artifact_name="QC report",
        )


def aggregate_qc_status(statuses: Iterable[QCStatus]) -> QCStatus:
    """Aggregate statuses with FAIL > WARN > PASS > NOT_CHECKED precedence."""
    supplied = tuple(statuses)
    if any(not isinstance(status, QCStatus) for status in supplied):
        raise InputValidationError(
            "QC status aggregation requires only QCStatus values."
        )
    for status in (QCStatus.FAIL, QCStatus.WARN, QCStatus.PASS):
        if status in supplied:
            return status
    return QCStatus.NOT_CHECKED


def run_input_dataset_qc(
    dataset: InputDataset,
    *,
    sample_id: str | None = None,
) -> QCResult:
    """Run low-cost filesystem and existing prefix validation for one dataset."""
    if not isinstance(dataset, InputDataset):
        raise InputValidationError(
            "run_input_dataset_qc requires an InputDataset instance."
        )
    if sample_id is not None:
        validate_sample_id(sample_id)
    _LOGGER.debug(
        "Running lightweight input QC: sample_id=%s, input_type=%s, file_count=%d",
        sample_id or "-",
        dataset.input_type.value,
        len(dataset.files),
    )

    file_sizes: list[int] = []
    file_metrics: list[QCMetric] = []
    for index, path in enumerate(dataset.files, start=1):
        if dataset.input_type is InputType.FASTQ:
            validated_path = validation.validate_fastq(path)
        else:
            validated_path = validation.validate_alignment_file(
                path,
                require_index=False,
            )
        try:
            size_bytes = validated_path.stat().st_size
        except OSError as error:
            raise InputValidationError(
                f"Unable to inspect input file '{validated_path}' during QC: {error}"
            ) from error
        file_sizes.append(size_bytes)
        file_metrics.extend(
            (
                QCMetric(f"file_{index}_path", str(validated_path.absolute())),
                QCMetric(f"file_{index}_size_bytes", size_bytes, unit="bytes"),
            )
        )

    metrics = [
        QCMetric("input_type", dataset.input_type.value),
        QCMetric("file_count", len(dataset.files)),
        QCMetric(
            "total_file_size_bytes",
            sum(file_sizes),
            unit="bytes",
        ),
    ]
    issues: list[QCIssue] = []
    status = QCStatus.PASS
    if dataset.input_type is InputType.FASTQ:
        metrics.append(QCMetric("compression", _fastq_compression(dataset.files)))
    else:
        index_present = _find_alignment_index(dataset.files[0]) is not None
        metrics.append(QCMetric("index_present", index_present))
        if not index_present:
            status = QCStatus.WARN
            issues.append(
                QCIssue(
                    code="ALIGNMENT_INDEX_MISSING",
                    level=QCIssueLevel.WARNING,
                    message=(
                        f"{dataset.input_type.value.upper()} input has no readable "
                        "conventional index; lightweight QC can continue."
                    ),
                )
            )
    metrics.extend(file_metrics)
    return QCResult(
        status=status,
        metrics=tuple(metrics),
        issues=tuple(issues),
        sample_id=sample_id,
    )


def run_input_qc(context: AnalysisContext) -> RunQCReport:
    """Run lightweight input QC for every sample in an AnalysisContext."""
    if not isinstance(context, AnalysisContext):
        raise InputValidationError("run_input_qc requires an AnalysisContext.")
    sample_results = tuple(
        run_input_dataset_qc(
            record.sample.input,
            sample_id=record.sample.sample_id,
        )
        for record in context.samples
    )
    report = RunQCReport(
        sample_results=sample_results,
        metrics=(
            QCMetric("reference_build", context.reference.build),
            QCMetric("reference_contig_count", len(context.reference.contigs)),
            QCMetric(
                "reference_checksum_available",
                context.reference.sha256 is not None,
            ),
        ),
    )
    counts = report.status_counts
    _LOGGER.info(
        "QC summary: %d PASS, %d WARN, %d FAIL, %d NOT_CHECKED",
        counts[QCStatus.PASS.value],
        counts[QCStatus.WARN.value],
        counts[QCStatus.FAIL.value],
        counts[QCStatus.NOT_CHECKED.value],
    )
    return report


def _validate_optional_text(value: str | None, label: str) -> None:
    if value is not None and (
        not isinstance(value, str) or not value.strip()
    ):
        raise InputValidationError(f"{label} must be non-empty text or None.")


def _validate_unique_metric_names(
    metrics: tuple[QCMetric, ...],
    label: str,
) -> None:
    names = [metric.name for metric in metrics]
    if len(set(names)) != len(names):
        raise InputValidationError(f"{label} contains duplicate metric names.")


def _validate_status_issue_consistency(
    status: QCStatus,
    issues: tuple[QCIssue, ...],
) -> None:
    if status in {QCStatus.PASS, QCStatus.NOT_CHECKED} and issues:
        raise InputValidationError(
            f"QC status '{status.value}' cannot contain warning/error issues."
        )
    if any(issue.level is QCIssueLevel.ERROR for issue in issues):
        if status is not QCStatus.FAIL:
            raise InputValidationError("Error-level QC issues require FAIL status.")
    if any(issue.level is QCIssueLevel.WARNING for issue in issues):
        if status not in {QCStatus.WARN, QCStatus.FAIL}:
            raise InputValidationError(
                "Warning-level QC issues require WARN or FAIL status."
            )
    if status is QCStatus.WARN and not issues:
        raise InputValidationError("WARN QC status requires at least one issue.")


def _fastq_compression(files: tuple[Path, ...]) -> str:
    compressed = [path.name.lower().endswith(".gz") for path in files]
    if all(compressed):
        return "gzip"
    if any(compressed):
        return "mixed"
    return "uncompressed"


def _find_alignment_index(path: Path) -> Path | None:
    index_suffix = ".bai" if path.suffix.lower() == ".bam" else ".crai"
    candidates = (Path(f"{path}{index_suffix}"), path.with_suffix(index_suffix))
    for candidate in candidates:
        if candidate.exists() and candidate.is_file() and os.access(candidate, os.R_OK):
            return candidate
    return None


__all__ = [
    "QCIssue",
    "QCIssueLevel",
    "QCMetric",
    "QCResult",
    "QCStatus",
    "RunQCReport",
    "aggregate_qc_status",
    "run_input_dataset_qc",
    "run_input_qc",
]
