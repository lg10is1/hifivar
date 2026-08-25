"""hap.py wrapper and robust small-variant benchmark summary parser."""

from __future__ import annotations

import csv
import gzip
import json
import math
import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hifivar.benchmark import BenchmarkMetric, BenchmarkRegion, BenchmarkVariantClass, TruthSet
from hifivar.command import CommandRunner, format_command
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar.reference import ReferenceGenome
from hifivar.validation import validate_file, validate_output_file

_VERSION = re.compile(r"(\d+(?:\.\d+)+)")
_CONFIGURED_VERSION = re.compile(r"^\d+(?:\.\d+)+(?:[-+._][A-Za-z0-9][A-Za-z0-9._-]*)?$")


@dataclass(frozen=True, slots=True)
class HappyRequest:
    benchmark_id: str
    sample_id: str
    reference: ReferenceGenome
    query_vcf: Path
    truth_set: TruthSet
    confident_regions: BenchmarkRegion
    output_prefix: Path
    engine: str = "xcmp"
    threads: int = 1
    summary_filter: str = "PASS"
    stratifications: tuple[BenchmarkRegion, ...] = ()
    overwrite: bool = False

    def __post_init__(self) -> None:
        for name in ("benchmark_id", "sample_id", "engine"):
            if not str(getattr(self, name)).strip():
                raise InputValidationError(f"hap.py {name} must be non-empty.")
        if isinstance(self.threads, bool) or self.threads < 1:
            raise InputValidationError("hap.py threads must be positive.")
        if self.reference.build and self.truth_set.reference_build != self.reference.build:
            raise InputValidationError("hap.py truth-set/reference build mismatch.")
        object.__setattr__(self, "query_vcf", Path(self.query_vcf))
        object.__setattr__(self, "output_prefix", Path(self.output_prefix))
        existing = tuple(self.output_prefix.parent.glob(f"{self.output_prefix.name}.*")) if self.output_prefix.parent.exists() else ()
        if existing and not self.overwrite:
            raise OutputValidationError(f"hap.py output prefix already has artifacts: '{self.output_prefix}'.")

    @property
    def summary_csv(self) -> Path:
        return Path(f"{self.output_prefix}.summary.csv")

    @property
    def metrics_json(self) -> Path:
        return Path(f"{self.output_prefix}.metrics.json")

    @property
    def metrics_json_gz(self) -> Path:
        return Path(f"{self.output_prefix}.metrics.json.gz")

    @property
    def metrics_candidates(self) -> tuple[Path, Path]:
        return (self.metrics_json_gz, self.metrics_json)

    @property
    def stratification_tsv(self) -> Path:
        return Path(f"{self.output_prefix}.stratification.tsv")


class HappyResultStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class HappyResult:
    request: HappyRequest
    status: HappyResultStatus
    command: tuple[str, ...]
    version: str | None = None
    version_source: str | None = None
    runtime_seconds: float = 0.0
    metrics: tuple[BenchmarkMetric, ...] = ()
    metrics_artifact: Path | None = None

    def to_dict(self) -> dict[str, object]:
        return {"benchmark_id": self.request.benchmark_id, "sample_id": self.request.sample_id,
                "status": self.status.value, "command": list(self.command),
                "display_command": format_command(self.command), "version": self.version,
                "version_source": self.version_source,
                "runtime_seconds": self.runtime_seconds,
                "outputs": [str(self.request.summary_csv), str(self.metrics_artifact)]
                if self.metrics_artifact is not None else [str(self.request.summary_csv)],
                "metrics": [metric.to_dict() for metric in self.metrics],
                "truth_set": self.request.truth_set.to_dict(),
                "confident_regions": self.request.confident_regions.to_dict()}


class HappyWrapper:
    def __init__(self, *, executable: str = "hap.py", configured_version: str | None = None,
                 runner: CommandRunner | None = None) -> None:
        if configured_version is not None and not _CONFIGURED_VERSION.fullmatch(configured_version.strip()):
            raise ToolVersionError(
                "Configured hap.py version must be an explicit numeric release, for example '0.3.15'."
            )
        self.executable = executable
        self.configured_version = configured_version.strip() if configured_version is not None else None
        self.runner = runner or CommandRunner()

    def plan_command(self, request: HappyRequest) -> tuple[str, ...]:
        command = [self.executable, str(request.truth_set.path.absolute()), str(request.query_vcf.absolute()),
                   "-f", str(request.confident_regions.path.absolute()), "-r", str(request.reference.fasta.absolute()),
                   "-o", str(request.output_prefix.absolute()), "--threads", str(request.threads), f"--engine={request.engine}"]
        if request.stratifications:
            command.extend(("--stratification", str(request.stratification_tsv.absolute())))
        return tuple(command)

    def detect_version(self) -> str:
        return self._resolve_version()[0]

    def _resolve_version(self) -> tuple[str, str]:
        self.runner.require_executable(self.executable)
        result = self.runner.run((self.executable, "--version"))
        match = _VERSION.search("\n".join(value for value in (result.stdout, result.stderr) if value))
        if match is not None:
            return match.group(1), "command"
        if self.configured_version is not None:
            return self.configured_version, "config"
        raise ToolVersionError(
            "hap.py completed '--version' but returned no parseable version. "
            "Set benchmark.small_variants.happy_version to the independently verified installed release."
        )

    def run(self, request: HappyRequest, *, dry_run: bool = False, stderr_path: Path | None = None) -> HappyResult:
        command = self.plan_command(request)
        if dry_run:
            self.runner.run(command, dry_run=True)
            return HappyResult(request, HappyResultStatus.PLANNED, command)
        for path in (request.reference.fasta, request.reference.fai, request.query_vcf,
                     request.truth_set.path, request.confident_regions.path):
            validate_file(path)
        for vcf in (request.query_vcf, request.truth_set.path):
            if str(vcf).endswith(".gz"):
                validate_file(Path(f"{vcf}.tbi"))
        for region in request.stratifications:
            validate_file(region.path)
        request.output_prefix.parent.mkdir(parents=True, exist_ok=True)
        if request.stratifications:
            _write_stratification_file(request.stratification_tsv, request.stratifications, overwrite=request.overwrite)
        version, version_source = self._resolve_version()
        execution = self.runner.run(command, stderr_path=stderr_path)
        validate_output_file(request.summary_csv)
        metrics_artifact = discover_happy_metrics(request.output_prefix)
        parse_happy_metrics_json(metrics_artifact)
        metrics = parse_happy_summary(request.summary_csv, summary_filter=request.summary_filter)
        return HappyResult(
            request, HappyResultStatus.COMPLETED, command, version, version_source,
            execution.duration_seconds, metrics, metrics_artifact,
        )


def discover_happy_metrics(output_prefix: Path) -> Path:
    """Return the sole hap.py metrics artifact, accepting plain or gzip JSON."""
    prefix = Path(output_prefix)
    candidates = (Path(f"{prefix}.metrics.json.gz"), Path(f"{prefix}.metrics.json"))
    present = tuple(path for path in candidates if path.exists())
    if not present:
        raise OutputValidationError(
            f"hap.py produced no metrics artifact; expected exactly one of: "
            f"'{candidates[0]}' or '{candidates[1]}'."
        )
    if len(present) != 1:
        raise OutputValidationError(
            f"hap.py metrics artifact is ambiguous; both compressed and uncompressed files exist: "
            f"'{present[0]}' and '{present[1]}'."
        )
    validate_output_file(present[0])
    return present[0]


def parse_happy_metrics_json(path: Path) -> object:
    """Parse a hap.py metrics JSON artifact regardless of gzip compression."""
    path = Path(path)
    try:
        if path.suffix == ".gz":
            handle = gzip.open(path, "rt", encoding="utf-8")
        else:
            handle = path.open("r", encoding="utf-8")
        with handle:
            return json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise OutputValidationError(f"Unable to parse hap.py metrics JSON '{path}': {error}") from error


def parse_happy_summary(path: Path, *, summary_filter: str = "PASS") -> tuple[BenchmarkMetric, ...]:
    """Parse SNP/INDEL aggregate rows by header names, never row positions."""
    path = Path(path)
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
    except OSError as error:
        raise OutputValidationError(f"Unable to read hap.py summary '{path}': {error}") from error
    metrics: list[BenchmarkMetric] = []
    for kind in ("SNP", "INDEL"):
        candidates = [row for row in rows if (row.get("Type") or "").upper() == kind and (row.get("Filter") or "").upper() == summary_filter.upper()]
        if not candidates:
            raise OutputValidationError(f"hap.py summary lacks {kind}/{summary_filter} aggregate row: '{path}'.")
        row = candidates[0]
        recall = _metric_value(row, "METRIC.Recall")
        precision = _metric_value(row, "METRIC.Precision")
        f1 = _optional_metric_value(row, "METRIC.F1_Score")
        if f1 is None and recall is not None and precision is not None and recall + precision:
            f1 = 2 * recall * precision / (recall + precision)
        for name, value, source in (("recall", recall, "METRIC.Recall"), ("precision", precision, "METRIC.Precision"), ("f1", f1, "METRIC.F1_Score")):
            metrics.append(BenchmarkMetric(name, value, BenchmarkVariantClass.SMALL_VARIANT, kind, source))
    return tuple(metrics)


def _metric_value(row: dict[str, str], name: str) -> float | None:
    if name not in row:
        raise OutputValidationError(f"hap.py summary is missing required column '{name}'.")
    return _coerce_metric(row.get(name))


def _optional_metric_value(row: dict[str, str], name: str) -> float | None:
    return _coerce_metric(row.get(name)) if name in row else None


def _coerce_metric(value: str | None) -> float | None:
    if value is None or value.strip().lower() in {"", ".", "na", "nan"}:
        return None
    try:
        number = float(value)
    except ValueError as error:
        raise OutputValidationError(f"Invalid hap.py metric value: {value!r}.") from error
    return None if math.isnan(number) else number


def _write_stratification_file(path: Path, regions: tuple[BenchmarkRegion, ...], *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise OutputValidationError(f"hap.py stratification file exists: '{path}'.")
    with path.open("w" if overwrite else "x", encoding="utf-8", newline="\n") as handle:
        for region in regions:
            handle.write(f"{region.name}\t{region.path.absolute()}\n")


__all__ = [
    "HappyRequest", "HappyResult", "HappyResultStatus", "HappyWrapper",
    "discover_happy_metrics", "parse_happy_metrics_json", "parse_happy_summary",
]
