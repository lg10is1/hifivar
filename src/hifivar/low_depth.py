"""Phase 15 deterministic low-depth experiment contracts.

The module plans and materializes read/template-level subsamples from an
existing coordinate-sorted BAM.  It does not invent an acceptable coverage
threshold: support is evaluated only against an explicit, preregistered metric
policy and replicated observations supplied by Phase 13 benchmark runs.
"""

from __future__ import annotations

import gzip
import hashlib
import re
from contextlib import ExitStack
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

from hifivar.command import CommandRunner
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar.serialization import utc_now_iso8601, write_json_atomic, write_yaml_atomic


_VERSION = re.compile(r"\bsamtools\s+([0-9]+(?:\.[0-9]+){1,2}(?:[-+._A-Za-z0-9]*)?)", re.I)
_SAFE_IDENTIFIER = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")


class LowDepthRunStatus(str, Enum):
    PLANNED = "PLANNED"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    NOT_RUN = "NOT_RUN"


@dataclass(frozen=True, slots=True)
class RawFastqSubset:
    """One deterministic molecule-level FASTQ subset."""

    target_coverage: float
    baseline_coverage: float
    sampling_fraction: float
    replicate: int
    seed: int
    output_fastq: Path

    def to_dict(self) -> dict[str, object]:
        return {
            "target_coverage": self.target_coverage,
            "baseline_coverage": self.baseline_coverage,
            "sampling_fraction": self.sampling_fraction,
            "replicate": self.replicate,
            "seed": self.seed,
            "output_fastq": str(self.output_fastq),
        }


@dataclass(frozen=True, slots=True)
class RawFastqSubsamplingPlan:
    """Shared raw-read subsets used by every parameter profile."""

    study_id: str
    sample_id: str
    source_fastqs: tuple[Path, ...]
    baseline_coverage: float
    subsets: tuple[RawFastqSubset, ...]
    created_at: str = field(default_factory=utc_now_iso8601)

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": "1.0",
            "study_id": self.study_id,
            "sample_id": self.sample_id,
            "source_fastqs": [str(path) for path in self.source_fastqs],
            "source_contract": "pacbio_hifi_fastq",
            "baseline_coverage": self.baseline_coverage,
            "subsets": [subset.to_dict() for subset in self.subsets],
            "created_at": self.created_at,
            "scientific_semantics": {
                "sampling_unit_is_pacbio_molecule": True,
                "subsamples_are_nested_across_depths": True,
                "subsamples_are_paired_across_profiles": True,
                "target_coverage_requires_post_alignment_measurement": True,
                "raw_sources_are_not_modified": True,
            },
        }

    def write(self, path: Path, *, overwrite: bool = False) -> Path:
        write_json_atomic(
            self.to_dict(), Path(path), overwrite=overwrite,
            artifact_name="Raw FASTQ subsampling plan",
        )
        return Path(path)


@dataclass(frozen=True, slots=True)
class RawFastqSubsetResult:
    plan: RawFastqSubsamplingPlan
    status: LowDepthRunStatus
    read_counts: Mapping[str, int]
    base_counts: Mapping[str, int]

    def to_dict(self) -> dict[str, object]:
        return {
            "study_id": self.plan.study_id,
            "status": self.status.value,
            "read_counts": dict(self.read_counts),
            "base_counts": dict(self.base_counts),
            "outputs": [str(item.output_fastq) for item in self.plan.subsets],
        }


@dataclass(frozen=True, slots=True)
class LowDepthBenchmarkRun:
    """Deterministic artifact contract for one raw-read benchmark run."""

    run_id: str
    sample_id: str
    target_coverage: float
    replicate: int
    input_fastq: Path
    alignment_bam: Path
    alignment_index: Path
    coverage_json: Path
    output_vcf: Path
    output_gvcf: Path
    happy_output_prefix: Path
    observation_json: Path

    @property
    def output_vcf_index(self) -> Path:
        return Path(f"{self.output_vcf}.tbi")

    @property
    def output_gvcf_index(self) -> Path:
        return Path(f"{self.output_gvcf}.tbi")

    @property
    def happy_summary(self) -> Path:
        return Path(f"{self.happy_output_prefix}.summary.csv")

    def happy_metrics(self, compression: str = "gzip") -> Path:
        suffix = ".metrics.json.gz" if compression == "gzip" else ".metrics.json"
        return Path(f"{self.happy_output_prefix}{suffix}")

    def to_dict(self, *, metrics_compression: str = "gzip") -> dict[str, object]:
        return {
            "run_id": self.run_id,
            "sample_id": self.sample_id,
            "target_coverage": self.target_coverage,
            "replicate": self.replicate,
            "input_fastq": str(self.input_fastq),
            "alignment_bam": str(self.alignment_bam),
            "alignment_index": str(self.alignment_index),
            "coverage_json": str(self.coverage_json),
            "output_vcf": str(self.output_vcf),
            "output_vcf_index": str(self.output_vcf_index),
            "output_gvcf": str(self.output_gvcf),
            "output_gvcf_index": str(self.output_gvcf_index),
            "happy_output_prefix": str(self.happy_output_prefix),
            "happy_summary": str(self.happy_summary),
            "happy_metrics": str(self.happy_metrics(metrics_compression)),
            "observation_json": str(self.observation_json),
        }


@dataclass(frozen=True, slots=True)
class LowDepthBenchmarkPlan:
    """Raw FASTQ to hap.py orchestration contract without scientific claims."""

    study_id: str
    sample_id: str
    runs: tuple[LowDepthBenchmarkRun, ...]
    metrics_compression: str = "gzip"

    def __post_init__(self) -> None:
        if self.metrics_compression not in {"gzip", "plain"}:
            raise InputValidationError(
                "Low-depth hap.py metrics_compression must be 'gzip' or 'plain'."
            )
        if not self.runs or len({run.run_id for run in self.runs}) != len(self.runs):
            raise InputValidationError(
                "Low-depth benchmark runs must be non-empty with unique run IDs."
            )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": "1.0",
            "study_id": self.study_id,
            "sample_id": self.sample_id,
            "metrics_compression": self.metrics_compression,
            "runs": [
                run.to_dict(metrics_compression=self.metrics_compression)
                for run in self.runs
            ],
            "scientific_semantics": {
                "each_subset_is_realigned_independently": True,
                "observed_coverage_is_measured_after_alignment": True,
                "hap_py_truth_inputs_are_explicit": True,
                "minimum_supported_depth_is_not_established_by_planning": True,
                "optimization_requires_independent_validation": True,
            },
        }


def build_low_depth_benchmark_plan(
    raw_plan: RawFastqSubsamplingPlan,
    output_root: Path,
    *,
    metrics_compression: str = "gzip",
) -> LowDepthBenchmarkPlan:
    """Map every raw subset to isolated alignment, calling and benchmark paths."""
    if not isinstance(raw_plan, RawFastqSubsamplingPlan):
        raise InputValidationError(
            "Low-depth end-to-end planning requires RawFastqSubsamplingPlan."
        )
    root = Path(output_root)
    runs = tuple(
        LowDepthBenchmarkRun(
            run_id=(
                f"{raw_plan.study_id}.{_coverage_label(item.target_coverage)}."
                f"r{item.replicate}"
            ),
            sample_id=raw_plan.sample_id,
            target_coverage=item.target_coverage,
            replicate=item.replicate,
            input_fastq=item.output_fastq,
            alignment_bam=(root / "alignment" / _coverage_label(item.target_coverage) / f"replicate-{item.replicate}.bam"),
            alignment_index=(root / "alignment" / _coverage_label(item.target_coverage) / f"replicate-{item.replicate}.bam.bai"),
            coverage_json=(root / "coverage" / _coverage_label(item.target_coverage) / f"replicate-{item.replicate}.coverage.json"),
            output_vcf=(root / "small" / _coverage_label(item.target_coverage) / f"replicate-{item.replicate}.small.vcf.gz"),
            output_gvcf=(root / "small" / _coverage_label(item.target_coverage) / f"replicate-{item.replicate}.g.vcf.gz"),
            happy_output_prefix=(root / "benchmark" / _coverage_label(item.target_coverage) / f"replicate-{item.replicate}.happy"),
            observation_json=(root / "observations" / _coverage_label(item.target_coverage) / f"replicate-{item.replicate}.json"),
        )
        for item in raw_plan.subsets
    )
    return LowDepthBenchmarkPlan(raw_plan.study_id, raw_plan.sample_id, runs, metrics_compression)


@dataclass(frozen=True, slots=True)
class AlignmentCoverageResult:
    """Observed reference-weighted mean depth from ``samtools coverage``."""

    sample_id: str
    alignment_bam: Path
    observed_coverage: float
    covered_bases: int
    reference_bases: int
    command: tuple[str, ...]
    tool_version: str | None

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": "1.0", "sample_id": self.sample_id,
            "alignment_bam": str(self.alignment_bam),
            "observed_coverage": self.observed_coverage,
            "covered_bases": self.covered_bases,
            "reference_bases": self.reference_bases,
            "command": list(self.command), "tool": "samtools",
            "tool_version": self.tool_version,
        }


class SamtoolsCoverageWrapper:
    """Measure alignment depth through CommandRunner without loading BAM data."""

    def __init__(self, *, runner: CommandRunner | None = None,
                 executable: str | Path = "samtools") -> None:
        self.runner = runner or CommandRunner()
        self.executable = str(executable)

    def plan_command(self, alignment_bam: Path) -> tuple[str, ...]:
        return (self.executable, "coverage", str(Path(alignment_bam)))

    def run(self, *, sample_id: str, alignment_bam: Path,
            alignment_index: Path, dry_run: bool = False,
            stderr_path: Path | None = None) -> AlignmentCoverageResult:
        bam = Path(alignment_bam)
        index = Path(alignment_index)
        command = self.plan_command(bam)
        if dry_run:
            self.runner.run(command, dry_run=True, stderr_path=stderr_path)
            return AlignmentCoverageResult(sample_id, bam, 0.0, 0, 0, command, None)
        for path in (bam, index):
            if not path.is_file() or path.stat().st_size == 0:
                raise InputValidationError(f"Coverage input is missing or empty: '{path}'.")
        version = SamtoolsDownsampleWrapper(
            runner=self.runner, executable=self.executable,
        ).detect_version()
        result = self.runner.run(command, stderr_path=stderr_path)
        coverage_text = result.stdout
        if coverage_text is None:
            if result.stdout_path is None:
                raise OutputValidationError(
                    "samtools coverage returned no captured output or persisted stdout path."
                )
            try:
                coverage_text = result.stdout_path.read_text(
                    encoding="utf-8", errors="replace",
                )
            except OSError as error:
                raise OutputValidationError(
                    "samtools coverage stdout log is missing or unreadable: "
                    f"'{result.stdout_path}'."
                ) from error
        reference_bases = 0
        covered_bases = 0
        depth_bases = 0.0
        for line in coverage_text.splitlines():
            if not line or line.startswith("#"):
                continue
            columns = line.split("\t")
            if len(columns) < 7:
                raise OutputValidationError("samtools coverage returned a malformed row.")
            try:
                length = int(columns[2]) - int(columns[1]) + 1
                row_covered = int(columns[4])
                mean_depth = float(columns[6])
            except ValueError as error:
                raise OutputValidationError(
                    "samtools coverage returned non-numeric coverage fields."
                ) from error
            if length <= 0 or row_covered < 0 or mean_depth < 0:
                raise OutputValidationError("samtools coverage returned invalid coverage values.")
            reference_bases += length
            covered_bases += row_covered
            depth_bases += mean_depth * length
        if reference_bases == 0:
            raise OutputValidationError("samtools coverage returned no reference rows.")
        return AlignmentCoverageResult(
            sample_id, bam, depth_bases / reference_bases, covered_bases,
            reference_bases, command, version,
        )


def build_raw_fastq_subsampling_plan(
    *, study_id: str, sample_id: str, source_fastqs: Sequence[Path],
    baseline_coverage: float, target_coverages: Sequence[float],
    replicates: int, base_seed: int, output_root: Path,
) -> RawFastqSubsamplingPlan:
    """Plan shared nested raw-read subsets without loading reads into memory."""
    if not isinstance(study_id, str) or _SAFE_IDENTIFIER.fullmatch(study_id) is None:
        raise InputValidationError("Raw FASTQ study_id must be a safe non-empty identifier.")
    if not isinstance(sample_id, str) or not sample_id.strip():
        raise InputValidationError("Raw FASTQ sample_id must be non-empty.")
    inputs = tuple(Path(path) for path in source_fastqs)
    if not inputs or len(set(inputs)) != len(inputs):
        raise InputValidationError("source_fastqs must contain unique FASTQ paths.")
    if any(not _is_fastq_path(path) for path in inputs):
        raise InputValidationError("Raw low-depth inputs must use .fastq, .fq, .fastq.gz, or .fq.gz.")
    if isinstance(baseline_coverage, bool) or not isinstance(baseline_coverage, (int, float)) or baseline_coverage <= 0:
        raise InputValidationError("baseline_coverage must be positive.")
    targets = tuple(float(value) for value in target_coverages)
    if (
        not targets
        or any(value <= 0 or value > float(baseline_coverage) for value in targets)
        or tuple(sorted(set(targets))) != targets
    ):
        raise InputValidationError(
            "Target coverages must be unique ascending positive values no greater than baseline coverage."
        )
    if not isinstance(replicates, int) or isinstance(replicates, bool) or replicates < 1:
        raise InputValidationError("replicates must be a positive integer.")
    if not isinstance(base_seed, int) or isinstance(base_seed, bool) or base_seed < 0:
        raise InputValidationError("base_seed must be a non-negative integer.")
    root = Path(output_root)
    subsets = tuple(
        RawFastqSubset(
            target_coverage=coverage,
            baseline_coverage=float(baseline_coverage),
            sampling_fraction=coverage / float(baseline_coverage),
            replicate=replicate,
            seed=base_seed + replicate,
            output_fastq=root / _coverage_label(coverage) / f"replicate-{replicate}.fastq.gz",
        )
        for replicate in range(1, replicates + 1)
        for coverage in targets
    )
    return RawFastqSubsamplingPlan(
        study_id, sample_id, inputs, float(baseline_coverage), subsets,
    )


def build_raw_fastq_subsampling_plan_from_config(
    section: Mapping[str, object], output_root: Path,
) -> RawFastqSubsamplingPlan:
    """Build a raw FASTQ plan from the validated effective config section."""
    try:
        source_fastqs = section["source_fastqs"]
        if not isinstance(source_fastqs, Sequence) or isinstance(source_fastqs, (str, bytes)):
            raise TypeError("source_fastqs must be a sequence")
        return build_raw_fastq_subsampling_plan(
            study_id=str(section["study_id"]),
            sample_id=str(section["sample_id"]),
            source_fastqs=tuple(Path(str(path)) for path in source_fastqs),
            baseline_coverage=float(section["baseline_coverage"]),
            target_coverages=section["target_coverages"],  # type: ignore[arg-type]
            replicates=int(section["replicates"]),
            base_seed=int(section["base_seed"]),
            output_root=Path(output_root),
        )
    except (KeyError, TypeError, ValueError) as error:
        raise InputValidationError(
            f"Invalid raw FASTQ low_depth configuration: {error}"
        ) from error


class FastqNestedDownsampler:
    """Stream FASTQ once per replicate and fan records into nested outputs."""

    def run(
        self, plan: RawFastqSubsamplingPlan, *, replicate: int,
        dry_run: bool = False, overwrite: bool = False,
    ) -> RawFastqSubsetResult:
        subsets = tuple(item for item in plan.subsets if item.replicate == replicate)
        if not subsets:
            raise InputValidationError(f"Unknown raw FASTQ replicate: {replicate}.")
        if dry_run:
            return RawFastqSubsetResult(plan, LowDepthRunStatus.PLANNED, {}, {})
        for source in plan.source_fastqs:
            if not source.is_file() or source.stat().st_size == 0:
                raise InputValidationError(f"Raw FASTQ input is missing or empty: '{source}'.")
        for subset in subsets:
            if subset.output_fastq.exists() and not overwrite:
                raise OutputValidationError(f"Raw FASTQ output exists: '{subset.output_fastq}'.")

        partials = {item: Path(f"{item.output_fastq}.partial") for item in subsets}
        for partial in partials.values():
            if partial.exists():
                raise OutputValidationError(f"Raw FASTQ partial output exists: '{partial}'.")
            partial.parent.mkdir(parents=True, exist_ok=True)
        read_counts = {str(item.output_fastq): 0 for item in subsets}
        base_counts = {str(item.output_fastq): 0 for item in subsets}
        try:
            with ExitStack() as stack:
                outputs = {
                    item: stack.enter_context(
                        gzip.GzipFile(filename=partials[item], mode="wb", mtime=0)
                    )
                    for item in subsets
                }
                for source in plan.source_fastqs:
                    with _open_fastq_binary(source) as handle:
                        for record_number, record in enumerate(_iter_fastq_records(handle), start=1):
                            header, sequence, plus, quality = record
                            _validate_fastq_record(source, record_number, header, sequence, plus, quality)
                            molecule = _pacbio_molecule_key(header)
                            score = _stable_sampling_score(molecule, subsets[0].seed)
                            payload = b"".join(record)
                            sequence_bases = len(sequence.rstrip(b"\r\n"))
                            for subset in subsets:
                                if score < subset.sampling_fraction:
                                    outputs[subset].write(payload)
                                    key = str(subset.output_fastq)
                                    read_counts[key] += 1
                                    base_counts[key] += sequence_bases
            for subset in subsets:
                partials[subset].replace(subset.output_fastq)
        except Exception:
            for partial in partials.values():
                if partial.is_file():
                    partial.unlink()
            raise
        for subset in subsets:
            if not subset.output_fastq.is_file() or subset.output_fastq.stat().st_size == 0:
                raise OutputValidationError(
                    f"Raw FASTQ subset is missing or empty: '{subset.output_fastq}'."
                )
        return RawFastqSubsetResult(
            plan, LowDepthRunStatus.COMPLETED, read_counts, base_counts,
        )


def _open_fastq_binary(path: Path):
    return gzip.open(path, "rb") if path.name.lower().endswith(".gz") else path.open("rb")


def _iter_fastq_records(handle):
    while True:
        header = handle.readline()
        if not header:
            return
        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()
        if not sequence or not plus or not quality:
            raise InputValidationError("FASTQ ended inside a four-line record.")
        yield header, sequence, plus, quality


def _validate_fastq_record(
    source: Path, record_number: int, header: bytes, sequence: bytes,
    plus: bytes, quality: bytes,
) -> None:
    sequence_value = sequence.rstrip(b"\r\n")
    quality_value = quality.rstrip(b"\r\n")
    if not header.startswith(b"@") or not plus.startswith(b"+"):
        raise InputValidationError(
            f"Malformed FASTQ record {record_number} in '{source}'."
        )
    if not sequence_value or len(sequence_value) != len(quality_value):
        raise InputValidationError(
            f"FASTQ sequence/quality length mismatch at record {record_number} in '{source}'."
        )


def _pacbio_molecule_key(header: bytes) -> bytes:
    read_name = header[1:].split(None, 1)[0]
    parts = read_name.split(b"/")
    return b"/".join(parts[:2]) if len(parts) >= 2 else read_name


def _stable_sampling_score(molecule_key: bytes, seed: int) -> float:
    digest = hashlib.blake2b(
        str(seed).encode("ascii") + b"\0" + molecule_key,
        digest_size=8,
        person=b"HiFiVarP15",
    ).digest()
    return int.from_bytes(digest, "big") / float(1 << 64)


def _is_fastq_path(path: Path) -> bool:
    name = path.name.lower()
    return name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))


@dataclass(frozen=True, slots=True)
class ParameterProfile:
    """Named, explicit caller/workflow settings; no profile is called optimal."""

    name: str
    settings: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or _SAFE_IDENTIFIER.fullmatch(self.name) is None:
            raise InputValidationError("Low-depth profile name must be a safe non-empty identifier.")
        if not isinstance(self.settings, Mapping):
            raise InputValidationError("Low-depth profile settings must be a mapping.")
        object.__setattr__(self, "settings", dict(self.settings))

    def to_dict(self) -> dict[str, object]:
        return {"name": self.name, "settings": dict(self.settings)}


@dataclass(frozen=True, slots=True)
class AccuracyThreshold:
    variant_class: str
    metric: str
    minimum: float

    def __post_init__(self) -> None:
        if not self.variant_class.strip() or not self.metric.strip():
            raise InputValidationError("Accuracy threshold class and metric must be non-empty.")
        if isinstance(self.minimum, bool) or not isinstance(self.minimum, (int, float)) or not 0 <= float(self.minimum) <= 1:
            raise InputValidationError("Accuracy threshold minimum must be between 0 and 1.")

    @property
    def key(self) -> str:
        return f"{self.variant_class}.{self.metric}"

    def to_dict(self) -> dict[str, object]:
        return {"variant_class": self.variant_class, "metric": self.metric, "minimum": float(self.minimum)}


@dataclass(frozen=True, slots=True)
class LowDepthRun:
    run_id: str
    sample_id: str
    profile: ParameterProfile
    target_coverage: float
    baseline_coverage: float
    sampling_fraction: float
    replicate: int
    seed: int
    output_bam: Path
    output_index: Path

    def to_dict(self) -> dict[str, object]:
        return {
            "run_id": self.run_id, "sample_id": self.sample_id,
            "profile": self.profile.to_dict(), "target_coverage": self.target_coverage,
            "baseline_coverage": self.baseline_coverage,
            "sampling_fraction": self.sampling_fraction, "replicate": self.replicate,
            "seed": self.seed, "output_bam": str(self.output_bam),
            "output_index": str(self.output_index),
        }


@dataclass(frozen=True, slots=True)
class LowDepthStudyPlan:
    study_id: str
    sample_id: str
    source_bam: Path
    baseline_coverage: float
    runs: tuple[LowDepthRun, ...]
    thresholds: tuple[AccuracyThreshold, ...]
    minimum_successful_replicates: int
    minimum_pass_fraction: float
    created_at: str = field(default_factory=utc_now_iso8601)

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": "1.0", "study_id": self.study_id,
            "sample_id": self.sample_id, "source_bam": str(self.source_bam),
            "source_contract": "coordinate_sorted_aligned_bam",
            "baseline_coverage": self.baseline_coverage,
            "runs": [run.to_dict() for run in self.runs],
            "thresholds": [item.to_dict() for item in self.thresholds],
            "minimum_successful_replicates": self.minimum_successful_replicates,
            "minimum_pass_fraction": self.minimum_pass_fraction,
            "created_at": self.created_at,
            "scientific_semantics": {
                "target_coverage_is_not_observed_coverage": True,
                "minimum_supported_depth_is_dataset_and_policy_specific": True,
                "parameter_profile_is_not_claimed_optimal": True,
                "subsamples_are_paired_across_profiles": True,
                "subsamples_are_nested_across_depths": True,
                "downsampling_does_not_modify_source": True,
            },
        }

    def write(self, path: Path, *, overwrite: bool = False) -> Path:
        write_json_atomic(self.to_dict(), Path(path), overwrite=overwrite, artifact_name="Low-depth plan")
        return Path(path)


def build_low_depth_plan(
    *, study_id: str, sample_id: str, source_bam: Path, baseline_coverage: float,
    target_coverages: Sequence[float], replicates: int, base_seed: int,
    profiles: Sequence[ParameterProfile], thresholds: Sequence[AccuracyThreshold],
    output_root: Path, minimum_successful_replicates: int,
    minimum_pass_fraction: float = 1.0,
) -> LowDepthStudyPlan:
    """Build a deterministic profile × coverage × replicate experiment matrix."""
    if not isinstance(study_id, str) or _SAFE_IDENTIFIER.fullmatch(study_id) is None:
        raise InputValidationError("Low-depth study_id must be a safe non-empty identifier.")
    if not isinstance(sample_id, str) or not sample_id.strip():
        raise InputValidationError("Low-depth sample_id must be non-empty.")
    if (
        isinstance(baseline_coverage, bool)
        or not isinstance(baseline_coverage, (int, float))
        or baseline_coverage <= 0
    ):
        raise InputValidationError("baseline_coverage must be positive.")
    targets = tuple(float(value) for value in target_coverages)
    if not targets or any(value <= 0 or value > baseline_coverage for value in targets):
        raise InputValidationError("Target coverages must be positive and no greater than baseline coverage.")
    if tuple(sorted(set(targets))) != targets:
        raise InputValidationError("Target coverages must be unique and ascending.")
    if not isinstance(replicates, int) or isinstance(replicates, bool) or replicates < 1:
        raise InputValidationError("replicates must be a positive integer.")
    if not isinstance(base_seed, int) or isinstance(base_seed, bool) or base_seed < 0:
        raise InputValidationError("base_seed must be a non-negative integer.")
    if (
        not isinstance(minimum_successful_replicates, int)
        or isinstance(minimum_successful_replicates, bool)
        or minimum_successful_replicates < 1
        or minimum_successful_replicates > replicates
    ):
        raise InputValidationError("minimum_successful_replicates must be between 1 and replicates.")
    if (
        isinstance(minimum_pass_fraction, bool)
        or not isinstance(minimum_pass_fraction, (int, float))
        or not 0 < minimum_pass_fraction <= 1
    ):
        raise InputValidationError("minimum_pass_fraction must be in (0, 1].")
    profiles = tuple(profiles)
    thresholds = tuple(thresholds)
    if not profiles or len({profile.name for profile in profiles}) != len(profiles):
        raise InputValidationError("At least one uniquely named low-depth profile is required.")
    if not thresholds or len({item.key for item in thresholds}) != len(thresholds):
        raise InputValidationError("At least one unique preregistered accuracy threshold is required.")
    root = Path(output_root)
    runs: list[LowDepthRun] = []
    for profile in profiles:
        for coverage in targets:
            label = _coverage_label(coverage)
            for replicate in range(1, replicates + 1):
                # samtools applies the fraction to a deterministic hash for a
                # seed. Reusing the replicate seed across depths makes the
                # lower-fraction sample a nested subset of the higher-fraction
                # sample, while separate replicates remain independent.
                seed = base_seed + replicate
                run_id = f"{profile.name}.{label}.r{replicate}"
                bam = root / label / f"replicate-{replicate}.bam"
                runs.append(LowDepthRun(run_id, sample_id, profile, coverage, float(baseline_coverage), coverage / float(baseline_coverage), replicate, seed, bam, Path(f"{bam}.bai")))
    return LowDepthStudyPlan(study_id, sample_id, Path(source_bam), float(baseline_coverage), tuple(runs), thresholds, minimum_successful_replicates, float(minimum_pass_fraction))


def build_low_depth_plan_from_config(
    section: Mapping[str, object],
    output_root: Path,
) -> LowDepthStudyPlan:
    """Build a validated plan from the effective ``low_depth`` config section."""
    try:
        profiles = tuple(
            ParameterProfile(str(item["name"]), item["settings"])
            for item in section["profiles"]  # type: ignore[union-attr]
        )
        thresholds = tuple(
            AccuracyThreshold(
                str(item["variant_class"]),
                str(item["metric"]),
                float(item["minimum"]),
            )
            for item in section["thresholds"]  # type: ignore[union-attr]
        )
        return build_low_depth_plan(
            study_id=str(section["study_id"]),
            sample_id=str(section["sample_id"]),
            source_bam=Path(str(section["source_bam"])),
            baseline_coverage=float(section["baseline_coverage"]),
            target_coverages=section["target_coverages"],  # type: ignore[arg-type]
            replicates=int(section["replicates"]),
            base_seed=int(section["base_seed"]),
            profiles=profiles,
            thresholds=thresholds,
            output_root=Path(output_root),
            minimum_successful_replicates=int(section["minimum_successful_replicates"]),
            minimum_pass_fraction=float(section["minimum_pass_fraction"]),
        )
    except (KeyError, TypeError, ValueError) as error:
        raise InputValidationError(
            f"Invalid effective low_depth configuration: {error}"
        ) from error


@dataclass(frozen=True, slots=True)
class LowDepthObservation:
    run_id: str
    profile: str
    target_coverage: float
    replicate: int
    status: LowDepthRunStatus
    metrics: Mapping[str, float]
    observed_coverage: float | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.status, LowDepthRunStatus):
            raise InputValidationError("Low-depth observation status is invalid.")
        if not isinstance(self.metrics, Mapping):
            raise InputValidationError("Low-depth observation metrics must be a mapping.")
        normalized: dict[str, float] = {}
        for key, value in self.metrics.items():
            if not isinstance(key, str) or not key.strip():
                raise InputValidationError("Low-depth metric names must be non-empty strings.")
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise InputValidationError(f"Low-depth metric '{key}' must be numeric.")
            normalized[key] = float(value)
        object.__setattr__(self, "metrics", normalized)
        if self.observed_coverage is not None and (
            isinstance(self.observed_coverage, bool)
            or not isinstance(self.observed_coverage, (int, float))
            or self.observed_coverage <= 0
        ):
            raise InputValidationError("observed_coverage must be positive or null.")

    def to_dict(self) -> dict[str, object]:
        return {
            "run_id": self.run_id,
            "profile": self.profile,
            "target_coverage": self.target_coverage,
            "observed_coverage": self.observed_coverage,
            "replicate": self.replicate,
            "status": self.status.value,
            "metrics": dict(self.metrics),
        }


@dataclass(frozen=True, slots=True)
class DepthAssessment:
    profile: str
    target_coverage: float
    successful_replicates: int
    passing_replicates: int
    pass_fraction: float
    supported: bool
    observed_coverage_range: tuple[float, float] | None
    metric_ranges: Mapping[str, tuple[float, float]]

    def to_dict(self) -> dict[str, object]:
        return {"profile": self.profile, "target_coverage": self.target_coverage,
                "successful_replicates": self.successful_replicates,
                "passing_replicates": self.passing_replicates,
                "pass_fraction": self.pass_fraction, "supported": self.supported,
                "observed_coverage_range": list(self.observed_coverage_range) if self.observed_coverage_range else None,
                "metric_ranges": {key: list(value) for key, value in self.metric_ranges.items()}}


@dataclass(frozen=True, slots=True)
class LowDepthStudyResult:
    plan: LowDepthStudyPlan
    observations: tuple[LowDepthObservation, ...]
    assessments: tuple[DepthAssessment, ...]
    minimum_supported_depth_by_profile: Mapping[str, float | None]

    @property
    def study_id(self) -> str:
        return self.plan.study_id

    def to_dict(self) -> dict[str, object]:
        return {"schema_version": "1.0", "study_id": self.study_id,
                "plan": self.plan.to_dict(),
                "observations": [item.to_dict() for item in self.observations],
                "assessments": [item.to_dict() for item in self.assessments],
                "minimum_supported_depth_by_profile": dict(self.minimum_supported_depth_by_profile),
                "scientific_semantics": {
                    "result_is_not_a_universal_minimum_depth": True,
                    "caller_optimization_requires_independent_validation": True,
                    "accuracy_is_not_clinical_utility": True,
                }}

    def write(self, root: Path, *, overwrite: bool = False) -> tuple[Path, Path]:
        root = Path(root); json_path = root / "low_depth_result.json"; yaml_path = root / "low_depth_result.yaml"
        write_json_atomic(self.to_dict(), json_path, overwrite=overwrite, artifact_name="Low-depth result")
        write_yaml_atomic(self.to_dict(), yaml_path, overwrite=overwrite, artifact_name="Low-depth result")
        return json_path, yaml_path


def evaluate_low_depth_study(plan: LowDepthStudyPlan, observations: Sequence[LowDepthObservation]) -> LowDepthStudyResult:
    """Evaluate replicated metrics without selecting a universal depth or profile."""
    run_map = {run.run_id: run for run in plan.runs}
    if len({item.run_id for item in observations}) != len(observations):
        raise InputValidationError("Low-depth observations contain duplicate run IDs.")
    for item in observations:
        run = run_map.get(item.run_id)
        if run is None or item.profile != run.profile.name or item.target_coverage != run.target_coverage or item.replicate != run.replicate:
            raise InputValidationError(f"Observation does not match the study plan: '{item.run_id}'.")
        if item.status is LowDepthRunStatus.COMPLETED:
            if item.observed_coverage is None:
                raise InputValidationError(
                    f"Completed observation '{item.run_id}' requires observed_coverage."
                )
            for threshold in plan.thresholds:
                if threshold.key not in item.metrics:
                    raise InputValidationError(
                        f"Completed observation '{item.run_id}' is missing metric '{threshold.key}'."
                    )
                value = item.metrics[threshold.key]
                if not 0 <= value <= 1:
                    raise InputValidationError(
                        f"Accuracy metric '{threshold.key}' for '{item.run_id}' must be between 0 and 1."
                    )
    grouped: dict[tuple[str, float], list[LowDepthObservation]] = {}
    for item in observations:
        grouped.setdefault((item.profile, item.target_coverage), []).append(item)
    assessments: list[DepthAssessment] = []
    for profile in dict.fromkeys(run.profile.name for run in plan.runs):
        for coverage in dict.fromkeys(run.target_coverage for run in plan.runs if run.profile.name == profile):
            successful = [item for item in grouped.get((profile, coverage), []) if item.status is LowDepthRunStatus.COMPLETED]
            passing = [item for item in successful if all(item.metrics.get(threshold.key, -1.0) >= threshold.minimum for threshold in plan.thresholds)]
            fraction = len(passing) / len(successful) if successful else 0.0
            ranges = {threshold.key: (min(values), max(values)) for threshold in plan.thresholds if (values := [item.metrics[threshold.key] for item in successful if threshold.key in item.metrics])}
            observed_values = [float(item.observed_coverage) for item in successful if item.observed_coverage is not None]
            observed_range = (min(observed_values), max(observed_values)) if observed_values else None
            supported = len(successful) >= plan.minimum_successful_replicates and fraction >= plan.minimum_pass_fraction and len(ranges) == len(plan.thresholds)
            assessments.append(DepthAssessment(profile, coverage, len(successful), len(passing), fraction, supported, observed_range, ranges))
    recommendations: dict[str, float | None] = {}
    for profile in dict.fromkeys(item.profile for item in assessments):
        ordered = [item for item in assessments if item.profile == profile]
        recommendations[profile] = next((item.target_coverage for index, item in enumerate(ordered) if item.supported and all(higher.supported for higher in ordered[index:])), None)
    return LowDepthStudyResult(plan, tuple(observations), tuple(assessments), recommendations)


@dataclass(frozen=True, slots=True)
class DownsampleResult:
    run: LowDepthRun
    status: LowDepthRunStatus
    commands: tuple[tuple[str, ...], ...]
    tool_version: str | None = None


class SamtoolsDownsampleWrapper:
    """Deterministically subsample and index one aligned BAM via CommandRunner."""

    def __init__(self, *, runner: CommandRunner | None = None, executable: str | Path = "samtools") -> None:
        self.runner = runner or CommandRunner(); self.executable = str(executable)

    def detect_version(self) -> str:
        self.runner.require_executable(self.executable)
        result = self.runner.run([self.executable, "--version"])
        match = _VERSION.search("\n".join(filter(None, (result.stdout, result.stderr))))
        if match is None:
            raise ToolVersionError("Unable to parse samtools version for low-depth downsampling.")
        return match.group(1)

    def plan_commands(self, source_bam: Path, run: LowDepthRun, *, threads: int = 1) -> tuple[tuple[str, ...], ...]:
        if threads < 1: raise InputValidationError("Downsampling threads must be positive.")
        view = (self.executable, "view", "-@", str(threads), "--subsample-seed", str(run.seed), "--subsample", format(run.sampling_fraction, ".12g"), "-b", "-o", str(run.output_bam), str(source_bam))
        index = (self.executable, "index", "-@", str(threads), "-b", str(run.output_bam), str(run.output_index))
        quickcheck = (self.executable, "quickcheck", "-v", str(run.output_bam))
        return view, index, quickcheck

    def validate_source(self, source_bam: Path) -> None:
        """Require a readable coordinate-sorted BAM before materialization."""
        result = self.runner.run([self.executable, "view", "-H", str(source_bam)])
        header_lines = result.stdout.splitlines()
        hd_fields = next(
            (line.split("\t") for line in header_lines if line.startswith("@HD\t")),
            (),
        )
        sort_order = next(
            (field.removeprefix("SO:") for field in hd_fields if field.startswith("SO:")),
            None,
        )
        if sort_order != "coordinate":
            raise InputValidationError(
                "Low-depth source BAM must declare @HD SO:coordinate; "
                f"observed {sort_order or 'missing'!r}."
            )

    def run(self, source_bam: Path, run: LowDepthRun, *, threads: int = 1, dry_run: bool = False, overwrite: bool = False, stderr_path: Path | None = None, redact_values: Collection[str] | None = None) -> DownsampleResult:
        source_bam = Path(source_bam); commands = self.plan_commands(source_bam, run, threads=threads)
        if dry_run:
            for command in commands: self.runner.run(command, dry_run=True, stderr_path=stderr_path, redact_values=redact_values)
            return DownsampleResult(run, LowDepthRunStatus.PLANNED, commands)
        if not source_bam.is_file() or source_bam.stat().st_size == 0:
            raise InputValidationError(f"Low-depth source BAM is missing or empty: '{source_bam}'.")
        for output in (run.output_bam, run.output_index):
            if output.exists() and not overwrite: raise OutputValidationError(f"Low-depth output exists: '{output}'.")
        version = self.detect_version()
        self.validate_source(source_bam)
        if overwrite:
            for output in (run.output_bam, run.output_index):
                if output.is_file():
                    output.unlink()
        run.output_bam.parent.mkdir(parents=True, exist_ok=True)
        for command in commands[:2]: self.runner.run(command, stderr_path=stderr_path, redact_values=redact_values)
        for output in (run.output_bam, run.output_index):
            if not output.is_file() or output.stat().st_size == 0: raise OutputValidationError(f"Low-depth output is missing or empty: '{output}'.")
        self.runner.run(commands[2], stderr_path=stderr_path, redact_values=redact_values)
        return DownsampleResult(run, LowDepthRunStatus.COMPLETED, commands, version)


def _coverage_label(value: float) -> str:
    return f"{value:g}".replace(".", "p") + "x"


__all__ = [
    "AccuracyThreshold", "AlignmentCoverageResult", "DepthAssessment", "DownsampleResult",
    "FastqNestedDownsampler", "LowDepthObservation", "LowDepthRun",
    "LowDepthBenchmarkPlan", "LowDepthBenchmarkRun",
    "LowDepthRunStatus", "LowDepthStudyPlan", "LowDepthStudyResult",
    "ParameterProfile", "RawFastqSubset", "RawFastqSubsetResult",
    "RawFastqSubsamplingPlan", "SamtoolsCoverageWrapper", "SamtoolsDownsampleWrapper",
    "build_low_depth_plan", "build_low_depth_plan_from_config",
    "build_low_depth_benchmark_plan",
    "build_raw_fastq_subsampling_plan",
    "build_raw_fastq_subsampling_plan_from_config",
    "evaluate_low_depth_study",
]
