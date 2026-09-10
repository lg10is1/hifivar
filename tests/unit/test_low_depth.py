from __future__ import annotations

from pathlib import Path
from unittest.mock import Mock

import pytest

from hifivar.command import CommandResult
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.low_depth import (
    AccuracyThreshold, FastqNestedDownsampler, LowDepthObservation,
    LowDepthRunStatus, ParameterProfile, SamtoolsDownsampleWrapper,
    SamtoolsCoverageWrapper,
    build_low_depth_plan, build_raw_fastq_subsampling_plan,
    build_low_depth_benchmark_plan,
    evaluate_low_depth_study,
)


def _plan(tmp_path: Path):
    return build_low_depth_plan(
        study_id="HG002-depth", sample_id="HG002", source_bam=tmp_path / "source.bam",
        baseline_coverage=30, target_coverages=(5, 10, 20, 30), replicates=2,
        base_seed=1701, profiles=(ParameterProfile("baseline"),),
        thresholds=(AccuracyThreshold("SNP", "f1", 0.98), AccuracyThreshold("INDEL", "f1", 0.90)),
        output_root=tmp_path / "out", minimum_successful_replicates=2,
    )


def test_plan_is_deterministic_complete_and_serializable(tmp_path: Path) -> None:
    first = _plan(tmp_path); second = _plan(tmp_path)
    assert [run.to_dict() for run in first.runs] == [run.to_dict() for run in second.runs]
    assert len(first.runs) == 8
    assert first.runs[0].sampling_fraction == pytest.approx(1 / 6)
    assert first.runs[-1].sampling_fraction == 1
    assert len({run.seed for run in first.runs}) == 2
    for replicate in (1, 2):
        assert len({run.seed for run in first.runs if run.replicate == replicate}) == 1
    assert first.to_dict()["scientific_semantics"]["minimum_supported_depth_is_dataset_and_policy_specific"] is True
    assert first.to_dict()["scientific_semantics"]["subsamples_are_nested_across_depths"] is True
    first.write(tmp_path / "plan.json")
    with pytest.raises(OutputValidationError): first.write(tmp_path / "plan.json")


def test_plan_rejects_upsampling_unsorted_targets_and_missing_thresholds(tmp_path: Path) -> None:
    common = dict(study_id="B", sample_id="S", source_bam=tmp_path/"x.bam", baseline_coverage=20,
                  replicates=2, base_seed=1, profiles=(ParameterProfile("p"),), output_root=tmp_path,
                  minimum_successful_replicates=2)
    with pytest.raises(InputValidationError): build_low_depth_plan(target_coverages=(10, 30), thresholds=(AccuracyThreshold("SNP","f1",.9),), **common)
    with pytest.raises(InputValidationError): build_low_depth_plan(target_coverages=(10, 5), thresholds=(AccuracyThreshold("SNP","f1",.9),), **common)
    with pytest.raises(InputValidationError): build_low_depth_plan(target_coverages=(5, 10), thresholds=(), **common)


def test_evaluation_requires_replicated_thresholds_and_monotonic_support(tmp_path: Path) -> None:
    plan = _plan(tmp_path); observations = []
    values = {5: (0.97, 0.88), 10: (0.981, 0.901), 20: (0.99, 0.92), 30: (0.995, 0.93)}
    for run in plan.runs:
        snp, indel = values[run.target_coverage]
        observations.append(LowDepthObservation(run.run_id, run.profile.name, run.target_coverage, run.replicate,
            LowDepthRunStatus.COMPLETED, {"SNP.f1": snp, "INDEL.f1": indel}, run.target_coverage * .98))
    result = evaluate_low_depth_study(plan, observations)
    assert result.minimum_supported_depth_by_profile == {"baseline": 10}
    assert result.to_dict()["scientific_semantics"]["result_is_not_a_universal_minimum_depth"] is True


def test_evaluation_rejects_duplicate_and_unknown_observations(tmp_path: Path) -> None:
    plan = _plan(tmp_path); run = plan.runs[0]
    observation = LowDepthObservation(run.run_id, "baseline", run.target_coverage, run.replicate, LowDepthRunStatus.COMPLETED, {})
    with pytest.raises(InputValidationError, match="duplicate"): evaluate_low_depth_study(plan, (observation, observation))
    bad = LowDepthObservation("unknown", "baseline", 5, 1, LowDepthRunStatus.COMPLETED, {})
    with pytest.raises(InputValidationError, match="does not match"): evaluate_low_depth_study(plan, (bad,))


def test_completed_observation_requires_measured_coverage_and_all_metrics(tmp_path: Path) -> None:
    plan = _plan(tmp_path)
    run = plan.runs[0]
    missing_coverage = LowDepthObservation(
        run.run_id, run.profile.name, run.target_coverage, run.replicate,
        LowDepthRunStatus.COMPLETED, {"SNP.f1": 0.99, "INDEL.f1": 0.91},
    )
    with pytest.raises(InputValidationError, match="observed_coverage"):
        evaluate_low_depth_study(plan, (missing_coverage,))
    missing_metric = LowDepthObservation(
        run.run_id, run.profile.name, run.target_coverage, run.replicate,
        LowDepthRunStatus.COMPLETED, {"SNP.f1": 0.99}, 4.9,
    )
    with pytest.raises(InputValidationError, match="missing metric"):
        evaluate_low_depth_study(plan, (missing_metric,))


def test_samtools_command_is_shell_free_deterministic_and_dry_run(tmp_path: Path) -> None:
    run = _plan(tmp_path).runs[0]; runner = Mock()
    runner.run.return_value = CommandResult(("samtools",), None, None, None, 0, None, False)
    wrapper = SamtoolsDownsampleWrapper(runner=runner)
    commands = wrapper.plan_commands(tmp_path/"source.bam", run, threads=4)
    assert commands[0][:4] == ("samtools", "view", "-@", "4")
    assert "--subsample-seed" in commands[0] and "--subsample" in commands[0]
    assert commands[1][1] == "index"
    assert commands[2][1:3] == ("quickcheck", "-v")
    result = wrapper.run(tmp_path/"missing.bam", run, dry_run=True)
    assert result.status is LowDepthRunStatus.PLANNED
    assert runner.run.call_count == 3


def test_samtools_execution_validates_outputs_and_preserves_source(tmp_path: Path) -> None:
    plan = _plan(tmp_path); run = plan.runs[0]; source = plan.source_bam
    source.write_bytes(b"source"); before = source.read_bytes(); runner = Mock()
    runner.require_executable.return_value = Path("samtools")
    def execute(command, **kwargs):
        if command[1] == "--version":
            return CommandResult(tuple(command), 0, "samtools 1.21\n", "", 0.1, None, True)
        if command[1:3] == ["view", "-H"]:
            return CommandResult(tuple(command), 0, "@HD\tVN:1.6\tSO:coordinate\n", "", 0.1, None, True)
        if command[1] == "quickcheck":
            return CommandResult(tuple(command), 0, "", "", 0.1, None, True)
        output = Path(command[command.index("-o") + 1]) if command[1] == "view" else Path(command[-1])
        output.write_bytes(b"output")
        return CommandResult(tuple(command), 0, "", "", 0.1, None, True)
    runner.run.side_effect = execute
    result = SamtoolsDownsampleWrapper(runner=runner).run(source, run, threads=2)
    assert result.status is LowDepthRunStatus.COMPLETED and result.tool_version == "1.21"
    assert source.read_bytes() == before
    with pytest.raises(OutputValidationError, match="exists"):
        SamtoolsDownsampleWrapper(runner=runner).run(source, run)


def test_samtools_execution_rejects_non_coordinate_source_before_overwrite(tmp_path: Path) -> None:
    run = _plan(tmp_path).runs[0]
    source = tmp_path / "source.bam"
    source.write_bytes(b"source")
    run.output_bam.parent.mkdir(parents=True)
    run.output_bam.write_bytes(b"existing")
    runner = Mock()
    runner.require_executable.return_value = Path("samtools")
    runner.run.side_effect = (
        CommandResult(("samtools", "--version"), 0, "samtools 1.21\n", "", 0.1, None, True),
        CommandResult(("samtools", "view", "-H"), 0, "@HD\tVN:1.6\tSO:queryname\n", "", 0.1, None, True),
    )
    with pytest.raises(InputValidationError, match="SO:coordinate"):
        SamtoolsDownsampleWrapper(runner=runner).run(source, run, overwrite=True)
    assert run.output_bam.read_bytes() == b"existing"


def _write_fastq(path: Path, count: int = 200, movie: str = "movie") -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for index in range(count):
            name = f"{movie}/{index}/ccs"
            sequence = "ACGT" * (1 + index % 3)
            handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def _fastq_names(path: Path) -> set[str]:
    import gzip

    with gzip.open(path, "rt", encoding="ascii") as handle:
        lines = handle.read().splitlines()
    return {lines[index][1:] for index in range(0, len(lines), 4)}


def test_raw_fastq_plan_and_materialization_are_nested_and_reproducible(tmp_path: Path) -> None:
    first_source = tmp_path / "part1.fastq"
    second_source = tmp_path / "part2.fastq"
    _write_fastq(first_source, 100)
    _write_fastq(second_source, 100, "movie2")
    plan = build_raw_fastq_subsampling_plan(
        study_id="HG002-raw", sample_id="HG002",
        source_fastqs=(first_source, second_source), baseline_coverage=20,
        target_coverages=(5, 10, 20), replicates=2, base_seed=1701,
        output_root=tmp_path / "raw",
    )
    assert len(plan.subsets) == 6
    assert {item.seed for item in plan.subsets if item.replicate == 1} == {1702}
    assert plan.to_dict()["source_contract"] == "pacbio_hifi_fastq"
    before = (first_source.read_bytes(), second_source.read_bytes())
    result = FastqNestedDownsampler().run(plan, replicate=1)
    by_depth = {
        item.target_coverage: _fastq_names(item.output_fastq)
        for item in plan.subsets if item.replicate == 1
    }
    assert by_depth[5] < by_depth[10] < by_depth[20]
    assert len(by_depth[20]) == 200
    assert result.status is LowDepthRunStatus.COMPLETED
    assert first_source.read_bytes() == before[0]
    assert second_source.read_bytes() == before[1]

    duplicate_plan = build_raw_fastq_subsampling_plan(
        study_id="HG002-copy", sample_id="HG002",
        source_fastqs=(first_source, second_source), baseline_coverage=20,
        target_coverages=(5, 10, 20), replicates=1, base_seed=1701,
        output_root=tmp_path / "copy",
    )
    FastqNestedDownsampler().run(duplicate_plan, replicate=1)
    for original, duplicate in zip(plan.subsets[:3], duplicate_plan.subsets, strict=True):
        assert _fastq_names(original.output_fastq) == _fastq_names(duplicate.output_fastq)
        assert original.output_fastq.read_bytes() == duplicate.output_fastq.read_bytes()


def test_raw_fastq_replicates_are_distinct_and_dry_run_is_non_mutating(tmp_path: Path) -> None:
    source = tmp_path / "reads.fastq"
    _write_fastq(source, 300)
    plan = build_raw_fastq_subsampling_plan(
        study_id="raw", sample_id="S", source_fastqs=(source,),
        baseline_coverage=20, target_coverages=(10,), replicates=2,
        base_seed=9, output_root=tmp_path / "out",
    )
    dry = FastqNestedDownsampler().run(plan, replicate=1, dry_run=True)
    assert dry.status is LowDepthRunStatus.PLANNED
    assert not plan.subsets[0].output_fastq.exists()
    FastqNestedDownsampler().run(plan, replicate=1)
    FastqNestedDownsampler().run(plan, replicate=2)
    assert _fastq_names(plan.subsets[0].output_fastq) != _fastq_names(plan.subsets[1].output_fastq)


def test_raw_fastq_rejects_malformed_input_and_protects_outputs(tmp_path: Path) -> None:
    bad = tmp_path / "bad.fastq"
    bad.write_text("@movie/1/ccs\nACGT\n+\nIII\n", encoding="ascii")
    plan = build_raw_fastq_subsampling_plan(
        study_id="raw", sample_id="S", source_fastqs=(bad,),
        baseline_coverage=10, target_coverages=(10,), replicates=1,
        base_seed=1, output_root=tmp_path / "out",
    )
    with pytest.raises(InputValidationError, match="length mismatch"):
        FastqNestedDownsampler().run(plan, replicate=1)
    assert not plan.subsets[0].output_fastq.exists()
    assert not Path(f"{plan.subsets[0].output_fastq}.partial").exists()

    _write_fastq(bad, 10)
    plan.subsets[0].output_fastq.parent.mkdir(parents=True, exist_ok=True)
    plan.subsets[0].output_fastq.write_bytes(b"keep")
    with pytest.raises(OutputValidationError, match="exists"):
        FastqNestedDownsampler().run(plan, replicate=1)
    assert plan.subsets[0].output_fastq.read_bytes() == b"keep"


def test_end_to_end_plan_has_isolated_deterministic_artifacts(tmp_path: Path) -> None:
    source = tmp_path / "reads.fastq"
    _write_fastq(source, 20)
    raw = build_raw_fastq_subsampling_plan(
        study_id="P15-e2e", sample_id="HG002", source_fastqs=(source,),
        baseline_coverage=20, target_coverages=(5, 10), replicates=2,
        base_seed=17, output_root=tmp_path / "raw",
    )
    plan = build_low_depth_benchmark_plan(raw, tmp_path / "results")
    assert len(plan.runs) == 4
    assert len({run.run_id for run in plan.runs}) == 4
    assert len({run.alignment_bam for run in plan.runs}) == 4
    assert plan.runs[0].alignment_bam.name == "replicate-1.bam"
    assert plan.runs[0].happy_metrics().name.endswith(".metrics.json.gz")
    assert plan.to_dict()["scientific_semantics"]["each_subset_is_realigned_independently"] is True


def test_samtools_coverage_is_reference_weighted_and_dry_run_safe(tmp_path: Path) -> None:
    bam = tmp_path / "sample.bam"
    index = tmp_path / "sample.bam.bai"
    bam.write_bytes(b"bam")
    index.write_bytes(b"index")
    runner = Mock()
    runner.require_executable.return_value = Path("samtools")
    runner.run.side_effect = (
        CommandResult(("samtools", "--version"), 0, "samtools 1.21\n", "", 0.1, None, True),
        CommandResult(
            ("samtools", "coverage", str(bam)), 0,
            "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
            "chr1\t1\t100\t1\t50\t50\t10\t60\t60\n"
            "chr2\t1\t300\t1\t30\t10\t2\t60\t60\n",
            "", 0.1, None, True,
        ),
    )
    result = SamtoolsCoverageWrapper(runner=runner).run(
        sample_id="HG002", alignment_bam=bam, alignment_index=index,
    )
    assert result.observed_coverage == pytest.approx(4.0)
    assert result.reference_bases == 400
    dry_runner = Mock()
    dry_runner.run.return_value = CommandResult(
        ("samtools", "coverage", "missing.bam"), None, None, None, 0, None, False,
    )
    dry = SamtoolsCoverageWrapper(runner=dry_runner).run(
        sample_id="HG002", alignment_bam=tmp_path / "missing.bam",
        alignment_index=tmp_path / "missing.bam.bai", dry_run=True,
    )
    assert dry.observed_coverage == 0


def test_samtools_coverage_reads_persisted_stdout_when_logging_stderr(
    tmp_path: Path,
) -> None:
    bam = tmp_path / "sample.bam"
    index = tmp_path / "sample.bam.bai"
    bam.write_bytes(b"bam")
    index.write_bytes(b"bai")
    stderr_path = tmp_path / "coverage.log"
    stdout_path = tmp_path / "coverage.stdout.log"
    stdout_path.write_text(
        "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
        "chr1\t1\t100\t1\t50\t50\t10\t60\t60\n",
        encoding="utf-8",
    )
    runner = Mock()
    runner.require_executable.return_value = Path("samtools")
    runner.run.side_effect = (
        CommandResult(("samtools", "--version"), 0, "samtools 1.21\n", "", 0.1, None, True),
        CommandResult(
            ("samtools", "coverage", str(bam)), 0, None, None, 0.1, None,
            True, stdout_path=stdout_path, stderr_path=stderr_path,
        ),
    )

    result = SamtoolsCoverageWrapper(runner=runner).run(
        sample_id="HG002", alignment_bam=bam, alignment_index=index,
        stderr_path=stderr_path,
    )

    assert result.observed_coverage == pytest.approx(10.0)
    assert result.reference_bases == 100
    assert result.covered_bases == 50


def test_samtools_coverage_rejects_missing_stdout_contract(tmp_path: Path) -> None:
    bam = tmp_path / "sample.bam"
    index = tmp_path / "sample.bam.bai"
    bam.write_bytes(b"bam")
    index.write_bytes(b"bai")
    runner = Mock()
    runner.require_executable.return_value = Path("samtools")
    runner.run.side_effect = (
        CommandResult(("samtools", "--version"), 0, "samtools 1.21\n", "", 0.1, None, True),
        CommandResult(("samtools", "coverage", str(bam)), 0, None, None, 0.1, None, True),
    )

    with pytest.raises(OutputValidationError, match="no captured output"):
        SamtoolsCoverageWrapper(runner=runner).run(
            sample_id="HG002", alignment_bam=bam, alignment_index=index,
        )
