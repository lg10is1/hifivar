"""Tests for the Phase 2.1 lightweight input QC framework."""

from __future__ import annotations

import gzip
import json
import re
from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest
import yaml

from hifivar import __version__
from hifivar.context import AnalysisContext
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.qc import (
    QCIssue,
    QCIssueLevel,
    QCMetric,
    QCResult,
    QCStatus,
    RunQCReport,
    aggregate_qc_status,
    run_input_dataset_qc,
    run_input_qc,
)
from hifivar.reference import ReferenceGenome
from hifivar.sample import InputDataset, InputType, Sample
from hifivar.sample_sheet import SampleRecord


FASTQ_TEXT = "@read1\nACGT\n+\nIIII\n"


def write_reference(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    fasta = root / "reference.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    return fasta


def write_fastq(root: Path, name: str = "sample.fastq") -> Path:
    root.mkdir(parents=True, exist_ok=True)
    path = root / name
    path.write_text(FASTQ_TEXT, encoding="utf-8")
    return path


def write_gzip_fastq(root: Path, name: str = "sample.fastq.gz") -> Path:
    root.mkdir(parents=True, exist_ok=True)
    path = root / name
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write(FASTQ_TEXT)
    return path


def write_alignment(root: Path, name: str) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    path = root / name
    path.write_bytes(b"alignment-placeholder")
    return path


def write_index(alignment: Path, *, adjacent: bool = True) -> Path:
    suffix = ".bai" if alignment.suffix.lower() == ".bam" else ".crai"
    path = Path(f"{alignment}{suffix}") if adjacent else alignment.with_suffix(suffix)
    path.write_bytes(b"index-placeholder")
    return path


def make_context(
    tmp_path: Path,
    records: tuple[SampleRecord, ...] | None = None,
) -> AnalysisContext:
    reference = ReferenceGenome.from_fasta(write_reference(tmp_path), build="GRCh38")
    if records is None:
        dataset = InputDataset.from_files((write_fastq(tmp_path),))
        records = (SampleRecord(Sample("S1", dataset)),)
    config = {
        "reference": {
            "fasta": str(reference.fasta.absolute()),
            "build": reference.build,
        },
        "samples": {"sheet": None},
    }
    return AnalysisContext(reference, records, config)


def warning_issue() -> QCIssue:
    return QCIssue(
        code="ALIGNMENT_INDEX_MISSING",
        level=QCIssueLevel.WARNING,
        message="Alignment index is missing.",
    )


def test_qc_status_values_are_stable() -> None:
    assert [status.value for status in QCStatus] == [
        "pass",
        "warn",
        "fail",
        "not_checked",
    ]


def test_qc_metric_serialization_uses_standard_types() -> None:
    metric = QCMetric("file_size_bytes", 42, unit="bytes", description="Size")
    assert metric.to_dict() == {
        "name": "file_size_bytes",
        "value": 42,
        "unit": "bytes",
        "description": "Size",
    }


@pytest.mark.parametrize("name", ("FileSize", "file-size", "_size", ""))
def test_qc_metric_rejects_unstable_names(name: str) -> None:
    with pytest.raises(InputValidationError, match="snake_case"):
        QCMetric(name, 1)


@pytest.mark.parametrize("value", (Path("data"), InputType.FASTQ, object(), [1]))
def test_qc_metric_rejects_non_scalar_values(value: object) -> None:
    with pytest.raises(InputValidationError, match="metric value"):
        QCMetric("value", value)  # type: ignore[arg-type]


@pytest.mark.parametrize("value", (float("nan"), float("inf")))
def test_qc_metric_rejects_non_finite_float(value: float) -> None:
    with pytest.raises(InputValidationError, match="finite"):
        QCMetric("value", value)


def test_qc_issue_serialization_is_machine_readable() -> None:
    assert warning_issue().to_dict() == {
        "code": "ALIGNMENT_INDEX_MISSING",
        "level": "warning",
        "message": "Alignment index is missing.",
    }


@pytest.mark.parametrize("code", ("alignment_missing", "BAD-CODE", ""))
def test_qc_issue_rejects_invalid_code(code: str) -> None:
    with pytest.raises(InputValidationError, match="uppercase"):
        QCIssue(code, QCIssueLevel.WARNING, "message")


def test_pass_result_has_no_issues() -> None:
    result = QCResult(QCStatus.PASS, metrics=(QCMetric("file_count", 1),))
    assert result.status is QCStatus.PASS
    assert result.issues == ()


def test_warn_result_contains_warning_issue() -> None:
    result = QCResult(QCStatus.WARN, issues=(warning_issue(),))
    assert result.to_dict()["status"] == "warn"


def test_pass_result_cannot_hide_warning_issue() -> None:
    with pytest.raises(InputValidationError, match="cannot contain"):
        QCResult(QCStatus.PASS, issues=(warning_issue(),))


def test_error_issue_requires_fail_status() -> None:
    issue = QCIssue("DATA_QUALITY_FAILURE", QCIssueLevel.ERROR, "Quality failed.")
    with pytest.raises(InputValidationError, match="require FAIL"):
        QCResult(QCStatus.WARN, issues=(issue,))


def test_warn_status_requires_an_issue() -> None:
    with pytest.raises(InputValidationError, match="requires at least one"):
        QCResult(QCStatus.WARN)


def test_duplicate_metric_names_are_rejected() -> None:
    with pytest.raises(InputValidationError, match="duplicate metric"):
        QCResult(
            QCStatus.PASS,
            metrics=(QCMetric("file_count", 1), QCMetric("file_count", 2)),
        )


def test_fail_status_has_highest_aggregation_precedence() -> None:
    assert aggregate_qc_status((QCStatus.PASS, QCStatus.FAIL)) is QCStatus.FAIL


def test_only_not_checked_aggregates_to_not_checked() -> None:
    assert aggregate_qc_status((QCStatus.NOT_CHECKED,)) is QCStatus.NOT_CHECKED


def test_empty_aggregation_is_not_checked() -> None:
    assert aggregate_qc_status(()) is QCStatus.NOT_CHECKED


@pytest.mark.parametrize(
    "statuses,expected",
    (
        ((QCStatus.PASS, QCStatus.WARN), QCStatus.WARN),
        ((QCStatus.WARN, QCStatus.FAIL), QCStatus.FAIL),
        ((QCStatus.PASS, QCStatus.NOT_CHECKED), QCStatus.PASS),
    ),
)
def test_status_aggregation_precedence(
    statuses: tuple[QCStatus, ...], expected: QCStatus
) -> None:
    assert aggregate_qc_status(statuses) is expected


def test_aggregation_rejects_plain_strings() -> None:
    with pytest.raises(InputValidationError, match="QCStatus"):
        aggregate_qc_status(("pass",))  # type: ignore[arg-type]


def test_single_fastq_lightweight_qc_passes(tmp_path: Path) -> None:
    fastq = write_fastq(tmp_path)
    result = run_input_dataset_qc(InputDataset.from_files((fastq,)))
    assert result.status is QCStatus.PASS
    assert result.get_metric("input_type").value == "fastq"
    assert result.get_metric("file_count").value == 1
    assert result.get_metric("total_file_size_bytes").value == fastq.stat().st_size


def test_multiple_fastq_qc_preserves_order_and_total_size(tmp_path: Path) -> None:
    first = write_fastq(tmp_path, "movie2.fastq")
    second = write_fastq(tmp_path, "movie1.fastq")
    dataset = InputDataset.from_files((first, second))
    before = dataset.files
    result = run_input_dataset_qc(dataset)
    assert result.get_metric("file_count").value == 2
    assert result.get_metric("total_file_size_bytes").value == (
        first.stat().st_size + second.stat().st_size
    )
    assert result.get_metric("file_1_path").value == str(first.absolute())
    assert result.get_metric("file_2_path").value == str(second.absolute())
    assert dataset.files == before


def test_gzip_fastq_reports_gzip_compression(tmp_path: Path) -> None:
    dataset = InputDataset.from_files((write_gzip_fastq(tmp_path),))
    assert run_input_dataset_qc(dataset).get_metric("compression").value == "gzip"


def test_mixed_fastq_compression_is_explicit(tmp_path: Path) -> None:
    dataset = InputDataset.from_files(
        (write_fastq(tmp_path), write_gzip_fastq(tmp_path))
    )
    assert run_input_dataset_qc(dataset).get_metric("compression").value == "mixed"


@pytest.mark.parametrize("adjacent", (True, False))
def test_bam_with_readable_index_passes(tmp_path: Path, adjacent: bool) -> None:
    bam = write_alignment(tmp_path, "sample.bam")
    write_index(bam, adjacent=adjacent)
    result = run_input_dataset_qc(InputDataset.from_files((bam,)))
    assert result.status is QCStatus.PASS
    assert result.get_metric("index_present").value is True


def test_bam_without_index_warns(tmp_path: Path) -> None:
    bam = write_alignment(tmp_path, "sample.bam")
    result = run_input_dataset_qc(InputDataset.from_files((bam,)))
    assert result.status is QCStatus.WARN
    assert result.get_metric("index_present").value is False
    assert result.issues[0].code == "ALIGNMENT_INDEX_MISSING"


def test_cram_with_index_passes(tmp_path: Path) -> None:
    cram = write_alignment(tmp_path, "sample.cram")
    write_index(cram)
    result = run_input_dataset_qc(InputDataset.from_files((cram,)))
    assert result.status is QCStatus.PASS
    assert result.get_metric("index_present").value is True


def test_cram_without_index_warns(tmp_path: Path) -> None:
    cram = write_alignment(tmp_path, "sample.cram")
    result = run_input_dataset_qc(InputDataset.from_files((cram,)))
    assert result.status is QCStatus.WARN
    assert result.issues[0].level is QCIssueLevel.WARNING


def test_fastq_has_no_index_metric(tmp_path: Path) -> None:
    result = run_input_dataset_qc(
        InputDataset.from_files((write_fastq(tmp_path),))
    )
    with pytest.raises(KeyError):
        result.get_metric("index_present")


def test_disappeared_input_remains_validation_error(tmp_path: Path) -> None:
    fastq = write_fastq(tmp_path)
    dataset = InputDataset.from_files((fastq,))
    fastq.unlink()
    with pytest.raises(InputValidationError, match="missing"):
        run_input_dataset_qc(dataset)


def test_input_dataset_is_not_mutated(tmp_path: Path) -> None:
    dataset = InputDataset.from_files((write_fastq(tmp_path),))
    original = InputDataset(dataset.input_type, dataset.files)
    run_input_dataset_qc(dataset)
    assert dataset == original


def test_sample_id_is_retained(tmp_path: Path) -> None:
    dataset = InputDataset.from_files((write_fastq(tmp_path),))
    assert run_input_dataset_qc(dataset, sample_id="HG002").sample_id == "HG002"


def test_dataset_qc_rejects_wrong_object() -> None:
    with pytest.raises(InputValidationError, match="InputDataset"):
        run_input_dataset_qc(object())  # type: ignore[arg-type]


def test_analysis_context_single_sample_qc(tmp_path: Path) -> None:
    report = run_input_qc(make_context(tmp_path))
    assert len(report.sample_results) == 1
    assert report.sample_results[0].sample_id == "S1"


def test_analysis_context_multi_sample_order_is_preserved(tmp_path: Path) -> None:
    records = (
        SampleRecord(
            Sample("S2", InputDataset.from_files((write_fastq(tmp_path, "S2.fastq"),)))
        ),
        SampleRecord(
            Sample("S1", InputDataset.from_files((write_fastq(tmp_path, "S1.fastq"),)))
        ),
    )
    report = run_input_qc(make_context(tmp_path, records))
    assert tuple(result.sample_id for result in report.sample_results) == ("S2", "S1")


def test_mixed_context_qc_aggregates_warn(tmp_path: Path) -> None:
    fastq_record = SampleRecord(
        Sample("FASTQ1", InputDataset.from_files((write_fastq(tmp_path),)))
    )
    bam_record = SampleRecord(
        Sample("BAM1", InputDataset.from_files((write_alignment(tmp_path, "S.bam"),)))
    )
    report = run_input_qc(make_context(tmp_path, (fastq_record, bam_record)))
    assert [result.status for result in report.sample_results] == [
        QCStatus.PASS,
        QCStatus.WARN,
    ]
    assert report.overall_status is QCStatus.WARN


def test_run_report_records_reference_metadata(tmp_path: Path) -> None:
    report = run_input_qc(make_context(tmp_path))
    assert report.get_metric("reference_build").value == "GRCh38"
    assert report.get_metric("reference_contig_count").value == 1
    assert report.get_metric("reference_checksum_available").value is False


def test_run_report_status_counts_include_zeroes() -> None:
    report = RunQCReport((QCResult(QCStatus.FAIL),))
    assert report.status_counts == {
        "pass": 0,
        "warn": 0,
        "fail": 1,
        "not_checked": 0,
    }


def test_empty_run_report_is_not_checked() -> None:
    assert RunQCReport(()).overall_status is QCStatus.NOT_CHECKED


def test_run_report_json_serialization(tmp_path: Path) -> None:
    report = run_input_qc(make_context(tmp_path))
    path = report.write_json(tmp_path / "reports" / "qc.json")
    assert json.loads(path.read_text(encoding="utf-8")) == report.to_dict()


def test_run_report_yaml_serialization(tmp_path: Path) -> None:
    report = run_input_qc(make_context(tmp_path))
    path = report.write_yaml(tmp_path / "reports" / "qc.yaml")
    assert yaml.safe_load(path.read_text(encoding="utf-8")) == report.to_dict()


def test_unicode_input_path_is_preserved_in_reports(tmp_path: Path) -> None:
    unicode_root = tmp_path / "测序数据"
    context = make_context(unicode_root)
    report = run_input_qc(context)
    json_path = report.write_json(tmp_path / "结果" / "qc.json")
    yaml_path = report.write_yaml(tmp_path / "结果" / "qc.yaml")
    assert "测序数据" in json_path.read_text(encoding="utf-8")
    assert "测序数据" in yaml_path.read_text(encoding="utf-8")


@pytest.mark.parametrize("method,suffix", (("write_json", "json"), ("write_yaml", "yaml")))
def test_report_refuses_overwrite_by_default(
    tmp_path: Path, method: str, suffix: str
) -> None:
    report = run_input_qc(make_context(tmp_path))
    path = tmp_path / f"qc.{suffix}"
    path.write_text("original", encoding="utf-8")
    with pytest.raises(OutputValidationError, match="already exists"):
        getattr(report, method)(path)
    assert path.read_text(encoding="utf-8") == "original"


@pytest.mark.parametrize("method,suffix", (("write_json", "json"), ("write_yaml", "yaml")))
def test_report_overwrite_true_replaces_existing(
    tmp_path: Path, method: str, suffix: str
) -> None:
    report = run_input_qc(make_context(tmp_path))
    path = tmp_path / f"qc.{suffix}"
    path.write_text("original", encoding="utf-8")
    getattr(report, method)(path, overwrite=True)
    assert path.read_text(encoding="utf-8") != "original"


def test_report_records_version_and_utc_timestamp(tmp_path: Path) -> None:
    report = run_input_qc(make_context(tmp_path))
    assert report.hifivar_version == __version__
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}T.+Z", report.created_at)


def test_lightweight_qc_never_computes_checksum(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def unexpected_checksum(path: object) -> str:
        raise AssertionError(f"unexpected checksum for {path}")

    monkeypatch.setattr("hifivar.qc.validation.compute_sha256", unexpected_checksum)
    run_input_qc(make_context(tmp_path))


def test_report_payload_contains_only_standard_types(tmp_path: Path) -> None:
    payload = run_input_qc(make_context(tmp_path)).to_dict()
    json.dumps(payload)
    yaml.safe_dump(payload)


def test_qc_models_are_frozen() -> None:
    metric = QCMetric("file_count", 1)
    with pytest.raises(FrozenInstanceError):
        metric.value = 2  # type: ignore[misc]


def test_run_input_qc_rejects_wrong_object() -> None:
    with pytest.raises(InputValidationError, match="AnalysisContext"):
        run_input_qc(object())  # type: ignore[arg-type]
