"""Tests for Phase 2.4 alignment artifacts, indexing plans, and QC."""

from __future__ import annotations

from pathlib import Path

import pytest

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
from hifivar.alignment_postprocess import (
    AlignmentArtifact,
    AlignmentIndexFormat,
    AlignmentIndexRequest,
    AlignmentSortOrder,
    AlignmentSource,
    choose_index_format,
    find_alignment_index,
    index_path_for,
    run_alignment_qc,
    validate_alignment_artifact,
)
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.qc import QCStatus
from hifivar.reference import ReferenceGenome
from hifivar.sample import InputDataset, Sample


def write_reference(root: Path, *, contig_length: int = 4) -> ReferenceGenome:
    root.mkdir(parents=True, exist_ok=True)
    fasta = root / "reference.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text(
        f"chr1\t{contig_length}\t6\t4\t5\n",
        encoding="utf-8",
    )
    return ReferenceGenome.from_fasta(fasta, build="GRCh38")


def write_fastq(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    path = root / "S1.fastq"
    path.write_text("@read1\nACGT\n+\nIIII\n", encoding="utf-8")
    return path


def write_alignment(root: Path, name: str = "S1.bam") -> Path:
    root.mkdir(parents=True, exist_ok=True)
    path = root / name
    path.write_bytes(b"alignment-placeholder")
    return path


def make_artifact(
    tmp_path: Path,
    *,
    name: str = "S1.bam",
    source: AlignmentSource = AlignmentSource.GENERATED,
    sort_order: AlignmentSortOrder = AlignmentSortOrder.COORDINATE,
    reference: ReferenceGenome | None = None,
    index_path: Path | None = None,
) -> AlignmentArtifact:
    path = write_alignment(tmp_path / "alignment", name)
    output_format = (
        AlignmentOutputFormat.CRAM
        if path.suffix.lower() == ".cram"
        else AlignmentOutputFormat.BAM
    )
    return AlignmentArtifact(
        sample_id="S1",
        path=path,
        output_format=output_format,
        reference=reference or write_reference(tmp_path / "reference"),
        source=source,
        sort_order=sort_order,
        index_path=index_path,
        tool=AlignmentTool.PBMM2 if source is AlignmentSource.GENERATED else None,
        tool_version="1.17.0" if source is AlignmentSource.GENERATED else None,
    )


def make_completed_alignment_result(tmp_path: Path) -> AlignmentResult:
    reference = write_reference(tmp_path / "reference")
    sample = Sample("S1", InputDataset.from_files((write_fastq(tmp_path / "reads"),)))
    output = write_alignment(tmp_path / "alignment")
    request = AlignmentRequest(
        sample=sample,
        reference=reference,
        output_path=output,
        tool=AlignmentTool.PBMM2,
        overwrite=True,
    )
    plan = AlignmentPlan(
        sample=sample,
        reference=reference,
        action=AlignmentAction.ALIGN,
        alignment_path=output,
        output_format=AlignmentOutputFormat.BAM,
        request=request,
    )
    return AlignmentResult(
        plan=plan,
        status=AlignmentResultStatus.COMPLETED,
        command=AlignmentCommandPlan(
            AlignmentTool.PBMM2,
            ("pbmm2", "align"),
        ),
        tool_version="1.17.0",
        duration_seconds=5.0,
    )


def test_completed_pbmm2_result_becomes_coordinate_sorted_artifact(
    tmp_path: Path,
) -> None:
    artifact = AlignmentArtifact.from_result(make_completed_alignment_result(tmp_path))
    assert artifact.source is AlignmentSource.GENERATED
    assert artifact.sort_order is AlignmentSortOrder.COORDINATE
    assert artifact.tool is AlignmentTool.PBMM2
    assert artifact.index_path is None


def test_reused_alignment_remains_unknown_sort_and_discovers_index(
    tmp_path: Path,
) -> None:
    reference = write_reference(tmp_path / "reference")
    bam = write_alignment(tmp_path / "alignment", "existing.bam")
    index = Path(f"{bam}.bai")
    index.write_bytes(b"index")
    sample = Sample("B1", InputDataset.from_files((bam,)))
    plan = AlignmentPlan(
        sample=sample,
        reference=reference,
        action=AlignmentAction.REUSE,
        alignment_path=bam,
        output_format=AlignmentOutputFormat.BAM,
    )
    artifact = AlignmentArtifact.from_result(
        AlignmentResult(plan, AlignmentResultStatus.REUSED)
    )
    assert artifact.source is AlignmentSource.EXISTING
    assert artifact.sort_order is AlignmentSortOrder.UNKNOWN
    assert artifact.index_path == index


def test_dry_run_result_cannot_claim_output_artifact(tmp_path: Path) -> None:
    completed = make_completed_alignment_result(tmp_path)
    planned = AlignmentResult(
        completed.plan,
        AlignmentResultStatus.PLANNED,
        completed.command,
    )
    with pytest.raises(OutputValidationError, match="dry-run"):
        AlignmentArtifact.from_result(planned)


def test_default_index_strategy_uses_bai_csi_and_crai(tmp_path: Path) -> None:
    normal = make_artifact(tmp_path / "normal")
    huge = make_artifact(
        tmp_path / "huge",
        reference=write_reference(
            tmp_path / "huge-reference",
            contig_length=2**29 + 1,
        ),
    )
    cram = make_artifact(tmp_path / "cram", name="S1.cram")
    assert choose_index_format(normal) is AlignmentIndexFormat.BAI
    assert choose_index_format(huge) is AlignmentIndexFormat.CSI
    assert choose_index_format(cram) is AlignmentIndexFormat.CRAI


def test_index_request_uses_conventional_adjacent_path(tmp_path: Path) -> None:
    artifact = make_artifact(tmp_path)
    request = AlignmentIndexRequest.create(artifact, threads=4)
    assert request.index_format is AlignmentIndexFormat.BAI
    assert request.output_path == Path(f"{artifact.path}.bai")
    assert request.to_dict()["threads"] == 4


def test_index_request_never_assumes_unknown_sort_is_safe(tmp_path: Path) -> None:
    artifact = make_artifact(
        tmp_path,
        source=AlignmentSource.EXISTING,
        sort_order=AlignmentSortOrder.UNKNOWN,
    )
    with pytest.raises(InputValidationError, match="UNKNOWN sort order"):
        AlignmentIndexRequest.create(artifact)


@pytest.mark.parametrize(
    "name,index_format,message",
    (
        ("S1.bam", AlignmentIndexFormat.CRAI, "BAM"),
        ("S1.cram", AlignmentIndexFormat.BAI, "CRAM"),
        ("S1.cram", AlignmentIndexFormat.CSI, "CRAM"),
    ),
)
def test_index_request_rejects_incompatible_formats(
    tmp_path: Path,
    name: str,
    index_format: AlignmentIndexFormat,
    message: str,
) -> None:
    with pytest.raises(InputValidationError, match=message):
        AlignmentIndexRequest.create(
            make_artifact(tmp_path, name=name),
            index_format=index_format,
        )


def test_existing_index_is_not_silently_replaced(tmp_path: Path) -> None:
    artifact = make_artifact(tmp_path)
    output = index_path_for(artifact.path, AlignmentIndexFormat.BAI)
    output.write_bytes(b"existing")
    with pytest.raises(OutputValidationError, match="already exists"):
        AlignmentIndexRequest.create(artifact)
    request = AlignmentIndexRequest.create(artifact, overwrite=True)
    assert request.overwrite is True
    assert output.read_bytes() == b"existing"


def test_find_alignment_index_supports_adjacent_and_replaced_names(
    tmp_path: Path,
) -> None:
    bam = write_alignment(tmp_path)
    replaced = bam.with_suffix(".bai")
    replaced.write_bytes(b"index")
    assert find_alignment_index(bam) == replaced


def test_validate_artifact_requires_index_only_when_requested(tmp_path: Path) -> None:
    artifact = make_artifact(tmp_path)
    assert validate_alignment_artifact(artifact) == artifact
    with pytest.raises(OutputValidationError, match="index missing"):
        validate_alignment_artifact(artifact, require_index=True)


def test_lightweight_qc_warns_for_existing_unknown_unindexed_alignment(
    tmp_path: Path,
) -> None:
    artifact = make_artifact(
        tmp_path,
        source=AlignmentSource.EXISTING,
        sort_order=AlignmentSortOrder.UNKNOWN,
    )
    result = run_alignment_qc(artifact)
    assert result.status is QCStatus.WARN
    assert tuple(issue.code for issue in result.issues) == (
        "ALIGNMENT_SORT_ORDER_UNKNOWN",
        "ALIGNMENT_INDEX_MISSING",
    )
    assert result.get_metric("index_present").value is False


def test_lightweight_qc_passes_coordinate_sorted_indexed_placeholder(
    tmp_path: Path,
) -> None:
    artifact = make_artifact(tmp_path)
    index = Path(f"{artifact.path}.bai")
    index.write_bytes(b"index-placeholder")
    result = run_alignment_qc(artifact.with_index(index))
    assert result.status is QCStatus.PASS
    assert result.issues == ()
    assert result.get_metric("sort_order").value == "coordinate"


def test_lightweight_qc_missing_alignment_is_validation_failure(
    tmp_path: Path,
) -> None:
    artifact = make_artifact(tmp_path)
    artifact.path.unlink()
    with pytest.raises(OutputValidationError, match="missing"):
        run_alignment_qc(artifact)
