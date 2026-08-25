"""Phase 2.3 result to Phase 2.4 indexing/QC integration test."""

from __future__ import annotations

from pathlib import Path

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
    AlignmentIndexRequest,
    AlignmentSortOrder,
    run_alignment_qc,
    validate_alignment_artifact,
)
from hifivar.command import CommandResult
from hifivar.qc import QCStatus
from hifivar.reference import ReferenceGenome
from hifivar.sample import InputDataset, Sample
from hifivar.samtools import IndexResultStatus, SamtoolsWrapper


class FakeSamtoolsRunner:
    """Create a tiny index while preserving the CommandRunner contract."""

    def require_executable(self, executable: str) -> Path:
        return Path("/opt/samtools")

    def run(self, command: object, **kwargs: object) -> CommandResult:
        args = tuple(str(arg) for arg in command)  # type: ignore[union-attr]
        if args[1:] == ("--version",):
            return CommandResult(
                args,
                0,
                "samtools 1.22.1\n",
                "",
                0.01,
                None,
                True,
            )
        Path(args[-1]).write_bytes(b"tiny-index")
        return CommandResult(args, 0, "", "", 0.1, None, True)


def test_completed_pbmm2_output_is_indexed_then_passes_lightweight_qc(
    tmp_path: Path,
) -> None:
    reference_dir = tmp_path / "reference"
    reference_dir.mkdir()
    fasta = reference_dir / "reference.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    reference = ReferenceGenome.from_fasta(fasta, build="GRCh38")
    reads_dir = tmp_path / "reads"
    reads_dir.mkdir()
    fastq = reads_dir / "sample.fastq"
    fastq.write_text("@read1\nACGT\n+\nIIII\n", encoding="utf-8")
    sample = Sample("HG002", InputDataset.from_files((fastq,)))
    output = tmp_path / "results" / "HG002.aligned.bam"
    output.parent.mkdir()
    output.write_bytes(b"tiny-coordinate-sorted-bam")
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
    result = AlignmentResult(
        plan=plan,
        status=AlignmentResultStatus.COMPLETED,
        command=AlignmentCommandPlan(
            AlignmentTool.PBMM2,
            ("pbmm2", "align"),
        ),
        tool_version="1.17.0",
        duration_seconds=1.0,
    )

    artifact = AlignmentArtifact.from_result(result)
    assert artifact.sort_order is AlignmentSortOrder.COORDINATE
    index_request = AlignmentIndexRequest.create(artifact, threads=2)
    index_result = SamtoolsWrapper(
        runner=FakeSamtoolsRunner(),  # type: ignore[arg-type]
    ).run_index(index_request)
    validated = validate_alignment_artifact(
        index_result.artifact,
        require_index=True,
    )
    qc = run_alignment_qc(validated)

    assert index_result.status is IndexResultStatus.COMPLETED
    assert validated.index_path == index_request.output_path
    assert qc.status is QCStatus.PASS
    assert qc.sample_id == "HG002"
