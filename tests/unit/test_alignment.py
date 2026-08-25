"""Tests for the Phase 2.2 tool-neutral alignment interface."""

from __future__ import annotations

import json
from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest
import yaml

from hifivar.alignment import (
    AlignmentAction,
    AlignmentBackend,
    AlignmentCommandPlan,
    AlignmentOutputFormat,
    AlignmentPlan,
    AlignmentRequest,
    AlignmentResources,
    AlignmentResult,
    AlignmentResultStatus,
    AlignmentTool,
    build_alignment_plans,
    build_alignment_requests,
    plan_alignment_command,
)
from hifivar.context import AnalysisContext
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.reference import ReferenceGenome
from hifivar.sample import InputDataset, Sample
from hifivar.sample_sheet import SampleRecord


def write_reference(root: Path) -> ReferenceGenome:
    root.mkdir(parents=True, exist_ok=True)
    fasta = root / "reference.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    return ReferenceGenome.from_fasta(fasta, build="GRCh38")


def write_fastq(root: Path, name: str = "sample.fastq") -> Path:
    root.mkdir(parents=True, exist_ok=True)
    path = root / name
    path.write_text("@read1\nACGT\n+\nIIII\n", encoding="utf-8")
    return path


def write_alignment(root: Path, name: str = "sample.bam") -> Path:
    root.mkdir(parents=True, exist_ok=True)
    path = root / name
    path.write_bytes(b"alignment-placeholder")
    return path


def make_sample(root: Path, sample_id: str = "S1") -> Sample:
    return Sample(
        sample_id,
        InputDataset.from_files((write_fastq(root, f"{sample_id}.fastq"),)),
    )


def make_context(
    root: Path,
    sample_ids: tuple[str, ...] = ("S1",),
) -> AnalysisContext:
    reference = write_reference(root / "reference")
    records = tuple(
        SampleRecord(make_sample(root / "reads", sample_id))
        for sample_id in sample_ids
    )
    return AnalysisContext(
        reference,
        records,
        {
            "reference": {
                "fasta": str(reference.fasta.absolute()),
                "build": reference.build,
            }
        },
    )


def make_request(tmp_path: Path, **overrides: object) -> AlignmentRequest:
    sample = make_sample(tmp_path / "reads")
    reference = write_reference(tmp_path / "reference")
    values: dict[str, object] = {
        "sample": sample,
        "reference": reference,
        "output_path": tmp_path / "results" / "S1.aligned.bam",
        "tool": AlignmentTool.PBMM2,
    }
    values.update(overrides)
    return AlignmentRequest(**values)  # type: ignore[arg-type]


def test_alignment_tool_values_are_stable() -> None:
    assert [tool.value for tool in AlignmentTool] == ["pbmm2", "minimap2"]


def test_alignment_output_formats_and_suffixes_are_stable() -> None:
    assert [(item.value, item.suffix) for item in AlignmentOutputFormat] == [
        ("bam", ".bam"),
        ("cram", ".cram"),
    ]


def test_valid_alignment_request_preserves_ordered_fastq_inputs(
    tmp_path: Path,
) -> None:
    first = write_fastq(tmp_path / "reads", "movie2.fastq")
    second = write_fastq(tmp_path / "reads", "movie1.fastq")
    sample = Sample("S1", InputDataset.from_files((first, second)))
    request = make_request(tmp_path, sample=sample, threads=8)

    assert request.input_paths == (first, second)
    assert request.threads == 8
    assert request.output_path == tmp_path / "results" / "S1.aligned.bam"


def test_cram_request_accepts_matching_case_insensitive_suffix(tmp_path: Path) -> None:
    request = make_request(
        tmp_path,
        output_path=tmp_path / "S1.CRAM",
        output_format=AlignmentOutputFormat.CRAM,
    )
    assert request.output_format is AlignmentOutputFormat.CRAM


def test_request_serialization_contains_only_standard_types(tmp_path: Path) -> None:
    request = make_request(tmp_path, tool=AlignmentTool.MINIMAP2, threads=16)
    payload = request.to_dict()

    assert payload["sample_id"] == "S1"
    assert payload["tool"] == "minimap2"
    assert payload["output"] == {
        "path": str(tmp_path / "results" / "S1.aligned.bam"),
        "format": "bam",
    }
    assert payload["threads"] == 16
    assert payload["resources"] == {
        "threads": 16,
        "memory_mb": None,
        "runtime_minutes": None,
    }
    assert payload["overwrite"] is False
    json.dumps(payload)
    yaml.safe_dump(payload)


@pytest.mark.parametrize("threads", (0, -1, True, 1.5, "8"))
def test_request_rejects_invalid_thread_count(
    tmp_path: Path, threads: object
) -> None:
    with pytest.raises(InputValidationError, match="positive integer"):
        make_request(tmp_path, threads=threads)


def test_request_rejects_bam_primary_input(tmp_path: Path) -> None:
    bam = write_alignment(tmp_path / "reads")
    sample = Sample("S1", InputDataset.from_files((bam,)))
    with pytest.raises(InputValidationError, match="requires FASTQ input"):
        make_request(tmp_path, sample=sample)


@pytest.mark.parametrize(
    "output_path,output_format",
    (
        ("S1.cram", AlignmentOutputFormat.BAM),
        ("S1.bam", AlignmentOutputFormat.CRAM),
        ("S1.sam", AlignmentOutputFormat.BAM),
    ),
)
def test_request_rejects_output_suffix_mismatch(
    tmp_path: Path,
    output_path: str,
    output_format: AlignmentOutputFormat,
) -> None:
    with pytest.raises(InputValidationError, match="expected suffix"):
        make_request(
            tmp_path,
            output_path=tmp_path / output_path,
            output_format=output_format,
        )


@pytest.mark.parametrize("output_path", ("", "   "))
def test_request_rejects_empty_output_path(
    tmp_path: Path, output_path: str
) -> None:
    with pytest.raises(InputValidationError, match="must not be empty"):
        make_request(tmp_path, output_path=output_path)


def test_request_rejects_output_overwriting_reference(tmp_path: Path) -> None:
    reference = write_reference(tmp_path / "reference")
    reference_named_bam = reference.fasta.with_suffix(".bam")
    fake_reference = ReferenceGenome(
        fasta=reference_named_bam,
        fai=reference.fai,
        build=reference.build,
        contigs=reference.contigs,
    )
    with pytest.raises(InputValidationError, match="conflicts"):
        make_request(
            tmp_path,
            reference=fake_reference,
            output_path=reference_named_bam,
        )


def test_request_rejects_plain_string_enums(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="AlignmentTool"):
        make_request(tmp_path, tool="pbmm2")
    with pytest.raises(InputValidationError, match="AlignmentOutputFormat"):
        make_request(tmp_path, output_format="bam")


def test_request_rejects_wrong_model_objects(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="Sample instance"):
        make_request(tmp_path, sample=object())
    with pytest.raises(InputValidationError, match="ReferenceGenome"):
        make_request(tmp_path, reference=object())


def test_request_is_frozen(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    with pytest.raises(FrozenInstanceError):
        request.threads = 2  # type: ignore[misc]


def test_build_requests_preserves_context_order_and_shared_reference(
    tmp_path: Path,
) -> None:
    context = make_context(tmp_path, ("S2", "S1"))
    output_root = tmp_path / "results"
    requests = build_alignment_requests(
        context,
        output_root,
        tool=AlignmentTool.PBMM2,
        threads=12,
    )

    assert tuple(item.sample.sample_id for item in requests) == ("S2", "S1")
    assert tuple(item.output_path.name for item in requests) == (
        "S2.aligned.bam",
        "S1.aligned.bam",
    )
    assert all(item.reference is context.reference for item in requests)
    assert all(item.threads == 12 for item in requests)
    assert not output_root.exists()


def test_build_requests_supports_cram_output_and_unicode_path(tmp_path: Path) -> None:
    context = make_context(tmp_path)
    output_root = tmp_path / "分析结果" / "比对"
    request = build_alignment_requests(
        context,
        output_root,
        tool=AlignmentTool.MINIMAP2,
        output_format=AlignmentOutputFormat.CRAM,
    )[0]

    assert request.output_path == output_root / "S1.aligned.cram"
    assert request.tool is AlignmentTool.MINIMAP2
    assert not output_root.exists()


def test_build_requests_rejects_mixed_input_context_atomically(
    tmp_path: Path,
) -> None:
    context = make_context(tmp_path)
    bam = write_alignment(tmp_path / "reads", "S2.bam")
    mixed = AnalysisContext(
        context.reference,
        context.samples
        + (SampleRecord(Sample("S2", InputDataset.from_files((bam,)))),),
        context.config,
    )
    output_root = tmp_path / "results"

    with pytest.raises(InputValidationError, match=r"S2 \(bam\)"):
        build_alignment_requests(
            mixed,
            output_root,
            tool=AlignmentTool.PBMM2,
        )
    assert not output_root.exists()


def test_build_requests_rejects_wrong_context_and_output_path(
    tmp_path: Path,
) -> None:
    with pytest.raises(InputValidationError, match="AnalysisContext"):
        build_alignment_requests(
            object(),  # type: ignore[arg-type]
            tmp_path,
            tool=AlignmentTool.PBMM2,
        )
    with pytest.raises(InputValidationError, match="must not be empty"):
        build_alignment_requests(
            make_context(tmp_path),
            "",
            tool=AlignmentTool.PBMM2,
        )


def test_build_requests_rejects_plain_string_output_format(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="AlignmentOutputFormat"):
        build_alignment_requests(
            make_context(tmp_path),
            tmp_path / "results",
            tool=AlignmentTool.PBMM2,
            output_format="bam",  # type: ignore[arg-type]
        )


def test_backend_protocol_exposes_tool_and_shell_free_command(
    tmp_path: Path,
) -> None:
    class FakeBackend:
        tool = AlignmentTool.PBMM2

        def build_command(self, request: AlignmentRequest) -> list[str]:
            return [self.tool.value, "align", str(request.reference.fasta)]

    backend: AlignmentBackend = FakeBackend()
    command = backend.build_command(make_request(tmp_path))
    assert backend.tool is AlignmentTool.PBMM2
    assert command == [
        "pbmm2",
        "align",
        str(tmp_path / "reference" / "reference.fa"),
    ]


def test_alignment_resources_are_explicit_and_serializable() -> None:
    resources = AlignmentResources(
        threads=24,
        memory_mb=48_000,
        runtime_minutes=720,
    )
    assert resources.to_dict() == {
        "threads": 24,
        "memory_mb": 48_000,
        "runtime_minutes": 720,
    }


@pytest.mark.parametrize(
    "kwargs",
    (
        {"threads": 0},
        {"threads": True},
        {"memory_mb": 0},
        {"runtime_minutes": -1},
    ),
)
def test_alignment_resources_reject_invalid_values(
    kwargs: dict[str, object],
) -> None:
    with pytest.raises(InputValidationError, match="positive integer"):
        AlignmentResources(**kwargs)  # type: ignore[arg-type]


def test_request_uses_canonical_resources(tmp_path: Path) -> None:
    resources = AlignmentResources(threads=12, memory_mb=32_000)
    request = make_request(tmp_path, resources=resources)
    assert request.resources is resources
    assert request.threads == 12


def test_request_rejects_conflicting_thread_shorthand(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="conflicts"):
        make_request(
            tmp_path,
            resources=AlignmentResources(threads=8),
            threads=4,
        )


def test_request_refuses_existing_output_without_explicit_overwrite(
    tmp_path: Path,
) -> None:
    output = tmp_path / "results" / "S1.aligned.bam"
    output.parent.mkdir(parents=True)
    output.write_bytes(b"existing")
    with pytest.raises(OutputValidationError, match="already exists"):
        make_request(tmp_path, output_path=output)

    request = make_request(tmp_path, output_path=output, overwrite=True)
    assert request.overwrite is True
    assert output.read_bytes() == b"existing"


def test_request_rejects_directory_as_output_even_with_overwrite(
    tmp_path: Path,
) -> None:
    output = tmp_path / "S1.bam"
    output.mkdir()
    with pytest.raises(OutputValidationError, match="directory"):
        make_request(tmp_path, output_path=output, overwrite=True)


def test_mixed_context_plans_fastq_and_reuses_existing_alignments(
    tmp_path: Path,
) -> None:
    context = make_context(tmp_path)
    bam = write_alignment(tmp_path / "reads", "B1.bam")
    cram = write_alignment(tmp_path / "reads", "C1.cram")
    mixed = AnalysisContext(
        context.reference,
        context.samples
        + (
            SampleRecord(Sample("B1", InputDataset.from_files((bam,)))),
            SampleRecord(Sample("C1", InputDataset.from_files((cram,)))),
        ),
        context.config,
    )
    resources = AlignmentResources(threads=16, memory_mb=64_000)
    output_root = tmp_path / "planned"

    plans = build_alignment_plans(
        mixed,
        output_root,
        tool=AlignmentTool.PBMM2,
        resources=resources,
    )

    assert tuple(plan.action for plan in plans) == (
        AlignmentAction.ALIGN,
        AlignmentAction.REUSE,
        AlignmentAction.REUSE,
    )
    assert plans[0].alignment_path == output_root / "S1.aligned.bam"
    assert plans[0].resources is resources
    assert plans[1].alignment_path == bam
    assert plans[1].output_format is AlignmentOutputFormat.BAM
    assert plans[1].resources is None
    assert plans[2].alignment_path == cram
    assert plans[2].output_format is AlignmentOutputFormat.CRAM
    assert not output_root.exists()


def test_alignment_plan_payload_makes_skip_behavior_explicit(tmp_path: Path) -> None:
    bam = write_alignment(tmp_path, "existing.bam")
    sample = Sample("B1", InputDataset.from_files((bam,)))
    reference = write_reference(tmp_path / "reference")
    plan = AlignmentPlan(
        sample=sample,
        reference=reference,
        action=AlignmentAction.REUSE,
        alignment_path=bam,
        output_format=AlignmentOutputFormat.BAM,
    )

    assert plan.requires_alignment is False
    assert plan.to_dict()["action"] == "reuse"
    assert plan.to_dict()["request"] is None


def test_generic_command_plan_is_shell_free_and_reproducible(tmp_path: Path) -> None:
    request = make_request(tmp_path)

    class FakeBackend:
        tool = AlignmentTool.PBMM2

        def build_command(self, request: AlignmentRequest) -> list[str]:
            return ["pbmm2", "align", str(request.reference.fasta)]

    command = plan_alignment_command(FakeBackend(), request)
    assert command.args[0] == "pbmm2"
    assert command.to_dict()["shell"] is False
    assert "pbmm2 align" in command.display


def test_generic_command_planner_rejects_backend_tool_mismatch(
    tmp_path: Path,
) -> None:
    request = make_request(tmp_path)

    class WrongBackend:
        tool = AlignmentTool.MINIMAP2

        def build_command(self, request: AlignmentRequest) -> list[str]:
            return ["minimap2"]

    with pytest.raises(InputValidationError, match="cannot handle"):
        plan_alignment_command(WrongBackend(), request)


def test_alignment_result_models_planned_and_reused_paths(tmp_path: Path) -> None:
    context = make_context(tmp_path)
    align_plan = build_alignment_plans(
        context,
        tmp_path / "results",
        tool=AlignmentTool.PBMM2,
    )[0]
    command = AlignmentCommandPlan(AlignmentTool.PBMM2, ("pbmm2", "align"))
    planned = AlignmentResult(
        plan=align_plan,
        status=AlignmentResultStatus.PLANNED,
        command=command,
    )
    assert planned.executed is False
    assert planned.to_dict()["status"] == "planned"

    bam = write_alignment(tmp_path, "existing.bam")
    sample = Sample("B1", InputDataset.from_files((bam,)))
    reuse_plan = AlignmentPlan(
        sample=sample,
        reference=context.reference,
        action=AlignmentAction.REUSE,
        alignment_path=bam,
        output_format=AlignmentOutputFormat.BAM,
    )
    reused = AlignmentResult(reuse_plan, AlignmentResultStatus.REUSED)
    assert reused.alignment_path == bam
    assert reused.to_dict()["executed"] is False
