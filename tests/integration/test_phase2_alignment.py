"""Phase 1 to Phase 2.2 alignment-planning integration test."""

from __future__ import annotations

import json
from pathlib import Path

import yaml

from hifivar.alignment import (
    AlignmentAction,
    AlignmentOutputFormat,
    AlignmentResources,
    AlignmentTool,
    build_alignment_plans,
    build_alignment_requests,
)
from hifivar.config import load_config
from hifivar.context import AnalysisContext


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = (
    PROJECT_ROOT / "src" / "hifivar" / "resources" / "configs" / "default.yaml"
)
STANDARD_PRESET = (
    PROJECT_ROOT
    / "src"
    / "hifivar"
    / "resources"
    / "configs"
    / "presets"
    / "standard.yaml"
)


def test_phase2_alignment_planning_end_to_end(tmp_path: Path) -> None:
    """Load a Unicode multi-sample context and create deterministic plans."""
    project = tmp_path / "科研项目" / "输入数据"
    project.mkdir(parents=True)
    fasta = project / "参考.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    for name in ("样本一.fastq", "样本二.fastq"):
        (project / name).write_text(
            "@read1\nACGT\n+\nIIII\n",
            encoding="utf-8",
        )
    sheet = project / "samples.tsv"
    sheet.write_text(
        "sample_id\tinput\nS2\t样本二.fastq\nS1\t样本一.fastq\n",
        encoding="utf-8",
    )
    user_config = project / "analysis.yaml"
    user_config.write_text(
        yaml.safe_dump(
            {
                "reference": {"fasta": "参考.fa", "build": "GRCh38"},
                "samples": {"sheet": "samples.tsv"},
            },
            allow_unicode=True,
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    effective = load_config(DEFAULT_CONFIG, STANDARD_PRESET, user_config)
    context = AnalysisContext.from_config(effective)
    output_root = tmp_path / "分析结果" / "alignment"
    requests = build_alignment_requests(
        context,
        output_root,
        tool=AlignmentTool.PBMM2,
        output_format=AlignmentOutputFormat.BAM,
        threads=24,
    )
    payload = [request.to_dict() for request in requests]

    assert tuple(request.sample.sample_id for request in requests) == ("S2", "S1")
    assert tuple(request.output_path.name for request in requests) == (
        "S2.aligned.bam",
        "S1.aligned.bam",
    )
    assert all(request.reference is context.reference for request in requests)
    assert all(request.threads == 24 for request in requests)
    assert not output_root.exists()
    assert "科研项目" in json.dumps(payload, ensure_ascii=False)
    assert yaml.safe_load(yaml.safe_dump(payload, allow_unicode=True)) == payload


def test_phase2_mixed_input_planning_skips_existing_alignment(
    tmp_path: Path,
) -> None:
    """FASTQ is planned for alignment while BAM is retained unchanged."""
    project = tmp_path / "inputs"
    project.mkdir()
    fasta = project / "reference.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    fastq = project / "reads.fastq"
    fastq.write_text("@read1\nACGT\n+\nIIII\n", encoding="utf-8")
    bam = project / "existing.bam"
    bam.write_bytes(b"alignment-placeholder")
    sheet = project / "samples.tsv"
    sheet.write_text(
        "sample_id\tinput\nFASTQ1\treads.fastq\nBAM1\texisting.bam\n",
        encoding="utf-8",
    )
    user_config = project / "analysis.yaml"
    user_config.write_text(
        yaml.safe_dump(
            {
                "reference": {"fasta": str(fasta), "build": "GRCh38"},
                "samples": {"sheet": str(sheet)},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    context = AnalysisContext.from_config(
        load_config(DEFAULT_CONFIG, STANDARD_PRESET, user_config)
    )

    plans = build_alignment_plans(
        context,
        tmp_path / "results",
        tool=AlignmentTool.PBMM2,
        resources=AlignmentResources(threads=8, memory_mb=16_000),
    )

    assert tuple(plan.action for plan in plans) == (
        AlignmentAction.ALIGN,
        AlignmentAction.REUSE,
    )
    assert plans[1].alignment_path == bam
    assert plans[1].request is None
