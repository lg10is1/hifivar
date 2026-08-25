"""Unit tests for Phase 3 orchestration and provenance."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import (
    AlignmentArtifact,
    AlignmentSortOrder,
    AlignmentSource,
)
from hifivar.command import CommandResult
from hifivar.context import AnalysisContext
from hifivar.deepvariant import DeepVariantWrapper
from hifivar.exceptions import InputValidationError, ReferenceError
from hifivar.phase3 import Phase3RunStatus, Phase3Settings, run_phase3
from hifivar.reference import ReferenceGenome
from hifivar.sample import InputDataset, Sample
from hifivar.small import SmallVariantResultStatus


class DryRunner:
    def __init__(self) -> None:
        self.commands: list[tuple[str, ...]] = []

    def require_executable(self, executable: str) -> Path:
        raise AssertionError("dry-run checked an executable")

    def run(self, command: object, **kwargs: object) -> CommandResult:
        args = tuple(str(value) for value in command)  # type: ignore[union-attr]
        self.commands.append(args)
        assert kwargs.get("dry_run") is True
        return CommandResult(args, None, None, None, 0.0, None, False)


def make_context(tmp_path: Path, *, fastq: bool = False) -> AnalysisContext:
    tmp_path.mkdir(parents=True, exist_ok=True)
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    if fastq:
        primary = tmp_path / "S1.fastq"
        primary.write_text("@r\nACG\n+\nIII\n", encoding="utf-8")
    else:
        primary = tmp_path / "S1.bam"
        primary.write_bytes(b"bam")
        Path(f"{primary}.bai").write_bytes(b"bai")
    config = {
        "reference": {"fasta": str(fasta), "build": "GRCh38"},
        "small": {
            "execution_mode": "native",
            "deepvariant_executable": "run_deepvariant",
            "deepvariant_image": None,
            "model_type": "PACBIO",
            "threads": 6,
            "memory_mb": 24000,
            "runtime_minutes": 300,
            "overwrite": False,
        },
    }
    return AnalysisContext.from_sample(
        ReferenceGenome.from_fasta(fasta, build="GRCh38"),
        Sample("S1", InputDataset.from_files((primary,))),
        config,
    )


def test_settings_load_resources_and_runtime(tmp_path: Path) -> None:
    settings = Phase3Settings.from_config(make_context(tmp_path).config)
    assert settings.resources.threads == 6
    assert settings.resources.memory_mb == 24000
    assert settings.runtime.mode.value == "native"
    assert settings.to_dict()["model_type"] == "PACBIO"


def test_existing_alignment_dry_run_is_ordered_and_serializable(tmp_path: Path) -> None:
    context = make_context(tmp_path)
    runner = DryRunner()
    report = run_phase3(
        context,
        tmp_path / "results",
        deepvariant_wrapper=DeepVariantWrapper(runner=runner),  # type: ignore[arg-type]
        dry_run=True,
    )
    assert report.status is Phase3RunStatus.PLANNED
    assert report.sample_results[0].call.status is SmallVariantResultStatus.PLANNED
    assert report.tool_versions == {"deepvariant": None}
    assert "S1.small.vcf.gz" in " ".join(runner.commands[0])
    assert not (tmp_path / "results").exists()
    json_path = report.write_json(tmp_path / "reports" / "phase3.json")
    yaml_path = report.write_yaml(tmp_path / "reports" / "phase3.yaml")
    assert json.loads(json_path.read_text(encoding="utf-8")) == yaml.safe_load(
        yaml_path.read_text(encoding="utf-8")
    )


def test_raw_fastq_requires_completed_alignment_handoff(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="raw FASTQ"):
        run_phase3(make_context(tmp_path, fastq=True), tmp_path / "results", dry_run=True)


def test_supplied_artifact_must_match_context_reference(tmp_path: Path) -> None:
    context = make_context(tmp_path / "context", fastq=True)
    other = tmp_path / "other"
    other.mkdir()
    fasta = other / "other.fa"
    fasta.write_text(">chr2\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr2\t4\t6\t4\t5\n", encoding="utf-8")
    bam = other / "S1.bam"
    bam.write_bytes(b"bam")
    bai = Path(f"{bam}.bai")
    bai.write_bytes(b"bai")
    artifact = AlignmentArtifact(
        sample_id="S1",
        path=bam,
        output_format=AlignmentOutputFormat.BAM,
        reference=ReferenceGenome.from_fasta(fasta, build="GRCh38"),
        source=AlignmentSource.EXISTING,
        sort_order=AlignmentSortOrder.COORDINATE,
        index_path=bai,
    )
    with pytest.raises((InputValidationError, ReferenceError)):
        run_phase3(
            context,
            tmp_path / "results",
            alignment_artifacts={"S1": artifact},
            dry_run=True,
        )
