"""Phase 1/2.2 to Phase 2.3 pbmm2 dry-run integration test."""

from __future__ import annotations

from pathlib import Path

import yaml

from hifivar.alignment import (
    AlignmentAction,
    AlignmentResources,
    AlignmentResultStatus,
    AlignmentTool,
    build_alignment_plans,
)
from hifivar.command import CommandResult
from hifivar.config import load_config
from hifivar.context import AnalysisContext
from hifivar.pbmm2 import Pbmm2Options, Pbmm2Wrapper


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


class DryRunOnlyRunner:
    """Prove that dry-run does not perform executable detection."""

    def __init__(self) -> None:
        self.command: tuple[str, ...] | None = None

    def require_executable(self, executable: str) -> Path:
        raise AssertionError(f"unexpected executable check for {executable}")

    def run(self, command: object, **kwargs: object) -> CommandResult:
        assert kwargs["dry_run"] is True
        self.command = tuple(str(arg) for arg in command)  # type: ignore[union-attr]
        return CommandResult(
            args=self.command,
            returncode=None,
            stdout=None,
            stderr=None,
            duration_seconds=0.0,
            cwd=None,
            executed=False,
        )


def test_phase2_pbmm2_dry_run_from_effective_config(tmp_path: Path) -> None:
    """Load real config/context models and preview a sorted pbmm2 command."""
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    fasta = inputs / "reference.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    fastq = inputs / "sample.fastq"
    fastq.write_text("@read1\nACGT\n+\nIIII\n", encoding="utf-8")
    sheet = inputs / "samples.tsv"
    sheet.write_text("sample_id\tinput\nHG002\tsample.fastq\n", encoding="utf-8")
    user_config = inputs / "analysis.yaml"
    user_config.write_text(
        yaml.safe_dump(
            {
                "reference": {"fasta": "reference.fa", "build": "GRCh38"},
                "samples": {"sheet": "samples.tsv"},
                "alignment": {"threads": 12, "pbmm2_log_level": "DEBUG"},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    effective = load_config(DEFAULT_CONFIG, STANDARD_PRESET, user_config)
    context = AnalysisContext.from_config(effective)
    alignment_config = effective["alignment"]
    assert isinstance(alignment_config, dict)
    resources = AlignmentResources(
        threads=alignment_config["threads"],  # type: ignore[arg-type]
        memory_mb=alignment_config["memory_mb"],  # type: ignore[arg-type]
        runtime_minutes=alignment_config["runtime_minutes"],  # type: ignore[arg-type]
    )
    plan = build_alignment_plans(
        context,
        tmp_path / "results" / "alignment",
        tool=AlignmentTool.PBMM2,
        resources=resources,
    )[0]
    assert plan.action is AlignmentAction.ALIGN
    assert plan.request is not None

    runner = DryRunOnlyRunner()
    result = Pbmm2Wrapper(
        runner=runner,  # type: ignore[arg-type]
        options=Pbmm2Options.from_config(effective),
    ).run(plan.request, dry_run=True)

    assert result.status is AlignmentResultStatus.PLANNED
    assert runner.command is not None
    assert runner.command[0:2] == ("pbmm2", "align")
    assert runner.command[runner.command.index("-j") + 1] == "12"
    assert runner.command[runner.command.index("--log-level") + 1] == "DEBUG"
    assert not plan.alignment_path.exists()
