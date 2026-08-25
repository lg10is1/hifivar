"""Tiny mock integration of the complete Phase 2 path."""

from __future__ import annotations

import json
from pathlib import Path

import yaml

from hifivar.command import CommandResult
from hifivar.config import load_config
from hifivar.context import AnalysisContext
from hifivar.pbmm2 import Pbmm2Options, Pbmm2Wrapper
from hifivar.phase2 import Phase2Settings, run_phase2
from hifivar.qc import QCStatus
from hifivar.samtools import SamtoolsWrapper


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


class FakeBioinformaticsRunner:
    """Mock pbmm2/samtools execution while creating tiny expected outputs."""

    def __init__(self, *, allow_execution: bool = True) -> None:
        self.allow_execution = allow_execution
        self.commands: list[tuple[str, ...]] = []
        self.require_calls: list[str] = []

    def require_executable(self, executable: str) -> Path:
        if not self.allow_execution:
            raise AssertionError(f"dry-run checked executable {executable}")
        self.require_calls.append(executable)
        return Path(f"/opt/{executable}")

    def run(self, command: object, **kwargs: object) -> CommandResult:
        args = tuple(str(arg) for arg in command)  # type: ignore[union-attr]
        self.commands.append(args)
        if kwargs.get("dry_run") is True:
            return CommandResult(args, None, None, None, 0.0, None, False)
        if args[1:] == ("--version",):
            version = "pbmm2 1.17.0\n" if args[0] == "pbmm2" else "samtools 1.22.1\n"
            return CommandResult(args, 0, version, "", 0.01, None, True)
        if args[0] == "pbmm2":
            Path(args[4]).write_bytes(b"tiny-coordinate-sorted-bam")
        elif args[0] == "samtools":
            Path(args[-1]).write_bytes(b"tiny-bai")
        else:  # pragma: no cover - integration guard
            raise AssertionError(f"unexpected executable: {args[0]}")
        return CommandResult(args, 0, "", "", 0.2, None, True)


def build_context(tmp_path: Path, *, mixed: bool) -> AnalysisContext:
    inputs = tmp_path / "输入数据"
    inputs.mkdir(parents=True)
    fasta = inputs / "参考.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    fastq = inputs / "样本.fastq"
    fastq.write_text("@read1\nACGT\n+\nIIII\n", encoding="utf-8")
    rows = ["FASTQ1\t样本.fastq"]
    if mixed:
        bam = inputs / "existing.bam"
        bam.write_bytes(b"existing-bam")
        Path(f"{bam}.bai").write_bytes(b"existing-bai")
        cram = inputs / "existing.cram"
        cram.write_bytes(b"existing-cram")
        rows.extend(("BAM1\texisting.bam", "CRAM1\texisting.cram"))
    sheet = inputs / "samples.tsv"
    sheet.write_text(
        "sample_id\tinput\n" + "\n".join(rows) + "\n",
        encoding="utf-8",
    )
    user_config = inputs / "analysis.yaml"
    user_config.write_text(
        yaml.safe_dump(
            {
                "reference": {"fasta": "参考.fa", "build": "GRCh38"},
                "samples": {"sheet": "samples.tsv"},
                "alignment": {
                    "threads": 12,
                    "memory_mb": 48000,
                    "runtime_minutes": 600,
                    "index_threads": 3,
                },
            },
            allow_unicode=True,
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    return AnalysisContext.from_config(
        load_config(DEFAULT_CONFIG, STANDARD_PRESET, user_config)
    )


def wrappers(
    runner: FakeBioinformaticsRunner,
) -> tuple[Pbmm2Wrapper, SamtoolsWrapper]:
    return (
        Pbmm2Wrapper(
            runner=runner,  # type: ignore[arg-type]
            options=Pbmm2Options(),
        ),
        SamtoolsWrapper(runner=runner),  # type: ignore[arg-type]
    )


def test_complete_phase2_mixed_path_with_mock_tools(tmp_path: Path) -> None:
    context = build_context(tmp_path, mixed=True)
    runner = FakeBioinformaticsRunner()
    pbmm2, samtools = wrappers(runner)
    report = run_phase2(
        context,
        tmp_path / "结果" / "alignment",
        settings=Phase2Settings.from_config(context.config),
        pbmm2_wrapper=pbmm2,
        samtools_wrapper=samtools,
    )

    assert tuple(item.sample_id for item in report.sample_results) == (
        "FASTQ1",
        "BAM1",
        "CRAM1",
    )
    fastq_result, bam_result, cram_result = report.sample_results
    assert fastq_result.plan.requires_alignment is True
    assert fastq_result.artifact is not None
    assert fastq_result.artifact.index_path is not None
    assert fastq_result.alignment_qc.status is QCStatus.PASS
    assert bam_result.plan.requires_alignment is False
    assert bam_result.alignment_result.status.value == "reused"
    assert bam_result.alignment_qc.status is QCStatus.WARN
    assert cram_result.plan.requires_alignment is False
    assert cram_result.alignment_qc.status is QCStatus.WARN
    assert report.overall_qc_status is QCStatus.WARN
    assert report.tool_versions == {"pbmm2": "1.17.0", "samtools": "1.22.1"}

    executed = [command for command in runner.commands if "--version" not in command]
    assert sum(command[0] == "pbmm2" for command in executed) == 1
    assert sum(command[0] == "samtools" for command in executed) == 1
    assert all(
        forbidden not in " ".join(command).lower()
        for command in executed
        for forbidden in ("deepvariant", "sniffles", "sawfish", "trgt")
    )

    json_path = report.write_json(tmp_path / "报告" / "phase2.json")
    yaml_path = report.write_yaml(tmp_path / "报告" / "phase2.yaml")
    json_payload = json.loads(json_path.read_text(encoding="utf-8"))
    yaml_payload = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
    assert json_payload == yaml_payload == report.to_dict()
    assert "输入数据" in json_path.read_text(encoding="utf-8")


def test_complete_phase2_dry_run_plans_both_tools_without_installation(
    tmp_path: Path,
) -> None:
    context = build_context(tmp_path, mixed=False)
    runner = FakeBioinformaticsRunner(allow_execution=False)
    pbmm2, samtools = wrappers(runner)
    output_root = tmp_path / "dry-results"
    report = run_phase2(
        context,
        output_root,
        pbmm2_wrapper=pbmm2,
        samtools_wrapper=samtools,
        dry_run=True,
    )

    sample = report.sample_results[0]
    assert sample.alignment_result.status.value == "planned"
    assert sample.index_result is not None
    assert sample.index_result.status.value == "planned"
    assert sample.artifact is None
    assert sample.alignment_qc.status is QCStatus.NOT_CHECKED
    assert report.tool_versions == {"pbmm2": None, "samtools": None}
    assert runner.require_calls == []
    assert tuple(command[0] for command in runner.commands) == (
        "pbmm2",
        "samtools",
    )
    assert not output_root.exists()


def test_existing_only_phase2_never_checks_or_runs_external_tools(
    tmp_path: Path,
) -> None:
    full_context = build_context(tmp_path, mixed=True)
    context = AnalysisContext(
        reference=full_context.reference,
        samples=full_context.samples[1:],
        config=full_context.config,
        source_sample_sheet=full_context.source_sample_sheet,
    )
    runner = FakeBioinformaticsRunner(allow_execution=False)
    pbmm2, samtools = wrappers(runner)
    output_root = tmp_path / "unused-output"

    report = run_phase2(
        context,
        output_root,
        pbmm2_wrapper=pbmm2,
        samtools_wrapper=samtools,
    )

    assert tuple(item.plan.action.value for item in report.sample_results) == (
        "reuse",
        "reuse",
    )
    assert report.tool_versions == {"pbmm2": None, "samtools": None}
    assert runner.require_calls == []
    assert runner.commands == []
    assert not output_root.exists()
