"""Launcher scope, input identity, resources and non-destructive resume contracts."""
from __future__ import annotations

import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from hifivar import launcher
from hifivar.cli import build_parser, main
from hifivar.config import HiFiVarConfig, load_config, validate_config
from hifivar.exceptions import ConfigurationError, InputValidationError, OutputValidationError, WorkflowError
from hifivar.launch_inputs import IndexedBamPreflight, unique_bam_index, validate_bam_header
from hifivar.reference import Contig, ReferenceGenome

HEADER = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:4\n@RG\tID:r\tSM:S1\tPL:PACBIO\n"
DEFAULT = Path(__file__).parents[2] / "src/hifivar/resources/configs/default.yaml"


class FakeRunner:
    """Deterministic process boundary, not a scientific output simulator."""
    def __init__(self, header=HEADER, fail_dag=False, interrupt=False):
        self.header, self.fail_dag, self.interrupt = header, fail_dag, interrupt
        self.commands = []

    def require_executable(self, executable):
        return Path(executable)

    def run(self, command, **kwargs):
        self.commands.append((command, kwargs))
        stdout = ""
        returncode = 0
        if command[0] == "samtools":
            if "--version" in command:
                stdout = "samtools 1.23.1\nUsing htslib 1.23.1\n"
            elif "-H" in command:
                stdout = self.header
            elif "-X" in command:
                stdout = "0\n"
        elif "--version" in command:
            stdout = "9.27.0\n"
        else:
            if self.interrupt:
                raise KeyboardInterrupt()
            returncode = 1 if self.fail_dag else 0
        for key, content in (("stdout_path", stdout), ("stderr_path", "")):
            if kwargs.get(key):
                Path(kwargs[key]).write_text(content, encoding="utf-8")
        return SimpleNamespace(args=tuple(command), stdout=stdout, returncode=returncode, duration_seconds=0.01)


@pytest.fixture
def config(tmp_path):
    (tmp_path / "ref.fa").write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    (tmp_path / "S1.bam").write_bytes(b"not-a-real-bam-test-process-boundary-only")
    (tmp_path / "S1.bam.bai").write_bytes(b"test-index")
    (tmp_path / "samples.tsv").write_text("sample_id\tinput\tsex\nS1\tS1.bam\tfemale\n", encoding="utf-8")
    (tmp_path / "user.yaml").write_text(yaml.safe_dump({
        "reference": {"fasta": "ref.fa", "build": "custom"}, "samples": {"sheet": "samples.tsv"},
        "runtime": {"threads": 2}, "small": {"enabled": True, "threads": 2, "memory_mb": 100},
        "workflow": {"memory_mb": 200},
    }), encoding="utf-8")
    return load_config(default_config=DEFAULT, user_config=tmp_path / "user.yaml")


@pytest.fixture
def linux(monkeypatch, tmp_path):
    monkeypatch.setattr(launcher.platform, "system", lambda: "Linux")
    wf = tmp_path / "workflow"
    wf.mkdir()
    (wf / "Snakefile").write_text("rule all:\n    input: []\n", encoding="utf-8")
    monkeypatch.setattr(launcher, "installed_workflow_root", lambda: wf)


def change(config, section, key, value):
    result = config.to_dict()
    result[section][key] = value
    return HiFiVarConfig(result, config.sources)


def test_prepare_resolves_inputs_and_preserves_original(config, tmp_path):
    original = config.to_dict()
    root = tmp_path / "run"
    effective, context, sheet, inputs = launcher.prepare_launch(config, root)
    assert effective["paths"]["outdir"] == str(root / "results")
    assert str(tmp_path / "S1.bam") in sheet and "female" in sheet
    assert len(inputs) == 6 and context.sample_ids == ("S1",)
    assert config.to_dict() == original and not root.exists()


@pytest.mark.parametrize("section,key,value,match", [
    ("small", "enabled", False, "at least one"),
    ("small", "threads", 3, "threads exceeds"),
    ("small", "memory_mb", 201, "memory_mb exceeds"),
    ("small", "overwrite", True, "overwrite=false"),
    ("paths", "outdir", "../escape", "owned run"),
    ("runtime", "tmpdir", "/tmp/escape", "owned run"),
    ("logging", "file", "outside.log", "owned run"),
])
def test_invalid_launch_options(config, tmp_path, section, key, value, match):
    with pytest.raises(ConfigurationError, match=match):
        launcher.prepare_launch(change(config, section, key, value), tmp_path / "run")


@pytest.mark.parametrize("name", launcher._UNSUPPORTED)
def test_unrouted_modules_fail_closed(config, name):
    values = config.to_dict()
    values[name]["enabled"] = True
    with pytest.raises(ConfigurationError, match="does not yet route"):
        launcher._enabled_tracks(values)


@pytest.mark.parametrize("value", [0, -1, True, "100"])
def test_workflow_memory_schema(value):
    with pytest.raises(ConfigurationError):
        validate_config({"workflow": {"memory_mb": value}})


@pytest.mark.parametrize("value", [None, "", "  ", 123])
def test_samtools_schema(value):
    with pytest.raises(ConfigurationError):
        validate_config({"workflow": {"samtools_executable": value}})


def test_command_has_explicit_all_and_no_shell(config, tmp_path):
    command = launcher.build_workflow_command(tmp_path, tmp_path, config.to_dict(), dry_run=True, rerun_incomplete=True)
    assert command[-2:] == ["--", "all"]
    assert "--dry-run" in command and "--rerun-incomplete" in command
    assert "mem_mb=200" in command and "deepvariant_slots=1" in command


def test_header_accepts_exact_identity(tmp_path):
    ref = ReferenceGenome(tmp_path / "r", tmp_path / "f", "custom", (Contig("chr1", 4),))
    validate_bam_header(HEADER, "S1", ref)


@pytest.mark.parametrize("header", [
    HEADER.replace("coordinate", "unsorted"), HEADER.replace("SM:S1", "SM:UnknownSampleName"),
    HEADER.replace("\tSM:S1", ""), HEADER.split("@RG")[0], HEADER.replace("LN:4", "LN:5"),
    HEADER.replace("SN:chr1", "SN:1"), HEADER + "@RG\tID:r\tSM:S1\n",
    HEADER + "@RG\tID:second\tSM:S2\n", HEADER + "@SQ\tSN:chr1\tLN:4\n",
    HEADER.replace("LN:4", "LN:bad"), HEADER.replace("ID:r", "ID:"),
    HEADER.replace("SM:S1", "SM:S1\tSM:S1"),
])
def test_bad_headers_rejected(header, tmp_path):
    ref = ReferenceGenome(tmp_path / "r", tmp_path / "f", "custom", (Contig("chr1", 4),))
    with pytest.raises(InputValidationError):
        validate_bam_header(header, "S1", ref)


def test_ambiguous_or_empty_indexes_rejected(config, tmp_path):
    bam = tmp_path / "S1.bam"
    (tmp_path / "S1.bai").write_bytes(b"other")
    with pytest.raises(InputValidationError):
        unique_bam_index(bam)
    (tmp_path / "S1.bai").unlink()
    (tmp_path / "S1.bam.bai").write_bytes(b"")
    with pytest.raises(InputValidationError):
        unique_bam_index(bam)


def test_hardlink_bam_reuse_rejected(config, tmp_path):
    os.link(tmp_path / "S1.bam", tmp_path / "alias.bam")
    (tmp_path / "alias.bam.bai").write_bytes(b"index")
    (tmp_path / "samples.tsv").write_text("sample_id\tinput\nS1\tS1.bam\nS2\talias.bam\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="same BAM"):
        launcher.prepare_launch(config, tmp_path / "run")


@pytest.mark.parametrize("extension,content", [("fastq", "@r\nACGT\n+\nIIII\n"), ("cram", "fake")])
def test_raw_and_cram_not_accepted(config, tmp_path, extension, content):
    (tmp_path / f"S1.{extension}").write_text(content, encoding="utf-8")
    (tmp_path / "samples.tsv").write_text(f"sample_id\tinput\nS1\tS1.{extension}\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="currently accepts"):
        launcher.prepare_launch(config, tmp_path / "run")


def test_dry_run_and_resume_preserve_inputs(config, tmp_path, linux, monkeypatch):
    monkeypatch.setenv("SNAKEMAKE_PROFILE", "UNTRUSTED_PROFILE")
    monkeypatch.setenv("PRIVATE_SECRET_TOKEN", "do-not-export")
    runner = FakeRunner()
    root = tmp_path / "run"
    assert launcher.WorkflowLauncher(runner).run(config, run_dir=root, dry_run=True) == "DRY_RUN_PASS"
    manifest = (root / "launch_manifest.json").read_bytes()
    assert launcher.WorkflowLauncher(runner).run(config, run_dir=root, resume=True) == "WORKFLOW_SUCCEEDED"
    assert (root / "launch_manifest.json").read_bytes() == manifest
    assert len(list((root / "attempts").iterdir())) == 2
    assert not (root / ".hifivar-launch.lock").exists()
    dag = [(c, k) for c, k in runner.commands if c[-2:] == ["--", "all"]]
    assert all(k["inherit_env"] is False and "PRIVATE_SECRET_TOKEN" not in k["env"] and "SNAKEMAKE_PROFILE" not in k["env"] for c, k in dag)
    assert all("shell" not in k for c, k in runner.commands)


@pytest.mark.parametrize("target", ["S1.bam", "ref.fa", "S1.bam.bai", "effective.yaml", "samples.tsv"])
def test_resume_identity_and_control_tamper(config, tmp_path, linux, target):
    obj = launcher.WorkflowLauncher(FakeRunner())
    root = tmp_path / "run"
    obj.run(config, run_dir=root, dry_run=True)
    file = root / target if target in {"effective.yaml", "samples.tsv"} else tmp_path / target
    with file.open("ab") as stream:
        stream.write(b"\nchanged")
    with pytest.raises(InputValidationError, match="RESUME_"):
        obj.run(config, run_dir=root, resume=True)


def test_resume_runtime_source_drift(config, tmp_path, linux):
    obj = launcher.WorkflowLauncher(FakeRunner())
    root = tmp_path / "run"
    obj.run(config, run_dir=root, dry_run=True)
    (tmp_path / "workflow/Snakefile").write_text("changed", encoding="utf-8")
    with pytest.raises(InputValidationError, match="RESUME_IDENTITY_MISMATCH"):
        obj.run(config, run_dir=root, resume=True)


def test_no_automatic_retry_and_failure_evidence(config, tmp_path, linux):
    runner = FakeRunner(fail_dag=True)
    root = tmp_path / "run"
    with pytest.raises(WorkflowError, match="No automatic retry"):
        launcher.WorkflowLauncher(runner).run(config, run_dir=root)
    status = json.loads(next((root / "attempts").glob("*/status.json")).read_text())
    assert status["status"] == "FAILED" and status["returncode"] == 1
    assert len([c for c, k in runner.commands if c[-2:] == ["--", "all"]]) == 1


def test_identity_failure_does_not_launch_dag(config, tmp_path, linux):
    runner = FakeRunner(header=HEADER.replace("SM:S1", "SM:other"))
    with pytest.raises(InputValidationError):
        launcher.WorkflowLauncher(runner).run(config, run_dir=tmp_path / "run")
    assert not any(c[-2:] == ["--", "all"] for c, k in runner.commands)


def test_interrupt_retains_lock(config, tmp_path, linux):
    root = tmp_path / "run"
    with pytest.raises(KeyboardInterrupt):
        launcher.WorkflowLauncher(FakeRunner(interrupt=True)).run(config, run_dir=root)
    assert (root / ".hifivar-launch.lock").is_file()
    status = json.loads(next((root / "attempts").glob("*/status.json")).read_text())
    assert status["status"] == "INTERRUPTED_LOCK_RETAINED"


def test_existing_or_locked_run_not_overwritten(config, tmp_path, linux):
    obj = launcher.WorkflowLauncher(FakeRunner())
    root = tmp_path / "run"
    obj.run(config, run_dir=root, dry_run=True)
    with pytest.raises(OutputValidationError):
        obj.run(config, run_dir=root)
    (root / ".hifivar-launch.lock").write_text("other", encoding="utf-8")
    with pytest.raises(OutputValidationError, match="locked"):
        obj.run(config, run_dir=root, resume=True)
    assert (root / ".hifivar-launch.lock").read_text() == "other"


def test_retry_requires_resume(config, tmp_path, linux):
    with pytest.raises(ConfigurationError, match="requires --resume"):
        launcher.WorkflowLauncher(FakeRunner()).run(config, run_dir=tmp_path / "run", rerun_incomplete=True)


def test_region_is_frozen_in_launcher_manifest(config, tmp_path, linux):
    obj = launcher.WorkflowLauncher(FakeRunner())
    root = tmp_path / "run"
    regional = change(config, "small", "region", "chr1:1-4")
    obj.run(regional, run_dir=root, dry_run=True)
    assert yaml.safe_load((root / "effective.yaml").read_text())["small"]["region"] == "chr1:1-4"
    with pytest.raises(InputValidationError, match="RESUME_IDENTITY_MISMATCH"):
        obj.run(change(regional, "small", "region", "chr1:2-4"), run_dir=root, resume=True)


def test_region_bounds_checked_before_launcher_creates_run(config, tmp_path):
    from hifivar.exceptions import ReferenceError
    root = tmp_path / "run"
    with pytest.raises(ReferenceError):
        launcher.prepare_launch(change(config, "small", "region", "chr1:1-5"), root)
    assert not root.exists()


def test_linux_only(config, tmp_path, monkeypatch):
    monkeypatch.setattr(launcher.platform, "system", lambda: "Windows")
    with pytest.raises(WorkflowError, match="requires Linux"):
        launcher.WorkflowLauncher(FakeRunner()).run(config, run_dir=tmp_path / "run")


def test_cli_parser_and_dispatch(config, tmp_path, monkeypatch, capsys):
    args = build_parser().parse_args(["--dry-run", "run", "--run-dir", "new-run", "--resume"])
    assert args.dry_run and args.resume
    calls = []
    monkeypatch.setattr(launcher.WorkflowLauncher, "run", lambda self, cfg, **kw: calls.append(kw) or "DRY_RUN_PASS")
    assert main(["--config", str(tmp_path / "user.yaml"), "--dry-run", "run", "--run-dir", str(tmp_path / "run")]) == 0
    assert calls[0]["dry_run"] and "DRY_RUN_PASS" in capsys.readouterr().out
