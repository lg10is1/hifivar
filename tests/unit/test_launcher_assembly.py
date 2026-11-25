"""Dual-context assembly launcher tests; simulated tools do not establish accuracy."""

from __future__ import annotations

import hashlib
import json
import os
import runpy
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from hifivar import launcher
from hifivar.command import CommandRunner
from hifivar.config import HiFiVarConfig, load_config, validate_config
from hifivar.exceptions import ConfigurationError, InputValidationError, WorkflowError
from hifivar.hifiasm import HifiasmWrapper
from hifivar.launch_alignment import prepare_aligned_samples


ROOT = Path(__file__).parents[2]
DEFAULT = ROOT / "src/hifivar/resources/configs/default.yaml"


@pytest.fixture
def setup(tmp_path, monkeypatch):
    (tmp_path / "ref.fa").write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    (tmp_path / "reads.fastq").write_text("@m/1/ccs\nACGT\n+\nIIII\n", encoding="utf-8")
    (tmp_path / "reads2.fastq").write_text("@m/2/ccs\nTGCA\n+\nIIII\n", encoding="utf-8")
    (tmp_path / "aligned.bam").write_bytes(b"fake-bam-input")
    (tmp_path / "aligned.bam.bai").write_bytes(b"fake-index")
    (tmp_path / "primary.tsv").write_text("sample_id\tinput\tsex\tgroup\nS1\treads.fastq\tfemale\tG\n", encoding="utf-8")
    (tmp_path / "assembly.tsv").write_text("sample_id\tinput\tsex\tgroup\nS1\treads.fastq\tfemale\tG\n", encoding="utf-8")
    (tmp_path / "PAV.Snakefile").write_text("# Explicit native PAV deployment fixture\n", encoding="utf-8")
    settings = {
        "reference": {"fasta": "ref.fa", "build": "custom"},
        "samples": {"sheet": "primary.tsv"}, "runtime": {"threads": 4},
        "workflow": {"memory_mb": 4000},
        "assembly": {"enabled": True, "threads": 2, "memory_mb": 1000},
    }
    path = tmp_path / "user.yaml"
    path.write_text(yaml.safe_dump(settings), encoding="utf-8")
    wf = tmp_path / "workflow"
    wf.mkdir()
    (wf / "Snakefile").write_text("rule all:\n    input: []\n", encoding="utf-8")
    monkeypatch.setattr(launcher.platform, "system", lambda: "Linux")
    monkeypatch.setattr(launcher, "installed_workflow_root", lambda: wf)
    return load_config(DEFAULT, user_config=path)


def change(config, updates):
    value = config.to_dict()
    for section, fields in updates.items():
        value[section].update(fields)
    return HiFiVarConfig(value, config.sources)


class SchedulerProcesses:
    """Fake scheduler boundary, not a fake scientific success claim."""
    def __init__(self, *, fail=False):
        self.calls, self.fail = [], fail

    def run(self, command, **kwargs):
        self.calls.append((list(command), kwargs))
        assert command[0] != "samtools", "assembly-only must not probe or align BAM"
        out = "9.26.1\n" if "--version" in command else ""
        code = 1 if self.fail and "--version" not in command else 0
        for key in ("stdout_path", "stderr_path"):
            if kwargs.get(key):
                Path(kwargs[key]).write_text(out if key == "stdout_path" else "", encoding="utf-8")
        return SimpleNamespace(stdout=out, returncode=code, duration_seconds=0.001)


def test_assembly_only_raw_fastq_needs_no_alignment_or_calling_enabled(setup, tmp_path):
    effective, context, _, inputs = launcher.prepare_launch(setup, tmp_path / "run")
    assert effective["assembly"]["sample_sheet"] == str((tmp_path / "primary.tsv").resolve())
    assert (tmp_path / "reads.fastq").resolve() in inputs
    assert context.sample_ids == ("S1",)
    runner = SchedulerProcesses()
    root = tmp_path / "run"
    assert launcher.WorkflowLauncher(runner).run(setup, run_dir=root, dry_run=True) == "DRY_RUN_PASS"
    binding = json.loads(next((root / "attempts").glob("*/assembly_input_context.json")).read_text())
    assert binding["implicit_bam_extraction"] is False
    assert binding["biological_identity"].startswith("USER_DECLARED")
    assert not (root / "results/assembly").exists()
    assert not any("pbmm2" in " ".join(command) for command, _ in runner.calls)


def test_bam_primary_with_explicit_assembly_fastq_context(setup, tmp_path):
    (tmp_path / "primary.tsv").write_text("sample_id\tinput\tsex\tgroup\nS1\taligned.bam\tfemale\tG\n", encoding="utf-8")
    config = change(setup, {"assembly": {"sample_sheet": "assembly.tsv"}, "small": {"enabled": True, "threads": 2, "memory_mb": 1000}})
    effective, primary, _, inputs = launcher.prepare_launch(config, tmp_path / "run")
    assert primary.samples[0].sample.input.files == (tmp_path / "aligned.bam",)
    assembly = launcher.build_assembly_context(config, primary)
    assert assembly.samples[0].sample.input.files == (tmp_path / "reads.fastq",)
    assert effective["assembly"]["sample_sheet"] != effective["samples"]["sheet"]
    assert (tmp_path / "reads.fastq").resolve() in inputs
    assert (tmp_path / "assembly.tsv").resolve() in inputs


@pytest.mark.parametrize("name", ["aligned.bam", "native_ubam.bam"])
def test_bam_or_ubam_never_implicitly_becomes_assembly_fastq(setup, tmp_path, name):
    (tmp_path / name).write_bytes(b"fake-bam")
    (tmp_path / "primary.tsv").write_text(f"sample_id\tinput\tsex\tgroup\nS1\t{name}\tfemale\tG\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="explicit assembly.sample_sheet"):
        launcher.prepare_launch(setup, tmp_path / "run")
    assert not (tmp_path / "run").exists()


def test_assembly_only_bam_primary_does_not_require_unused_bam_index(setup, tmp_path):
    (tmp_path / "aligned.bam.bai").unlink()
    (tmp_path / "primary.tsv").write_text("sample_id\tinput\tsex\tgroup\nS1\taligned.bam\tfemale\tG\n", encoding="utf-8")
    config = change(setup, {"assembly": {"sample_sheet": "assembly.tsv"}})
    launcher.prepare_launch(config, tmp_path / "run")
    assert launcher.WorkflowLauncher(SchedulerProcesses()).run(config, run_dir=tmp_path / "run") == "WORKFLOW_SUCCEEDED"


@pytest.mark.parametrize("row,match", [
    ("S2\treads.fastq\tfemale\tG", "ASSEMBLY_SAMPLE_SET_MISMATCH"),
    ("S1\treads.fastq\tmale\tG", "ASSEMBLY_SAMPLE_METADATA_MISMATCH"),
    ("S1\treads.fastq\tfemale\tOTHER", "ASSEMBLY_SAMPLE_METADATA_MISMATCH"),
    ("S1\taligned.bam\tfemale\tG", "implicit BAM"),
])
def test_explicit_sheet_must_match_identity_metadata_and_fastq_type(setup, tmp_path, row, match):
    (tmp_path / "assembly.tsv").write_text("sample_id\tinput\tsex\tgroup\n" + row + "\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match=match):
        launcher.prepare_launch(change(setup, {"assembly": {"sample_sheet": "assembly.tsv"}}), tmp_path / "run")


def test_matching_sample_set_accepts_declared_order_difference_by_name(setup, tmp_path):
    (tmp_path / "primary.tsv").write_text("sample_id\tinput\tsex\tgroup\nS1\treads.fastq\tfemale\tG\nS2\treads2.fastq\tmale\tG\n", encoding="utf-8")
    (tmp_path / "assembly.tsv").write_text("sample_id\tinput\tsex\tgroup\nS2\treads2.fastq\tmale\tG\nS1\treads.fastq\tfemale\tG\n", encoding="utf-8")
    config = change(setup, {"assembly": {"sample_sheet": "assembly.tsv"}})
    context = launcher.build_assembly_context(config)
    assert context.sample_ids == ("S2", "S1")
    launcher.prepare_launch(config, tmp_path / "run")


def test_multiple_fastq_order_is_preserved(setup, tmp_path):
    (tmp_path / "primary.tsv").write_text("sample_id\tinput\tsex\tgroup\nS1\treads2.fastq;reads.fastq\tfemale\tG\n", encoding="utf-8")
    context = launcher.build_assembly_context(setup)
    assert [path.name for path in context.samples[0].sample.input.files] == ["reads2.fastq", "reads.fastq"]


def test_fastq_hardlink_aliases_are_not_assembled_twice(setup, tmp_path):
    os.link(tmp_path / "reads.fastq", tmp_path / "alias.fastq")
    (tmp_path / "primary.tsv").write_text("sample_id\tinput\tsex\tgroup\nS1\treads.fastq;alias.fastq\tfemale\tG\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="assigned more than once"):
        launcher.prepare_launch(setup, tmp_path / "run")


@pytest.mark.parametrize("updates,match", [
    ({"assembly": {"threads": 5}}, "threads exceeds"),
    ({"assembly": {"memory_mb": 4001}}, "memory_mb exceeds"),
    ({"assembly": {"overwrite": True}}, "overwrite=false"),
    ({"assembly": {"enabled": False}, "assembly_sv": {"enabled": True}}, "requires assembly.enabled"),
    ({"assembly_sv": {"enabled": True, "overwrite": True}}, "overwrite=false"),
])
def test_assembly_launch_resource_and_overwrite_guards(setup, tmp_path, updates, match):
    with pytest.raises(ConfigurationError, match=match):
        launcher.prepare_launch(change(setup, updates), tmp_path / "run")


def _pav_config(config, tmp_path):
    result = config.to_dict()
    result["assembly_sv"]["enabled"] = True
    result["assembly_sv"]["pav"].update(snakefile="PAV.Snakefile", version="2.4.6", threads=2, memory_mb=1000)
    result["assembly_sv"]["svim_asm"].update(threads=2, memory_mb=1000)
    return HiFiVarConfig(result, config.sources)


def test_pav_snakefile_is_explicit_hashed_dependency_and_resume_guard(setup, tmp_path):
    config = _pav_config(setup, tmp_path)
    effective, _, _, inputs = launcher.prepare_launch(config, tmp_path / "run")
    assert effective["assembly_sv"]["pav"]["snakefile"] == str((tmp_path / "PAV.Snakefile").resolve())
    assert (tmp_path / "PAV.Snakefile").resolve() in inputs
    root = tmp_path / "run"
    launcher.WorkflowLauncher(SchedulerProcesses()).run(config, run_dir=root, dry_run=True)
    manifest = json.loads((root / "launch_manifest.json").read_text())
    records = {item["path"]: item for item in manifest["identity"]["inputs"]}
    assert records[str((tmp_path / "PAV.Snakefile").resolve())]["sha256"] == hashlib.sha256((tmp_path / "PAV.Snakefile").read_bytes()).hexdigest()
    with (tmp_path / "PAV.Snakefile").open("a") as handle:
        handle.write("#changed\n")
    with pytest.raises(InputValidationError, match="RESUME_IDENTITY_MISMATCH"):
        launcher.WorkflowLauncher(SchedulerProcesses()).run(config, run_dir=root, resume=True)


def test_missing_pav_snakefile_stops_before_creating_run(setup, tmp_path):
    config = _pav_config(setup, tmp_path)
    (tmp_path / "PAV.Snakefile").unlink()
    with pytest.raises(InputValidationError, match="missing or empty"):
        launcher.prepare_launch(config, tmp_path / "run")
    assert not (tmp_path / "run").exists()


@pytest.mark.parametrize("value", [123, [], {}, ""])
def test_new_assembly_sample_sheet_schema_rejects_invalid_types(value):
    with pytest.raises(ConfigurationError):
        validate_config({"assembly": {"sample_sheet": value}})


def test_native_pav_and_svim_are_both_selected_but_no_container_support_is_invented(setup, tmp_path):
    config = _pav_config(setup, tmp_path)
    names = [name for name, _ in launcher._enabled_tracks(config.to_dict())]
    assert names == ["assembly", "assembly_sv.pav", "assembly_sv.svim_asm"]
    assert "image" not in config["assembly_sv"]["pav"]


class PythonExecutableFixtureRunner(CommandRunner):
    """Run a Python fake binary via the real shell-free CommandRunner."""
    def __init__(self, script):
        super().__init__()
        self.script = script

    def require_executable(self, executable, **kwargs):
        return super().require_executable(sys.executable, **kwargs)

    def run(self, command, **kwargs):
        return super().run([sys.executable, "-B", str(self.script), *map(str, command[1:])], **kwargs)


def test_hifiasm_bridge_uses_declared_fastq_with_real_fake_executable(setup, tmp_path, monkeypatch):
    import hifivar.hifiasm as module
    script = tmp_path / "fake_hifiasm.py"
    script.write_text(
        "import pathlib,sys\n"
        "if '--version' in sys.argv: print('hifiasm 0.25.0-r726'); raise SystemExit(0)\n"
        "prefix=sys.argv[sys.argv.index('-o')+1]\n"
        "for suffix in ('.bp.p_ctg.gfa','.bp.hap1.p_ctg.gfa','.bp.hap2.p_ctg.gfa'):\n"
        " pathlib.Path(prefix+suffix).write_text('H\\tVN:Z:1.0\\nS\\tctg1\\tACGT\\n')\n",
        encoding="utf-8",
    )
    runner = PythonExecutableFixtureRunner(script)
    real_wrapper = HifiasmWrapper
    monkeypatch.setattr(module, "HifiasmWrapper", lambda **kw: real_wrapper(runner=runner, **kw))
    context = launcher.build_assembly_context(setup)
    snakemake = SimpleNamespace(
        wildcards=SimpleNamespace(sample="S1"), input=SimpleNamespace(reads=list(context.samples[0].sample.input.files)),
        params=SimpleNamespace(output_prefix=tmp_path / "work/S1.asm", assembly_directory=tmp_path / "results"),
        threads=2, resources=SimpleNamespace(mem_mb=1000, runtime_min=1),
        config=setup.to_dict(), log=[tmp_path / "hifiasm.stderr.log"],
    )
    before = (tmp_path / "reads.fastq").read_bytes()
    runpy.run_path(str(ROOT / "workflow/scripts/run_hifiasm.py"), init_globals={"snakemake": snakemake})
    assert len(list((tmp_path / "results").glob("S1.*.fa"))) == 3
    assert (tmp_path / "reads.fastq").read_bytes() == before


def test_raw_alignment_calling_config_preserves_independent_assembly_context(setup, tmp_path, monkeypatch):
    config = change(setup, {"alignment": {"enabled": True, "threads": 2, "index_threads": 1, "memory_mb": 1000}, "small": {"enabled": True, "threads": 2, "memory_mb": 1000}})
    effective, context, _, _ = launcher.prepare_launch(config, tmp_path / "run")
    root = tmp_path / "run"
    root.mkdir()
    attempt = root / "attempt"
    attempt.mkdir()
    # This synthetic test isolates independent assembly-sheet preservation, not
    # native BAM preparation or a launcher execution receipt. Those contracts
    # have genuine subprocess tests and a separate real-tool Linux gate.
    monkeypatch.setattr("hifivar.launch_alignment.classify_launch_inputs", lambda *args: [{"action": "REUSE", "kind": "ALIGNED_BAM"}])
    monkeypatch.setattr("hifivar.launch_alignment.IndexedBamPreflight.run", lambda *args: {})
    monkeypatch.setattr("hifivar.launch_alignment.publish_preparation", lambda *args, **kwargs: None)
    monkeypatch.setattr("hifivar.launch_alignment.seal_handoff", lambda *args: None)
    # Supply a BAM primary context with the original FASTQ sheet declared for assembly.
    (tmp_path / "primary.tsv").write_text("sample_id\tinput\tsex\tgroup\nS1\taligned.bam\tfemale\tG\n", encoding="utf-8")
    config = change(config, {"assembly": {"sample_sheet": "assembly.tsv"}})
    effective, context, _, _ = launcher.prepare_launch(config, root)
    calling = prepare_aligned_samples(context, effective, root, attempt, SchedulerProcesses(), dry_run=False, rerun_incomplete=False)
    generated = yaml.safe_load(calling.read_text(encoding="utf-8"))
    assert generated["samples"]["sheet"] == str(root / "calling.samples.tsv")
    assert generated["assembly"]["sample_sheet"] == str((tmp_path / "assembly.tsv").resolve())
    assembly = launcher.build_assembly_context(generated)
    assert assembly.samples[0].sample.input.files == (tmp_path / "reads.fastq",)
