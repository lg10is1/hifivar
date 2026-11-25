"""Raw launcher routing and failure safety, with explicit simulated processes."""
from pathlib import Path
from types import SimpleNamespace
import csv
import json

import pytest
import yaml

from hifivar import launcher
from hifivar.alignment import AlignmentRequest, AlignmentResources, AlignmentTool
from hifivar.config import load_config, validate_config
from hifivar.exceptions import ConfigurationError, InputValidationError, OutputValidationError, WorkflowError
from hifivar.launch_alignment import ScopedAlignmentRunner, validate_alignment_budget
from hifivar.launch_inputs import unique_bam_index
from hifivar.pbmm2 import Pbmm2Options, Pbmm2Wrapper
from hifivar.ubam import Pbmm2UbamWrapper, validate_ubam_header

DEFAULT = Path(__file__).parents[2] / "src/hifivar/resources/configs/default.yaml"
UBAM_HEADER = "@HD\tVN:1.6\tSO:unsorted\n@RG\tID:raw\tSM:S1\tPL:PACBIO\n"
BAM_HEADER = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:4\n@RG\tID:raw\tSM:S1\tPL:PACBIO\n"


class Processes:
    """Fake binary I/O for orchestration tests; not biological validation."""

    def __init__(self, *, fail_align=False, header=UBAM_HEADER, mapped=0, invalid=0):
        self.commands = []
        self.fail_align, self.header, self.mapped, self.invalid = fail_align, header, mapped, invalid

    def require_executable(self, executable):
        return Path(executable)

    def run(self, command, **kw):
        self.commands.append((command, kw))
        out = ""
        if "--version" in command:
            out = {"samtools": "samtools 1.23.1\n", "pbmm2": "pbmm2 1.17.0\n"}.get(command[0], "9.27.0\n")
        elif command[0] == "pbmm2":
            if self.fail_align:
                raise WorkflowError("simulated alignment failure")
            Path(command[4]).write_bytes(b"fake-coordinate-bam")
        elif command[0] == "samtools":
            if command[1] == "index":
                Path(command[-1]).write_bytes(b"fake-index")
            elif "-H" in command:
                out = self.header if "raw.bam" in command[-1] else BAM_HEADER
            elif "-c" in command:
                out = str(self.mapped if "-F" in command else self.invalid if "--rf" in command else 0 if "-X" in command else 3)
        for key, content in (("stdout_path", out), ("stderr_path", "")):
            if kw.get(key):
                Path(kw[key]).write_text(content, encoding="utf-8")
        return SimpleNamespace(args=tuple(command), stdout=out, stderr="", returncode=0, duration_seconds=0.001,
                               executed=True, stdout_path=kw.get("stdout_path"), stderr_path=kw.get("stderr_path"))


@pytest.fixture
def setup(tmp_path, monkeypatch):
    # This fixture tests simulated orchestration only. Strict observed receipt
    # contracts are covered independently by test_alignment_preparation.py.
    from hifivar import alignment_preparation, launch_alignment
    def simulated_receipt(root, attempt, route, bam, index, *args, **kwargs):
        path = alignment_preparation.preparation_path(attempt, route["sample_id"])
        path.parent.mkdir(exist_ok=True)
        path.write_text('{"fixture_scope":"SIMULATED_NOT_NATIVE"}')
    def simulated_handoff(root, attempt):
        with (root / "calling.samples.tsv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        items = [launcher.file_identity(path) for row in rows for path in
                 (Path(row["input"]), unique_bam_index(Path(row["input"])))]
        (attempt / "alignment_handoff.json").write_text(json.dumps({"selected_artifacts": items}))
    def simulated_verify(root, attempt, binding, **kwargs):
        return json.loads((attempt / "alignment_handoff.json").read_text())
    monkeypatch.setattr(launch_alignment, "publish_preparation", simulated_receipt)
    monkeypatch.setattr(launch_alignment, "seal_handoff", simulated_handoff)
    monkeypatch.setattr(alignment_preparation, "verify_handoff", simulated_verify)
    (tmp_path / "ref.fa").write_text(">chr1\nACGT\n")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n")
    (tmp_path / "raw.fastq").write_text("@m/1/ccs\nACGT\n+\nIIII\n")
    (tmp_path / "raw.bam").write_bytes(b"fake-unaligned-bam")
    (tmp_path / "samples.tsv").write_text("sample_id\tinput\tsex\nS1\traw.fastq\tfemale\n")
    settings = {"reference": {"fasta": "ref.fa", "build": "custom"}, "samples": {"sheet": "samples.tsv"},
                "runtime": {"threads": 4}, "alignment": {"enabled": True, "threads": 2, "index_threads": 1},
                "small": {"enabled": True, "threads": 2}}
    (tmp_path / "user.yaml").write_text(yaml.safe_dump(settings))
    wf = tmp_path / "workflow"
    wf.mkdir()
    (wf / "Snakefile").write_text("rule all:\n    input: []\n")
    monkeypatch.setattr(launcher.platform, "system", lambda: "Linux")
    monkeypatch.setattr(launcher, "installed_workflow_root", lambda: wf)
    return load_config(default_config=DEFAULT, user_config=tmp_path / "user.yaml")


@pytest.mark.parametrize("key,value", [("enabled", "true"), ("pbmm2_executable", ""), ("minimap2_executable", None), ("pbmm2_sort_threads", 0), ("pbmm2_sort_threads", True)])
def test_raw_schema_rejects_invalid(key, value):
    with pytest.raises(ConfigurationError):
        validate_config({"alignment": {key: value}})


@pytest.mark.parametrize("key,value", [("threads", 5), ("memory_mb", 64001), ("index_threads", 4),
                                       ("pbmm2_sort_threads", 2), ("overwrite", True), ("output_format", "cram")])
def test_resource_and_output_boundary(setup, key, value):
    cfg = setup.to_dict()
    cfg["alignment"][key] = value
    with pytest.raises(ConfigurationError):
        validate_alignment_budget(cfg)


@pytest.mark.parametrize("header", [UBAM_HEADER.replace("SM:S1", "SM:S2"), UBAM_HEADER.replace("\tSM:S1", ""),
                                     UBAM_HEADER.replace("PL:PACBIO", "PL:ILLUMINA"), UBAM_HEADER.replace("unsorted", "coordinate"),
                                     UBAM_HEADER + "@SQ\tSN:chr1\tLN:4\n", UBAM_HEADER + "@RG\tID:raw\tSM:S1\tPL:PACBIO\n",
                                     UBAM_HEADER.replace("ID:raw", "ID:raw\tID:other"), ""])
def test_ubam_header_rejects_ambiguity(header):
    with pytest.raises(InputValidationError):
        validate_ubam_header(header, "S1")


def test_ubam_native_command_and_no_conversion(setup, tmp_path):
    context = launcher.AnalysisContext.from_config(setup)
    runner = Processes()
    result = Pbmm2UbamWrapper(runner, "pbmm2", "samtools", Pbmm2Options(sort_threads=1)).run(
        tmp_path / "raw.bam", "S1", context.reference, tmp_path / "new.bam", threads=1, timeout=60)
    cmd = result["command"]
    assert cmd[3].endswith("raw.bam") and "--rg" not in cmd and "--sample" not in cmd
    assert result["conversion"] == "NONE_NATIVE_UBAM" and "-J" in cmd
    assert not any("fastq" in c for c, k in runner.commands)


@pytest.mark.parametrize("kwargs", [{"mapped": 1}, {"invalid": 1}])
def test_nonprimary_or_mapped_ubam_rejected(setup, tmp_path, kwargs):
    wrapper = Pbmm2UbamWrapper(Processes(**kwargs), "pbmm2", "samtools", Pbmm2Options(sort_threads=1))
    with pytest.raises(InputValidationError, match="unpaired primary unmapped"):
        wrapper.inspect(tmp_path / "raw.bam", "S1", timeout=60)


def test_pbmm2_explicit_sort_setting_preserves_legacy_default(setup, tmp_path):
    ctx = launcher.AnalysisContext.from_config(setup)
    req = AlignmentRequest(sample=ctx.samples[0].sample, reference=ctx.reference, output_path=tmp_path / "out.bam",
                           resources=AlignmentResources(threads=1), tool=AlignmentTool.PBMM2)
    assert "-J" not in Pbmm2Wrapper().build_command(req)
    command = Pbmm2Wrapper(options=Pbmm2Options(sort_threads=2)).build_command(req)
    assert command[-2:] == ["-J", "2"]


def test_raw_plan_has_no_alignment_or_fake_bam(setup, tmp_path):
    runner = Processes()
    root = tmp_path / "run"
    assert launcher.WorkflowLauncher(runner).run(setup, run_dir=root, dry_run=True) == "RAW_INPUT_PLAN_ONLY"
    assert not list(root.rglob("*.bam")) and not (root / "calling.yaml").exists()
    assert not any(c[0] == "pbmm2" or c[-2:] == ["--", "all"] for c, k in runner.commands)
    status = json.loads(next((root / "attempts").glob("*/status.json")).read_text())
    assert status["calling_dag"] == "NOT_RUN"


@pytest.mark.parametrize("kind", ["fastq", "bam"])
def test_run_then_resume_reuses_sealed_alignment(setup, tmp_path, kind):
    (tmp_path / "samples.tsv").write_text(f"sample_id\tinput\tsex\nS1\traw.{kind}\tfemale\n")
    runner = Processes()
    root = tmp_path / "run"
    job = launcher.WorkflowLauncher(runner)
    before = (tmp_path / f"raw.{kind}").read_bytes()
    assert job.run(setup, run_dir=root) == "WORKFLOW_SUCCEEDED"
    receipt = root / "results/alignment/S1/completed.json"
    saved = receipt.read_bytes()
    assert job.run(setup, run_dir=root, resume=True) == "WORKFLOW_SUCCEEDED"
    assert receipt.read_bytes() == saved and (tmp_path / f"raw.{kind}").read_bytes() == before
    assert len([c for c, k in runner.commands if c[0] == "pbmm2" and "align" in c]) == 1
    assert "aligned.bam" in (root / "calling.samples.tsv").read_text()
    assert all(str(root / "calling.yaml") in c for c, k in runner.commands if c[-2:] == ["--", "all"])


def test_failed_alignment_requires_explicit_retry_and_preserves_evidence(setup, tmp_path):
    runner = Processes(fail_align=True)
    root = tmp_path / "run"
    job = launcher.WorkflowLauncher(runner)
    with pytest.raises(WorkflowError, match="simulated"):
        job.run(setup, run_dir=root)
    assert not (root / "results/alignment/S1/completed.json").exists()
    assert not any(c[-2:] == ["--", "all"] for c, k in runner.commands)
    with pytest.raises(OutputValidationError, match="rerun-incomplete"):
        job.run(setup, run_dir=root, resume=True)
    runner.fail_align = False
    assert job.run(setup, run_dir=root, resume=True, rerun_incomplete=True) == "WORKFLOW_SUCCEEDED"
    assert len(list((root / "results/alignment/S1").glob("attempt_*"))) == 2


@pytest.mark.parametrize("target", ["bam", "index", "calling.yaml", "calling.samples.tsv"])
def test_derived_alignment_or_handoff_tampering_blocks_resume(setup, tmp_path, target):
    root = tmp_path / "run"
    job = launcher.WorkflowLauncher(Processes())
    job.run(setup, run_dir=root)
    receipt = json.loads((root / "results/alignment/S1/completed.json").read_text())
    path = Path(receipt["artifacts"][0 if target == "bam" else 1]["path"]) if target in {"bam", "index"} else root / target
    with path.open("ab") as out:
        out.write(b"tamper")
    with pytest.raises(InputValidationError):
        job.run(setup, run_dir=root, resume=True)


def test_scoped_environment_does_not_export_secrets(tmp_path, monkeypatch):
    monkeypatch.setenv("API_TOKEN", "sensitive")
    monkeypatch.setenv("PYTHONPATH", "outside")
    runner = Processes()
    ScopedAlignmentRunner(runner, tmp_path).run(["samtools", "--version"])
    kw = runner.commands[0][1]
    assert kw["inherit_env"] is False and "API_TOKEN" not in kw["env"] and "PYTHONPATH" not in kw["env"]
    assert kw["env"]["TMPDIR"] == str(tmp_path)
    assert (tmp_path / "command_0001.json").exists()
