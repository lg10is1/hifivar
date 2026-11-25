"""Two-stage launcher scheduling with genuine tiny-process producer receipts.

The scheduler and samtools metadata boundary are simulated; no native Snakemake,
caller, BAM validation or annotation/cohort computation is claimed. Producer
receipts observe real isolated Python children writing structural BGZF/TBI test
files. The independent package/Linux tests must validate the actual workflow.
"""

from __future__ import annotations

import base64
import copy
import csv
import io
import json
import os
import struct
import sys
import zlib
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

import hifivar.attempt_receipts as receipts
import hifivar.launch_inputs as launch_inputs
import hifivar.launcher as launcher
import hifivar.package_resources as resources
from hifivar.command import CommandRunner
from hifivar.config import load_config
from hifivar.exceptions import InputValidationError, OutputValidationError, WorkflowError


DEFAULT = Path(__file__).parents[2] / "src/hifivar/resources/configs/default.yaml"


def _bgzf(data: bytes) -> bytes:
    compressor = zlib.compressobj(wbits=-15)
    payload = compressor.compress(data) + compressor.flush()
    size = 18 + len(payload) + 8
    header = (b"\x1f\x8b\x08\x04" + b"\0" * 4 + b"\0\xff\x06\0BC\x02\0"
              + struct.pack("<H", size - 1))
    return header + payload + struct.pack("<II", zlib.crc32(data), len(data))


def _vcf(sample: str) -> bytes:
    return _bgzf((
        "##fileformat=VCFv4.3\n##reference=GRCh38\n##contig=<ID=chr1,length=1000>\n"
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
        "chr1\t100\tv1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL\tGT\t0/1\n"
    ).encode("utf-8"))


@pytest.fixture
def generated_config(tmp_path: Path, monkeypatch):
    (tmp_path / "ref.fa").write_text(">chr1\n" + "A" * 1000 + "\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t1000\t6\t1000\t1001\n", encoding="utf-8")
    for sample in ("S1", "S2"):
        (tmp_path / f"{sample}.bam").write_bytes(b"synthetic-BAM-process-boundary-only")
        (tmp_path / f"{sample}.bam.bai").write_bytes(b"synthetic-BAI-process-boundary-only")
    (tmp_path / "samples.tsv").write_text("sample_id\tinput\nS2\tS2.bam\nS1\tS1.bam\n", encoding="utf-8")
    (tmp_path / "regions.bed").write_text("chr1\t90\t170\tGENE1\n", encoding="utf-8")
    rows = [{"sample_id": sample, "track": "sv", "origin": "GENERATED",
             "producer": "sv.sawfish", "output_role": "vcf", "reference_build": "GRCh38"}
            for sample in ("S1", "S2")]
    annotation = {**rows[0], "variant_ids": ["v1"]}
    (tmp_path / "sources.yaml").write_text(yaml.safe_dump({
        "schema_version": 1, "cohort": rows, "annotation": [annotation],
    }), encoding="utf-8")
    user = {
        "reference": {"fasta": "ref.fa", "build": "GRCh38"},
        "samples": {"sheet": "samples.tsv"}, "runtime": {"threads": 2},
        "workflow": {"memory_mb": 4096, "downstream_source_plan": "sources.yaml"},
        "sv": {"enabled": True, "sawfish": {"enabled": True, "threads": 1, "memory_mb": 1000},
               "sniffles2": {"enabled": False}, "pbsv": {"enabled": False}, "cutesv": {"enabled": False}},
        "annotation": {"enabled": True, "overlap_enabled": True, "gene_bed": "regions.bed",
                       "gene_version": "unit-database-1", "threads": 1, "memory_mb": 1000},
        "cohort": {"enabled": True, "cohort_id": "C1", "sv": {"enabled": True, "threads": 1, "memory_mb": 1000}},
    }
    (tmp_path / "user.yaml").write_text(yaml.safe_dump(user), encoding="utf-8")
    workflow = tmp_path / "installed-workflow"
    workflow.mkdir()
    (workflow / "Snakefile").write_text("rule all:\n    input: []\n", encoding="utf-8")
    runtime = {"hifivar_version": "synthetic-test-runtime", "payload_sha256": "c" * 64}
    monkeypatch.setattr(launcher.platform, "system", lambda: "Linux")
    monkeypatch.setattr(launcher, "runtime_identity", lambda _: runtime)
    monkeypatch.setattr(launcher, "installed_workflow_root", lambda: workflow)
    monkeypatch.setattr(resources, "installed_workflow_root", lambda: workflow)
    # Only synthetic keys are inspected. Never enumerate/read the real env.
    environment = {"PATH": "synthetic-safe-path", "LANG": "C", "SystemRoot": r"C:\Windows",
                   "UNIT_TEST_API_KEY": "SYNTHETIC_NOT_A_REAL_SECRET",
                   "SNAKEMAKE_PROFILE": "synthetic-untrusted-profile", "PYTHONPATH": "synthetic-untrusted-path"}
    monkeypatch.setattr(launcher, "os", SimpleNamespace(
        environ=environment, getpid=os.getpid, walk=os.walk,
    ))
    monkeypatch.setattr(launch_inputs, "os", SimpleNamespace(environ=environment))
    receipt_proxy = SimpleNamespace(environ={}, fdopen=os.fdopen, fsync=os.fsync, link=os.link)
    monkeypatch.setattr(receipts, "os", receipt_proxy)
    config = load_config(DEFAULT, user_config=tmp_path / "user.yaml")
    return SimpleNamespace(config=config, root=tmp_path / "run", base=tmp_path,
                           receipt_proxy=receipt_proxy)


class PortableScheduler:
    """Simulate scheduling only; genuine producer command witnesses are separate."""

    def __init__(self, fixture, *, producer_exit=0, downstream_exit=0,
                 missing_receipts=False, missing_sample=False, tamper=None):
        self.fixture = fixture
        self.producer_exit, self.downstream_exit = producer_exit, downstream_exit
        self.missing_receipts, self.missing_sample = missing_receipts, missing_sample
        self.tamper = tamper
        self.calls, self.dags, self.native_outputs = [], [], []

    def require_executable(self, executable, **kwargs):
        return Path(executable)

    def _produce(self, root, environment):
        self.fixture.receipt_proxy.environ = {receipts.ATTEMPT_CONTEXT_ENV: environment[receipts.ATTEMPT_CONTEXT_ENV]}
        if self.missing_receipts:
            return
        for sample in (("S1",) if self.missing_sample else ("S1", "S2")):
            output = root / "results" / sample / "sawfish.vcf.gz"
            index = Path(f"{output}.tbi")
            outputs = {"vcf": output, "index": index}
            payloads = {str(output): base64.b64encode(_vcf(sample)).decode("ascii"),
                        str(index): base64.b64encode(_bgzf(b"TBI\x01structural-index-fixture")).decode("ascii")}
            with receipts.attempt_step("sv.sawfish", (sample,), outputs=outputs) as receipt:
                CommandRunner().run([
                    sys.executable, "-I", "-B", "-c",
                    "import json,base64,sys;from pathlib import Path;"
                    "[(Path(p).parent.mkdir(parents=True,exist_ok=True),"
                    "Path(p).write_bytes(base64.b64decode(b))) for p,b in json.loads(sys.argv[1]).items()]",
                    json.dumps(payloads),
                ], inherit_env=False, env={"SystemRoot": r"C:\Windows"} if os.name == "nt" else {}, timeout=20)
                receipt.complete(
                    outputs=outputs, tool_versions={"sawfish": "unit-process-1"},
                    validation={"validated": True, "scope": "SYNTHETIC_UNIT_BYTES_NOT_NATIVE_CALLER"},
                    metadata={"reference_build": "GRCh38", "tool_version_sources": {"sawfish": "UNIT_COMMAND_RECORDED"}},
                )
            self.native_outputs.extend(outputs.values())

    def run(self, command, **kwargs):
        command = list(command)
        self.calls.append((command, dict(kwargs)))
        stdout, code = "", 0
        if command[0] == "samtools":
            if "--version" in command:
                stdout = "samtools 1.23.1\nUsing htslib 1.23.1\n"
            elif "-H" in command:
                sample = Path(command[-1]).stem
                stdout = ("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"
                          f"@RG\tID:r\tSM:{sample}\tPL:PACBIO\n")
            elif "-X" in command:
                stdout = "0\n"
        elif "--version" in command:
            stdout = "9.26.1\n"
        else:
            assert command[-2:] == ["--", "all"]
            config_path = Path(command[command.index("--configfile") + 1])
            config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
            stage = "producer" if not self.dags else "downstream"
            self.dags.append({"stage": stage, "command": command, "config": config, "path": config_path})
            root = Path(kwargs["cwd"])
            if stage == "producer":
                if "--dry-run" not in command:
                    self._produce(root, kwargs["env"])
                code = self.producer_exit
                if self.tamper == "producer_config":
                    with config_path.open("a", encoding="utf-8") as handle:
                        handle.write("\n# altered during producer scheduler boundary\n")
                if self.tamper == "frozen_control":
                    with (root / "effective.yaml").open("a", encoding="utf-8") as handle:
                        handle.write("\n# altered frozen control\n")
                if self.tamper == "input":
                    with (self.fixture.base / "ref.fa").open("a", encoding="utf-8") as handle:
                        handle.write("\n")
            else:
                code = self.downstream_exit
                if self.tamper == "materialization":
                    target = config_path.parent / "materialization.json"
                    value = json.loads(target.read_text(encoding="utf-8"))
                    value["accuracy"] = "incorrectly-rewritten"
                    target.write_text(json.dumps(value), encoding="utf-8")
                if self.tamper == "downstream_config":
                    with config_path.open("a", encoding="utf-8") as handle:
                        handle.write("\n# altered downstream control\n")
                if self.tamper == "producer_receipt":
                    attempt = config_path.parent.parent
                    target = attempt / "receipts" / receipts.receipt_name("sv.sawfish", ("S1",))
                    with target.open("a", encoding="utf-8") as handle:
                        handle.write("\n ")
                if self.tamper == "producer_output":
                    with self.native_outputs[0].open("ab") as handle:
                        handle.write(b"tamper")
        for key, text in (("stdout_path", stdout), ("stderr_path", "")):
            if kwargs.get(key):
                Path(kwargs[key]).write_text(text, encoding="utf-8")
        return SimpleNamespace(args=tuple(command), stdout=stdout, returncode=code, duration_seconds=0.01)


def _attempt(fixture):
    attempts = list((fixture.root / "attempts").iterdir())
    assert len(attempts) == 1
    return attempts[0]


def _status(fixture):
    return json.loads((_attempt(fixture) / "status.json").read_text(encoding="utf-8"))


def test_explicit_generated_plan_runs_producers_once_then_only_consumers(generated_config):
    fixture = generated_config
    original = copy.deepcopy(fixture.config.to_dict())
    runner = PortableScheduler(fixture)
    assert launcher.WorkflowLauncher(runner).run(fixture.config, run_dir=fixture.root) == "WORKFLOW_SUCCEEDED"
    assert len(runner.dags) == 2
    producer, consumer = runner.dags
    assert producer["config"]["sv"] == original["sv"]
    assert producer["config"]["sv"]["enabled"]
    assert not producer["config"]["annotation"]["enabled"] and not producer["config"]["cohort"]["enabled"]
    assert producer["config"]["workflow"]["downstream_source_plan"] is None
    assert consumer["config"]["annotation"]["enabled"] and consumer["config"]["cohort"]["enabled"]
    assert not any(consumer["config"][name]["enabled"] for name in (
        "alignment", "small", "sv", "tr", "phasing", "assembly", "assembly_sv"))
    assert consumer["config"]["workflow"]["downstream_source_plan"] is None
    assert "--rerun-incomplete" not in consumer["command"]
    assert fixture.config.to_dict() == original
    materialization = json.loads((_attempt(fixture) / "downstream/materialization.json").read_text(encoding="utf-8"))
    assert len(materialization["source_receipts"]) == 3  # S1 annotation and S1/S2 cohort.
    assert materialization["accuracy"] == "NOT_ASSESSED"
    rows = list(csv.DictReader(io.StringIO((_attempt(fixture) / "downstream/cohort.inputs.tsv").read_text()), delimiter="\t"))
    assert [row["sample"] for row in rows] == ["S2", "S1"]
    assert _status(fixture)["producer_returncode"] == 0
    assert all(path.is_file() for path in runner.native_outputs)


def test_generated_plan_dry_run_is_plan_only_never_called_manifest(generated_config):
    fixture = generated_config
    runner = PortableScheduler(fixture)
    assert launcher.WorkflowLauncher(runner).run(fixture.config, run_dir=fixture.root, dry_run=True) == "GENERATED_SOURCE_PLAN_ONLY"
    assert len(runner.dags) == 1 and "--dry-run" in runner.dags[0]["command"]
    assert not (_attempt(fixture) / "downstream").exists()
    assert not (_attempt(fixture) / "receipts").exists()
    assert not runner.native_outputs
    status = _status(fixture)
    assert status["downstream_dag"] == "NOT_RUN_WAITING_FOR_GENERATED_RECEIPTS"
    assert status["dry_run"] is True


@pytest.mark.parametrize("failure", ["producer_exit", "missing_receipts", "missing_sample"])
def test_producer_failure_or_incomplete_receipts_stop_before_consumers(generated_config, failure):
    fixture = generated_config
    kwargs = {failure: 7 if failure == "producer_exit" else True}
    runner = PortableScheduler(fixture, **kwargs)
    with pytest.raises((InputValidationError, OutputValidationError, WorkflowError)):
        launcher.WorkflowLauncher(runner).run(fixture.config, run_dir=fixture.root)
    assert len(runner.dags) == 1
    assert not (_attempt(fixture) / "downstream").exists()
    assert _status(fixture)["status"] == "FAILED"
    assert all(path.is_file() for path in runner.native_outputs)


def test_downstream_failure_preserves_successful_producer_artifacts(generated_config):
    fixture = generated_config
    runner = PortableScheduler(fixture, downstream_exit=9)
    with pytest.raises(WorkflowError, match="Downstream Snakemake failed"):
        launcher.WorkflowLauncher(runner).run(fixture.config, run_dir=fixture.root)
    status = _status(fixture)
    assert len(runner.dags) == 2 and status["status"] == "FAILED" and status["returncode"] == 9
    assert status["producer_returncode"] == 0
    assert (_attempt(fixture) / "downstream/materialization.json").is_file()
    assert all(path.is_file() for path in runner.native_outputs)
    assert len(list((_attempt(fixture) / "receipts").iterdir())) == 2


@pytest.mark.parametrize("tamper", ["producer_config", "frozen_control", "input", "materialization",
                                  "downstream_config", "producer_receipt", "producer_output"])
def test_input_controls_materialization_and_source_tamper_fail_closed(generated_config, tamper):
    fixture = generated_config
    runner = PortableScheduler(fixture, tamper=tamper)
    with pytest.raises((InputValidationError, OutputValidationError, WorkflowError)):
        launcher.WorkflowLauncher(runner).run(fixture.config, run_dir=fixture.root)
    assert _status(fixture)["status"] == "FAILED"
    assert all(path.is_file() for path in runner.native_outputs)
    if tamper in {"producer_config", "frozen_control", "input"}:
        assert len(runner.dags) == 1, "Changed producer inputs must fail before consumer execution"


def test_scheduler_and_bam_preflight_do_not_inherit_synthetic_secret_environment(generated_config):
    fixture = generated_config
    runner = PortableScheduler(fixture)
    launcher.WorkflowLauncher(runner).run(fixture.config, run_dir=fixture.root)
    assert runner.calls
    for _, kwargs in runner.calls:
        assert kwargs.get("inherit_env") is False
        environment = kwargs["env"]
        assert "UNIT_TEST_API_KEY" not in environment
        assert "SNAKEMAKE_PROFILE" not in environment and "PYTHONPATH" not in environment
        assert "SYNTHETIC_NOT_A_REAL_SECRET" not in json.dumps(environment)
    for _, kwargs in runner.calls:
        if receipts.ATTEMPT_CONTEXT_ENV in kwargs["env"]:
            assert kwargs["env"]["TMPDIR"] == str(fixture.root / "tmp")
            assert kwargs["env"]["XDG_CACHE_HOME"] == str(fixture.root / "cache")
