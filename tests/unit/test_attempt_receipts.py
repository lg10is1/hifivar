"""Attempt receipts use real tiny processes, never fabricated tool success."""

from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

import pytest
import yaml

from hifivar import launcher
from hifivar.attempt_receipts import (
    ATTEMPT_CONTEXT_ENV, artifact_identity, attempt_step, create_attempt_context,
    read_step_receipt, receipt_name,
)
from hifivar.command import CommandRunner
from hifivar.command_observation import observe_commands
from hifivar.exceptions import CommandExecutionError, InputValidationError, OutputValidationError


@pytest.fixture
def owned(tmp_path, monkeypatch):
    root = tmp_path / "run"
    root.mkdir()
    source = tmp_path / "source.vcf"
    source.write_text("synthetic-source-bytes\n", encoding="utf-8")
    source_manifest = root / "annotation.inputs.tsv"
    source_manifest.write_text("sample\tvariant_category\tsource_vcf\tsource_tool\tsource_variant_ids\n"
        + "".join(f"S1\t{category}\t{source}\tunit-source\tv1\n" for category in ("sv", "tr")), encoding="utf-8")
    config = {
        "small": {"enabled": True},
        "annotation": {"enabled": True, "overlap_enabled": True, "vep_enabled": True,
                       "input_manifest": str(source_manifest)},
        "cohort": {"enabled": True, "sv": {"enabled": True}},
    }
    sheet = "sample_id\tinput\nS1\t/external/S1.bam\nS2\t/external/S2.bam\n"
    (root / "effective.yaml").write_text(yaml.safe_dump(config), encoding="utf-8")
    (root / "samples.tsv").write_text(sheet, encoding="utf-8")
    runtime = {"payload_sha256": "a" * 64, "hifivar_version": "test-runtime"}
    monkeypatch.setattr(launcher, "runtime_identity", lambda workflow: runtime)
    # Isolated synthetic metadata fixture; not an installed-package gate.
    import hifivar.package_resources
    monkeypatch.setattr(hifivar.package_resources, "installed_workflow_root", lambda: root)
    manifest = {
        "identity": {"schema": 2, "run_dir": str(root), "config": config, "sample_sheet": sheet, "runtime": runtime,
                     "inputs": [artifact_identity(source)]},
        "control_sha256": {name: artifact_identity(root / name)["sha256"] for name in ("effective.yaml", "samples.tsv", "annotation.inputs.tsv")},
    }
    (root / "launch_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    attempt = root / "attempts" / "attempt_a"
    attempt.mkdir(parents=True)
    (attempt / "status.json").write_text(json.dumps({"status": "DAG_RUNNING", "dry_run": False}), encoding="utf-8")
    context = create_attempt_context(root, attempt, dry_run=False)
    monkeypatch.setenv(ATTEMPT_CONTEXT_ENV, str(context))
    return root, attempt, context


def role_paths(path):
    return {"vcf": path, "index": Path(str(path) + ".tbi"),
            "gvcf": Path(str(path) + ".g.vcf.gz"),
            "gvcf_index": Path(str(path) + ".g.vcf.gz.tbi")}


def native_write(path: Path, *, code=0, secret=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    # Only this synthetic key can reach the child; no inherited credentials.
    environment = {}
    if os.name == "nt":
        system_root = os.environ.get("SystemRoot") or os.environ.get("SYSTEMROOT")
        if system_root:
            environment["SystemRoot"] = system_root
    if secret:
        environment["TEST_API_KEY"] = secret
    return CommandRunner().run([
        sys.executable, "-I", "-B", "-c",
        "import pathlib,sys;[(pathlib.Path(sys.argv[1]+s).write_bytes(b'validated-test-output')) for s in ('','.tbi','.g.vcf.gz','.g.vcf.gz.tbi')];sys.exit(int(sys.argv[2]))",
        str(path), str(code),
    ], env=environment, inherit_env=False)


def finish(receipt, output, **kwargs):
    receipt.complete(outputs=role_paths(output), tool_versions={"synthetic-process": "1"}, validation={"validated": True, "scope": "TEST_BYTES_ONLY"}, **kwargs)


def test_success_is_linked_to_exact_attempt_and_streamed_bytes(owned):
    root, attempt, _ = owned
    output = root / "results" / "small.vcf.gz"
    with attempt_step("small", ("S1",), outputs=role_paths(output)) as receipt:
        native_write(output)
        finish(receipt, output)
    value = read_step_receipt(root, attempt.name, "small", ("S1",))
    assert value["attempt_id"] == value["producer_attempt_id"] == "attempt_a"
    assert value["outputs"]["vcf"] == artifact_identity(output)
    assert value["commands"][0]["executed"] is True
    assert value["commands"][0]["returncode"] == 0
    assert value["scientific_accuracy"] == "NOT_ASSESSED"
    assert value["provenance_class"] == "THIS_ATTEMPT_EXECUTED_NOT_REUSED"
    assert "stdout" not in value["commands"][0] and "stderr" not in value["commands"][0]


def test_without_launcher_context_no_metadata_is_created(tmp_path, monkeypatch):
    monkeypatch.delenv(ATTEMPT_CONTEXT_ENV, raising=False)
    output = tmp_path / "independent-output"
    output.write_text("pre-existing standalone workflow artifact", encoding="utf-8")
    with attempt_step("small", ("S1",)) as receipt:
        finish(receipt, output)
    assert set(tmp_path.iterdir()) == {output}


def test_no_commands_cannot_certify_a_native_step(owned):
    root, attempt, _ = owned
    output = root / "results" / "small.vcf.gz"
    with pytest.raises(OutputValidationError, match="observed CommandRunner"):
        with attempt_step("small", ("S1",), outputs=role_paths(output)) as receipt:
            output.parent.mkdir()
            output.write_bytes(b"not-native-proof")
            finish(receipt, output)
    assert not (attempt / "receipts").exists()


def test_version_probe_alone_is_not_an_analysis(owned):
    root, attempt, _ = owned
    output = root / "results" / "small.vcf.gz"
    with pytest.raises(OutputValidationError, match="probes alone"):
        with attempt_step("small", ("S1",), outputs=role_paths(output)) as receipt:
            CommandRunner().run([sys.executable, "--version"], inherit_env=False)
            output.parent.mkdir()
            output.write_bytes(b"old-or-synthetic")
            finish(receipt, output)
    assert not (attempt / "receipts").exists()


def test_failed_command_preserves_error_and_cannot_emit_receipt(owned):
    root, attempt, _ = owned
    output = root / "results" / "small.vcf.gz"
    with pytest.raises(CommandExecutionError):
        with attempt_step("small", ("S1",), outputs=role_paths(output)) as receipt:
            native_write(output, code=7)
            finish(receipt, output)
    assert output.exists()
    assert not (attempt / "receipts").exists()


def test_unsealed_owned_bridge_is_not_a_success(owned):
    root, _, _ = owned
    with pytest.raises(OutputValidationError, match="did not seal"):
        with attempt_step("small", ("S1",), outputs={"vcf": root / "results" / "x"}):
            pass


def test_preexisting_output_is_rejected_before_execution(owned):
    root, _, _ = owned
    output = root / "results" / "x"
    output.parent.mkdir()
    output.write_bytes(b"old")
    with pytest.raises(OutputValidationError, match="must be fresh"):
        with attempt_step("small", ("S1",), outputs=role_paths(output)):
            pytest.fail("must not execute")
    assert output.read_bytes() == b"old"


def test_missing_validation_and_changed_output_set_fail(owned):
    root, _, _ = owned
    output = root / "results" / "x"
    with pytest.raises(OutputValidationError, match="explicit outputs"):
        with attempt_step("small", ("S1",), outputs=role_paths(output)) as receipt:
            native_write(output)
            receipt.complete(outputs=role_paths(output), tool_versions={"x": "1"}, validation={"validated": False})
    other = root / "results" / "y"
    with pytest.raises(OutputValidationError, match="fresh pre-execution"):
        with attempt_step("small", ("S1",), outputs=role_paths(other)) as receipt:
            native_write(other)
            finish(receipt, output)


@pytest.mark.parametrize("change", ["effective.yaml", "samples.tsv", "launch_manifest.json"])
def test_frozen_control_changes_are_rejected(owned, change):
    root, _, _ = owned
    with (root / change).open("a", encoding="utf-8") as stream:
        stream.write(" ")
    with pytest.raises(InputValidationError, match="changed"):
        with attempt_step("small", ("S1",), outputs={"vcf": root / "results" / "x"}):
            pytest.fail("changed controls")


def test_runtime_drift_is_rejected_before_native_execution(owned, monkeypatch):
    root, _, _ = owned
    monkeypatch.setattr(launcher, "runtime_identity", lambda _: {"payload_sha256": "different"})
    with pytest.raises(InputValidationError, match="runtime differs"):
        with attempt_step("small", ("S1",), outputs={"vcf": root / "results" / "x"}):
            pytest.fail("drift")


@pytest.mark.parametrize("state", ["FAILED", "WORKFLOW_SUCCEEDED", "PREFLIGHT_RUNNING"])
def test_inactive_attempts_cannot_emit_receipts(owned, state):
    root, attempt, _ = owned
    (attempt / "status.json").write_text(json.dumps({"status": state, "dry_run": False}), encoding="utf-8")
    with pytest.raises(InputValidationError, match="active native"):
        with attempt_step("small", ("S1",), outputs={"vcf": root / "results" / "x"}):
            pytest.fail("inactive")


@pytest.mark.parametrize("track,samples", [("tr", ("S1",)), ("other", ("S1",)), ("small", ("unknown",)), ("small", ("S1", "S1")), ("cohort.sv", ("S2", "S1"))])
def test_disabled_unknown_duplicate_or_wrong_cohort_order_is_rejected(owned, track, samples):
    root, _, _ = owned
    with pytest.raises(InputValidationError):
        with attempt_step(track, samples, outputs={"vcf": root / "results" / "x"}):
            pytest.fail("invalid unit")


def test_native_output_mutation_after_completion_is_detected(owned):
    root, attempt, _ = owned
    output = root / "results" / "x"
    with attempt_step("small", ("S1",), outputs=role_paths(output)) as receipt:
        native_write(output)
        finish(receipt, output)
    output.write_bytes(b"changed")
    with pytest.raises(OutputValidationError, match="bytes differ"):
        read_step_receipt(root, attempt.name, "small", ("S1",))


def test_wrong_attempt_or_sample_does_not_find_an_old_receipt(owned):
    root, attempt, _ = owned
    output = root / "results" / "x"
    with attempt_step("small", ("S1",), outputs=role_paths(output)) as receipt:
        native_write(output)
        finish(receipt, output)
    assert read_step_receipt(root, "attempt_other", "small", ("S1",)) is None
    assert read_step_receipt(root, attempt.name, "small", ("S2",)) is None


def test_annotation_categories_have_distinct_explicit_units(owned):
    root, attempt, _ = owned
    for category in ("sv", "tr"):
        output = root / "results" / (category + ".overlap.tsv")
        with attempt_step("annotation.overlap", ("S1",), outputs={"tsv": output}, unit_id=category,
                          inputs={"S1": {"vcf": root.parent / "source.vcf", "index": None}},
                          input_manifest=root / "annotation.inputs.tsv") as receipt:
            output.parent.mkdir(exist_ok=True)
            output.write_text("valid-header-only\n", encoding="utf-8")
            receipt.complete(outputs={"tsv": output}, tool_versions={"hifivar": "test"}, validation={"validated": True}, execution_kind="IN_PROCESS")
    assert receipt_name("annotation.overlap", ("S1",), "sv") != receipt_name("annotation.overlap", ("S1",), "tr")
    assert read_step_receipt(root, attempt.name, "annotation.overlap", ("S1",), unit_id="sv")["unit_id"] == "sv"


def test_command_observation_is_nested_and_has_no_raw_capture(owned):
    root, _, _ = owned
    with observe_commands() as outer:
        native_write(root / "results" / "a")
        with observe_commands() as inner:
            native_write(root / "results" / "b")
        native_write(root / "results" / "c")
    assert len(outer) == 2 and len(inner) == 1
    assert all("stdout" not in record.to_dict() for record in outer)
    assert not hasattr(outer, "append")


def test_sensitive_metadata_is_redacted_without_changing_native_artifacts(owned):
    root, attempt, _ = owned
    output = root / "results" / "x"
    with attempt_step("small", ("S1",), outputs=role_paths(output)) as receipt:
        native_write(output)
        finish(receipt, output, metadata={"api_key": "synthetic-never-real-secret", "safe": "retained"})
    value = read_step_receipt(root, attempt.name, "small", ("S1",))
    assert value["metadata"] == {"api_key": "***", "safe": "retained"}
    assert "synthetic-never-real-secret" not in json.dumps(value)


def test_duplicate_receipt_never_overwrites_original(owned):
    root, attempt, _ = owned
    output = root / "results" / "x"
    with attempt_step("small", ("S1",), outputs=role_paths(output)) as receipt:
        native_write(output)
        finish(receipt, output)
    path = attempt / "receipts" / receipt_name("small", ("S1",))
    before = path.read_bytes()
    other = root / "results" / "other"
    with pytest.raises(OutputValidationError, match="already exists"):
        with attempt_step("small", ("S1",), outputs=role_paths(other)) as receipt:
            native_write(other)
            finish(receipt, other)
    assert path.read_bytes() == before


def test_dry_run_context_cannot_become_native_completion(owned, monkeypatch):
    root, attempt, _ = owned
    second = root / "attempts" / "attempt_dry"
    second.mkdir()
    (second / "status.json").write_text(json.dumps({"status": "DAG_RUNNING", "dry_run": True}), encoding="utf-8")
    context = create_attempt_context(root, second, dry_run=True)
    monkeypatch.setenv(ATTEMPT_CONTEXT_ENV, str(context))
    with pytest.raises(InputValidationError, match="active native"):
        with attempt_step("small", ("S1",), outputs={"vcf": root / "results" / "x"}):
            pytest.fail("dry-run")


def test_context_creation_does_not_overwrite(owned):
    root, attempt, context = owned
    before = context.read_bytes()
    with pytest.raises(OutputValidationError, match="already exists"):
        create_attempt_context(root, attempt, dry_run=False)
    assert context.read_bytes() == before


def test_fake_dictionary_is_not_a_command_witness(owned):
    root, _, _ = owned
    output = root / "results" / "x"
    with pytest.raises(OutputValidationError, match="originate"):
        with attempt_step("small", ("S1",), outputs=role_paths(output)) as receipt:
            native_write(output)
            receipt.commands = [{"args": ["native", "run"], "executed": True, "returncode": 0, "duration_seconds": 1}]
            finish(receipt, output)


def test_native_caller_cannot_be_marked_in_process(owned):
    root, _, _ = owned
    output = root / "results" / "x"
    with pytest.raises(OutputValidationError, match="in-process"):
        with attempt_step("small", ("S1",), outputs=role_paths(output)) as receipt:
            output.parent.mkdir()
            output.write_bytes(b"fake")
            finish(receipt, output, execution_kind="IN_PROCESS")
