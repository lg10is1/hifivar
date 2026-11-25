"""Independent adversarial checks of receipt linkage, not a Python sandbox.

All native processes are tiny isolated Python children with an explicit,
synthetic-only environment. No server, biological tool or real credential is
read. Fixtures deliberately substitute runtime identity, so these are unit
tests rather than installed-package or scientific validation.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

import hifivar.attempt_receipts as receipts
import hifivar.launcher as launcher
import hifivar.package_resources as resources
from hifivar.command import CommandRunner
from hifivar.command_observation import observe_commands
from hifivar.exceptions import InputValidationError, OutputValidationError


@pytest.fixture
def isolated_attempt(tmp_path: Path, monkeypatch):
    """Construct only owned, synthetic launcher metadata and no inherited env."""
    root = tmp_path / "run"
    root.mkdir()
    config = {
        "small": {"enabled": True},
        "annotation": {
            "enabled": True, "annovar_enabled": True,
            "vep_enabled": True, "overlap_enabled": True,
        },
        "cohort": {
            "enabled": True, "small_variants": {"enabled": True},
            "sv": {"enabled": True}, "tr": {"enabled": True},
        },
    }
    sample_sheet = "sample_id\tinput\nS1\t/external/S1.bam\nS2\t/external/S2.bam\n"
    (root / "effective.yaml").write_text(yaml.safe_dump(config), encoding="utf-8")
    (root / "samples.tsv").write_text(sample_sheet, encoding="utf-8")
    runtime = {"payload_sha256": "a" * 64, "hifivar_version": "isolated-unit-fixture"}
    monkeypatch.setattr(launcher, "runtime_identity", lambda _: runtime)
    monkeypatch.setattr(resources, "installed_workflow_root", lambda: root)
    manifest = {
        "identity": {
            "schema": 2, "run_dir": str(root), "config": config,
            "sample_sheet": sample_sheet, "runtime": runtime,
        },
        "control_sha256": {
            name: receipts.artifact_identity(root / name)["sha256"]
            for name in ("effective.yaml", "samples.tsv")
        },
    }
    (root / "launch_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    attempt = root / "attempts" / "attempt_a"
    attempt.mkdir(parents=True)
    (attempt / "status.json").write_text(
        json.dumps({"status": "DAG_RUNNING", "dry_run": False}), encoding="utf-8",
    )
    context = receipts.create_attempt_context(root, attempt, dry_run=False)
    # Proxy only this module's environment lookup, without reading os.environ.
    proxy = SimpleNamespace(
        environ={receipts.ATTEMPT_CONTEXT_ENV: str(context)},
        fdopen=os.fdopen, fsync=os.fsync, link=os.link,
    )
    monkeypatch.setattr(receipts, "os", proxy)
    return root, attempt, context, proxy


def _outputs(path: Path) -> dict[str, Path]:
    """Declare the complete small-variant output roles; bytes stay synthetic."""
    return {
        "vcf": path, "index": path.with_name(path.name + ".tbi"),
        "gvcf": path.with_name(path.name + ".gvcf"),
        "gvcf_index": path.with_name(path.name + ".gvcf.tbi"),
    }


def _manual_outputs(path: Path) -> None:
    """Only negative tests manually write fake files, never positive tool proof."""
    path.parent.mkdir(parents=True, exist_ok=True)
    for value in _outputs(path).values():
        value.write_bytes(b"manually-created-not-native-proof")


def _native_write(path: Path, *, returncode: int = 0, dry_run: bool = False,
                  secret: str | None = None):
    path.parent.mkdir(parents=True, exist_ok=True)
    environment = {"SystemRoot": r"C:\Windows"} if os.name == "nt" else {}
    if secret is not None:
        environment["UNIT_TEST_API_KEY"] = secret
    return CommandRunner().run(
        [sys.executable, "-I", "-B", "-c",
         "from pathlib import Path;import sys;"
         "p=Path(sys.argv[1]);p.write_bytes(b'synthetic-output');"
         "[p.with_name(p.name+s).write_bytes(b'synthetic-companion') "
         "for s in ('.tbi','.gvcf','.gvcf.tbi')];"
         "print('synthetic-stdout');print('synthetic-stderr',file=sys.stderr);"
         "sys.exit(int(sys.argv[2]))", str(path), str(returncode)],
        inherit_env=False, env=environment, check=False, dry_run=dry_run,
        timeout=20,
    )


def _complete(receipt, output: Path, **kwargs):
    receipt.complete(
        outputs=_outputs(output), tool_versions={"unit-process": "1.0"},
        validation={"validated": True, "scope": "SYNTHETIC_BYTES_ONLY"}, **kwargs,
    )


@pytest.mark.parametrize("area", ["effective.yaml", "logs/fake", "attempts/attempt_a/fake"])
def test_control_log_or_attempt_files_cannot_be_declared_as_outputs(isolated_attempt, area):
    root, attempt, _, _ = isolated_attempt
    with pytest.raises(InputValidationError, match="results/work"):
        with receipts.attempt_step("small", ("S1",), outputs=_outputs(root / area)):
            pytest.fail("An unowned output reached execution")
    assert not (attempt / "receipts").exists()


@pytest.mark.parametrize("execution_kind,returncode,dry_run", [
    ("NATIVE", 9, False), ("NATIVE", 0, True), ("IN_PROCESS", 0, False),
])
def test_nonzero_planned_or_inprocess_caller_cannot_certify_native_success(
    isolated_attempt, execution_kind, returncode, dry_run,
):
    root, attempt, _, _ = isolated_attempt
    output = root / "results" / "failure.bin"
    with pytest.raises(OutputValidationError):
        with receipts.attempt_step("small", ("S1",), outputs=_outputs(output)) as receipt:
            _native_write(output, returncode=returncode, dry_run=dry_run)
            if dry_run:
                _manual_outputs(output)
            _complete(receipt, output, execution_kind=execution_kind)
    assert not (attempt / "receipts").exists()


@pytest.mark.parametrize("probe", ["--version", "--help", "-h"])
def test_only_probe_command_cannot_certify_new_manually_created_artifact(isolated_attempt, probe):
    root, attempt, _, _ = isolated_attempt
    output = root / "results" / "probe.bin"
    with pytest.raises(OutputValidationError, match="probes alone"):
        with receipts.attempt_step("small", ("S1",), outputs=_outputs(output)) as receipt:
            CommandRunner().run(
                [sys.executable, probe], inherit_env=False,
                env={"SystemRoot": r"C:\Windows"} if os.name == "nt" else {},
            )
            _manual_outputs(output)
            _complete(receipt, output)
    assert not (attempt / "receipts").exists()


def test_command_witness_copy_is_readonly_and_does_not_retain_capture(isolated_attempt):
    root, _, _, _ = isolated_attempt
    with observe_commands() as observed:
        _native_write(root / "work" / "observer.bin")
    witness = tuple(observed)[0]
    display = witness.to_dict()
    assert "stdout" not in display and "stderr" not in display and "env" not in display
    display["returncode"] = 19
    assert witness.to_dict()["returncode"] == 0
    assert not hasattr(observed, "append")
    with pytest.raises((AttributeError, TypeError)):
        witness._record = {}


def test_semantic_secrets_are_redacted_across_all_receipt_fields(isolated_attempt):
    root, attempt, _, _ = isolated_attempt
    secret = "UNIT_TEST_ONLY_FAKE_SECRET_5bd6d0"
    output = root / "results" / "redacted.bin"
    with receipts.attempt_step("small", ("S1",), outputs=_outputs(output)) as receipt:
        _native_write(output, secret=secret)
        receipt.complete(
            outputs=_outputs(output),
            tool_versions={"unit-process": "1.0", "api_key": secret},
            validation={"validated": True, "nested": {"access_token": secret}},
            metadata={"password": secret, "safe": "retained"},
        )
    result = receipts.read_step_receipt(root, attempt.name, "small", ("S1",))
    serialized = json.dumps(result)
    assert secret not in serialized
    assert result["validation"]["nested"]["access_token"] == "***"
    assert result["tool_versions"]["api_key"] == "***"
    assert result["metadata"]["safe"] == "retained"


def test_output_roles_cannot_alias_one_inode(isolated_attempt):
    root, attempt, _, _ = isolated_attempt
    first, second = root / "results" / "a.bin", root / "results" / "b.bin"
    outputs = _outputs(first)
    outputs["index"] = second
    with pytest.raises(OutputValidationError, match="hardlink"):
        with receipts.attempt_step("small", ("S1",), outputs=outputs) as receipt:
            _native_write(first)
            try:
                os.link(first, second)
            except OSError as error:
                pytest.skip(f"Hardlink fixture unavailable: {type(error).__name__}")
            receipt.complete(
                outputs=outputs, tool_versions={"unit": "1"},
                validation={"validated": True},
            )
    assert not (attempt / "receipts").exists()


def test_changed_sample_control_mid_step_is_rejected_before_sealing(isolated_attempt):
    root, attempt, _, _ = isolated_attempt
    output = root / "results" / "control.bin"
    with pytest.raises(InputValidationError, match="changed"):
        with receipts.attempt_step("small", ("S1",), outputs=_outputs(output)) as receipt:
            _native_write(output)
            (root / "samples.tsv").write_text("sample_id\tinput\nS2\tx\n", encoding="utf-8")
            _complete(receipt, output)
    assert not (attempt / "receipts").exists()


def test_attempt_failure_mid_step_does_not_become_completion(isolated_attempt):
    root, attempt, _, _ = isolated_attempt
    output = root / "results" / "failed-attempt.bin"
    with pytest.raises(InputValidationError, match="active native"):
        with receipts.attempt_step("small", ("S1",), outputs=_outputs(output)) as receipt:
            _native_write(output)
            (attempt / "status.json").write_text(
                json.dumps({"status": "FAILED", "dry_run": False}), encoding="utf-8",
            )
            _complete(receipt, output)
    assert not (attempt / "receipts").exists()


def test_copy_of_valid_receipt_cannot_link_to_another_attempt(isolated_attempt):
    root, attempt, _, _ = isolated_attempt
    output = root / "results" / "original.bin"
    with receipts.attempt_step("small", ("S1",), outputs=_outputs(output)) as receipt:
        _native_write(output)
        _complete(receipt, output)
    second = root / "attempts" / "attempt_b"
    second.mkdir()
    receipts.create_attempt_context(root, second, dry_run=False)
    target = second / "receipts" / receipts.receipt_name("small", ("S1",))
    target.parent.mkdir()
    target.write_bytes((attempt / "receipts" / target.name).read_bytes())
    with pytest.raises(InputValidationError, match="selected native attempt"):
        receipts.read_step_receipt(root, second.name, "small", ("S1",))


def test_completion_twice_preserves_first_receipt_bytes(isolated_attempt):
    root, attempt, _, _ = isolated_attempt
    output = root / "results" / "once.bin"
    with receipts.attempt_step("small", ("S1",), outputs=_outputs(output)) as receipt:
        _native_write(output)
        _complete(receipt, output)
        target = attempt / "receipts" / receipts.receipt_name("small", ("S1",))
        before = target.read_bytes()
        with pytest.raises(OutputValidationError, match="twice"):
            _complete(receipt, output)
        assert target.read_bytes() == before


def test_no_context_has_no_launcher_success_claim(isolated_attempt):
    root, attempt, _, proxy = isolated_attempt
    proxy.environ = {}
    output = root / "results" / "standalone.bin"
    with receipts.attempt_step("small", ("S1",)) as receipt:
        output.parent.mkdir(parents=True)
        output.write_bytes(b"independent-module-output")
        _complete(receipt, output)
    assert not (attempt / "receipts").exists()


def test_old_observation_cannot_replace_this_steps_native_execution(isolated_attempt):
    """Wrongly wiring a real but unrelated observation must not prove this step."""
    root, attempt, _, _ = isolated_attempt
    with observe_commands() as old_observation:
        _native_write(root / "work" / "previous-step.bin")
    output = root / "results" / "unexecuted-step.bin"
    with pytest.raises(OutputValidationError):
        with receipts.attempt_step("small", ("S1",), outputs=_outputs(output)) as receipt:
            _manual_outputs(output)
            receipt.commands = old_observation
            _complete(receipt, output)
    assert not (attempt / "receipts").exists()


@pytest.mark.parametrize("corruption", ["returncode_bool", "not_executed", "duration_nan", "validation_false"])
def test_corrupt_persisted_execution_record_is_not_accepted(isolated_attempt, corruption):
    root, attempt, _, _ = isolated_attempt
    output = root / "results" / "persisted.bin"
    with receipts.attempt_step("small", ("S1",), outputs=_outputs(output)) as receipt:
        _native_write(output)
        _complete(receipt, output)
    path = attempt / "receipts" / receipts.receipt_name("small", ("S1",))
    value = json.loads(path.read_text(encoding="utf-8"))
    if corruption == "returncode_bool":
        value["commands"][0]["returncode"] = True
    elif corruption == "not_executed":
        value["commands"][0]["executed"] = False
    elif corruption == "duration_nan":
        value["commands"][0]["duration_seconds"] = float("nan")
    else:
        value["validation"]["validated"] = False
    path.write_text(json.dumps(value), encoding="utf-8")
    with pytest.raises((InputValidationError, OutputValidationError)):
        receipts.read_step_receipt(root, attempt.name, "small", ("S1",))


def test_reading_persisted_receipt_rechecks_owned_output_area(isolated_attempt):
    """Reject a valid hash whose stored output path points into frozen controls."""
    root, attempt, _, _ = isolated_attempt
    output = root / "results" / "owned.bin"
    with receipts.attempt_step("small", ("S1",), outputs=_outputs(output)) as receipt:
        _native_write(output)
        _complete(receipt, output)
    path = attempt / "receipts" / receipts.receipt_name("small", ("S1",))
    value = json.loads(path.read_text(encoding="utf-8"))
    value["outputs"]["vcf"] = receipts.artifact_identity(root / "effective.yaml")
    path.write_text(json.dumps(value), encoding="utf-8")
    with pytest.raises((InputValidationError, OutputValidationError)):
        receipts.read_step_receipt(root, attempt.name, "small", ("S1",))


def test_failure_after_one_sealed_step_does_not_rewrite_attempt_status(isolated_attempt):
    """Partial native success remains scoped, never overriding overall failure."""
    root, attempt, _, _ = isolated_attempt
    output = root / "results" / "partial.bin"
    with receipts.attempt_step("small", ("S1",), outputs=_outputs(output)) as receipt:
        _native_write(output)
        _complete(receipt, output)
    failed = {"status": "FAILED", "dry_run": False, "exit_code": 9}
    (attempt / "status.json").write_text(json.dumps(failed), encoding="utf-8")
    value = receipts.read_step_receipt(root, attempt.name, "small", ("S1",))
    assert value["status"] == "COMPLETE"
    assert value["scientific_accuracy"] == "NOT_ASSESSED"
    assert json.loads((attempt / "status.json").read_text(encoding="utf-8")) == failed


@pytest.mark.parametrize("missing_role", ["vcf", "index", "gvcf", "gvcf_index"])
def test_missing_required_small_variant_output_role_cannot_seal_completion(isolated_attempt, missing_role):
    root, attempt, _, _ = isolated_attempt
    outputs = _outputs(root / "results" / "incomplete.bin")
    del outputs[missing_role]
    with pytest.raises((InputValidationError, OutputValidationError)):
        with receipts.attempt_step("small", ("S1",), outputs=outputs) as receipt:
            _native_write(root / "results" / "incomplete.bin")
            receipt.complete(
                outputs=outputs, tool_versions={"unit-process": "1"},
                validation={"validated": True, "scope": "SYNTHETIC_BYTES_ONLY"},
            )
    assert not (attempt / "receipts").exists()


def test_persisted_receipt_missing_required_role_is_rejected(isolated_attempt):
    root, attempt, _, _ = isolated_attempt
    output = root / "results" / "complete.bin"
    with receipts.attempt_step("small", ("S1",), outputs=_outputs(output)) as receipt:
        _native_write(output)
        _complete(receipt, output)
    path = attempt / "receipts" / receipts.receipt_name("small", ("S1",))
    value = json.loads(path.read_text(encoding="utf-8"))
    del value["outputs"]["gvcf_index"]
    path.write_text(json.dumps(value), encoding="utf-8")
    with pytest.raises((InputValidationError, OutputValidationError)):
        receipts.read_step_receipt(root, attempt.name, "small", ("S1",))
