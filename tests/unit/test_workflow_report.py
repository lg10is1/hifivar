"""Honest local report adaptation from explicit launcher attempts only."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest
import yaml

from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.workflow_report import build_launcher_report


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False), encoding="utf-8")


def _status(state="WORKFLOW_SUCCEEDED", **overrides) -> dict:
    value = {
        "started": "2026-10-16T00:00:00Z", "finished": "2026-10-16T01:00:00Z",
        "status": state, "dry_run": False, "returncode": 0,
        "input_immutability": "PASS", "snakemake_version": "9.26.1",
    }
    value.update(overrides)
    return value


def _owned_run(tmp_path: Path, *, state="WORKFLOW_SUCCEEDED", overrides=None, config=None, attempt="attempt_a") -> Path:
    root = tmp_path / "原始运行"
    root.mkdir()
    cfg = config or {
        "reference": {"build": "GRCh38", "fasta": "/data/private/ref.fa"},
        "small": {"enabled": True}, "sv": {"enabled": False},
        "service_secret": "not-for-rendering",
    }
    sheet = "sample_id\tinput\tinput_type\nHG002\t/data/private/HG002.bam\tbam\n"
    effective = root / "effective.yaml"
    effective.write_text(yaml.safe_dump(cfg, allow_unicode=True), encoding="utf-8")
    (root / "samples.tsv").write_text(sheet, encoding="utf-8")
    identity = {
        "schema": 2, "kind": "hifivar_two_stage_launcher", "run_dir": str(root.absolute()),
        "config": cfg, "sample_sheet": sheet,
        "runtime": {"hifivar_version": "0.1.0rc5", "payload_sha256": "a" * 64},
        "inputs": [{"path": "/data/private/ref.fa", "sha256": "b" * 64}],
    }
    _write_json(root / "launch_manifest.json", {
        "identity": identity,
        "control_sha256": {name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in ("effective.yaml", "samples.tsv")},
    })
    if attempt is not None:
        if state is None:
            (root / "attempts" / attempt).mkdir(parents=True)
        else:
            _write_json(root / "attempts" / attempt / "status.json", _status(state, **(overrides or {})))
    return root


def _payload(result: dict) -> dict:
    return json.loads(Path(result["output_paths"]["json"]).read_text(encoding="utf-8"))


def _tracks(payload: dict) -> dict:
    return {track["name"]: track for track in payload["tracks"]}


def _fingerprint(root: Path) -> dict[str, str]:
    return {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest() for path in root.rglob("*") if path.is_file()}


def test_successful_scheduler_is_distinct_from_missing_per_track_receipts(tmp_path):
    run = _owned_run(tmp_path)
    before = _fingerprint(run)
    result = build_launcher_report(run, tmp_path / "独立报告")
    payload = _payload(result)
    tracks = _tracks(payload)
    assert result["execution_status"] == "WORKFLOW_SUCCEEDED"
    assert result["status"] == "PARTIAL"
    assert tracks["launcher_execution"]["status"] == "COMPLETE"
    assert tracks["small"]["status"] == "PARTIAL"
    assert tracks["assembly"]["status"] == "DISABLED"
    assert tracks["benchmark"]["status"] == "DISABLED"
    assert payload["provenance"]["scientific_accuracy"] == "NOT_ASSESSED"
    assert payload["reference"]["recorded_sha256"] == "b" * 64
    assert payload["reference"]["checksum_verification"] == "NOT_RUN_METADATA_ONLY"
    assert payload["config"]["service_secret"] == "***"
    assert payload["git_sha"] is None
    assert payload["project_version"] == "0.1.0rc5"
    assert payload["provenance"]["publication"] == "LOCAL_ONLY_NOT_PUBLISHED"
    assert _fingerprint(run) == before
    assert set(result["output_paths"]) == {"json", "yaml", "md", "html"}
    assert all(Path(path).is_file() for path in result["output_paths"].values())
    assert yaml.safe_load(Path(result["output_paths"]["yaml"]).read_text(encoding="utf-8")) == payload
    assert "Accuracy is NOT_ASSESSED" in Path(result["output_paths"]["md"]).read_text(encoding="utf-8")
    html = Path(result["output_paths"]["html"]).read_text(encoding="utf-8")
    assert "Accuracy is NOT_ASSESSED" in html and "<script" not in html


def test_new_failed_attempt_cannot_be_hidden_by_older_success_or_vcfs(tmp_path):
    run = _owned_run(tmp_path, attempt="attempt_old")
    _write_json(run / "attempts" / "attempt_new" / "status.json", _status("FAILED", returncode=1, error="tool failure"))
    stale = run / "results" / "small" / "HG002.small.vcf.gz"
    stale.parent.mkdir(parents=True)
    stale.write_bytes(b"old-success-output")
    before = _fingerprint(run)
    result = build_launcher_report(run, tmp_path / "failed-report", "attempt_new")
    payload = _payload(result)
    assert result["status"] == "FAILED"
    assert result["execution_status"] == "FAILED"
    assert payload["provenance"]["attempt_id"] == "attempt_new"
    assert _tracks(payload)["small"]["artifacts"] == []
    assert _fingerprint(run) == before


def test_multiple_attempts_require_explicit_identity_not_newest_name_or_mtime(tmp_path):
    run = _owned_run(tmp_path, attempt="z_latest_name")
    _write_json(run / "attempts" / "a_newest_time" / "status.json", _status("FAILED", returncode=1))
    output = tmp_path / "not-created"
    with pytest.raises(InputValidationError, match="MULTIPLE_ATTEMPTS"):
        build_launcher_report(run, output)
    assert not output.exists()
    result = build_launcher_report(run, output, "z_latest_name")
    assert result["execution_status"] == "WORKFLOW_SUCCEEDED"


@pytest.mark.parametrize("state", ["DRY_RUN_PASS", "RAW_INPUT_PLAN_ONLY", "GENERATED_SOURCE_PLAN_ONLY"])
def test_dry_run_and_raw_planning_are_not_analysis_completion(tmp_path, state):
    run = _owned_run(tmp_path, state=state, overrides={"dry_run": True})
    result = build_launcher_report(run, tmp_path / "planned-report")
    assert result["status"] == "NOT_RUN"
    assert result["execution_status"] == state
    assert _tracks(_payload(result))["small"]["status"] == "NOT_RUN"


@pytest.mark.parametrize("attempt", [None, "attempt_without_status"])
def test_no_attempt_or_no_status_is_not_run_even_with_existing_outputs(tmp_path, attempt):
    run = _owned_run(tmp_path, state=None, attempt=attempt)
    (run / "old_vcf.vcf.gz").write_bytes(b"not-authority")
    result = build_launcher_report(run, tmp_path / "missing-status-report")
    assert result["status"] == "NOT_RUN"
    assert result["execution_status"] == "NO_STATUS_RECORD"
    assert _tracks(_payload(result))["small"]["status"] == "NOT_RUN"


@pytest.mark.parametrize("state", ["DAG_RUNNING", "PREFLIGHT_RUNNING", "ALIGNMENT_PREPARATION", "INTERRUPTED_LOCK_RETAINED",
                                    "PRODUCER_DAG_RUNNING", "DOWNSTREAM_PREPARATION", "DOWNSTREAM_DAG_RUNNING"])
def test_unfinished_and_interrupted_attempts_are_partial_not_live_claims(tmp_path, state):
    run = _owned_run(tmp_path, state=state, overrides={"finished": None, "returncode": None})
    result = build_launcher_report(run, tmp_path / "partial-report")
    assert result["status"] == "PARTIAL"
    assert "no process liveness check" in " ".join(_payload(result)["warnings"])


def test_nonzero_exit_wins_over_inconsistent_success_label(tmp_path):
    run = _owned_run(tmp_path, overrides={"returncode": 9})
    result = build_launcher_report(run, tmp_path / "bad-exit-report")
    assert result["status"] == "FAILED"


@pytest.mark.parametrize("changes", [
    {"returncode": None}, {"dry_run": True}, {"finished": None}, {"input_immutability": "NOT_CHECKED"},
    {"returncode": "0"}, {"returncode": True}, {"dry_run": "false"},
    {"status": "PASS"}, {"status": ["FAILED"]},
])
def test_incoherent_or_unknown_status_never_becomes_success(tmp_path, changes):
    run = _owned_run(tmp_path, overrides=changes)
    output = tmp_path / "invalid-status-report"
    with pytest.raises(InputValidationError):
        build_launcher_report(run, output)
    assert not output.exists()


@pytest.mark.parametrize("changes", [
    {"started": None}, {"started": "yesterday"}, {"started": "2026-10-16T00:00:00"},
    {"finished": "2026-10-15T23:00:00Z"}, {"finished": "2026-10-16T01:00:00"},
])
def test_attempt_must_have_coherent_traceable_times(tmp_path, changes):
    run = _owned_run(tmp_path, overrides=changes)
    with pytest.raises(InputValidationError):
        build_launcher_report(run, tmp_path / "bad-time-report")


@pytest.mark.parametrize("name", ["../attempt_a", "a/b", "a\\b", "", ".", "..", "/tmp/attempt_a"])
def test_attempt_id_cannot_escape_the_attempts_directory(tmp_path, name):
    run = _owned_run(tmp_path)
    with pytest.raises(InputValidationError):
        build_launcher_report(run, tmp_path / "no-traversal", name)


def test_nonexistent_explicit_attempt_fails(tmp_path):
    run = _owned_run(tmp_path)
    with pytest.raises(InputValidationError, match="does not exist"):
        build_launcher_report(run, tmp_path / "no-attempt", "attempt_missing")


def test_missing_manifest_is_not_a_launcher_run(tmp_path):
    run = tmp_path / "unowned"
    run.mkdir()
    with pytest.raises(InputValidationError, match="missing or empty"):
        build_launcher_report(run, tmp_path / "no-unowned-report")


@pytest.mark.parametrize("name", ["effective.yaml", "samples.tsv"])
def test_control_file_changes_are_rejected_without_output(tmp_path, name):
    run = _owned_run(tmp_path)
    with (run / name).open("a", encoding="utf-8") as handle:
        handle.write("\n# changed\n")
    with pytest.raises(InputValidationError, match="LAUNCHER_CONTROL_CHANGED|LAUNCHER_CONFIG_MISMATCH"):
        build_launcher_report(run, tmp_path / "no-changed-report")


@pytest.mark.parametrize("change", ["schema", "kind", "run_dir", "config", "sample_sheet", "control_sha256"])
def test_mismatched_manifest_identity_is_rejected(tmp_path, change):
    run = _owned_run(tmp_path)
    path = run / "launch_manifest.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if change == "control_sha256":
        manifest.pop(change)
    else:
        manifest["identity"][change] = "wrong"
    _write_json(path, manifest)
    with pytest.raises(InputValidationError):
        build_launcher_report(run, tmp_path / "no-identity-report")


def test_report_destination_must_be_fresh_and_outside_original_run(tmp_path):
    run = _owned_run(tmp_path)
    before = _fingerprint(run)
    with pytest.raises(OutputValidationError, match="independent"):
        build_launcher_report(run, run / "reports")
    output = tmp_path / "report"
    build_launcher_report(run, output)
    frozen = _fingerprint(output)
    with pytest.raises(OutputValidationError, match="already exists"):
        build_launcher_report(run, output)
    assert _fingerprint(output) == frozen
    assert _fingerprint(run) == before


def test_dotdot_in_output_or_run_path_is_rejected(tmp_path):
    run = _owned_run(tmp_path)
    with pytest.raises(InputValidationError, match="must not contain"):
        build_launcher_report(run, tmp_path / "child" / ".." / "report")
    with pytest.raises(InputValidationError, match="must not contain"):
        build_launcher_report(run / ".." / run.name, tmp_path / "report")


@pytest.mark.parametrize("target", ["run", "output", "attempt", "status", "manifest", "attempts"])
def test_redirected_paths_are_rejected_without_following_symlinks(tmp_path, monkeypatch, target):
    run = _owned_run(tmp_path)
    output = tmp_path / "report"
    redirect = {
        "run": run, "output": output, "attempt": run / "attempts" / "attempt_a",
        "status": run / "attempts" / "attempt_a" / "status.json",
        "manifest": run / "launch_manifest.json", "attempts": run / "attempts",
    }[target]
    original = Path.is_symlink
    monkeypatch.setattr(Path, "is_symlink", lambda self: self == redirect or original(self))
    with pytest.raises(InputValidationError, match="Redirected"):
        build_launcher_report(run, output)
    assert not output.exists()


def test_per_caller_enablement_is_not_inferred_from_present_artifacts(tmp_path):
    config = {
        "reference": {"build": "GRCh38"}, "small": {"enabled": False},
        "sv": {"enabled": True, "sawfish": {"enabled": True}, "sniffles2": {"enabled": False}, "harmonization": {"enabled": True}},
    }
    run = _owned_run(tmp_path, config=config)
    result = build_launcher_report(run, tmp_path / "sv-report")
    tracks = _tracks(_payload(result))
    assert tracks["small"]["status"] == "DISABLED"
    assert tracks["sv.sawfish"]["status"] == "PARTIAL"
    assert tracks["sv.sniffles2"]["status"] == "DISABLED"
    assert tracks["sv.harmonization"]["status"] == "PARTIAL"
    assert tracks["sv.harmonization"]["phase"] == "9"


def test_invalid_enabled_type_is_not_silently_coerced(tmp_path):
    run = _owned_run(tmp_path, config={"reference": {}, "small": {"enabled": "false"}})
    with pytest.raises(InputValidationError, match="boolean enabled"):
        build_launcher_report(run, tmp_path / "no-coercion")


def test_malformed_json_and_yaml_fail_before_creating_destination(tmp_path):
    run = _owned_run(tmp_path)
    (run / "attempts" / "attempt_a" / "status.json").write_text("{bad", encoding="utf-8")
    with pytest.raises(InputValidationError, match="Invalid launcher metadata"):
        build_launcher_report(run, tmp_path / "bad-json")
    (run / "effective.yaml").write_text("[bad", encoding="utf-8")
    with pytest.raises(InputValidationError, match="Invalid launcher metadata"):
        build_launcher_report(run, tmp_path / "bad-yaml")


def test_failure_report_preserves_machine_error_and_private_paths_only_locally(tmp_path):
    run = _owned_run(tmp_path, state="FAILED", overrides={"returncode": 1, "error_type": "WorkflowError", "error": "tool stopped"})
    result = build_launcher_report(run, tmp_path / "local-failure")
    payload = _payload(result)
    assert payload["provenance"]["selected_attempt_status"]["error_type"] == "WorkflowError"
    assert payload["samples"][0]["input"] == "/data/private/HG002.bam"
    assert result["output_publication"] == "LOCAL_ONLY_NOT_PUBLISHED"
    assert payload["tracks"][0]["qc"]["accuracy_status"] == "NOT_ASSESSED"
