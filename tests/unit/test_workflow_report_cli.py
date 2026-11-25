"""CLI adaptation preserves selected attempt state without rerunning analyses."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest
import yaml

from hifivar.cli import build_parser, main


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


def _status(state: str) -> dict:
    return {"started": "2026-10-16T00:00:00Z", "finished": "2026-10-16T01:00:00Z",
            "status": state, "dry_run": False, "returncode": 1 if state == "FAILED" else 0,
            "input_immutability": "PASS", "snakemake_version": "9.26.1"}


@pytest.fixture
def existing_run(tmp_path: Path) -> Path:
    root = tmp_path / "existing-run"
    root.mkdir()
    config = {"reference": {"build": "GRCh38", "fasta": "ref.fa"}, "small": {"enabled": True}}
    sheet = "sample_id\tinput\tinput_type\nHG002\tHG002.bam\tbam\n"
    (root / "effective.yaml").write_text(yaml.safe_dump(config), encoding="utf-8")
    (root / "samples.tsv").write_text(sheet, encoding="utf-8")
    _write_json(root / "launch_manifest.json", {
        "identity": {"schema": 2, "kind": "hifivar_two_stage_launcher", "run_dir": str(root.absolute()),
                     "config": config, "sample_sheet": sheet, "runtime": {"hifivar_version": "0.1.0rc5"}, "inputs": []},
        "control_sha256": {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                           for name in ("effective.yaml", "samples.tsv")},
    })
    _write_json(root / "attempts" / "attempt_a" / "status.json", _status("WORKFLOW_SUCCEEDED"))
    return root


def _fingerprint(root: Path) -> dict[str, str]:
    return {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in root.rglob("*") if path.is_file()}


def _args(run: Path, output: Path, *extra: str) -> list[str]:
    return ["report", "--run-dir", str(run), "--output-dir", str(output), *extra]


def test_report_parser_registers_explicit_attempt(existing_run: Path, tmp_path: Path) -> None:
    args = build_parser().parse_args(_args(existing_run, tmp_path / "report", "--attempt-id", "attempt_a"))
    assert args.command == "report" and args.attempt_id == "attempt_a"
    assert args.run_dir == existing_run


def test_partial_report_generation_succeeds_but_does_not_claim_calling_failure(
    existing_run: Path, tmp_path: Path, capsys,
) -> None:
    before = _fingerprint(existing_run)
    output = tmp_path / "report"
    assert main(_args(existing_run, output)) == 0
    captured = capsys.readouterr()
    result = json.loads(captured.out)
    assert result["status"] == "PARTIAL" and result["execution_status"] == "WORKFLOW_SUCCEEDED"
    assert result["scientific_accuracy"] == "NOT_ASSESSED"
    assert "does not by itself mean calling failed" in captured.err
    assert set(result["output_paths"]) == {"json", "yaml", "md", "html"}
    assert all(Path(path).is_file() for path in result["output_paths"].values())
    payload = json.loads((output / "final_report.json").read_text(encoding="utf-8"))
    assert payload["status"] == "PARTIAL"
    assert payload["provenance"]["attempt_id"] == "attempt_a"
    assert _fingerprint(existing_run) == before


def test_failed_attempt_remains_failed_in_all_report_formats_and_cli(
    existing_run: Path, tmp_path: Path, capsys,
) -> None:
    _write_json(existing_run / "attempts" / "attempt_b" / "status.json", _status("FAILED"))
    before = _fingerprint(existing_run)
    output = tmp_path / "failed-report"
    assert main(_args(existing_run, output, "--attempt-id", "attempt_b")) == 0
    captured = capsys.readouterr()
    result = json.loads(captured.out)
    assert result["status"] == result["execution_status"] == "FAILED"
    assert "selected attempt remains FAILED" in captured.err
    assert result["attempt_id"] == "attempt_b"
    for path in result["output_paths"].values():
        assert "FAILED" in Path(path).read_text(encoding="utf-8")
    assert _fingerprint(existing_run) == before


def test_global_dry_run_rejected_without_backend_call_or_output(
    existing_run: Path, tmp_path: Path, monkeypatch, capsys,
) -> None:
    def forbidden(*args, **kwargs):
        pytest.fail("report dry-run must not call the writing adapter")
    monkeypatch.setattr("hifivar.workflow_report.build_launcher_report", forbidden)
    output = tmp_path / "not-created"
    assert main(["--dry-run", *_args(existing_run, output)]) == 1
    assert "not supported" in capsys.readouterr().err
    assert not output.exists()


def test_existing_output_is_preserved_and_cli_returns_one(existing_run: Path, tmp_path: Path) -> None:
    output = tmp_path / "owned-output"
    output.mkdir()
    sentinel = output / "owner.txt"
    sentinel.write_text("preserve", encoding="utf-8")
    assert main(_args(existing_run, output)) == 1
    assert sentinel.read_text() == "preserve"
    assert list(output.iterdir()) == [sentinel]


@pytest.mark.parametrize("attempt", [None, "../attempt_a", "missing"])
def test_multiple_or_invalid_attempts_fail_without_output(
    existing_run: Path, tmp_path: Path, attempt: str | None,
) -> None:
    _write_json(existing_run / "attempts" / "attempt_b" / "status.json", _status("FAILED"))
    output = tmp_path / "invalid-report"
    extra = ("--attempt-id", attempt) if attempt is not None else ()
    assert main(_args(existing_run, output, *extra)) == 1
    assert not output.exists()


def test_mutated_control_is_error_not_partial_success(existing_run: Path, tmp_path: Path) -> None:
    with (existing_run / "samples.tsv").open("a", encoding="utf-8") as stream:
        stream.write("EXTRA\tmore.bam\tbam\n")
    output = tmp_path / "invalid-report"
    assert main(_args(existing_run, output)) == 1
    assert not output.exists()


def test_report_cannot_write_inside_original_run(existing_run: Path) -> None:
    before = _fingerprint(existing_run)
    output = existing_run / "report"
    assert main(_args(existing_run, output)) == 1
    assert not output.exists() and _fingerprint(existing_run) == before


def test_logging_override_does_not_create_unowned_file(existing_run: Path, tmp_path: Path) -> None:
    log = tmp_path / "outside" / "cli.log"
    assert main(["--log-file", str(log), *_args(existing_run, tmp_path / "report")]) == 0
    assert not log.parent.exists()


def test_existing_cli_subcommands_remain_registered() -> None:
    parser = build_parser()
    assert parser.parse_args(["benchmark", "sv", "--manifest", "job.yaml", "--output-dir", "new"]).command == "benchmark"
    assert parser.parse_args(["review", "prepare", "--output-dir", "new"]).command == "review"
    assert parser.parse_args(["run", "--run-dir", "new"]).command == "run"
