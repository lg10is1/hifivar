"""CLI dispatch only for support-vs-depth analysis; no biological execution."""
from __future__ import annotations

import importlib
import json
import logging
from pathlib import Path
import subprocess

import pytest

from hifivar.cli import build_parser, main
from hifivar.command import CommandRunner
from hifivar.logging_utils import HIFIVAR_LOGGER_NAME


@pytest.fixture(autouse=True)
def isolate_logging_and_forbid_external_execution(monkeypatch):
    logger = logging.getLogger(HIFIVAR_LOGGER_NAME)
    handlers = list(logger.handlers)
    original = (logger.level, logger.propagate, logger.disabled)
    for handler in handlers:
        logger.removeHandler(handler)
    logger.setLevel(logging.NOTSET);logger.propagate = True;logger.disabled = False
    def forbidden(*args, **kwargs):
        raise AssertionError("CLI dispatch tests must not run callers, Jasmine or benchmarks")
    monkeypatch.setattr(subprocess, "Popen", forbidden)
    monkeypatch.setattr(subprocess, "run", forbidden)
    monkeypatch.setattr(CommandRunner, "run", forbidden)
    monkeypatch.setattr(CommandRunner, "require_executable", forbidden)
    yield
    for handler in list(logger.handlers):
        logger.removeHandler(handler);handler.close()
    for handler in handlers:
        logger.addHandler(handler)
    logger.setLevel(original[0]);logger.propagate = original[1];logger.disabled = original[2]


@pytest.fixture
def support_api():
    # Import the actual module after the no-execution guard is active. Dispatch
    # is replaced below; this is not an invented sys.modules substitute.
    return importlib.import_module("hifivar.sv_depth_support")


def args_for(manifest, output):
    return ["benchmark", "sv-depth-support", "--manifest", str(manifest),
            "--output-dir", str(output)]


def test_support_cli_help_lists_only_supported_analysis_options(capsys):
    with pytest.raises(SystemExit) as caught:
        main(["benchmark", "sv-depth-support", "--help"])
    assert caught.value.code == 0
    help_text = capsys.readouterr().out
    assert "--manifest" in help_text and "--output-dir" in help_text and "--unit-id" in help_text
    assert "--caller" not in help_text and "--resume" not in help_text


def test_support_parser_keeps_path_objects_and_no_caller_resume_namespace(tmp_path):
    manifest = tmp_path / "plan with spaces.yaml";output = tmp_path / "output with spaces"
    parsed = build_parser().parse_args(args_for(manifest, output))
    assert parsed.command == "benchmark" and parsed.benchmark_command == "sv-depth-support"
    assert parsed.manifest == manifest and isinstance(parsed.manifest, Path)
    assert parsed.output_dir == output and isinstance(parsed.output_dir, Path)
    assert parsed.unit_ids is None and parsed.dry_run is False
    assert not hasattr(parsed, "callers") and not hasattr(parsed, "resume")


@pytest.mark.parametrize("dry_run", [False, True])
def test_support_dispatch_passes_exact_paths_and_global_dry_run(tmp_path, monkeypatch, capsys, support_api, dry_run):
    manifest = tmp_path / "manifest.yaml";output = tmp_path / "not_created"
    seen = []
    payload = {"status": "PLANNED" if dry_run else "COMPLETE", "SYNTHETIC_DISPATCH_ONLY": True}
    def fake(plan, destination, **kwargs):
        seen.append((plan, destination, kwargs))
        return payload
    monkeypatch.setattr(support_api, "run_sv_depth_support", fake)
    argv = (["--dry-run"] if dry_run else []) + args_for(manifest, output)
    assert main(argv) == 0
    assert seen == [(manifest, output, {"dry_run": dry_run, "unit_ids": None})]
    assert json.loads(capsys.readouterr().out) == payload
    assert not output.exists()


def test_support_dispatch_preserves_repeated_unit_ids_without_new_options(tmp_path, monkeypatch, support_api):
    manifest = tmp_path / "manifest.yaml";output = tmp_path / "not_created"
    seen = []
    def fake(plan, destination, **kwargs):
        seen.append((plan, destination, kwargs))
        return {"status": "PLANNED", "SYNTHETIC_DISPATCH_ONLY": True}
    monkeypatch.setattr(support_api, "run_sv_depth_support", fake)
    argv = ["--dry-run", *args_for(manifest, output),
            "--unit-id", "depth10_r2", "--unit-id", "depth5_r1"]
    assert main(argv) == 0
    assert seen == [(manifest, output, {"dry_run": True, "unit_ids": ["depth10_r2", "depth5_r1"]})]
    assert not output.exists()


@pytest.mark.parametrize("payload,expected", [
    ({"status": "PARTIAL", "selection_status": "COMPLETE"}, 0),
    ({"status": "COMPLETE", "selection_status": "FAILED"}, 1),
    ({"status": "COMPLETE", "selection_status": "PARTIAL"}, 1),
    ({"status": "FAILED", "selection_status": "PLANNED"}, 0),
    ({"status": "COMPLETE"}, 0),
    ({"status": "PLANNED"}, 0),
    ({"status": "FAILED"}, 1),
    ({"status": "RUNNING"}, 1),
    ({}, 1),
])
def test_support_cli_uses_selection_status_then_global_status(tmp_path, monkeypatch, capsys, support_api, payload, expected):
    monkeypatch.setattr(support_api, "run_sv_depth_support", lambda *a, **k: dict(payload))
    output = tmp_path / "not_created"
    assert main(args_for(tmp_path / "manifest.yaml", output)) == expected
    assert json.loads(capsys.readouterr().out) == payload
    assert not output.exists()


@pytest.mark.parametrize("unsupported", [["--caller", "sniffles2"], ["--resume"]])
def test_support_cli_rejects_calling_and_resume_options_before_dispatch(tmp_path, monkeypatch, support_api, unsupported):
    calls = []
    def forbidden(*args, **kwargs):
        calls.append(True)
        raise AssertionError("Unsupported parser options must not reach analysis")
    monkeypatch.setattr(support_api, "run_sv_depth_support", forbidden)
    with pytest.raises(SystemExit) as caught:
        main([*args_for(tmp_path / "plan.yaml", tmp_path / "out"), *unsupported])
    assert caught.value.code == 2 and calls == []


@pytest.mark.parametrize("missing", ["manifest", "output"])
def test_support_cli_requires_manifest_and_output_before_dispatch(tmp_path, monkeypatch, support_api, missing):
    calls = []
    monkeypatch.setattr(support_api, "run_sv_depth_support", lambda *a, **k: calls.append(True))
    argv = ["benchmark", "sv-depth-support"]
    argv += ["--output-dir", str(tmp_path / "out")] if missing == "manifest" else ["--manifest", str(tmp_path / "plan.yaml")]
    with pytest.raises(SystemExit) as caught:
        main(argv)
    assert caught.value.code == 2 and calls == []
