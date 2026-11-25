"""CLI dispatch and zero-write plans; fixtures do not establish native TRGT success."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from hifivar.cli import build_parser, main
from hifivar.exceptions import InputValidationError


def test_review_parser_has_three_explicit_actions():
    args = build_parser().parse_args(["review", "tr-plot", "--plan", "selected.yaml"])
    assert args.review_command == "tr-plot" and args.plan == Path("selected.yaml")
    with pytest.raises(SystemExit):
        build_parser().parse_args(["review", "tr-plot"])


@pytest.mark.parametrize("status,exit_code", [("COMPLETE_TECHNICAL_SVG", 0), ("PLANNED", 0), ("FAILED", 1), ("COMPLETE", 1)])
def test_tr_plot_dispatch_is_explicit_and_returns_technical_status(monkeypatch, capsys, tmp_path, status, exit_code):
    observed = []
    def run(path, *, dry_run):
        observed.append((path, dry_run))
        return {"status": status, "manual_review": "NOT_REVIEWED", "accuracy": "NOT_ESTABLISHED"}
    monkeypatch.setattr("hifivar.trgt_plot_plan.run_trgt_plot_plan", run)
    monkeypatch.setattr("hifivar.review_cli.import_review_decisions", lambda *a, **kw: pytest.fail("wrong action"))
    plan = tmp_path / "plot.yaml"
    assert main(["--dry-run", "review", "tr-plot", "--plan", str(plan)]) == exit_code
    assert observed == [(plan, True)]
    assert json.loads(capsys.readouterr().out)["manual_review"] == "NOT_REVIEWED"


def test_invalid_plan_uses_existing_cli_error_contract(monkeypatch, capsys, tmp_path):
    def fail(*args, **kwargs):
        raise InputValidationError("TR plot fixture plan is invalid")
    monkeypatch.setattr("hifivar.trgt_plot_plan.run_trgt_plot_plan", fail)
    assert main(["review", "tr-plot", "--plan", str(tmp_path / "missing.yaml")]) == 1
    captured = capsys.readouterr()
    assert not captured.out and "ERROR:" in captured.err and "Traceback" not in captured.err


def test_tr_plot_dry_run_does_not_open_global_log(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr("hifivar.trgt_plot_plan.run_trgt_plot_plan", lambda *a, **kw: {"status": "PLANNED"})
    log = tmp_path / "must_not_be_created.log"
    assert main(["--log-file", str(log), "--dry-run", "review", "tr-plot", "--plan", str(tmp_path / "plot.yaml")]) == 0
    assert not log.exists() and json.loads(capsys.readouterr().out)["status"] == "PLANNED"


@pytest.mark.parametrize("action", ["prepare", "import-decisions"])
def test_prior_review_actions_do_not_enter_tr_plot(monkeypatch, capsys, tmp_path, action):
    monkeypatch.setattr("hifivar.trgt_plot_plan.run_trgt_plot_plan", lambda *a, **kw: pytest.fail("wrong TR plot action"))
    observed = []
    if action == "prepare":
        def prepare(config, output, *, dry_run):
            observed.append((output, dry_run))
            return {"status": "REVIEW_PREPARED"}
        monkeypatch.setattr("hifivar.review_cli.prepare_review", prepare)
        args = ["review", "prepare", "--output-dir", str(tmp_path / "old_prepare")]
    else:
        def decisions(bundle, table, output, *, dry_run):
            observed.append((output, dry_run))
            return {"status": "DECISIONS_IMPORTED"}
        monkeypatch.setattr("hifivar.review_cli.import_review_decisions", decisions)
        args = ["review", "import-decisions", "--bundle-dir", str(tmp_path / "bundle"),
                "--decisions", str(tmp_path / "review.tsv"), "--output-dir", str(tmp_path / "old_decisions")]
    assert main(args) == 0 and observed[0][1] is False
    assert json.loads(capsys.readouterr().out)["status"] in {"REVIEW_PREPARED", "DECISIONS_IMPORTED"}
