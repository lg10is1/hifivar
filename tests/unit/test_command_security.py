"""Synthetic-only regressions for safe command diagnostic surfaces."""

from __future__ import annotations

import json
import logging
import subprocess
import sys
import traceback
from pathlib import Path

import pytest

from hifivar.command import CommandResult, CommandRunner, _build_environment
from hifivar.exceptions import CommandExecutionError, InputValidationError, ToolNotFoundError


CANARY = "SYNTHETIC_NOT_A_REAL_CREDENTIAL_20261016"


def _python(code: str, *arguments: str) -> list[str]:
    return [sys.executable, "-c", code, *arguments]


@pytest.mark.parametrize(
    "key",
    (
        "TOKEN", "PASSWORD", "SERVICE_API_KEY", "service-api-key",
        "ALPHA_GENOME_API_KEY", "SERVICE_SECRET", "SERVICE_CREDENTIALS",
        "SERVICE_AUTHORIZATION", "SERVICE_PRIVATE_KEY", "SERVICE_PASSWD",
        "AWS_SECRET_ACCESS_KEY", "AWS_ACCESS_KEY_ID",
    ),
)
def test_semantic_environment_values_redact_failure_logs_and_messages(
    key: str, caplog: pytest.LogCaptureFixture,
) -> None:
    command = _python(
        "import os,sys; print(os.environ[sys.argv[1]], file=sys.stderr); sys.exit(8)",
        key,
    )
    with caplog.at_level(logging.DEBUG, logger="hifivar"):
        with pytest.raises(CommandExecutionError) as captured:
            CommandRunner().run(command, env={key: CANARY})

    assert CANARY not in str(captured.value)
    assert CANARY not in repr(captured.value)
    assert CANARY not in caplog.text
    assert "return code 8" in str(captured.value)
    assert "***" in str(captured.value)


def test_inherited_secret_is_preserved_for_child_but_not_diagnostics(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture,
) -> None:
    monkeypatch.setenv("HIFIVAR_TEST_API_KEY", CANARY)
    monkeypatch.setenv("HIFIVAR_NONSENSITIVE_MODE", "reproducible-mode")
    with caplog.at_level(logging.DEBUG, logger="hifivar"):
        result = CommandRunner().run(_python(
            "import os; print(os.environ['HIFIVAR_TEST_API_KEY']); "
            "print(os.environ['HIFIVAR_NONSENSITIVE_MODE'])"
        ))

    assert result.stdout == CANARY + "\nreproducible-mode\n"
    assert CANARY not in caplog.text
    assert CANARY not in repr(result)
    assert CANARY not in str(result)
    safe = result.to_dict()
    assert CANARY not in json.dumps(safe)
    assert safe["stdout"] == "***\nreproducible-mode\n"
    assert safe["returncode"] == 0
    assert safe["executed"] is True
    assert "environment" not in safe


def test_environment_repr_does_not_display_any_values() -> None:
    environment = _build_environment(
        {"HIFIVAR_TEST_API_KEY": CANARY, "HIFIVAR_NONSENSITIVE_MODE": "visible-only-to-child"},
        inherit_env=False,
    )
    assert environment["HIFIVAR_TEST_API_KEY"] == CANARY
    assert environment["HIFIVAR_NONSENSITIVE_MODE"] == "visible-only-to-child"
    assert CANARY not in repr(environment)
    assert "visible-only-to-child" not in str(environment)
    assert "2 variables" in repr(environment)


@pytest.mark.parametrize(
    "argument",
    (
        "--token=" + CANARY,
        "https://user:" + CANARY + "@example.invalid/resource",
        "https://example.invalid/resource?token=" + CANARY,
    ),
)
def test_explicit_redactions_cover_embedded_argv_in_safe_result(
    argument: str, caplog: pytest.LogCaptureFixture,
) -> None:
    with caplog.at_level(logging.DEBUG, logger="hifivar"):
        result = CommandRunner().run(
            _python("import sys; print(sys.argv[1])", argument),
            redact_values=[CANARY],
        )

    assert result.args[-1] == argument
    assert result.stdout == argument + "\n"
    assert CANARY not in repr(result)
    assert CANARY not in json.dumps(result.to_dict())
    assert CANARY not in caplog.text


def test_check_false_failure_result_has_safe_display() -> None:
    result = CommandRunner().run(
        _python("import os,sys; print(os.environ['TEST_TOKEN']); sys.exit(5)"),
        env={"TEST_TOKEN": CANARY}, check=False,
    )
    assert result.stdout == CANARY + "\n"
    assert result.returncode == 5
    assert CANARY not in repr(result)
    assert result.to_dict()["stdout"] == "***\n"


@pytest.mark.parametrize("kind", ("timeout", "missing", "permission"))
def test_sensitive_low_level_causes_do_not_reappear_in_tracebacks(
    kind: str, monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    def fake_run(*args: object, **kwargs: object) -> None:
        if kind == "timeout":
            raise subprocess.TimeoutExpired([sys.executable, CANARY], 1, stderr=CANARY)
        if kind == "missing":
            raise FileNotFoundError("tool vanished " + CANARY)
        raise PermissionError("execution denied " + CANARY)

    monkeypatch.setattr("hifivar.command.subprocess.run", fake_run)
    expected = ToolNotFoundError if kind == "missing" else CommandExecutionError
    with caplog.at_level(logging.DEBUG, logger="hifivar"):
        with pytest.raises(expected) as captured:
            CommandRunner().run(
                _python("pass", CANARY), timeout=1, redact_values=[CANARY],
            )

    formatted = "".join(traceback.format_exception(captured.value))
    assert CANARY not in str(captured.value)
    assert CANARY not in repr(captured.value)
    assert CANARY not in formatted
    assert CANARY not in caplog.text
    assert captured.value.__suppress_context__ is True
    assert "***" in formatted


def test_missing_executable_error_redacts_explicit_value() -> None:
    with pytest.raises(ToolNotFoundError) as captured:
        CommandRunner().run(
            ["missing_executable_" + CANARY], redact_values=[CANARY],
        )
    assert CANARY not in str(captured.value)
    assert "***" in str(captured.value)


def test_invalid_cwd_error_redacts_explicit_value(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError) as captured:
        CommandRunner().run(
            _python("pass"), cwd=tmp_path / CANARY,
            redact_values=[CANARY],
        )
    assert CANARY not in str(captured.value)
    assert "does not exist" in str(captured.value)


def test_output_open_failure_does_not_reveal_sensitive_path(
    tmp_path: Path, caplog: pytest.LogCaptureFixture,
) -> None:
    blocker = tmp_path / CANARY
    blocker.write_text("block", encoding="utf-8")
    with caplog.at_level(logging.DEBUG, logger="hifivar"):
        with pytest.raises(Exception) as captured:
            CommandRunner().run(
                _python("pass"), stdout_path=blocker / "output.log",
                redact_values=[CANARY],
            )
    assert type(captured.value).__name__ == "OutputValidationError"
    assert CANARY not in str(captured.value)
    assert CANARY not in "".join(traceback.format_exception(captured.value))
    assert CANARY not in caplog.text


def test_safe_dry_run_does_not_write_secret_or_output(tmp_path: Path) -> None:
    output = tmp_path / CANARY / "out.log"
    result = CommandRunner().run(
        ["uninstalled-tool", "--token", CANARY],
        env={"TEST_API_KEY": CANARY}, stdout_path=output, dry_run=True,
    )
    assert result.args[-1] == CANARY
    assert result.executed is False
    assert CANARY not in repr(result)
    assert CANARY not in json.dumps(result.to_dict())
    assert not output.parent.exists()


def test_safe_record_is_independent_copy() -> None:
    result = CommandRunner().run(_python("print('safe marker')"))
    first = result.to_dict()
    first["args"] = [CANARY]
    assert CANARY not in repr(result)
    assert CANARY not in json.dumps(result.to_dict())
    assert result.to_dict()["stdout"] == "safe marker\n"


def test_manually_constructed_result_omits_unverified_raw_fields() -> None:
    result = CommandResult(
        ("tool", CANARY), 1, CANARY, CANARY, 0.1, Path(CANARY), True,
        Path(CANARY), Path(CANARY),
    )
    assert result.args[-1] == CANARY
    assert CANARY not in repr(result)
    assert CANARY not in json.dumps(result.to_dict())
    assert result.to_dict()["returncode"] == 1


def test_original_tool_file_bytes_are_explicitly_not_redacted(tmp_path: Path) -> None:
    """Binary/file output compatibility is not confused with diagnostic safety."""
    output = tmp_path / "private-tool-output.log"
    result = CommandRunner().run(
        _python("import os; print(os.environ['TEST_API_KEY'])"),
        env={"TEST_API_KEY": CANARY}, stdout_path=output,
    )
    assert result.stdout is None
    assert output.read_text(encoding="utf-8") == CANARY + "\n"
    assert CANARY not in repr(result)
    assert CANARY not in json.dumps(result.to_dict())
