"""Tests for the Phase 0.7 HiFiVar command-line foundation."""

from __future__ import annotations

import logging
import subprocess
import sys
from collections.abc import Iterator
from importlib.metadata import version
from pathlib import Path

import pytest
import yaml

from hifivar import __version__
from hifivar.cli import build_parser, main
from hifivar.config import load_yaml
from hifivar.logging_utils import HIFIVAR_LOGGER_NAME


@pytest.fixture(autouse=True)
def isolate_cli_logging() -> Iterator[None]:
    """Restore the HiFiVar logger after each CLI invocation test."""
    namespace_logger = logging.getLogger(HIFIVAR_LOGGER_NAME)
    original_handlers = list(namespace_logger.handlers)
    original_level = namespace_logger.level
    original_propagate = namespace_logger.propagate
    original_disabled = namespace_logger.disabled

    for handler in original_handlers:
        namespace_logger.removeHandler(handler)
    namespace_logger.setLevel(logging.NOTSET)
    namespace_logger.propagate = True
    namespace_logger.disabled = False

    yield

    for handler in list(namespace_logger.handlers):
        namespace_logger.removeHandler(handler)
        handler.close()
    for handler in original_handlers:
        namespace_logger.addHandler(handler)
    namespace_logger.setLevel(original_level)
    namespace_logger.propagate = original_propagate
    namespace_logger.disabled = original_disabled


def parse_stdout_yaml(output: str) -> dict[str, object]:
    """Parse CLI YAML output and assert its root shape."""
    parsed = yaml.safe_load(output)
    assert isinstance(parsed, dict)
    return parsed


def test_build_parser_returns_argparse_parser() -> None:
    """Parser construction should remain available to tests and extensions."""
    assert build_parser().prog == "hifivar"


def test_help_succeeds(capsys: pytest.CaptureFixture[str]) -> None:
    """Top-level help should use argparse's successful exit behavior."""
    with pytest.raises(SystemExit) as exit_info:
        main(["--help"])

    assert exit_info.value.code == 0
    assert "usage: hifivar" in capsys.readouterr().out


def test_version_contains_package_version(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """CLI version output should reuse the package version source."""
    with pytest.raises(SystemExit) as exit_info:
        main(["--version"])

    assert exit_info.value.code == 0
    assert capsys.readouterr().out.strip() == f"hifivar {__version__}"
    assert version("hifivar") == __version__


def test_python_module_version_matches_cli(tmp_path: Path) -> None:
    """``python -m hifivar`` should route through the same CLI main function."""
    completed = subprocess.run(
        [sys.executable, "-m", "hifivar", "--version"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        encoding="utf-8",
        check=False,
    )

    assert completed.returncode == 0
    assert completed.stdout.strip() == f"hifivar {__version__}"


def test_no_arguments_prints_help_and_returns_zero(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """A bare command should offer guidance rather than fail mysteriously."""
    assert main([]) == 0
    assert "usage: hifivar" in capsys.readouterr().out


def test_unknown_argument_uses_argparse_exit_code() -> None:
    """Parser errors should retain argparse's conventional exit code 2."""
    with pytest.raises(SystemExit) as exit_info:
        main(["--this-does-not-exist"])

    assert exit_info.value.code == 2


def test_config_show_outputs_valid_yaml(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Default effective config should be emitted as YAML, not a repr."""
    assert main(["config", "show"]) == 0

    shown = parse_stdout_yaml(capsys.readouterr().out)
    assert shown["project"]["name"] == "hifivar"  # type: ignore[index]
    assert shown["workflow"]["preset"] == "standard"  # type: ignore[index]


@pytest.mark.parametrize("preset", ("fast", "comprehensive"))
def test_preset_changes_effective_workflow(
    preset: str,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Packaged presets should override the default workflow preset."""
    assert main(["--preset", preset, "config", "show"]) == 0

    shown = parse_stdout_yaml(capsys.readouterr().out)
    assert shown["workflow"]["preset"] == preset  # type: ignore[index]


def test_unknown_preset_is_parser_error() -> None:
    """Preset names remain a fixed Phase 0.7 argparse choice."""
    with pytest.raises(SystemExit) as exit_info:
        main(["--preset", "unknown", "config", "show"])

    assert exit_info.value.code == 2


def test_user_config_overrides_packaged_config(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """User YAML should retain its higher precedence over a preset."""
    user_config = tmp_path / "user.yaml"
    user_config.write_text("runtime:\n  threads: 8\n", encoding="utf-8")

    assert main(["--config", str(user_config), "config", "show"]) == 0

    shown = parse_stdout_yaml(capsys.readouterr().out)
    assert shown["runtime"]["threads"] == 8  # type: ignore[index]


def test_cli_log_level_overrides_user_yaml(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Explicit CLI logging options should have the highest precedence."""
    user_config = tmp_path / "user.yaml"
    user_config.write_text("logging:\n  level: INFO\n", encoding="utf-8")

    assert (
        main(
            [
                "--config",
                str(user_config),
                "--log-level",
                "debug",
                "config",
                "show",
            ]
        )
        == 0
    )

    shown = parse_stdout_yaml(capsys.readouterr().out)
    assert shown["logging"]["level"] == "DEBUG"  # type: ignore[index]


def test_invalid_log_level_is_parser_error() -> None:
    """CLI logging validation should reuse the shared level parser."""
    with pytest.raises(SystemExit) as exit_info:
        main(["--log-level", "INVALID_LEVEL", "config", "show"])

    assert exit_info.value.code == 2


def test_missing_user_config_returns_runtime_error(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Missing YAML should produce a concise error and exit code 1."""
    missing = tmp_path / "missing.yaml"

    assert main(["--config", str(missing), "config", "validate"]) == 1

    error_output = capsys.readouterr().err
    assert "ERROR:" in error_output
    assert "missing.yaml" in error_output


def test_invalid_yaml_returns_runtime_error(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Parser details should be wrapped without a normal-user traceback."""
    invalid = tmp_path / "invalid.yaml"
    invalid.write_text("runtime:\n  threads: [\n", encoding="utf-8")

    assert main(["--config", str(invalid), "config", "validate"]) == 1

    error_output = capsys.readouterr().err
    assert "ERROR:" in error_output
    assert "Invalid YAML" in error_output
    assert "Traceback" not in error_output


def test_config_validate_succeeds(capsys: pytest.CaptureFixture[str]) -> None:
    """The validate action should confirm the already merged effective config."""
    assert main(["config", "validate"]) == 0
    assert capsys.readouterr().out.strip() == "Configuration valid."


def test_config_validate_rejects_invalid_value(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Config validation errors should cross the CLI boundary as exit 1."""
    invalid = tmp_path / "invalid.yaml"
    invalid.write_text("runtime:\n  threads: 0\n", encoding="utf-8")

    assert main(["--config", str(invalid), "config", "validate"]) == 1

    assert "positive integer" in capsys.readouterr().err


def test_dump_effective_creates_reloadable_yaml(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The CLI should reuse the Phase 0.5 effective-config writer."""
    output = tmp_path / "nested" / "effective.yaml"

    assert (
        main(
            [
                "--preset",
                "trio",
                "config",
                "dump-effective",
                "--output",
                str(output),
            ]
        )
        == 0
    )

    effective = load_yaml(output)
    assert effective["workflow"]["preset"] == "trio"  # type: ignore[index]
    assert "written" in capsys.readouterr().out


def test_utf8_user_config_is_shown(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Unicode config values should survive CLI display."""
    user_config = tmp_path / "unicode.yaml"
    user_config.write_text(
        'project:\n  name: "HiFiVar 测试"\n',
        encoding="utf-8",
    )

    assert main(["--config", str(user_config), "config", "show"]) == 0

    assert "HiFiVar 测试" in capsys.readouterr().out


def test_log_file_override_configures_file_logging(
    tmp_path: Path,
) -> None:
    """CLI should apply its log-file override before dispatch."""
    log_file = tmp_path / "nested" / "hifivar.log"

    assert main(["--log-file", str(log_file), "doctor"]) == 0

    log_text = log_file.read_text(encoding="utf-8")
    assert "HiFiVar CLI initialized" in log_text


def test_log_file_os_error_returns_runtime_error(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Expected logging filesystem failures should not show a traceback."""
    def fail_logging(*args: object, **kwargs: object) -> None:
        raise OSError("log destination denied")

    monkeypatch.setattr("hifivar.cli.configure_logging", fail_logging)

    assert main(["doctor"]) == 1

    error_output = capsys.readouterr().err
    assert "log destination denied" in error_output
    assert "Traceback" not in error_output


def test_doctor_reports_foundation_environment(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Doctor should inspect HiFiVar itself, not future external tools."""
    assert main(["doctor"]) == 0

    output = capsys.readouterr().out
    assert "HiFiVar doctor" in output
    assert f"Version: {__version__}" in output
    assert "Python:" in output
    assert "Platform:" in output
    assert "Configuration: OK" in output


def test_packaged_configs_work_outside_repository_cwd(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Default and preset lookup must not depend on repository cwd."""
    monkeypatch.chdir(tmp_path)

    assert main(["--preset", "cohort", "config", "show"]) == 0

    shown = parse_stdout_yaml(capsys.readouterr().out)
    assert shown["workflow"]["preset"] == "cohort"  # type: ignore[index]


def test_global_options_after_subcommand_are_not_supported() -> None:
    """Phase 0.7 intentionally requires global options before subcommands."""
    with pytest.raises(SystemExit) as exit_info:
        main(["config", "show", "--preset", "fast"])

    assert exit_info.value.code == 2
