"""Cross-platform tests for the HiFiVar command execution layer."""

from __future__ import annotations

import logging
import os
import subprocess
import sys
from pathlib import Path

import pytest

from hifivar.command import CommandRunner, format_command
from hifivar.exceptions import (
    CommandExecutionError,
    InputValidationError,
    ToolNotFoundError,
)


PYTHON = Path(sys.executable)
MISSING_EXECUTABLE = "hifivar_tool_that_does_not_exist_123456"


def python_command(code: str, *args: str | Path) -> list[str | Path]:
    """Build a platform-independent tiny external command."""
    return [PYTHON, "-c", code, *args]


def test_successful_command_returns_captured_result() -> None:
    """A successful process should expose output and execution metadata."""
    result = CommandRunner().run(python_command("print('hello')"))

    assert result.executed is True
    assert result.returncode == 0
    assert result.stdout is not None
    assert "hello" in result.stdout
    assert result.stderr == ""


def test_stderr_is_captured() -> None:
    """Small stderr output should be retained in memory by default."""
    result = CommandRunner().run(
        python_command("import sys; print('warning', file=sys.stderr)")
    )

    assert result.stderr is not None
    assert "warning" in result.stderr


def test_nonzero_exit_raises_by_default() -> None:
    """A non-zero return code should become a clear HiFiVar error."""
    with pytest.raises(CommandExecutionError, match="return code 7"):
        CommandRunner().run(python_command("import sys; sys.exit(7)"))


def test_check_false_returns_nonzero_result() -> None:
    """Explicit ``check=False`` should preserve a non-zero result."""
    result = CommandRunner().run(
        python_command("import sys; sys.exit(7)"),
        check=False,
    )

    assert result.executed is True
    assert result.returncode == 7


def test_missing_executable_raises_tool_not_found_error() -> None:
    """Missing tools should not leak a low-level FileNotFoundError."""
    with pytest.raises(ToolNotFoundError, match=MISSING_EXECUTABLE):
        CommandRunner().run([MISSING_EXECUTABLE, "--version"])


def test_executable_race_is_wrapped_as_tool_not_found(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A tool disappearing after availability checks should stay explicit."""
    def fake_run(*args: object, **kwargs: object) -> None:
        raise FileNotFoundError("executable vanished")

    monkeypatch.setattr("hifivar.command.subprocess.run", fake_run)

    with pytest.raises(ToolNotFoundError, match="could not be started"):
        CommandRunner().run(python_command("pass"))


def test_process_creation_os_error_is_wrapped(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Permission and related process-creation errors use HiFiVar errors."""
    def fake_run(*args: object, **kwargs: object) -> None:
        raise PermissionError("execution denied")

    monkeypatch.setattr("hifivar.command.subprocess.run", fake_run)

    with pytest.raises(CommandExecutionError, match="execution denied"):
        CommandRunner().run(python_command("pass"))


def test_find_and_require_executable() -> None:
    """Availability helpers should locate the current Python executable."""
    runner = CommandRunner()

    assert runner.find_executable(PYTHON) is not None
    assert runner.require_executable(PYTHON).is_file()


def test_empty_command_is_rejected() -> None:
    """An empty command cannot be executed or logged meaningfully."""
    with pytest.raises(InputValidationError, match="must not be empty"):
        CommandRunner().run([])


def test_shell_string_command_is_rejected() -> None:
    """A whole shell command string must not be accepted as an arg sequence."""
    with pytest.raises(InputValidationError, match="sequence"):
        CommandRunner().run("python -c pass")  # type: ignore[arg-type]


def test_non_string_command_element_is_rejected() -> None:
    """Unexpected objects should not be silently stringified."""
    with pytest.raises(InputValidationError, match="index 1"):
        CommandRunner().run([PYTHON, 42])  # type: ignore[list-item]


def test_subprocess_receives_args_list_and_shell_false(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Execution must never reinterpret the display as a shell command."""
    observed: dict[str, object] = {}

    def fake_run(args: object, **kwargs: object) -> subprocess.CompletedProcess[str]:
        observed["args"] = args
        observed["shell"] = kwargs["shell"]
        return subprocess.CompletedProcess(args, 0, "", "")

    monkeypatch.setattr("hifivar.command.subprocess.run", fake_run)

    CommandRunner().run([PYTHON, "-c", "pass"])

    assert observed["args"] == (str(PYTHON), "-c", "pass")
    assert observed["shell"] is False


def test_path_command_argument_is_converted_safely(tmp_path: Path) -> None:
    """Path arguments should reach the process as one unquoted argument."""
    path_argument = tmp_path / "sample with spaces.bam"
    result = CommandRunner().run(
        python_command("import sys; print(sys.argv[1])", path_argument)
    )

    assert result.stdout is not None
    assert result.stdout.strip() == str(path_argument)
    assert result.args[-1] == str(path_argument)


def test_command_runs_in_requested_cwd(tmp_path: Path) -> None:
    """The child process should observe the requested working directory."""
    result = CommandRunner().run(
        python_command("from pathlib import Path; print(Path.cwd())"),
        cwd=tmp_path,
    )

    assert result.stdout is not None
    assert Path(result.stdout.strip()).resolve() == tmp_path.resolve()
    assert result.cwd == tmp_path


def test_missing_cwd_is_rejected(tmp_path: Path) -> None:
    """CommandRunner must not create a requested working directory."""
    missing_cwd = tmp_path / "missing"

    with pytest.raises(InputValidationError, match="does not exist"):
        CommandRunner().run(python_command("pass"), cwd=missing_cwd)

    assert not missing_cwd.exists()


def test_file_cannot_be_used_as_cwd(tmp_path: Path) -> None:
    """An existing non-directory cwd should fail before process creation."""
    file_path = tmp_path / "not-a-directory"
    file_path.write_text("data", encoding="utf-8")

    with pytest.raises(InputValidationError, match="not a directory"):
        CommandRunner().run(python_command("pass"), cwd=file_path)


def test_environment_override_preserves_existing_environment() -> None:
    """Overrides should be layered over, rather than replace, os.environ."""
    code = (
        "import os; "
        "print(os.environ['HIFIVAR_TEST_VALUE']); "
        "print(bool(os.environ.get('PATH')))"
    )
    result = CommandRunner().run(
        python_command(code),
        env={"HIFIVAR_TEST_VALUE": "abc123"},
    )

    assert result.stdout is not None
    assert result.stdout.splitlines() == ["abc123", "True"]


def test_invalid_environment_override_is_rejected() -> None:
    """Environment keys and values must remain strings for subprocess."""
    with pytest.raises(InputValidationError, match="strings to strings"):
        CommandRunner().run(
            python_command("pass"),
            env={"VALUE": 8},  # type: ignore[dict-item]
        )


def test_timeout_is_wrapped_as_command_execution_error() -> None:
    """TimeoutExpired should not leak through the public API."""
    with pytest.raises(CommandExecutionError, match=r"timed out after 0\.05"):
        CommandRunner().run(
            python_command("import time; time.sleep(2)"),
            timeout=0.05,
        )


@pytest.mark.parametrize("timeout", (0, -1, True, "10"))
def test_invalid_timeout_is_rejected(timeout: object) -> None:
    """Timeouts must be positive numeric seconds."""
    with pytest.raises(InputValidationError, match="timeout"):
        CommandRunner().run(
            python_command("pass"),
            timeout=timeout,  # type: ignore[arg-type]
        )


def test_dry_run_does_not_execute_or_create_side_effect(tmp_path: Path) -> None:
    """Dry-run metadata must not imply that outputs were produced."""
    side_effect = tmp_path / "created.txt"
    result = CommandRunner().run(
        python_command(
            "from pathlib import Path; Path(__import__('sys').argv[1]).touch()",
            side_effect,
        ),
        dry_run=True,
    )

    assert result.executed is False
    assert result.returncode is None
    assert result.duration_seconds == 0.0
    assert not side_effect.exists()


def test_dry_run_can_preview_an_unavailable_tool() -> None:
    """A preview should not require tools installed only on the target HPC."""
    result = CommandRunner().run(
        [MISSING_EXECUTABLE, "--version"],
        dry_run=True,
    )

    assert result.executed is False
    assert result.args[0] == MISSING_EXECUTABLE


def test_duration_is_nonnegative() -> None:
    """Execution duration should use a monotonic high-resolution clock."""
    result = CommandRunner().run(python_command("pass"))

    assert result.duration_seconds >= 0


def test_utf8_stdout_is_captured() -> None:
    """UTF-8 sample text should decode consistently across platforms."""
    code = "import sys; sys.stdout.buffer.write('样本 HG002'.encode('utf-8'))"
    result = CommandRunner().run(python_command(code))

    assert result.stdout == "样本 HG002"


def test_invalid_utf8_output_uses_replacement_character() -> None:
    """Malformed tool bytes should not crash text capture."""
    result = CommandRunner().run(
        python_command("import sys; sys.stdout.buffer.write(b'bad:\\xff')")
    )

    assert result.stdout == "bad:�"


def test_redaction_hides_secret_but_preserves_execution_args(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Secrets should reach the child while remaining absent from logs."""
    secret = "ABC_SECRET_123"
    command = python_command("import sys; print(sys.argv[1])", secret)

    with caplog.at_level(logging.DEBUG, logger="hifivar"):
        result = CommandRunner().run(command, redact_values={secret})

    assert result.stdout is not None
    assert result.stdout.strip() == secret
    assert result.args[-1] == secret
    assert secret not in caplog.text
    assert "***" in caplog.text


def test_format_command_redacts_embedded_secret() -> None:
    """Redaction should also cover values embedded in one display argument."""
    secret = "TOKEN_VALUE"

    display = format_command(
        ["tool", f"--api-key={secret}"],
        redact_values={secret},
    )

    assert secret not in display
    assert display == "tool --api-key=***"


def test_same_runner_supports_repeated_runs() -> None:
    """Runner instances should not retain per-command subprocess state."""
    runner = CommandRunner()

    first = runner.run(python_command("print('first')"))
    second = runner.run(python_command("print('second')"))

    assert first.stdout is not None and first.stdout.strip() == "first"
    assert second.stdout is not None and second.stdout.strip() == "second"


def test_stdout_path_streams_output_to_file(tmp_path: Path) -> None:
    """Explicit stdout files should avoid retaining output in memory."""
    stdout_path = tmp_path / "stdout.log"
    result = CommandRunner().run(
        python_command("print('file stdout')"),
        stdout_path=stdout_path,
    )

    assert result.stdout is None
    assert result.stdout_path == stdout_path
    assert stdout_path.read_text(encoding="utf-8").strip() == "file stdout"


def test_stderr_path_also_persists_stdout_to_sibling_log(tmp_path: Path) -> None:
    """Tool logs retain both streams when callers provide only stderr_path."""
    stderr_path = tmp_path / "stderr.log"
    result = CommandRunner().run(
        python_command(
            "import sys; print('file stdout'); "
            "print('file stderr', file=sys.stderr)"
        ),
        stderr_path=stderr_path,
    )

    stdout_path = tmp_path / "stderr.stdout.log"
    assert result.stdout is None
    assert result.stderr is None
    assert result.stdout_path == stdout_path
    assert result.stderr_path == stderr_path
    assert stdout_path.read_text(encoding="utf-8").strip() == "file stdout"
    assert stderr_path.read_text(encoding="utf-8").strip() == "file stderr"


def test_output_parent_directories_are_created(tmp_path: Path) -> None:
    """Only explicitly requested output-file parents should be created."""
    stdout_path = tmp_path / "nested" / "logs" / "stdout.log"

    CommandRunner().run(
        python_command("print('nested')"),
        stdout_path=stdout_path,
    )

    assert stdout_path.read_text(encoding="utf-8").strip() == "nested"


def test_same_stdout_and_stderr_path_is_rejected(tmp_path: Path) -> None:
    """Opening the same output twice is intentionally unsupported."""
    output_path = tmp_path / "combined.log"

    with pytest.raises(InputValidationError, match="different files"):
        CommandRunner().run(
            python_command("pass"),
            stdout_path=output_path,
            stderr_path=output_path,
        )


def test_original_command_list_is_not_mutated(tmp_path: Path) -> None:
    """Normalization must copy rather than rewrite caller-owned arguments."""
    path_argument = tmp_path / "sample.bam"
    command: list[str | Path] = python_command(
        "import sys; print(sys.argv[1])",
        path_argument,
    )
    original = list(command)

    CommandRunner().run(command)

    assert command == original
    assert isinstance(command[-1], Path)


def test_success_log_contains_executing_command(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Every real execution should have a reproducibility log entry."""
    with caplog.at_level(logging.INFO, logger="hifivar"):
        CommandRunner().run(python_command("pass"))

    assert "Executing command" in caplog.text


def test_failure_log_contains_return_code(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Failure logs should expose the return code without brittle formatting."""
    with caplog.at_level(logging.ERROR, logger="hifivar"):
        with pytest.raises(CommandExecutionError):
            CommandRunner().run(python_command("import sys; sys.exit(9)"))

    assert "return code 9" in caplog.text


def test_failure_exception_contains_stderr_summary() -> None:
    """Captured stderr should aid diagnosis without exposing raw exceptions."""
    code = (
        "import sys; "
        "print('diagnostic failure', file=sys.stderr); "
        "sys.exit(4)"
    )

    with pytest.raises(CommandExecutionError, match="diagnostic failure"):
        CommandRunner().run(python_command(code))


def test_output_paths_are_not_created_during_dry_run(tmp_path: Path) -> None:
    """Dry-run should not create even explicitly named log outputs."""
    stdout_path = tmp_path / "nested" / "stdout.log"

    result = CommandRunner().run(
        python_command("print('unused')"),
        stdout_path=stdout_path,
        dry_run=True,
    )

    assert result.stdout_path == stdout_path
    assert not stdout_path.exists()


def test_environment_override_does_not_mutate_os_environ() -> None:
    """Per-command environment values must not leak into the parent process."""
    key = "HIFIVAR_COMMAND_PARENT_ENV_TEST"
    original = os.environ.get(key)

    CommandRunner().run(
        python_command("pass"),
        env={key: "child-only"},
    )

    assert os.environ.get(key) == original
