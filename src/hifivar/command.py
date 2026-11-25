"""Safe, reproducible external command execution for HiFiVar."""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from collections.abc import Collection, Mapping, Sequence
from contextlib import ExitStack
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import BinaryIO

from hifivar.exceptions import (
    CommandExecutionError,
    InputValidationError,
    OutputValidationError,
    ToolNotFoundError,
)
from hifivar.logging_utils import get_logger
from hifivar.command_observation import _record_command_metadata


CommandArg = str | Path
CommandPath = str | Path

_LOGGER = get_logger(__name__)
_REDACTION_PLACEHOLDER = "***"
_LOG_OUTPUT_LIMIT = 2_000
_ERROR_OUTPUT_LIMIT = 500
_SENSITIVE_ENV_KEYS = frozenset(
    {
        "token", "password", "passwd", "secret", "api_key", "apikey",
        "access_key", "access_key_id", "credential", "credentials",
        "authorization", "private_key",
    }
)


@dataclass(frozen=True, slots=True, repr=False)
class CommandResult:
    """Stable result returned by :class:`CommandRunner`.

    Captured output is stored as text. When an output path is requested, the
    corresponding value is ``None`` because bytes are streamed directly to the
    file instead of being retained in memory. Dry runs have ``executed=False``
    and no return code or output. Raw attributes remain available to tool
    adapters; use ``to_dict()`` rather than ``dataclasses.asdict`` for diagnostics
    or provenance. ``str``/``repr`` and ``to_dict`` expose only a redacted copy.
    """

    args: tuple[str, ...]
    returncode: int | None
    stdout: str | None
    stderr: str | None
    duration_seconds: float
    cwd: Path | None
    executed: bool
    stdout_path: Path | None = None
    stderr_path: Path | None = None
    _safe_record: dict[str, object] | None = field(default=None, compare=False)

    def to_dict(self) -> dict[str, object]:
        """Return display-safe metadata without serializing an environment."""
        if self._safe_record is not None:
            return deepcopy(self._safe_record)
        # Manually constructed fake results have no verified redaction context.
        return {
            "args": ["<raw arguments omitted>"],
            "returncode": self.returncode,
            "stdout": "<captured output omitted>" if self.stdout is not None else None,
            "stderr": "<captured output omitted>" if self.stderr is not None else None,
            "duration_seconds": self.duration_seconds,
            "cwd": "<path omitted>" if self.cwd is not None else None,
            "executed": self.executed,
            "stdout_path": "<path omitted>" if self.stdout_path is not None else None,
            "stderr_path": "<path omitted>" if self.stderr_path is not None else None,
        }

    def __repr__(self) -> str:
        """Avoid dataclass repr exposing raw tool arguments or captured output."""
        fields = ", ".join(f"{key}={value!r}" for key, value in self.to_dict().items())
        return f"CommandResult({fields})"


class _CommandEnvironment(dict[str, str]):
    """Preserve child inheritance while omitting values from diagnostic repr."""

    def __repr__(self) -> str:
        return f"<CommandEnvironment: {len(self)} variables; values not displayed>"

    __str__ = __repr__


class CommandRunner:
    """Execute external programs without invoking a command shell."""

    def find_executable(
        self,
        executable: CommandArg,
        *,
        env: Mapping[str, str] | None = None,
        inherit_env: bool = True,
        cwd: CommandPath | None = None,
    ) -> Path | None:
        """Return the available executable path, or ``None`` when absent."""
        executable_text = _normalize_executable(executable)
        environment = _build_environment(env, inherit_env=inherit_env)
        working_directory = _validate_cwd(cwd)
        return _find_executable(
            executable_text,
            environment,
            working_directory,
        )

    def require_executable(
        self,
        executable: CommandArg,
        *,
        env: Mapping[str, str] | None = None,
        inherit_env: bool = True,
        cwd: CommandPath | None = None,
    ) -> Path:
        """Return an executable path or raise :class:`ToolNotFoundError`."""
        executable_text = _normalize_executable(executable)
        found = self.find_executable(
            executable_text, env=env, inherit_env=inherit_env, cwd=cwd
        )
        if found is None:
            redactions = _environment_redactions(
                _build_environment(env, inherit_env=inherit_env)
            )
            raise ToolNotFoundError(
                "Required executable was not found: "
                f"{_redact_text(executable_text, redactions)}"
            )
        return found

    def run(
        self,
        command: Sequence[CommandArg],
        *,
        cwd: CommandPath | None = None,
        env: Mapping[str, str] | None = None,
        inherit_env: bool = True,
        timeout: float | None = None,
        check: bool = True,
        dry_run: bool = False,
        redact_values: Collection[str] | None = None,
        stdout_path: CommandPath | None = None,
        stderr_path: CommandPath | None = None,
    ) -> CommandResult:
        """Run one command using an argument sequence and ``shell=False``.

        Environment values override a copy of the current process environment.
        With ``inherit_env=False``, only the explicit ``env`` mapping is passed
        to the child; use this for tools that serialize their environment.
        Captured output uses UTF-8 with invalid bytes replaced. Explicit output
        paths receive the process's original bytes directly, preventing large
        tool output from accumulating in Python memory. Output paths are opened
        relative to the caller's process, not relative to ``cwd``.

        When a caller supplies only ``stderr_path``, stdout is persisted to a
        sibling ``*.stdout.*`` log by default. An explicit ``stdout_path``
        continues to take precedence, including binary pipeline outputs.

        Values listed in ``redact_values`` and non-empty values of semantic
        secret environment keys are replaced in logs, errors and result display
        metadata. Raw arguments and captured output remain available to tool
        adapters. Explicit output files retain the original tool bytes and are
        not safe publication artifacts merely because diagnostic redaction is on.
        """
        args = _normalize_command(command)
        environment = _build_environment(env, inherit_env=inherit_env)
        redactions = _normalize_redactions(redact_values) | _environment_redactions(
            environment
        )
        display_command = format_command(args, redact_values=redactions)
        try:
            working_directory = _validate_cwd(cwd)
            timeout_seconds = _validate_timeout(timeout)
            normalized_stdout_path = _normalize_output_path(stdout_path)
            normalized_stderr_path = _normalize_output_path(stderr_path)
            if normalized_stdout_path is None and normalized_stderr_path is not None:
                normalized_stdout_path = _default_stdout_log_path(
                    normalized_stderr_path
                )
            _validate_distinct_output_paths(
                normalized_stdout_path,
                normalized_stderr_path,
            )
        except InputValidationError as error:
            raise InputValidationError(_redact_text(str(error), redactions)) from None

        _LOGGER.info("Executing command: %s", display_command)

        if dry_run:
            _LOGGER.debug("Dry run; command was not executed: %s", display_command)
            _record_command_metadata(_result_display_record(
                args, None, None, None, 0.0, working_directory, False,
                normalized_stdout_path, normalized_stderr_path, redactions,
            ))
            return CommandResult(
                args=args,
                returncode=None,
                stdout=None,
                stderr=None,
                duration_seconds=0.0,
                cwd=working_directory,
                executed=False,
                stdout_path=normalized_stdout_path,
                stderr_path=normalized_stderr_path,
                _safe_record=_result_display_record(
                    args, None, None, None, 0.0, working_directory, False,
                    normalized_stdout_path, normalized_stderr_path, redactions,
                ),
            )

        try:
            self.require_executable(
                args[0],
                env=environment,
                inherit_env=False,
                cwd=working_directory,
            )
        except ToolNotFoundError as error:
            raise ToolNotFoundError(_redact_text(str(error), redactions)) from None

        start_time = time.perf_counter()
        try:
            with ExitStack() as output_stack:
                stdout_target = _prepare_output_target(
                    output_stack,
                    normalized_stdout_path,
                )
                stderr_target = _prepare_output_target(
                    output_stack,
                    normalized_stderr_path,
                )
                completed = subprocess.run(
                    args,
                    cwd=working_directory,
                    env=environment,
                    stdout=stdout_target,
                    stderr=stderr_target,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    timeout=timeout_seconds,
                    check=False,
                    shell=False,
                )
        except subprocess.TimeoutExpired as error:
            duration = max(time.perf_counter() - start_time, 0.0)
            stderr = _coerce_timeout_output(error.stderr)
            summary = _summarize_output(
                stderr,
                redactions,
                limit=_ERROR_OUTPUT_LIMIT,
            )
            message = (
                f"Command timed out after {timeout_seconds:g} seconds: "
                f"{display_command}"
            )
            if summary:
                message += f". stderr: {summary}"
            elif normalized_stderr_path is not None:
                message += (
                    f". stderr was written to: "
                    f"{_redact_text(str(normalized_stderr_path), redactions)}"
                )
            if normalized_stdout_path is not None:
                message += (
                    f". stdout was written to: "
                    f"{_redact_text(str(normalized_stdout_path), redactions)}"
                )
            _LOGGER.error("%s (duration %.3f seconds)", message, duration)
            raise CommandExecutionError(message) from None
        except FileNotFoundError as error:
            message = (
                f"Executable disappeared or could not be started while running "
                f"command: {display_command}. {type(error).__name__}: "
                f"{_redact_text(str(error), redactions)}"
            )
            _LOGGER.error("%s", message)
            raise ToolNotFoundError(message) from None
        except OutputValidationError as error:
            message = _redact_text(str(error), redactions)
            _LOGGER.error("%s", message)
            raise OutputValidationError(message) from None
        except OSError as error:
            message = (
                f"Unable to execute command {display_command}: "
                f"{type(error).__name__}: {_redact_text(str(error), redactions)}"
            )
            _LOGGER.error("%s", message)
            raise CommandExecutionError(message) from None

        duration = max(time.perf_counter() - start_time, 0.0)
        stdout = completed.stdout
        stderr = completed.stderr
        _log_command_output(
            stdout,
            stderr,
            normalized_stdout_path,
            normalized_stderr_path,
            redactions,
        )

        _record_command_metadata(_result_display_record(
            args, completed.returncode, None, None, duration,
            working_directory, True, normalized_stdout_path,
            normalized_stderr_path, redactions,
        ))
        if completed.returncode != 0:
            message = _failure_message(
                display_command,
                completed.returncode,
                stdout,
                normalized_stdout_path,
                stderr,
                normalized_stderr_path,
                redactions,
            )
            if check:
                _LOGGER.error("%s", message)
                raise CommandExecutionError(message)
            _LOGGER.warning("%s", message)
        else:
            _LOGGER.debug(
                "Command completed with return code 0 in %.3f seconds: %s",
                duration,
                display_command,
            )

        return CommandResult(
            args=args,
            returncode=completed.returncode,
            stdout=stdout,
            stderr=stderr,
            duration_seconds=duration,
            cwd=working_directory,
            executed=True,
            stdout_path=normalized_stdout_path,
            stderr_path=normalized_stderr_path,
            _safe_record=_result_display_record(
                args, completed.returncode, stdout, stderr, duration,
                working_directory, True, normalized_stdout_path,
                normalized_stderr_path, redactions,
            ),
        )


def format_command(
    command: Sequence[CommandArg],
    *,
    redact_values: Collection[str] | None = None,
) -> str:
    """Create a redacted human-readable display, never an execution string."""
    args = _normalize_command(command)
    redactions = _normalize_redactions(redact_values)
    displayed_args = (
        _format_display_arg(_redact_text(arg, redactions)) for arg in args
    )
    return " ".join(displayed_args)


def _normalize_command(command: Sequence[CommandArg]) -> tuple[str, ...]:
    """Validate and copy command elements into an immutable string tuple."""
    if isinstance(command, (str, Path)):
        raise InputValidationError(
            "Command must be a non-empty sequence of string or Path arguments."
        )
    if not command:
        raise InputValidationError("Command must not be empty.")

    normalized: list[str] = []
    for index, value in enumerate(command):
        if not isinstance(value, (str, Path)):
            raise InputValidationError(
                f"Command argument at index {index} must be str or Path, "
                f"not {type(value).__name__}."
            )
        normalized.append(str(value))

    if not normalized[0].strip():
        raise InputValidationError("The first command argument must be an executable.")
    return tuple(normalized)


def _normalize_executable(executable: CommandArg) -> str:
    """Validate one executable argument."""
    if not isinstance(executable, (str, Path)):
        raise InputValidationError("Executable must be a string or Path.")
    normalized = str(executable)
    if not normalized.strip():
        raise InputValidationError("Executable must not be empty.")
    return normalized


def _build_environment(
    overrides: Mapping[str, str] | None,
    *,
    inherit_env: bool = True,
) -> dict[str, str]:
    """Return the inherited or isolated environment with validated overrides."""
    if not isinstance(inherit_env, bool):
        raise InputValidationError("inherit_env must be boolean.")
    environment = _CommandEnvironment(os.environ if inherit_env else {})
    if overrides is None:
        return environment

    for key, value in overrides.items():
        if not isinstance(key, str) or not isinstance(value, str):
            raise InputValidationError(
                "Environment overrides must map strings to strings."
            )
        environment[key] = value
    return environment


def _environment_redactions(environment: Mapping[str, str]) -> frozenset[str]:
    """Collect semantic-secret values without retaining or logging the mapping."""
    return frozenset(
        value
        for key, value in environment.items()
        if value and (
            (normalized := key.casefold().replace("-", "_")) in _SENSITIVE_ENV_KEYS
            or any(normalized.endswith(f"_{suffix}") for suffix in _SENSITIVE_ENV_KEYS)
        )
    )


def _result_display_record(
    args: tuple[str, ...],
    returncode: int | None,
    stdout: str | None,
    stderr: str | None,
    duration: float,
    cwd: Path | None,
    executed: bool,
    stdout_path: Path | None,
    stderr_path: Path | None,
    redactions: Collection[str],
) -> dict[str, object]:
    """Build a safe copy once; never store exact redaction values in a result."""
    return {
        "args": [_redact_text(value, redactions) for value in args],
        "returncode": returncode,
        "stdout": _redact_text(stdout, redactions) if stdout is not None else None,
        "stderr": _redact_text(stderr, redactions) if stderr is not None else None,
        "duration_seconds": duration,
        "cwd": _redact_text(str(cwd), redactions) if cwd is not None else None,
        "executed": executed,
        "stdout_path": _redact_text(str(stdout_path), redactions) if stdout_path is not None else None,
        "stderr_path": _redact_text(str(stderr_path), redactions) if stderr_path is not None else None,
    }


def _validate_cwd(cwd: CommandPath | None) -> Path | None:
    """Validate an existing working directory without creating it."""
    if cwd is None:
        return None
    if not isinstance(cwd, (str, Path)):
        raise InputValidationError("cwd must be a string, Path, or None.")

    working_directory = Path(cwd).expanduser()
    if not working_directory.exists():
        raise InputValidationError(
            f"Command working directory does not exist: {working_directory}"
        )
    if not working_directory.is_dir():
        raise InputValidationError(
            f"Command working directory is not a directory: {working_directory}"
        )
    return working_directory


def _validate_timeout(timeout: float | None) -> float | None:
    """Require timeout seconds to be positive when provided."""
    if timeout is None:
        return None
    if isinstance(timeout, bool) or not isinstance(timeout, (int, float)):
        raise InputValidationError("timeout must be a positive number of seconds.")
    if timeout <= 0:
        raise InputValidationError("timeout must be greater than zero seconds.")
    return float(timeout)


def _normalize_redactions(
    redact_values: Collection[str] | None,
) -> frozenset[str]:
    """Validate exact secret values used only for display redaction."""
    if redact_values is None:
        return frozenset()

    normalized: set[str] = set()
    for value in redact_values:
        if not isinstance(value, str):
            raise InputValidationError("Redacted values must be strings.")
        if value:
            normalized.add(value)
    return frozenset(normalized)


def _normalize_output_path(path: CommandPath | None) -> Path | None:
    """Normalize an optional process-output destination."""
    if path is None:
        return None
    if not isinstance(path, (str, Path)):
        raise InputValidationError("Output path must be a string, Path, or None.")
    return Path(path).expanduser()


def _default_stdout_log_path(stderr_path: Path) -> Path:
    """Derive a sibling stdout log when a caller supplies only stderr_path."""
    suffix = stderr_path.suffix
    if suffix:
        return stderr_path.with_name(
            f"{stderr_path.stem}.stdout{suffix}"
        )
    return stderr_path.with_name(f"{stderr_path.name}.stdout.log")


def _validate_distinct_output_paths(
    stdout_path: Path | None,
    stderr_path: Path | None,
) -> None:
    """Reject ambiguous attempts to open one file twice for process output."""
    if stdout_path is None or stderr_path is None:
        return
    stdout_identity = os.path.normcase(str(stdout_path.absolute()))
    stderr_identity = os.path.normcase(str(stderr_path.absolute()))
    if stdout_identity == stderr_identity:
        raise InputValidationError(
            "stdout_path and stderr_path must refer to different files."
        )


def _find_executable(
    executable: str,
    environment: Mapping[str, str],
    cwd: Path | None,
) -> Path | None:
    """Locate a command, accounting for explicit paths relative to ``cwd``."""
    candidate = Path(executable)
    lookup = executable
    if candidate.name != executable and not candidate.is_absolute() and cwd:
        lookup = str(cwd / candidate)

    found = shutil.which(lookup, path=environment.get("PATH"))
    return Path(found) if found is not None else None


def _prepare_output_target(
    stack: ExitStack,
    output_path: Path | None,
) -> int | BinaryIO:
    """Return a PIPE or a binary file receiving original process bytes."""
    if output_path is None:
        return subprocess.PIPE
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        return stack.enter_context(output_path.open("wb"))
    except OSError as error:
        raise OutputValidationError(
            f"Unable to open command output file '{output_path}': {error}"
        ) from error


def _coerce_timeout_output(output: str | bytes | None) -> str | None:
    """Normalize TimeoutExpired output across Python implementations."""
    if isinstance(output, bytes):
        return output.decode("utf-8", errors="replace")
    return output


def _failure_message(
    display_command: str,
    returncode: int,
    stdout: str | None,
    stdout_path: Path | None,
    stderr: str | None,
    stderr_path: Path | None,
    redactions: Collection[str],
) -> str:
    """Build a concise, redacted non-zero-return diagnostic."""
    message = f"Command failed with return code {returncode}: {display_command}"
    stderr_summary = _summarize_output(
        stderr,
        redactions,
        limit=_ERROR_OUTPUT_LIMIT,
    )
    if stderr_summary:
        message += f". stderr: {stderr_summary}"
    elif stderr_path is not None:
        safe_stderr_path = _redact_text(str(stderr_path), redactions)
        message += f". stderr was written to: {safe_stderr_path}"
    stdout_summary = _summarize_output(
        stdout,
        redactions,
        limit=_ERROR_OUTPUT_LIMIT,
    )
    if stdout_summary:
        message += f". stdout: {stdout_summary}"
    elif stdout_path is not None:
        safe_stdout_path = _redact_text(str(stdout_path), redactions)
        message += f". stdout was written to: {safe_stdout_path}"
    return message


def _log_command_output(
    stdout: str | None,
    stderr: str | None,
    stdout_path: Path | None,
    stderr_path: Path | None,
    redactions: Collection[str],
) -> None:
    """Log bounded output summaries without changing returned content."""
    if stdout_path is not None:
        _LOGGER.debug(
            "Command stdout was written to %s",
            _redact_text(str(stdout_path), redactions),
        )
    elif stdout:
        _LOGGER.debug(
            "Command stdout: %s",
            _summarize_output(stdout, redactions, limit=_LOG_OUTPUT_LIMIT),
        )

    if stderr_path is not None:
        _LOGGER.debug(
            "Command stderr was written to %s",
            _redact_text(str(stderr_path), redactions),
        )
    elif stderr:
        _LOGGER.debug(
            "Command stderr: %s",
            _summarize_output(stderr, redactions, limit=_LOG_OUTPUT_LIMIT),
        )


def _summarize_output(
    output: str | None,
    redactions: Collection[str],
    *,
    limit: int,
) -> str:
    """Redact, flatten, and bound output for safe diagnostic logging."""
    if not output:
        return ""
    summary = _redact_text(output, redactions).strip()
    summary = summary.replace("\r", "\\r").replace("\n", "\\n")
    if len(summary) > limit:
        return f"{summary[:limit]}... [truncated]"
    return summary


def _redact_text(text: str, redactions: Collection[str]) -> str:
    """Replace sensitive substrings in display text only."""
    redacted = text
    for value in sorted(redactions, key=len, reverse=True):
        redacted = redacted.replace(value, _REDACTION_PLACEHOLDER)
    return redacted


def _format_display_arg(argument: str) -> str:
    """Keep simple args readable while visibly delimiting ambiguous values."""
    if not argument or any(character.isspace() for character in argument):
        return repr(argument)
    if "'" in argument or '"' in argument:
        return repr(argument)
    return argument


__all__ = [
    "CommandArg",
    "CommandResult",
    "CommandRunner",
    "format_command",
]
