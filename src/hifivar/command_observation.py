"""Opt-in, task-local observation of safe CommandRunner result metadata.

No environment, captured stdout or stderr is retained. Observation does not
change subprocess execution or turn a planned command into an executed one.
"""

from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from copy import deepcopy
from dataclasses import dataclass
from typing import Iterator

from hifivar.exceptions import OutputValidationError


_WITNESS = object()


@dataclass(frozen=True, slots=True, init=False)
class ObservedCommand:
    """Opaque internal runner witness, not a user-supplied success dictionary.

    This prevents accidental fabrication/mutation, not hostile Python running
    inside the same interpreter. It is not a signature or a sandbox boundary.
    """

    _record: dict
    _witness: object

    def to_dict(self) -> dict:
        if self._witness is not _WITNESS:
            raise OutputValidationError("Unobserved command result cannot certify execution.")
        return deepcopy(self._record)


class CommandObservation:
    """Read-only consumer view of task-local native result witnesses."""

    def __init__(self) -> None:
        self._records: list[ObservedCommand] = []

    def __iter__(self):
        return iter(tuple(self._records))

    def __len__(self) -> int:
        return len(self._records)


_CURRENT: ContextVar[CommandObservation | None] = ContextVar("hifivar_commands", default=None)


@contextmanager
def observe_commands() -> Iterator[CommandObservation]:
    """Collect one native step, without sharing mutable state across tasks."""
    records = CommandObservation()
    token = _CURRENT.set(records)
    try:
        yield records
    finally:
        _CURRENT.reset(token)


def _record_command_metadata(safe_record: dict) -> None:
    """Accept only an already-redacted CommandRunner display record."""
    records = _CURRENT.get()
    if records is None:
        return
    if len(records) >= 1024:
        raise OutputValidationError("Attempt command observation exceeds 1024 commands.")
    keys = (
        "args", "returncode", "duration_seconds", "cwd", "executed",
        "stdout_path", "stderr_path",
    )
    witness = object.__new__(ObservedCommand)
    object.__setattr__(witness, "_record", deepcopy({key: safe_record.get(key) for key in keys}))
    object.__setattr__(witness, "_witness", _WITNESS)
    records._records.append(witness)


def is_current_observation(records: CommandObservation) -> bool:
    """Refuse accidental reuse of witnesses from another step or task."""
    return _CURRENT.get() is records


__all__ = ["observe_commands"]
