"""Shared standard-type conversion and atomic UTF-8 report writers."""

from __future__ import annotations

import json
import os
import tempfile
from collections.abc import Mapping
from copy import deepcopy
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any

import yaml

from hifivar.exceptions import InputValidationError, OutputValidationError


_REDACTED = "***"
_SENSITIVE_KEYS = frozenset({"token", "password", "secret", "api_key", "apikey"})


def redact_sensitive_data(value: object) -> Any:
    """Recursively copy data while redacting semantic secret-key variants."""
    if isinstance(value, Mapping):
        redacted: dict[str, object] = {}
        for key, item in value.items():
            key_text = str(key)
            normalized = key_text.casefold().replace("-", "_")
            sensitive = normalized in _SENSITIVE_KEYS or normalized.endswith(
                ("_token", "_password", "_secret", "_api_key", "_apikey")
            )
            redacted[key_text] = _REDACTED if sensitive else redact_sensitive_data(item)
        return redacted
    if isinstance(value, (list, tuple)):
        return [redact_sensitive_data(item) for item in value]
    return deepcopy(value)


def standardize_data(value: object, *, context: str = "Serialized value") -> Any:
    """Convert nested values to JSON/YAML containers and primitive scalars."""
    if isinstance(value, Enum):
        return standardize_data(value.value, context=context)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {
            str(key): standardize_data(item, context=context)
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [standardize_data(item, context=context) for item in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return deepcopy(value)
    raise InputValidationError(
        f"{context} has unsupported type: {type(value).__name__}."
    )


def utc_now_iso8601() -> str:
    """Return the current UTC time as an ISO 8601 string ending in ``Z``."""
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def write_json_atomic(
    payload: object,
    path: str | Path,
    *,
    overwrite: bool = False,
    artifact_name: str = "Output",
) -> Path:
    """Serialize standard data and atomically write indented UTF-8 JSON."""
    try:
        serialized = json.dumps(
            standardize_data(payload, context=f"{artifact_name} value"),
            ensure_ascii=False,
            indent=2,
        ) + "\n"
    except (TypeError, ValueError) as error:
        raise OutputValidationError(
            f"Unable to serialize {artifact_name.lower()} as JSON: {error}"
        ) from error
    return _atomic_write(
        path,
        serialized,
        overwrite=overwrite,
        artifact_name=artifact_name,
    )


def write_yaml_atomic(
    payload: object,
    path: str | Path,
    *,
    overwrite: bool = False,
    artifact_name: str = "Output",
) -> Path:
    """Serialize standard data and atomically write UTF-8 YAML."""
    try:
        serialized = yaml.safe_dump(
            standardize_data(payload, context=f"{artifact_name} value"),
            allow_unicode=True,
            default_flow_style=False,
            sort_keys=False,
        )
    except yaml.YAMLError as error:
        raise OutputValidationError(
            f"Unable to serialize {artifact_name.lower()} as YAML: {error}"
        ) from error
    return _atomic_write(
        path,
        serialized,
        overwrite=overwrite,
        artifact_name=artifact_name,
    )


def _atomic_write(
    path: str | Path,
    content: str,
    *,
    overwrite: bool,
    artifact_name: str,
) -> Path:
    """Write through an owned sibling temporary file and replace atomically."""
    destination = Path(path).expanduser()
    if destination.exists() and not overwrite:
        raise OutputValidationError(
            f"{artifact_name} output already exists: '{destination}'. "
            "Pass overwrite=True to replace it."
        )
    temporary_path: Path | None = None
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            mode="wt",
            encoding="utf-8",
            newline="\n",
            prefix=f".{destination.name}.",
            suffix=".tmp",
            dir=destination.parent,
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, destination)
    except OSError as error:
        if temporary_path is not None:
            try:
                temporary_path.unlink(missing_ok=True)
            except OSError:
                pass
        raise OutputValidationError(
            f"Unable to write {artifact_name.lower()} '{destination}': {error}"
        ) from error
    return destination


__all__ = [
    "redact_sensitive_data",
    "standardize_data",
    "utc_now_iso8601",
    "write_json_atomic",
    "write_yaml_atomic",
]
