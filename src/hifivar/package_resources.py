"""Locate installed HiFiVar runtime resources without source-tree assumptions."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, distribution
from pathlib import Path

from hifivar.exceptions import ConfigurationError


def installed_workflow_root() -> Path:
    """Return the packaged workflow root for wheel or editable installs."""
    try:
        files = distribution("hifivar").files or ()
    except PackageNotFoundError as error:
        raise ConfigurationError("HiFiVar distribution metadata is unavailable.") from error
    for item in files:
        normalized = str(item).replace("\\", "/")
        if normalized.endswith("hifivar/workflow/Snakefile") or normalized == "workflow/Snakefile":
            candidate = Path(item.locate()).parent
            if candidate.joinpath("Snakefile").is_file():
                return candidate
    source_candidate = Path(__file__).resolve().parents[2] / "workflow"
    if source_candidate.joinpath("Snakefile").is_file():
        return source_candidate
    raise ConfigurationError("Packaged HiFiVar workflow resources are missing.")


__all__ = ["installed_workflow_root"]
