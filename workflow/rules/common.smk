"""Shared Phase 0 workflow configuration and resource conventions."""

from pathlib import Path

from hifivar.exceptions import WorkflowError


DEFAULT_WORKDIR = "work"
DEFAULT_OUTDIR = "results"
DEFAULT_MEM_MB = 64
DEFAULT_RUNTIME_MIN = 1


def _required_mapping(section_name):
    """Return one required effective-config section."""
    section = config.get(section_name)
    if not isinstance(section, dict):
        raise WorkflowError(
            f"Effective config requires a mapping section: {section_name}."
        )
    return section


def _required_text(section, key, field_name):
    """Return one required non-empty string value."""
    value = section.get(key)
    if not isinstance(value, str) or not value.strip():
        raise WorkflowError(
            f"Effective config requires a non-empty string: {field_name}."
        )
    return value


def _configured_path(value, default, field_name):
    """Use a relative workflow default only when effective config has null."""
    if value is None:
        return Path(default)
    if not isinstance(value, str) or not value.strip():
        raise WorkflowError(
            f"Effective config path must be a non-empty string or null: "
            f"{field_name}."
        )
    return Path(value).expanduser()


PROJECT_CONFIG = _required_mapping("project")
RUNTIME_CONFIG = _required_mapping("runtime")
PATH_CONFIG = _required_mapping("paths")
WORKFLOW_CONFIG = _required_mapping("workflow")

PROJECT_NAME = _required_text(PROJECT_CONFIG, "name", "project.name")
WORKFLOW_PRESET = _required_text(WORKFLOW_CONFIG, "preset", "workflow.preset")

RUNTIME_THREADS = RUNTIME_CONFIG.get("threads")
if (
    not isinstance(RUNTIME_THREADS, int)
    or isinstance(RUNTIME_THREADS, bool)
    or RUNTIME_THREADS <= 0
):
    raise WorkflowError(
        "Effective config requires runtime.threads to be a positive integer."
    )

WORK_ROOT = _configured_path(
    PATH_CONFIG.get("workdir"),
    DEFAULT_WORKDIR,
    "paths.workdir",
)
OUTPUT_ROOT = _configured_path(
    PATH_CONFIG.get("outdir"),
    DEFAULT_OUTDIR,
    "paths.outdir",
)
LOG_ROOT = Path("logs")

PHASE0_PREPARE_MARKER = (WORK_ROOT / "phase0" / "config_ready.done").as_posix()
PHASE0_SMOKE_MARKER = (
    OUTPUT_ROOT / "phase0" / "snakemake_smoke.done"
).as_posix()
PHASE0_PREPARE_LOG = (LOG_ROOT / "phase0" / "phase0_prepare.log").as_posix()
PHASE0_SMOKE_LOG = (LOG_ROOT / "phase0" / "phase0_smoke.log").as_posix()
