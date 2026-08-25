"""Integration tests for the Phase 0.9 Snakemake infrastructure."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from hifivar.config import load_config, write_effective_config


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = REPOSITORY_ROOT / "workflow" / "Snakefile"
COMMON_RULES = REPOSITORY_ROOT / "workflow" / "rules" / "common.smk"
SMOKE_RULES = REPOSITORY_ROOT / "workflow" / "rules" / "smoke.smk"
SMALL_RULES = REPOSITORY_ROOT / "workflow" / "rules" / "small.smk"
MARKER_SCRIPT = (
    REPOSITORY_ROOT / "workflow" / "scripts" / "write_phase0_marker.py"
)
DEFAULT_CONFIG = (
    REPOSITORY_ROOT / "src" / "hifivar" / "resources" / "configs" / "default.yaml"
)
STANDARD_PRESET = (
    REPOSITORY_ROOT
    / "src"
    / "hifivar"
    / "resources"
    / "configs"
    / "presets"
    / "standard.yaml"
)


def _snakemake_executable() -> Path:
    """Locate Snakemake beside the active Python or on PATH."""
    executable_name = "snakemake.exe" if sys.platform == "win32" else "snakemake"
    adjacent = Path(sys.executable).with_name(executable_name)
    if adjacent.is_file():
        return adjacent

    discovered = shutil.which("snakemake")
    if discovered is None:
        pytest.skip("Snakemake is not installed; install the workflow extra.")
    return Path(discovered)


def _write_effective_config(
    root: Path,
    *,
    workdir: Path | None = None,
    outdir: Path | None = None,
) -> Path:
    """Generate an effective config through the production Config API."""
    user_config: Path | None = None
    if workdir is not None or outdir is not None:
        user_config = root / "user_config.yaml"
        user_config.write_text(
            yaml.safe_dump(
                {
                    "paths": {
                        "workdir": str(workdir) if workdir is not None else None,
                        "outdir": str(outdir) if outdir is not None else None,
                    }
                },
                allow_unicode=True,
                sort_keys=False,
            ),
            encoding="utf-8",
        )

    effective = load_config(
        DEFAULT_CONFIG,
        STANDARD_PRESET,
        user_config,
    )
    destination = root / "effective_config.yaml"
    write_effective_config(effective, destination)
    return destination


def _run_snakemake(
    cwd: Path,
    config_path: Path,
    *extra_arguments: str,
) -> subprocess.CompletedProcess[str]:
    """Invoke the public Snakemake CLI from an isolated working directory."""
    command = [
        str(_snakemake_executable()),
        "--snakefile",
        str(SNAKEFILE),
        "--configfile",
        str(config_path),
        "--cores",
        "1",
        "--shared-fs-usage",
        "input-output",
        "persistence",
        "software-deployment",
        "software-deployment-cache",
        "sources",
        "storage-local-copies",
        *extra_arguments,
    ]
    environment = os.environ.copy()
    cache_root = cwd / ".cache"
    environment["XDG_CACHE_HOME"] = str(cache_root)
    if sys.platform == "win32":
        environment["LOCALAPPDATA"] = str(cwd / ".local-app-data")
        environment["APPDATA"] = str(cwd / ".app-data")
    return subprocess.run(
        command,
        cwd=cwd,
        env=environment,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )


def _combined_output(result: subprocess.CompletedProcess[str]) -> str:
    """Combine CLI streams for stable cross-version diagnostics."""
    return f"{result.stdout}\n{result.stderr}"


def test_phase0_workflow_files_exist() -> None:
    """The entry point, included rules, env location, and config source exist."""
    assert SNAKEFILE.is_file()
    assert COMMON_RULES.is_file()
    assert SMOKE_RULES.is_file()
    assert SMALL_RULES.is_file()
    assert MARKER_SCRIPT.is_file()
    assert (REPOSITORY_ROOT / "workflow" / "envs").is_dir()
    assert DEFAULT_CONFIG.is_file()


def test_snakefile_keeps_modular_phase0_and_phase3_rule_modules() -> None:
    """The entry point should stay small and use paths relative to itself."""
    snakefile_text = SNAKEFILE.read_text(encoding="utf-8")

    assert 'include: "rules/common.smk"' in snakefile_text
    assert 'include: "rules/smoke.smk"' in snakefile_text
    assert 'include: "rules/small.smk"' in snakefile_text
    assert "rule all:" in snakefile_text
    assert "rule phase0" not in snakefile_text
    assert "sys.path" not in snakefile_text


def test_workflow_contains_no_out_of_scope_callers() -> None:
    """Phase 3 may add DeepVariant but no Phase 4+ callers."""
    workflow_text = "\n".join(
        path.read_text(encoding="utf-8")
        for path in (SNAKEFILE, COMMON_RULES, SMOKE_RULES, SMALL_RULES)
    ).lower()
    forbidden = (
        "sawfish_discover",
        "sniffles",
        "pbsv",
        "cutesv",
        "trgt_genotype",
        "hifiasm_assemble",
    )

    for text in forbidden:
        assert text not in workflow_text


def test_snakemake_dry_run_builds_phase0_dag(tmp_path: Path) -> None:
    """A real Snakemake dry-run should resolve includes and the smoke DAG."""
    effective_config = _write_effective_config(tmp_path)

    result = _run_snakemake(tmp_path, effective_config, "--dry-run")

    assert result.returncode == 0, _combined_output(result)
    output = _combined_output(result)
    assert "phase0_prepare" in output
    assert "phase0_smoke" in output
    assert "snakemake_smoke.done" in output
    assert "deepvariant_small" not in output
    assert not (tmp_path / "results" / "phase0" / "snakemake_smoke.done").exists()


def test_snakemake_executes_marker_with_custom_paths(tmp_path: Path) -> None:
    """Execution should honor workdir/outdir and write deterministic content."""
    configured_workdir = Path("custom-work")
    configured_outdir = Path("custom-results")
    effective_config = _write_effective_config(
        tmp_path,
        workdir=configured_workdir,
        outdir=configured_outdir,
    )

    result = _run_snakemake(tmp_path, effective_config)

    assert result.returncode == 0, _combined_output(result)
    custom_workdir = tmp_path / configured_workdir
    custom_outdir = tmp_path / configured_outdir
    prepare_marker = custom_workdir / "phase0" / "config_ready.done"
    smoke_marker = custom_outdir / "phase0" / "snakemake_smoke.done"
    assert prepare_marker.read_text(encoding="utf-8") == (
        "HiFiVar effective config accepted\n"
    )
    assert smoke_marker.read_text(encoding="utf-8") == (
        "HiFiVar Snakemake infrastructure OK\n"
        "project: hifivar\n"
        "preset: standard\n"
    )
    assert (tmp_path / "logs" / "phase0" / "phase0_prepare.log").is_file()
    assert (tmp_path / "logs" / "phase0" / "phase0_smoke.log").is_file()


def test_snakemake_dry_run_from_non_repository_cwd(tmp_path: Path) -> None:
    """Absolute Snakefile use must not make includes depend on repository cwd."""
    assert REPOSITORY_ROOT not in (tmp_path, *tmp_path.parents)
    effective_config = _write_effective_config(tmp_path)

    result = _run_snakemake(tmp_path, effective_config, "--dry-run")

    assert result.returncode == 0, _combined_output(result)
    assert "phase0_smoke" in _combined_output(result)
