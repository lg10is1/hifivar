"""End-to-end release-candidate checks for the complete Phase 0 foundation."""

from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
import sys
from collections.abc import Iterator
from contextlib import ExitStack
from importlib.metadata import version
from importlib.resources import as_file, files
from pathlib import Path

import pytest
import yaml

import hifivar
from hifivar.command import CommandRunner
from hifivar.config import load_config, write_effective_config
from hifivar.exceptions import CommandExecutionError, ReferenceError
from hifivar.logging_utils import (
    HIFIVAR_LOGGER_NAME,
    configure_logging,
    get_logger,
)
from hifivar.validation import (
    compute_sha256,
    validate_bed,
    validate_contig_compatibility,
    validate_fasta,
    validate_fastq,
    validate_vcf,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = REPOSITORY_ROOT / "workflow" / "Snakefile"


@pytest.fixture(autouse=True)
def restore_hifivar_logging() -> Iterator[None]:
    """Keep integration logging independent of all other test ordering."""
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


def _installed_executable(name: str, *, optional: bool = False) -> Path:
    """Locate an installed console script beside Python or on PATH."""
    executable_name = f"{name}.exe" if sys.platform == "win32" else name
    adjacent = Path(sys.executable).with_name(executable_name)
    if adjacent.is_file():
        return adjacent

    discovered = shutil.which(name)
    if discovered is not None:
        return Path(discovered)
    if optional:
        pytest.skip(f"{name} is not installed; install the workflow extra.")
    pytest.fail(f"Required installed console script was not found: {name}")


def _run_cli(
    command: list[str],
    *,
    cwd: Path,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run one installed CLI command with deterministic UTF-8 capture."""
    return subprocess.run(
        command,
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )


def _run_snakemake(
    cwd: Path,
    config_path: Path,
    *,
    dry_run: bool,
) -> subprocess.CompletedProcess[str]:
    """Run Snakemake without writing shared user cache during tests."""
    command = [
        str(_installed_executable("snakemake", optional=True)),
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
    ]
    if dry_run:
        command.append("--dry-run")

    environment = os.environ.copy()
    environment["XDG_CACHE_HOME"] = str(cwd / ".cache")
    if sys.platform == "win32":
        environment["LOCALAPPDATA"] = str(cwd / ".local-app-data")
        environment["APPDATA"] = str(cwd / ".app-data")
    return _run_cli(command, cwd=cwd, env=environment)


def _assert_success(result: subprocess.CompletedProcess[str]) -> None:
    """Expose both subprocess streams when an integration command fails."""
    assert result.returncode == 0, f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"


def test_phase0_end_to_end_release_candidate(tmp_path: Path) -> None:
    """Exercise the installed Phase 0 foundation as one isolated user flow."""
    assert hifivar.__version__ == version("hifivar")
    assert re.fullmatch(
        r"\d+\.\d+\.\d+(?:(?:a|b|rc)\d+|\.dev\d+)?",
        hifivar.__version__,
    )

    unicode_root = tmp_path / "科研测试"
    outside_repo = unicode_root / "outside_repo"
    input_root = unicode_root / "样本测试"
    outside_repo.mkdir(parents=True)
    input_root.mkdir()
    assert REPOSITORY_ROOT not in (outside_repo, *outside_repo.parents)

    fasta = input_root / "tiny.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text(
        "chr1\t4\t6\t4\t5\n",
        encoding="utf-8",
    )
    fastq = input_root / "tiny.fastq"
    fastq.write_text("@read1\nACGT\n+\nIIII\n", encoding="utf-8")
    vcf = input_root / "tiny.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "##contig=<ID=chr1,length=4>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t2\t.\tC\tT\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    bed = input_root / "tiny.bed"
    bed.write_text("chr1\t0\t4\n", encoding="utf-8")

    log_path = unicode_root / "logs" / "样本测试.log"
    user_config = unicode_root / "user_config.yaml"
    user_config.write_text(
        yaml.safe_dump(
            {
                "project": {"name": "phase0-e2e"},
                "runtime": {"threads": 2},
                "paths": {
                    "workdir": "工作目录",
                    "outdir": "结果目录",
                },
                "logging": {
                    "level": "DEBUG",
                    "file": str(log_path),
                },
            },
            allow_unicode=True,
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    hifivar_executable = _installed_executable("hifivar")
    version_result = _run_cli(
        [str(hifivar_executable), "--version"],
        cwd=outside_repo,
    )
    _assert_success(version_result)
    assert version_result.stdout.strip() == f"hifivar {hifivar.__version__}"

    module_result = _run_cli(
        [sys.executable, "-m", "hifivar", "--version"],
        cwd=outside_repo,
    )
    _assert_success(module_result)
    assert module_result.stdout.strip() == f"hifivar {hifivar.__version__}"

    show_result = _run_cli(
        [
            str(hifivar_executable),
            "--config",
            str(user_config),
            "--preset",
            "standard",
            "config",
            "show",
        ],
        cwd=outside_repo,
    )
    _assert_success(show_result)
    shown_config = yaml.safe_load(show_result.stdout)
    assert shown_config["project"]["name"] == "phase0-e2e"
    assert shown_config["runtime"]["threads"] == 2
    assert shown_config["workflow"]["preset"] == "standard"

    config_root = files("hifivar").joinpath("resources", "configs")
    with ExitStack() as stack:
        default_config = Path(
            stack.enter_context(as_file(config_root.joinpath("default.yaml")))
        )
        preset_config = Path(
            stack.enter_context(
                as_file(config_root.joinpath("presets", "standard.yaml"))
            )
        )
        effective = load_config(default_config, preset_config, user_config)

    assert effective["project"]["name"] == "phase0-e2e"  # type: ignore[index]
    assert effective["runtime"]["threads"] == 2  # type: ignore[index]
    assert effective["workflow"]["preset"] == "standard"  # type: ignore[index]
    effective_config = outside_repo / "effective_config.yaml"
    write_effective_config(effective, effective_config)
    assert effective_config.is_file()

    logging_config = effective["logging"]
    assert isinstance(logging_config, dict)
    configure_logging(
        level=str(logging_config["level"]),
        log_file=str(logging_config["file"]),
    )
    get_logger("phase0_e2e").info("Phase 0 integration test – 样本测试")
    for handler in logging.getLogger(HIFIVAR_LOGGER_NAME).handlers:
        handler.flush()
    log_text = log_path.read_text(encoding="utf-8")
    assert "Phase 0 integration test – 样本测试" in log_text

    assert validate_fasta(fasta, require_fai=True) == fasta
    assert validate_fastq(fastq) == fastq
    assert validate_vcf(vcf) == vcf
    assert validate_bed(bed) == bed
    validate_contig_compatibility(["chr1"], ["chr1"])
    with pytest.raises(ReferenceError, match="REFERENCE_CONTIG_MISMATCH"):
        validate_contig_compatibility(["chr1"], ["1"])

    checksum = compute_sha256(fasta)
    assert re.fullmatch(r"[0-9a-f]{64}", checksum)

    runner = CommandRunner()
    command_result = runner.run(
        [sys.executable, "-c", "print('HiFiVar Phase 0 CommandRunner OK')"],
        cwd=outside_repo,
    )
    assert command_result.executed is True
    assert command_result.returncode == 0
    assert command_result.duration_seconds >= 0
    assert command_result.stdout is not None
    assert command_result.stdout.strip() == "HiFiVar Phase 0 CommandRunner OK"
    with pytest.raises(CommandExecutionError, match="return code 3"):
        runner.run([sys.executable, "-c", "import sys; sys.exit(3)"])

    dry_run = _run_snakemake(outside_repo, effective_config, dry_run=True)
    _assert_success(dry_run)
    dry_run_output = f"{dry_run.stdout}\n{dry_run.stderr}"
    assert "phase0_prepare" in dry_run_output
    assert "phase0_smoke" in dry_run_output

    smoke_run = _run_snakemake(outside_repo, effective_config, dry_run=False)
    _assert_success(smoke_run)
    smoke_marker = outside_repo / "结果目录" / "phase0" / "snakemake_smoke.done"
    assert smoke_marker.read_text(encoding="utf-8") == (
        "HiFiVar Snakemake infrastructure OK\n"
        "project: phase0-e2e\n"
        "preset: standard\n"
    )
    assert (outside_repo / "工作目录" / "phase0" / "config_ready.done").is_file()
