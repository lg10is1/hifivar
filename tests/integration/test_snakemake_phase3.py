"""Snakemake dry-run coverage for the Phase 3 DeepVariant rule."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from hifivar.config import load_config, write_effective_config


ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = ROOT / "workflow" / "Snakefile"
DEFAULT_CONFIG = ROOT / "src" / "hifivar" / "resources" / "configs" / "default.yaml"


def _snakemake() -> Path:
    name = "snakemake.exe" if sys.platform == "win32" else "snakemake"
    adjacent = Path(sys.executable).with_name(name)
    if adjacent.is_file():
        return adjacent
    found = shutil.which("snakemake")
    if found is None:
        pytest.skip("Snakemake is not installed; install the workflow extra.")
    return Path(found)


def test_deepvariant_rule_dry_run_has_deterministic_handoff(tmp_path: Path) -> None:
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    bam = tmp_path / "S1.bam"
    bam.write_bytes(b"alignment")
    Path(f"{bam}.bai").write_bytes(b"index")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        f"sample_id\tinput\tinput_type\nS1\t{bam}\tbam\n",
        encoding="utf-8",
    )
    user = tmp_path / "phase3.yaml"
    user.write_text(
        yaml.safe_dump(
            {
                "reference": {"fasta": str(fasta), "build": "GRCh38"},
                "samples": {"sheet": str(sheet)},
                "paths": {
                    "workdir": str(tmp_path / "work"),
                    "outdir": str(tmp_path / "results"),
                },
                "small": {"enabled": True, "threads": 7, "memory_mb": 12345},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    config = load_config(DEFAULT_CONFIG, user_config=user)
    effective = tmp_path / "effective.yaml"
    write_effective_config(config, effective)
    environment = os.environ.copy()
    environment["XDG_CACHE_HOME"] = str(tmp_path / ".cache")
    if sys.platform == "win32":
        environment["LOCALAPPDATA"] = str(tmp_path / ".local")
        environment["APPDATA"] = str(tmp_path / ".app")
    result = subprocess.run(
        [
            str(_snakemake()),
            "--snakefile",
            str(SNAKEFILE),
            "--configfile",
            str(effective),
            "--cores",
            "7",
            "--dry-run",
            "--printshellcmds",
            "--shared-fs-usage",
            "input-output",
            "persistence",
            "software-deployment",
            "software-deployment-cache",
            "sources",
            "storage-local-copies",
        ],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
        timeout=60,
    )
    output = f"{result.stdout}\n{result.stderr}"
    assert result.returncode == 0, output
    assert "deepvariant_small" in output
    assert "S1.small.vcf.gz" in output
    assert "S1.g.vcf.gz" in output
    assert "threads: 7" in output
    assert "deepvariant_slots=1" in output
    assert not (tmp_path / "results" / "small" / "S1.small.vcf.gz").exists()


def test_phase3_rule_and_bridge_keep_command_ownership_in_wrapper() -> None:
    rule_text = (ROOT / "workflow" / "rules" / "small.smk").read_text(encoding="utf-8")
    bridge_text = (ROOT / "workflow" / "scripts" / "run_deepvariant.py").read_text(encoding="utf-8")
    assert "/opt/deepvariant" not in rule_text
    assert "--model_type" not in rule_text
    assert "DeepVariantWrapper" in bridge_text
    assert "temporary_directory" in bridge_text
    assert "deepvariant_slots=1" in rule_text
    assert "subprocess" not in bridge_text
    assert "os.system" not in bridge_text
