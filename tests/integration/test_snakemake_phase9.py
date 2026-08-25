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
DEFAULT = ROOT / "src" / "hifivar" / "resources" / "configs" / "default.yaml"


def snakemake():
    name = "snakemake.exe" if sys.platform == "win32" else "snakemake"
    adjacent = Path(sys.executable).with_name(name)
    found = adjacent if adjacent.is_file() else shutil.which("snakemake")
    if found is None:
        pytest.skip("Snakemake is not installed.")
    return Path(found)


def test_phase9_read_only_dry_run_is_deterministic(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\nACGT\n")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n")
    bam = tmp_path / "S1.bam"
    bam.write_bytes(b"BAM")
    Path(f"{bam}.bai").write_bytes(b"BAI")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(f"sample_id\tinput\tinput_type\nS1\t{bam}\tbam\n")
    user = tmp_path / "phase9.yaml"
    user.write_text(
        yaml.safe_dump(
            {
                "reference": {"fasta": str(fasta), "build": "GRCh38"},
                "samples": {"sheet": str(sheet)},
                "paths": {
                    "workdir": str(tmp_path / "work"),
                    "outdir": str(tmp_path / "results"),
                },
                "sv": {
                    "enabled": True,
                    "sawfish": {"enabled": True},
                    "sniffles2": {"enabled": False},
                    "pbsv": {"enabled": False},
                    "cutesv": {"enabled": False},
                    "harmonization": {"enabled": True},
                },
            },
            sort_keys=False,
        )
    )
    effective = tmp_path / "effective.yaml"
    loaded = load_config(DEFAULT, user_config=user)
    assert loaded["sv"]["harmonization"]["backend"] == "jasmine"
    write_effective_config(loaded, effective)
    env = os.environ.copy()
    env["XDG_CACHE_HOME"] = str(tmp_path / ".cache")
    if sys.platform == "win32":
        env["LOCALAPPDATA"] = str(tmp_path / ".local")
        env["APPDATA"] = str(tmp_path / ".app")
    result = subprocess.run(
        [
            str(snakemake()), "--snakefile", str(SNAKEFILE),
            "--configfile", str(effective), "--cores", "8", "--dry-run",
            "--printshellcmds", "--shared-fs-usage", "input-output",
            "persistence", "software-deployment", "software-deployment-cache",
            "sources", "storage-local-copies",
        ],
        cwd=tmp_path, env=env, capture_output=True, text=True,
        encoding="utf-8", errors="replace", check=False, timeout=60,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "read_based_sv" in output and "harmonize_sv" in output
    assert "S1.sawfish.sv.vcf.gz" in output
    assert "S1.harmonized.sv.vcf.gz" in output
    assert "pav_assembly_sv" not in output and "svim_asm_call" not in output


def test_phase9_rule_delegates_to_wrappers():
    bridge = (ROOT / "workflow" / "scripts" / "run_harmonization.py").read_text()
    assert "JasmineWrapper" in bridge and "TruvariWrapper" in bridge
    assert "subprocess" not in bridge and "os.system" not in bridge
