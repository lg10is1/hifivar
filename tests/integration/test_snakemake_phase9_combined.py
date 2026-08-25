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


def test_phase9_combines_external_read_vcf_with_assembly_branches(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\nACGT\n")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n")
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r\nACGT\n+\n!!!!\n")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(f"sample_id\tinput\tinput_type\nS1\t{reads}\tfastq\n")
    external = tmp_path / "S1.sawfish.sv.vcf.gz"
    external.write_bytes(b"VCF")
    Path(f"{external}.tbi").write_bytes(b"TBI")
    pav_snake = tmp_path / "pav_site" / "Snakefile"
    pav_snake.parent.mkdir()
    pav_snake.write_text("# PAV\n")
    user = tmp_path / "phase9-combined.yaml"
    user.write_text(
        yaml.safe_dump(
            {
                "reference": {"fasta": str(fasta), "build": "GRCh38"},
                "samples": {"sheet": str(sheet)},
                "paths": {
                    "workdir": str(tmp_path / "work"),
                    "outdir": str(tmp_path / "results"),
                },
                "assembly": {"enabled": True},
                "assembly_sv": {
                    "enabled": True,
                    "pav": {"snakefile": str(pav_snake)},
                },
                "sv": {
                    "harmonization": {
                        "enabled": True,
                        "input_vcfs": {"S1": {"sawfish": str(external)}},
                    }
                },
            },
            sort_keys=False,
        )
    )
    effective = tmp_path / "effective.yaml"
    write_effective_config(load_config(DEFAULT, user_config=user), effective)
    env = os.environ.copy()
    env["XDG_CACHE_HOME"] = str(tmp_path / ".cache")
    if sys.platform == "win32":
        env["LOCALAPPDATA"] = str(tmp_path / ".local")
        env["APPDATA"] = str(tmp_path / ".app")
    result = subprocess.run(
        [
            str(snakemake()), "--snakefile", str(SNAKEFILE),
            "--configfile", str(effective), "--cores", "32", "--dry-run",
            "--printshellcmds", "--shared-fs-usage", "input-output",
            "persistence", "software-deployment", "software-deployment-cache",
            "sources", "storage-local-copies",
        ],
        cwd=tmp_path, env=env, capture_output=True, text=True,
        encoding="utf-8", errors="replace", check=False, timeout=60,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    for rule in ("hifiasm_assemble", "pav_assembly_sv", "svim_asm_call", "harmonize_sv"):
        assert rule in output
    assert str(external).replace("\\", "/") in output.replace("\\", "/")
    assert "S1.pav.assembly.sv.vcf.gz" in output
    assert "S1.svim_asm.assembly.sv.vcf.gz" in output
