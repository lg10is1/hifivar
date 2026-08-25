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
    if found is None: pytest.skip("Snakemake is not installed.")
    return Path(found)

def test_phase8_dry_run_has_two_independent_assembly_sv_branches(tmp_path):
    fasta = tmp_path / "ref.fa"; fasta.write_text(">chr1\nACGT\n")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n")
    reads = tmp_path / "reads.fastq"; reads.write_text("@r\nACGT\n+\n!!!!\n")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(f"sample_id\tinput\tinput_type\nS1\t{reads}\tfastq\n")
    pav_snake = tmp_path / "pav_site" / "Snakefile"
    pav_snake.parent.mkdir(); pav_snake.write_text("# PAV\n")
    user = tmp_path / "phase8.yaml"
    user.write_text(yaml.safe_dump({
        "reference": {"fasta": str(fasta), "build": "GRCh38"},
        "samples": {"sheet": str(sheet)},
        "paths": {"workdir": str(tmp_path / "work"), "outdir": str(tmp_path / "results")},
        "assembly": {"enabled": True},
        "assembly_sv": {"enabled": True, "pav": {"snakefile": str(pav_snake)}},
    }, sort_keys=False))
    effective = tmp_path / "effective.yaml"
    write_effective_config(load_config(DEFAULT, user_config=user), effective)
    env = os.environ.copy(); env["XDG_CACHE_HOME"] = str(tmp_path / ".cache")
    if sys.platform == "win32":
        env["LOCALAPPDATA"] = str(tmp_path / ".local"); env["APPDATA"] = str(tmp_path / ".app")
    result = subprocess.run([
        str(snakemake()), "--snakefile", str(SNAKEFILE), "--configfile", str(effective),
        "--cores", "32", "--dry-run", "--printshellcmds", "--shared-fs-usage",
        "input-output", "persistence", "software-deployment", "software-deployment-cache",
        "sources", "storage-local-copies",
    ], cwd=tmp_path, env=env, capture_output=True, text=True, encoding="utf-8",
       errors="replace", check=False, timeout=60)
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    for rule in ("hifiasm_assemble", "pav_assembly_sv", "svim_asm_call"):
        assert rule in output
    assert "S1.pav.assembly.sv.vcf.gz" in output
    assert "S1.svim_asm.assembly.sv.vcf.gz" in output
    assert "jasmine" not in output.casefold() and "truvari" not in output.casefold()

def test_phase8_rules_delegate_to_dedicated_wrappers():
    pav = (ROOT / "workflow" / "scripts" / "run_pav.py").read_text()
    svim = (ROOT / "workflow" / "scripts" / "run_svim_asm.py").read_text()
    assert "PavWrapper" in pav and "SvimAsmWrapper" in svim
    assert "subprocess" not in pav + svim and "os.system" not in pav + svim
