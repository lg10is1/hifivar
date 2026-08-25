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


def test_phase5_trgt_dry_run_has_independent_tr_outputs(tmp_path: Path) -> None:
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t8\t6\t8\t9\n", encoding="utf-8")
    bam = tmp_path / "S1.bam"
    bam.write_bytes(b"BAM")
    Path(f"{bam}.bai").write_bytes(b"BAI")
    catalog = tmp_path / "catalog.bed"
    catalog.write_text("chr1\t1\t7\tID=TR1;MOTIFS=CG;STRUC=<TR>\n", encoding="utf-8")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(f"sample_id\tinput\tinput_type\tsex\nS1\t{bam}\tbam\tfemale\n", encoding="utf-8")
    user = tmp_path / "phase5.yaml"
    user.write_text(yaml.safe_dump({
        "reference": {"fasta": str(fasta), "build": "GRCh38"},
        "samples": {"sheet": str(sheet)},
        "paths": {"workdir": str(tmp_path / "work"), "outdir": str(tmp_path / "results")},
        "tr": {"enabled": True, "catalog": str(catalog), "catalog_reference_build": "GRCh38"},
    }, sort_keys=False), encoding="utf-8")
    effective = tmp_path / "effective.yaml"
    write_effective_config(load_config(DEFAULT_CONFIG, user_config=user), effective)
    environment = os.environ.copy()
    environment["XDG_CACHE_HOME"] = str(tmp_path / ".cache")
    if sys.platform == "win32":
        environment["LOCALAPPDATA"] = str(tmp_path / ".local")
        environment["APPDATA"] = str(tmp_path / ".app")
    result = subprocess.run([
        str(_snakemake()), "--snakefile", str(SNAKEFILE), "--configfile", str(effective),
        "--cores", "8", "--dry-run", "--printshellcmds", "--shared-fs-usage",
        "input-output", "persistence", "software-deployment", "software-deployment-cache",
        "sources", "storage-local-copies",
    ], cwd=tmp_path, env=environment, capture_output=True, text=True, encoding="utf-8", errors="replace", check=False, timeout=60)
    output = f"{result.stdout}\n{result.stderr}"
    assert result.returncode == 0, output
    assert "tandem_repeat" in output
    for expected in ("S1.tr.vcf.gz", "S1.tr.vcf.gz.tbi", "S1.tr.spanning.bam", "S1.tr.spanning.bam.bai"):
        assert expected in output
    assert not (tmp_path / "results" / "tr").exists()


def test_phase5_rule_delegates_to_trgt_wrapper() -> None:
    rule = (ROOT / "workflow" / "rules" / "tr.smk").read_text(encoding="utf-8")
    bridge = (ROOT / "workflow" / "scripts" / "run_trgt.py").read_text(encoding="utf-8")
    assert "trgt genotype" not in rule and "bcftools sort" not in rule
    assert "TrgtWrapper" in bridge
    assert "subprocess" not in bridge and "os.system" not in bridge
