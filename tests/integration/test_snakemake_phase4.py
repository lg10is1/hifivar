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


def test_phase4_four_caller_dry_run_has_independent_outputs(tmp_path: Path) -> None:
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    bam = tmp_path / "S1.bam"
    bam.write_bytes(b"alignment")
    Path(f"{bam}.bai").write_bytes(b"index")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(f"sample_id\tinput\tinput_type\nS1\t{bam}\tbam\n", encoding="utf-8")
    user = tmp_path / "phase4.yaml"
    user.write_text(yaml.safe_dump({
        "reference": {"fasta": str(fasta), "build": "GRCh38"},
        "samples": {"sheet": str(sheet)},
        "paths": {"workdir": str(tmp_path / "work"), "outdir": str(tmp_path / "results")},
        "sv": {"enabled": True, "sawfish": {"threads": 7}},
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
        "--cores", "32", "--dry-run", "--printshellcmds", "--shared-fs-usage",
        "input-output", "persistence", "software-deployment", "software-deployment-cache",
        "sources", "storage-local-copies",
    ], cwd=tmp_path, env=environment, capture_output=True, text=True, encoding="utf-8", errors="replace", check=False, timeout=60)
    output = f"{result.stdout}\n{result.stderr}"
    assert result.returncode == 0, output
    assert output.count("read_based_sv") >= 4
    for caller in ("sawfish", "sniffles2", "pbsv", "cutesv"):
        assert f"S1.{caller}.sv.vcf.gz" in output
    assert "S1.sv.vcf.gz" not in output
    assert not (tmp_path / "results" / "sv").exists()


def test_phase4_rule_delegates_all_commands_to_wrappers() -> None:
    rule = (ROOT / "workflow" / "rules" / "sv.smk").read_text(encoding="utf-8")
    bridge = (ROOT / "workflow" / "scripts" / "run_sv_caller.py").read_text(encoding="utf-8")
    assert "pbsv discover" not in rule and "--max_cluster_bias_INS" not in rule
    for wrapper in ("SawfishWrapper", "Sniffles2Wrapper", "PbsvWrapper", "CuteSvWrapper"):
        assert wrapper in bridge
    assert "subprocess" not in bridge and "os.system" not in bridge
