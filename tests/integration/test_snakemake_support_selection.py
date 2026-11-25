"""DAG contracts for optional caller-support bundles; no real calling here."""
from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml

import hifivar
from hifivar.config import load_config, write_effective_config


ROOT = Path(hifivar.__file__).resolve().parents[2]
DEFAULT = ROOT / "src/hifivar/resources/configs/default.yaml"


def _snakemake():
    name = "snakemake.exe" if sys.platform == "win32" else "snakemake"
    adjacent = Path(sys.executable).with_name(name)
    executable = adjacent if adjacent.is_file() else shutil.which("snakemake")
    if executable is None:
        pytest.skip("Snakemake is not installed.")
    return str(executable)


@pytest.mark.parametrize("minimum,preexisting", [(1, False), (2, False), (2, True), (3, False)])
def test_support_bundle_dag_and_guards(tmp_path, minimum, preexisting):
    executable = _snakemake()
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(str(fasta) + ".fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    bam = tmp_path / "S1.bam"
    bam.write_bytes(b"SYNTHETIC_NOT_BAM")
    Path(str(bam) + ".bai").write_bytes(b"SYNTHETIC_NOT_INDEX")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(f"sample_id\tinput\tinput_type\nS1\t{bam}\tbam\n", encoding="utf-8")
    user = tmp_path / "config.yaml"
    user.write_text(yaml.safe_dump({
        "reference": {"fasta": str(fasta), "build": "GRCh38"}, "samples": {"sheet": str(sheet)},
        "paths": {"workdir": str(tmp_path / "work"), "outdir": str(tmp_path / "results")},
        "sv": {"enabled": True, "sawfish": {"enabled": True}, "sniffles2": {"enabled": True},
            "pbsv": {"enabled": False}, "cutesv": {"enabled": False},
            "harmonization": {"enabled": True, "minimum_caller_support": minimum, "overwrite": preexisting}}}),
        encoding="utf-8")
    effective = tmp_path / "effective.yaml"
    write_effective_config(load_config(DEFAULT, user_config=user), effective)
    protected = tmp_path / "results/sv_harmonized/S1/support-ge-2/protected.txt"
    if preexisting:
        protected.parent.mkdir(parents=True)
        protected.write_bytes(b"PRESERVE_BEFORE_SCHEDULING")
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(ROOT / "src") + (os.pathsep + environment["PYTHONPATH"] if environment.get("PYTHONPATH") else "")
    environment["XDG_CACHE_HOME"] = str(tmp_path / ".cache")
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    if sys.platform == "win32":
        environment["LOCALAPPDATA"] = str(tmp_path / ".local")
        environment["APPDATA"] = str(tmp_path / ".app")
    command = [executable, "--snakefile", str(ROOT / "workflow/Snakefile"), "--configfile", str(effective),
        "--cores", "8", "--dry-run", "--printshellcmds", "--shared-fs-usage", "input-output",
        "persistence", "software-deployment", "software-deployment-cache", "sources", "storage-local-copies"]
    result = subprocess.run(command, cwd=tmp_path, env=environment, capture_output=True, text=True,
        encoding="utf-8", errors="replace", check=False, timeout=60)
    output = result.stdout + result.stderr
    if preexisting:
        assert result.returncode != 0
        assert "fresh per-sample bundle before scheduling" in output
        assert protected.read_bytes() == b"PRESERVE_BEFORE_SCHEDULING"
    elif minimum > 2:
        assert result.returncode != 0
        assert "minimum_caller_support exceeds the actual input callers" in output
        assert not (tmp_path / "results").exists()
    else:
        assert result.returncode == 0, output
        assert "harmonize_sv" in output
        assert "S1.harmonized.sv.vcf.gz" in output
        assert "S1.sawfish.sv.vcf.gz" in output and "S1.sniffles2.sv.vcf.gz" in output
        assert ("support-ge-2" in output) is (minimum == 2)
        assert not (tmp_path / "results").exists()
    assert not (tmp_path / "work").exists()
