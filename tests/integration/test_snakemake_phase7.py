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
        pytest.skip("Snakemake is not installed.")
    return Path(found)


def test_phase7_hifiasm_dry_run_is_independent_and_deterministic(tmp_path: Path) -> None:
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t8\t6\t8\t9\n", encoding="utf-8")
    first = tmp_path / "reads1.fastq"
    second = tmp_path / "reads2.fastq"
    first.write_text("@one\nACGT\n+\n!!!!\n", encoding="utf-8")
    second.write_text("@two\nTGCA\n+\n!!!!\n", encoding="utf-8")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        f"sample_id\tinput\tinput_type\nS1\t{first};{second}\tfastq\n",
        encoding="utf-8",
    )
    user = tmp_path / "phase7.yaml"
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
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    effective = tmp_path / "effective.yaml"
    write_effective_config(load_config(DEFAULT_CONFIG, user_config=user), effective)
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
            "32",
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
    assert "hifiasm_assemble" in output
    for expected in (
        "S1.asm.bp.p_ctg.gfa",
        "S1.asm.bp.hap1.p_ctg.gfa",
        "S1.asm.bp.hap2.p_ctg.gfa",
        "S1.primary.fa",
        "S1.hap1.fa",
        "S1.hap2.fa",
    ):
        assert expected in output
    assert "pav" not in output.casefold()
    assert "svim-asm" not in output.casefold()
    assert "dipcall" not in output.casefold()
    assert not (tmp_path / "results" / "assembly").exists()


def test_phase7_rule_delegates_to_hifiasm_wrapper() -> None:
    rule = (ROOT / "workflow" / "rules" / "assembly.smk").read_text(
        encoding="utf-8"
    )
    bridge = (ROOT / "workflow" / "scripts" / "run_hifiasm.py").read_text(
        encoding="utf-8"
    )
    assert "hifiasm -o" not in rule
    assert "HifiasmWrapper" in bridge
    assert "subprocess" not in bridge and "os.system" not in bridge
    for forbidden in ("PAV", "SVIM", "dipcall"):
        assert forbidden not in bridge
