from __future__ import annotations

import os
from pathlib import Path
import subprocess

import pytest
import yaml

from hifivar.config import load_config, write_effective_config
from test_snakemake_phase6 import _snakemake


ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "src/hifivar/resources/configs/default.yaml"


@pytest.mark.parametrize("caller,include_tr", [("pbsv", True), ("pbsv", False), (None, True)])
def test_two_sample_joint_phasing_dag_keeps_category_outputs(tmp_path, caller, include_tr):
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
    Path(f"{reference}.fai").write_text("chr1\t8\t6\t8\t9\n", encoding="utf-8")
    sheet = tmp_path / "samples.tsv"
    rows = ["sample_id\tinput\tinput_type"]
    for sample in ("S1", "S2"):
        bam = tmp_path / f"{sample}.bam"
        bam.write_bytes(b"BAM")
        Path(f"{bam}.bai").write_bytes(b"BAI")
        rows.append(f"{sample}\t{bam}\tbam")
    sheet.write_text("\n".join(rows) + "\n", encoding="utf-8")
    catalog = tmp_path / "catalog.bed"
    catalog.write_text("chr1\t0\t4\tID=L1;MOTIFS=AC;STRUC=(AC)n\n", encoding="utf-8")
    user = tmp_path / "joint.yaml"
    user.write_text(yaml.safe_dump({"reference": {"fasta": str(reference), "build": "GRCh38"},
        "samples": {"sheet": str(sheet)}, "paths": {"workdir": str(tmp_path / "work"), "outdir": str(tmp_path / "results")},
        "small": {"enabled": True},
        "sv": {"enabled": caller is not None, "pbsv": {"enabled": caller == "pbsv"},
               "sawfish": {"enabled": False}, "sniffles2": {"enabled": False}, "cutesv": {"enabled": False}},
        "tr": {"enabled": include_tr, "catalog": str(catalog), "catalog_reference_build": "GRCh38", "karyotype": "XX"},
        "phasing": {"enabled": True, "sv_caller": caller, "include_tr": include_tr}}), encoding="utf-8")
    effective = tmp_path / "effective.yaml"
    write_effective_config(load_config(DEFAULT, user_config=user), effective)
    environment = dict(os.environ, XDG_CACHE_HOME=str(tmp_path / "cache"),
                       LOCALAPPDATA=str(tmp_path / "local"), APPDATA=str(tmp_path / "app"))
    result = subprocess.run([str(_snakemake()), "--snakefile", str(ROOT / "workflow/Snakefile"),
        "--configfile", str(effective), "--cores", "16", "--dry-run", "--printshellcmds",
        "--shared-fs-usage", "input-output", "persistence", "software-deployment", "software-deployment-cache",
        "sources", "storage-local-copies"], cwd=tmp_path, env=environment, capture_output=True,
        text=True, encoding="utf-8", errors="replace", timeout=60, check=False)
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "deepvariant_small" in output and "hiphase_phasing" in output
    assert ("read_based_sv" in output) == bool(caller)
    assert ("tandem_repeat" in output) == include_tr
    for sample in ("S1", "S2"):
        assert f"{sample}.phased.vcf.gz" in output
        assert (f"{sample}.sv.phased.vcf.gz" in output) == bool(caller)
        assert (f"{sample}.tr.phased.vcf.gz" in output) == include_tr
        assert f"{sample}.sniffles2.sv.vcf.gz" not in output
    assert not (tmp_path / "results/phasing").exists()
