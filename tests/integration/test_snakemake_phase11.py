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


def snakemake() -> Path:
    name = "snakemake.exe" if sys.platform == "win32" else "snakemake"
    adjacent = Path(sys.executable).with_name(name)
    found = adjacent if adjacent.is_file() else shutil.which("snakemake")
    if found is None:
        pytest.skip("Snakemake is not installed.")
    return Path(found)


def test_phase11_annotation_is_optional_independent_dry_run(tmp_path):
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{reference}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    bam = tmp_path / "S1.bam"
    bam.write_bytes(b"BAM")
    Path(f"{bam}.bai").write_bytes(b"BAI")
    samples = tmp_path / "samples.tsv"
    samples.write_text(f"sample_id\tinput\tinput_type\nS1\t{bam}\tbam\n", encoding="utf-8")
    source = tmp_path / "S1.small.vcf"
    source.write_text(
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t2\tv1\tA\tT\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    manifest = tmp_path / "annotation_inputs.tsv"
    manifest.write_text(
        "sample\tvariant_category\tsource_vcf\tsource_tool\tsource_variant_ids\n"
        f"S1\tsmall\t{source}\tdeepvariant\tv1\n",
        encoding="utf-8",
    )
    user = tmp_path / "phase11.yaml"
    user.write_text(
        yaml.safe_dump(
            {
                "reference": {"fasta": str(reference), "build": "GRCh38"},
                "samples": {"sheet": str(samples)},
                "paths": {"workdir": str(tmp_path / "work"), "outdir": str(tmp_path / "results")},
                "annotation": {
                    "enabled": True,
                    "input_manifest": str(manifest),
                    "vep_enabled": True,
                    "vep_cache_directory": str(tmp_path / "vep-cache"),
                    "vep_cache_version": "115",
                    "vep_species": "homo_sapiens",
                    "vep_assembly": "GRCh38",
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    effective = tmp_path / "effective.yaml"
    write_effective_config(load_config(DEFAULT, user_config=user), effective)
    environment = os.environ.copy()
    environment["XDG_CACHE_HOME"] = str(tmp_path / ".cache")
    if sys.platform == "win32":
        environment["LOCALAPPDATA"] = str(tmp_path / ".local")
        environment["APPDATA"] = str(tmp_path / ".app")
    result = subprocess.run(
        [
            str(snakemake()), "--snakefile", str(SNAKEFILE), "--configfile", str(effective),
            "--cores", "1", "--dry-run", "--printshellcmds", "--shared-fs-usage",
            "input-output", "persistence", "software-deployment",
            "software-deployment-cache", "sources", "storage-local-copies",
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
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "annotate_variants" in output
    assert "results/annotation/S1/small" in output.replace("\\", "/")
    assert "deepvariant_small" not in output and "manual_review" not in output


def test_phase11_bridge_uses_dedicated_wrappers():
    bridge = (ROOT / "workflow" / "scripts" / "run_annotation.py").read_text(encoding="utf-8")
    assert "AnnovarWrapper" in bridge and "VepWrapper" in bridge and "run_phase11" in bridge
    assert "subprocess" not in bridge and "os.system" not in bridge
