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
    if found is None: pytest.skip("Snakemake is not installed.")
    return Path(found)


def test_phase12_tracks_are_optional_and_independent_dry_run(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"; reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{reference}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    samples = tmp_path / "samples.tsv"
    bam1, bam2 = tmp_path / "S1.bam", tmp_path / "S2.bam"
    bam1.write_bytes(b"BAM"); bam2.write_bytes(b"BAM")
    samples.write_text(f"sample_id\tinput\tinput_type\nS1\t{bam1}\tbam\nS2\t{bam2}\tbam\n", encoding="utf-8")
    source = tmp_path / "S1.sv.vcf.gz"; source.write_bytes(b"vcf"); Path(f"{source}.tbi").write_bytes(b"index")
    manifest = tmp_path / "cohort.tsv"
    manifest.write_text(
        "sample\ttrack\tstate\tsource_path\tindex_path\tsource_tool\tsource_version\treference_build\tcatalog_id\n"
        f"S1\tsv\tCALLED\t{source}\t{source}.tbi\tjasmine\t1.1.5\tGRCh38\t\n"
        "S2\tsv\tFAILED\t\t\tjasmine\t1.1.5\tGRCh38\t\n", encoding="utf-8")
    user = tmp_path / "phase12.yaml"
    user.write_text(yaml.safe_dump({
        "reference": {"fasta": str(reference), "build": "GRCh38"},
        "samples": {"sheet": str(samples)},
        "paths": {"workdir": str(tmp_path / "work"), "outdir": str(tmp_path / "results")},
        "cohort": {"enabled": True, "cohort_id": "C1", "input_manifest": str(manifest), "sv": {"enabled": True}},
    }, sort_keys=False), encoding="utf-8")
    effective = tmp_path / "effective.yaml"; write_effective_config(load_config(DEFAULT, user_config=user), effective)
    environment = os.environ.copy(); environment["XDG_CACHE_HOME"] = str(tmp_path / ".cache")
    if sys.platform == "win32": environment["LOCALAPPDATA"] = str(tmp_path / ".local"); environment["APPDATA"] = str(tmp_path / ".app")
    result = subprocess.run([str(snakemake()), "--snakefile", str(SNAKEFILE), "--configfile", str(effective), "--cores", "1", "--dry-run", "--printshellcmds", "--shared-fs-usage", "input-output", "persistence", "software-deployment", "software-deployment-cache", "sources", "storage-local-copies"], cwd=tmp_path, env=environment, capture_output=True, text=True, encoding="utf-8", errors="replace", check=False, timeout=60)
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "cohort_sv" in output and "cohort_manifest" in output
    assert "cohort_small_variants" not in output and "cohort_tr" not in output


def test_phase12_workflow_bridges_do_not_call_subprocess() -> None:
    for name in ("run_cohort_small.py", "run_cohort_sv.py", "run_cohort_tr.py"):
        text = (ROOT / "workflow" / "scripts" / name).read_text(encoding="utf-8")
        assert "subprocess" not in text and "os.system" not in text


def test_small_cohort_requires_and_exposes_explicit_memory(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{reference}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    samples = tmp_path / "samples.tsv"
    sample_rows = ["sample_id\tinput\tinput_type"]
    manifest_rows = [
        "sample\ttrack\tstate\tsource_path\tindex_path\tsource_tool\tsource_version\treference_build\tcatalog_id"
    ]
    for sample in ("S1", "S2"):
        bam = tmp_path / f"{sample}.bam"
        bam.write_bytes(b"BAM")
        gvcf = tmp_path / f"{sample}.g.vcf.gz"
        index = Path(f"{gvcf}.tbi")
        gvcf.write_bytes(b"gVCF")
        index.write_bytes(b"index")
        sample_rows.append(f"{sample}\t{bam}\tbam")
        manifest_rows.append(
            f"{sample}\tsmall_variants\tCALLED\t{gvcf}\t{index}\tdeepvariant\t1.10.0\tGRCh38\t"
        )
    samples.write_text("\n".join(sample_rows) + "\n", encoding="utf-8")
    manifest = tmp_path / "cohort.tsv"
    manifest.write_text("\n".join(manifest_rows) + "\n", encoding="utf-8")
    user = tmp_path / "phase12-small.yaml"
    user.write_text(
        yaml.safe_dump(
            {
                "reference": {"fasta": str(reference), "build": "GRCh38"},
                "samples": {"sheet": str(samples)},
                "paths": {
                    "workdir": str(tmp_path / "work"),
                    "outdir": str(tmp_path / "results"),
                },
                "cohort": {
                    "enabled": True,
                    "cohort_id": "C1",
                    "input_manifest": str(manifest),
                    "small_variants": {"enabled": True, "memory_gb": 192},
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
            str(snakemake()),
            "--snakefile",
            str(SNAKEFILE),
            "--configfile",
            str(effective),
            "--cores",
            "1",
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
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "cohort_small_variants" in output
    assert "mem_mb=196608" in output
