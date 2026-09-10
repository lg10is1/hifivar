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


def _snakemake() -> Path:
    executable = "snakemake.exe" if sys.platform == "win32" else "snakemake"
    adjacent = Path(sys.executable).with_name(executable)
    found = adjacent if adjacent.is_file() else shutil.which("snakemake")
    if found is None:
        pytest.skip("Snakemake is not installed.")
    return Path(found)


def test_phase15_low_depth_branch_is_optional_and_deterministic(tmp_path: Path) -> None:
    source = tmp_path / "source.bam"
    source.write_bytes(b"tiny-placeholder")
    user = tmp_path / "config.yaml"
    user.write_text(
        yaml.safe_dump(
            {
                "paths": {"workdir": str(tmp_path / "work"), "outdir": str(tmp_path / "results")},
                "low_depth": {
                    "enabled": True,
                    "study_id": "P15",
                    "sample_id": "HG002",
                    "source_bam": str(source),
                    "baseline_coverage": 20,
                    "target_coverages": [5, 10],
                    "replicates": 2,
                    "profiles": [{"name": "baseline", "settings": {}}],
                    "thresholds": [{"variant_class": "SNP", "metric": "f1", "minimum": 0.98}],
                    "minimum_successful_replicates": 2,
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    effective = tmp_path / "effective.yaml"
    write_effective_config(load_config(DEFAULT, user_config=user), effective)
    env = os.environ.copy()
    env.update({"XDG_CACHE_HOME": str(tmp_path / "cache"), "LOCALAPPDATA": str(tmp_path / "local"), "APPDATA": str(tmp_path / "app")})
    result = subprocess.run(
        [
            str(_snakemake()), "--snakefile", str(SNAKEFILE), "--configfile", str(effective),
            "--cores", "1", "--dry-run", "--shared-fs-usage", "input-output", "persistence",
            "software-deployment", "software-deployment-cache", "sources", "storage-local-copies",
        ],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=60,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "low_depth_plan" in output
    assert "low_depth_downsample" in output
    assert "replicate-1.bam" in output
    assert "replicate-2.bam.bai" in output


def test_phase15_bridges_never_invoke_subprocess_directly() -> None:
    for name in (
        "write_low_depth_plan.py", "run_low_depth_downsample.py",
        "run_low_depth_raw_fastq.py",
        "run_low_depth_alignment.py", "run_low_depth_coverage.py",
        "run_low_depth_deepvariant.py", "run_low_depth_happy.py",
        "write_low_depth_observation.py",
    ):
        text = (ROOT / "workflow" / "scripts" / name).read_text(encoding="utf-8")
        assert "subprocess" not in text
        assert "os.system" not in text


def test_phase15_raw_fastq_branch_is_one_nested_materialization_job(tmp_path: Path) -> None:
    source = tmp_path / "reads.fastq.gz"
    source.write_bytes(b"placeholder")
    user = tmp_path / "raw.yaml"
    user.write_text(
        yaml.safe_dump(
            {
                "paths": {"workdir": str(tmp_path / "work"), "outdir": str(tmp_path / "results")},
                "low_depth": {
                    "enabled": True,
                    "study_id": "P15-raw",
                    "sample_id": "HG002",
                    "input_mode": "raw_fastq",
                    "source_fastqs": [str(source)],
                    "baseline_coverage": 20,
                    "target_coverages": [5, 10],
                    "replicates": 2,
                    "profiles": [{"name": "baseline", "settings": {}}],
                    "thresholds": [{"variant_class": "SNP", "metric": "f1", "minimum": 0.98}],
                    "minimum_successful_replicates": 2,
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    effective = tmp_path / "effective.yaml"
    write_effective_config(load_config(DEFAULT, user_config=user), effective)
    env = os.environ.copy()
    env.update({"XDG_CACHE_HOME": str(tmp_path / "cache"), "LOCALAPPDATA": str(tmp_path / "local"), "APPDATA": str(tmp_path / "app")})
    result = subprocess.run(
        [
            str(_snakemake()), "--snakefile", str(SNAKEFILE), "--configfile", str(effective),
            "--cores", "1", "--dry-run", "--shared-fs-usage", "input-output", "persistence",
            "software-deployment", "software-deployment-cache", "sources", "storage-local-copies",
        ],
        cwd=tmp_path, env=env, capture_output=True, text=True, encoding="utf-8",
        errors="replace", timeout=60,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "low_depth_raw_fastq_nested_subsamples" in output
    assert "replicate-1.fastq.gz" in output
    assert "replicate-2.fastq.gz" in output
    assert "low_depth_downsample" not in output


def test_phase15_end_to_end_branch_connects_wrappers_and_benchmark(tmp_path: Path) -> None:
    source = tmp_path / "reads.fastq.gz"
    reference = tmp_path / "GRCh38.fa"
    truth = tmp_path / "truth.vcf.gz"
    confident = tmp_path / "confident.bed"
    for path in (source, reference, Path(f"{reference}.fai"), truth, Path(f"{truth}.tbi"), confident):
        path.write_bytes(b"placeholder")
    user = tmp_path / "e2e.yaml"
    user.write_text(yaml.safe_dump({
        "reference": {"fasta": str(reference), "build": "GRCh38"},
        "paths": {"workdir": str(tmp_path / "work"), "outdir": str(tmp_path / "results")},
        "low_depth": {
            "enabled": True, "study_id": "P15-e2e", "sample_id": "HG002",
            "input_mode": "raw_fastq", "source_fastqs": [str(source)],
            "baseline_coverage": 20, "target_coverages": [5, 10],
            "replicates": 1, "profiles": [{"name": "baseline", "settings": {}}],
            "thresholds": [{"variant_class": "SNP", "metric": "f1", "minimum": 0.98}],
            "minimum_successful_replicates": 1,
            "end_to_end": {
                "enabled": True, "truth_vcf": str(truth),
                "truth_version": "GIAB-4.2.1", "truth_source": "GIAB",
                "confident_bed": str(confident),
                "confident_bed_version": "GIAB-4.2.1",
            },
        },
    }, sort_keys=False), encoding="utf-8")
    effective = tmp_path / "effective.yaml"
    write_effective_config(load_config(DEFAULT, user_config=user), effective)
    env = os.environ.copy()
    env.update({"XDG_CACHE_HOME": str(tmp_path / "cache"), "LOCALAPPDATA": str(tmp_path / "local"), "APPDATA": str(tmp_path / "app")})
    result = subprocess.run([
        str(_snakemake()), "--snakefile", str(SNAKEFILE), "--configfile", str(effective),
        "--cores", "1", "--dry-run", "--shared-fs-usage", "input-output", "persistence",
        "software-deployment", "software-deployment-cache", "sources", "storage-local-copies",
    ], cwd=tmp_path, env=env, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=60)
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    for rule in (
        "low_depth_align", "low_depth_measure_coverage",
        "low_depth_deepvariant", "low_depth_happy", "low_depth_observation",
    ):
        assert rule in output
    assert "replicate-1.small.vcf.gz" in output
    assert "replicate-1.coverage.json" in output
