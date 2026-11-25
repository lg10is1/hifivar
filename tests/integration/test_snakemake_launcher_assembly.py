"""Dry-run acceptance of the launcher calling-BAM / assembly-FASTQ boundary."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from hifivar.config import load_config
from hifivar.launcher import prepare_launch


ROOT = Path(__file__).parents[2]
DEFAULT = ROOT / "src/hifivar/resources/configs/default.yaml"


@pytest.mark.parametrize("calling", [False, True])
def test_launcher_dual_context_assembly_dag_dry_run(tmp_path, calling):
    pytest.importorskip("snakemake")
    (tmp_path / "ref.fa").write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    (tmp_path / "reads.fastq").write_text("@m/1/ccs\nACGT\n+\nIIII\n", encoding="utf-8")
    (tmp_path / "S1.bam").write_bytes(b"dryrun-metadata-only-not-a-real-bam")
    (tmp_path / "S1.bam.bai").write_bytes(b"dryrun-index")
    header = "sample_id\tinput\tsex\n"
    (tmp_path / "primary.tsv").write_text(header + f"S1\t{'S1.bam' if calling else 'reads.fastq'}\tfemale\n", encoding="utf-8")
    (tmp_path / "assembly.tsv").write_text(header + "S1\treads.fastq\tfemale\n", encoding="utf-8")
    (tmp_path / "PAV.Snakefile").write_text("# native deployment dependency\n", encoding="utf-8")
    settings = {
        "reference": {"fasta": "ref.fa", "build": "custom"}, "samples": {"sheet": "primary.tsv"},
        "runtime": {"threads": 4}, "workflow": {"memory_mb": 4000},
        "small": {"enabled": calling, "threads": 2, "memory_mb": 1000},
        "assembly": {"enabled": True, "sample_sheet": "assembly.tsv" if calling else None, "threads": 2, "memory_mb": 1000},
        "assembly_sv": {"enabled": True,
            "pav": {"snakefile": "PAV.Snakefile", "version": "2.4.6", "threads": 2, "memory_mb": 1000},
            "svim_asm": {"threads": 2, "memory_mb": 1000}},
    }
    user = tmp_path / "user.yaml"
    user.write_text(yaml.safe_dump(settings), encoding="utf-8")
    run = tmp_path / "run"
    effective, _, sheet, _ = prepare_launch(load_config(DEFAULT, user_config=user), run)
    run.mkdir()
    (run / "effective.yaml").write_text(yaml.safe_dump(effective), encoding="utf-8")
    (run / "samples.tsv").write_text(sheet, encoding="utf-8")
    environment = os.environ.copy()
    environment.update(PYTHONPATH=str(ROOT / "src"), TMP=str(tmp_path), TEMP=str(tmp_path), TMPDIR=str(tmp_path),
                       XDG_CACHE_HOME=str(tmp_path / "cache"), XDG_CONFIG_HOME=str(tmp_path / "config"),
                       LOCALAPPDATA=str(tmp_path / "local"), APPDATA=str(tmp_path / "app"),
                       PYTHONDONTWRITEBYTECODE="1")
    result = subprocess.run(
        [sys.executable, "-B", "-m", "snakemake", "--snakefile", str(ROOT / "workflow/Snakefile"),
         "--configfile", str(run / "effective.yaml"), "--cores", "4", "--dry-run", "--printshellcmds",
         "--shared-fs-usage", "input-output", "persistence", "software-deployment", "software-deployment-cache",
         "sources", "storage-local-copies", "--", "all"],
        cwd=run, env=environment, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=60,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    for name in ("hifiasm_assemble", "pav_assembly_sv", "svim_asm_call"):
        assert name in output
    assert "S1.pav.assembly.sv.vcf.gz" in output
    assert "S1.svim_asm.assembly.sv.vcf.gz" in output
    assert ("deepvariant_small" in output) is calling
    assert "dipcall" not in output.casefold()
    assert not (run / "results").exists()
    assert not (run / "work/hifiasm").exists()
