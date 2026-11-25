"""Real DAG planning from sealed ready-input manifests, without caller execution."""

from __future__ import annotations

import csv
import os
import struct
import subprocess
import sys
import zlib
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from hifivar import launcher
from hifivar.config import load_config


ROOT = Path(__file__).parents[2]
DEFAULT = ROOT / "src/hifivar/resources/configs/default.yaml"


def _bgzf(data):
    compressor = zlib.compressobj(wbits=-15)
    payload = compressor.compress(data) + compressor.flush()
    size = 18 + len(payload) + 8
    return (b"\x1f\x8b\x08\x04" + b"\0" * 4 + b"\0\xff\x06\0BC\x02\0" + struct.pack("<H", size - 1)
            + payload + struct.pack("<II", zlib.crc32(data), len(data)))


def _source(path, sample, category):
    marker = {"small_variants": ('##ALT=<ID=NON_REF,Description="Other">\n'
                                 '##INFO=<ID=END,Number=1,Type=Integer,Description="Reference-block end">\n'
                                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                                 '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Genotype likelihoods">\n'),
              "sv": '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n',
              "tr": '##INFO=<ID=TRID,Number=1,Type=String,Description="Locus">\n##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motifs">\n##INFO=<ID=STRUC,Number=1,Type=String,Description="Structure">\n'}[category]
    record = ("chr1\t100\t.\tA\t<NON_REF>\t.\t.\tEND=150\tGT:PL\t0/0:0,30,300\n"
              if category == "small_variants" else "chr1\t100\tv1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=150\tGT\t0/1\n")
    data = ("##fileformat=VCFv4.3\n##reference=GRCh38\n##contig=<ID=chr1,length=1000>\n" + marker
            + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n" + record)
    path.write_bytes(_bgzf(data.encode()))
    # Structural fixture only; this test plans the DAG and never tabix-queries it.
    Path(f"{path}.tbi").write_bytes(_bgzf(b"TBI\x01fixture"))


class _SchedulerBoundary:
    def run(self, command, **kwargs):
        out = "9.26.1\n" if "--version" in command else ""
        for key in ("stdout_path", "stderr_path"):
            if kwargs.get(key):
                Path(kwargs[key]).write_text(out if key == "stdout_path" else "", encoding="utf-8")
        return SimpleNamespace(stdout=out, returncode=0, duration_seconds=0.001)


@pytest.mark.parametrize("mode", ["annotation", "small_variants", "sv", "tr", "combined"])
def test_ready_downstream_launcher_dag(tmp_path, monkeypatch, mode):
    pytest.importorskip("snakemake")
    (tmp_path / "ref.fa").write_text(">chr1\n" + "A" * 1000 + "\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t1000\t6\t1000\t1001\n", encoding="utf-8")
    for sample in ("S2", "S1"):
        (tmp_path / f"{sample}.bam").write_bytes(b"not-a-real-bam-unused-downstream-context")
        for category in ("small_variants", "sv", "tr"):
            _source(tmp_path / f"{sample}.{category}.vcf.gz", sample, category)
    (tmp_path / "samples.tsv").write_text("sample_id\tinput\nS2\tS2.bam\nS1\tS1.bam\n", encoding="utf-8")
    (tmp_path / "regions.bed").write_text("chr1\t90\t160\tFEATURE\n", encoding="utf-8")
    (tmp_path / "annotation.tsv").write_text("sample\tvariant_category\tsource_vcf\tsource_tool\tsource_variant_ids\nS1\tsv\tS1.sv.vcf.gz\tcaller\tv1\n", encoding="utf-8")
    fields = ("sample", "track", "state", "source_path", "index_path", "source_tool", "source_version", "reference_build", "catalog_id")
    with (tmp_path / "cohort.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for category in ("small_variants", "sv", "tr"):
            for sample in ("S1", "S2"):
                filename = f"{sample}.{category}.vcf.gz"
                writer.writerow((sample, category, "CALLED", filename, filename + ".tbi", "declared-tool", "1.0", "GRCh38", "CAT1" if category == "tr" else ""))
    settings = {
        "reference": {"fasta": "ref.fa", "build": "GRCh38"}, "samples": {"sheet": "samples.tsv"},
        "runtime": {"threads": 4}, "workflow": {"memory_mb": 4096},
        "annotation": {"enabled": mode in {"annotation", "combined"}, "input_manifest": "annotation.tsv", "overlap_enabled": True,
                       "gene_bed": "regions.bed", "gene_version": "1", "threads": 1, "memory_mb": 1000},
        "cohort": {"enabled": mode != "annotation", "cohort_id": "C1", "input_manifest": "cohort.tsv",
                   "small_variants": {"enabled": mode in {"small_variants", "combined"}, "threads": 1, "memory_gb": 1},
                   "sv": {"enabled": mode in {"sv", "combined"}, "threads": 1, "memory_mb": 1000},
                   "tr": {"enabled": mode in {"tr", "combined"}, "threads": 1, "memory_mb": 1000}},
    }
    user = tmp_path / "user.yaml"
    user.write_text(yaml.safe_dump(settings), encoding="utf-8")
    run = tmp_path / "run"
    monkeypatch.setattr(launcher.platform, "system", lambda: "Linux")
    monkeypatch.setattr(launcher, "installed_workflow_root", lambda: ROOT / "workflow")
    assert launcher.WorkflowLauncher(_SchedulerBoundary()).run(load_config(DEFAULT, user_config=user), run_dir=run, dry_run=True) == "DRY_RUN_PASS"
    environment = os.environ.copy()
    environment.update(PYTHONPATH=str(ROOT / "src"), TMP=str(tmp_path), TEMP=str(tmp_path), TMPDIR=str(tmp_path),
                       XDG_CACHE_HOME=str(tmp_path / "cache"), XDG_CONFIG_HOME=str(tmp_path / "config"),
                       LOCALAPPDATA=str(tmp_path / "local"), APPDATA=str(tmp_path / "app"), PYTHONDONTWRITEBYTECODE="1")
    result = subprocess.run([sys.executable, "-B", "-m", "snakemake", "--snakefile", str(ROOT / "workflow/Snakefile"),
        "--configfile", str(run / "effective.yaml"), "--cores", "4", "--dry-run", "--printshellcmds",
        "--shared-fs-usage", "input-output", "persistence", "software-deployment", "software-deployment-cache", "sources", "storage-local-copies", "--", "all"],
        cwd=run, env=environment, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=60)
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert ("annotate_variants" in output) is (mode in {"annotation", "combined"})
    for category, rule in (("small_variants", "cohort_small_variants"), ("sv", "cohort_sv"), ("tr", "cohort_tr")):
        assert (rule in output) is (mode in {category, "combined"})
    assert "deepvariant_small" not in output and "pbmm2" not in output
    assert not (run / "results/annotation").exists() and not (run / "results/cohort").exists()
