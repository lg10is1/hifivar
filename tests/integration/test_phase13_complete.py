"""Tiny Phase 13 integration without real WGS or external benchmark tools."""
from pathlib import Path
import json
from hifivar.benchmark import (BenchmarkManifest, BenchmarkResult, BenchmarkStatus,
    BenchmarkVariantClass, TruthSet, compare_tr_vcfs)
from hifivar.happy import parse_happy_summary


def _tr(path: Path, allele: str) -> None:
    path.write_text("##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        f"chr1\t10\tL1\tA\t<STR>\t.\tPASS\tTRID=L1;END=12\tGT:AL\t0/1:{allele}\n",encoding="utf-8")


def test_phase13_tiny_multitrack_manifest_keeps_queries_immutable(tmp_path: Path) -> None:
    summary=tmp_path/"happy.summary.csv"
    summary.write_text("Type,Filter,METRIC.Recall,METRIC.Precision\nSNP,PASS,1,1\nINDEL,PASS,0.5,1\n",encoding="utf-8")
    small_query=tmp_path/"S1.small.vcf.gz"; small_query.write_bytes(b"immutable-small")
    small_truth=TruthSet(tmp_path/"GIAB.vcf.gz","v4.2.1","GRCh38","GIAB")
    small=BenchmarkResult("B1","S1",BenchmarkVariantClass.SMALL_VARIANT,BenchmarkStatus.PASS,
        small_query,small_truth,"hap.py","0.3.15",(summary,),parse_happy_summary(summary))
    tr_query,tr_truth=tmp_path/"S1.tr.vcf",tmp_path/"truth.tr.vcf"; _tr(tr_query,"10,12"); _tr(tr_truth,"10,12")
    tr=compare_tr_vcfs(benchmark_id="B1",sample_id="S1",query_vcf=tr_query,
        truth_set=TruthSet(tr_truth,"truth-v1","GRCh38","curated","catalog-v1"),
        query_catalog_id="catalog-v1",output_tsv=tmp_path/"tr.tsv")
    outputs=BenchmarkManifest("B1","GRCh38",(small,tr),{"tracks":["small_variant","tr"]}).write(tmp_path/"results")
    payload=json.loads(outputs[0].read_text(encoding="utf-8"))
    assert [item["status"] for item in payload["results"]] == ["PASS","PASS"]
    assert small_query.read_bytes()==b"immutable-small" and tr_query.read_text(encoding="utf-8").startswith("##fileformat")
    assert all(item["scientific_semantics"]["benchmark_performance_is_not_pathogenicity"] for item in payload["results"])
