from __future__ import annotations

import json
from pathlib import Path
import pytest

from hifivar.benchmark import (BenchmarkManifest, BenchmarkResult, BenchmarkStatus,
    BenchmarkVariantClass, TruthSet, compare_tr_vcfs)
from hifivar.exceptions import InputValidationError, OutputValidationError


def truth(path: Path, *, catalog: str | None = None) -> TruthSet:
    return TruthSet(path, "v1.0", "GRCh38", "synthetic-test", catalog)


def test_truth_version_must_be_explicit_and_pinned(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="latest"):
        TruthSet(tmp_path / "truth.vcf.gz", "latest", "GRCh38", "GIAB")


def test_manifest_serialization_preserves_scientific_boundary(tmp_path: Path) -> None:
    query = tmp_path / "query.vcf.gz"
    result = BenchmarkResult("B1", "S1", BenchmarkVariantClass.SV, BenchmarkStatus.NOT_RUN,
        query, truth(tmp_path / "truth.vcf.gz"), "truvari", None,
        notes=("Missing truth resource is not zero performance.",))
    outputs = BenchmarkManifest("B1", "GRCh38", (result,), {"enabled": True}).write(tmp_path / "out")
    payload = json.loads(outputs[0].read_text(encoding="utf-8"))
    assert payload["results"][0]["status"] == "NOT_RUN"
    assert payload["results"][0]["scientific_semantics"]["benchmark_performance_is_not_pathogenicity"] is True
    assert payload["results"][0]["scientific_semantics"]["f1_is_not_clinical_utility"] is True
    assert payload["results"][0]["scientific_semantics"]["caller_count_is_not_benchmark_truth"] is True
    assert not query.exists()


def _tr_vcf(path: Path, allele: str, gt: str = "0/1") -> None:
    path.write_text("##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        f"chr1\t10\tL1\tA\t<STR10>\t.\tPASS\tTRID=L1;END=12\tGT:AL\t{gt}:{allele}\n", encoding="utf-8")


def test_tr_exact_catalog_comparison_and_no_call(tmp_path: Path) -> None:
    query, expected = tmp_path / "query.vcf", tmp_path / "truth.vcf"
    _tr_vcf(query, "10,12"); _tr_vcf(expected, "10,12")
    result = compare_tr_vcfs(benchmark_id="B1", sample_id="S1", query_vcf=query,
        truth_set=truth(expected, catalog="catalog-sha256"), query_catalog_id="catalog-sha256",
        output_tsv=tmp_path / "comparison.tsv")
    assert result.status is BenchmarkStatus.PASS
    assert {m.name: m.value for m in result.metrics}["exact_allele_agreement"] == 1.0
    assert query.read_text(encoding="utf-8").startswith("##fileformat")


def test_tr_catalog_mismatch_is_not_silently_compared(tmp_path: Path) -> None:
    query, expected = tmp_path / "query.vcf", tmp_path / "truth.vcf"
    _tr_vcf(query, "10,12"); _tr_vcf(expected, "10,12")
    with pytest.raises(InputValidationError, match="catalog"):
        compare_tr_vcfs(benchmark_id="B", sample_id="S1", query_vcf=query,
            truth_set=truth(expected, catalog="A"), query_catalog_id="B", output_tsv=tmp_path / "out.tsv")


def test_manifest_refuses_overwrite(tmp_path: Path) -> None:
    result = BenchmarkResult("B", "S", BenchmarkVariantClass.TR, BenchmarkStatus.UNSUPPORTED,
        tmp_path / "q", truth(tmp_path / "t"), "internal", "1")
    manifest = BenchmarkManifest("B", "GRCh38", (result,)); manifest.write(tmp_path / "out")
    with pytest.raises(OutputValidationError): manifest.write(tmp_path / "out")

def test_mixed_statuses_remain_distinct_and_missing_truth_has_no_zero_metric(tmp_path: Path) -> None:
    source=truth(tmp_path/"truth")
    results=tuple(BenchmarkResult("B", "S", kind, status, tmp_path/f"{kind.value}.vcf", source, "tool", None)
        for kind,status in ((BenchmarkVariantClass.SMALL_VARIANT,BenchmarkStatus.PASS),(BenchmarkVariantClass.SV,BenchmarkStatus.PARTIAL),(BenchmarkVariantClass.TR,BenchmarkStatus.NOT_RUN)))
    payload=BenchmarkManifest("B","GRCh38",results).to_dict()["results"]
    assert [item["status"] for item in payload] == ["PASS","PARTIAL","NOT_RUN"]
    assert payload[-1]["metrics"] == []
