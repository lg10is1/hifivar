"""Snakemake bridge for the Phase 13 hap.py track."""
from pathlib import Path
from hifivar.benchmark import BenchmarkRegion, TruthSet
from hifivar.context import AnalysisContext
from hifivar.happy import HappyRequest, HappyWrapper
from hifivar.serialization import write_json_atomic

section = snakemake.config["benchmark"]  # type: ignore[name-defined]
track = section["small_variants"]
context = AnalysisContext.from_config(snakemake.config)  # type: ignore[name-defined]
truth = TruthSet(Path(track["truth_vcf"]), track["truth_version"], context.reference.build, track["truth_source"])
region = BenchmarkRegion("confident", Path(track["confident_bed"]), track["confident_bed_version"])
strata = tuple(BenchmarkRegion(item["name"], Path(item["path"]), item["version"], item["region_class"]) for item in track.get("stratifications", []))
request = HappyRequest(section["benchmark_id"], section["sample_id"], context.reference, Path(track["query_vcf"]), truth, region, Path(str(snakemake.output.summary)).with_suffix("").with_suffix(""), track["engine"], int(snakemake.threads), track["summary_filter"], strata, bool(section["overwrite"]))  # type: ignore[name-defined]
result = HappyWrapper(
    executable=track["happy_executable"], configured_version=track.get("happy_version")
).run(request, stderr_path=Path(str(snakemake.log[0])))  # type: ignore[name-defined]
if result.metrics_artifact != Path(str(snakemake.output.metrics)):  # type: ignore[name-defined]
    raise RuntimeError(
        f"hap.py metrics output did not match configured workflow contract: "
        f"expected '{snakemake.output.metrics}', discovered '{result.metrics_artifact}'."  # type: ignore[name-defined]
    )
payload = result.to_dict(); payload["execution_status"] = payload.pop("status"); payload["status"] = "PASS"; payload["variant_class"] = "small_variant"; payload["query_path"] = str(request.query_vcf); payload["tool"] = "hap.py"
write_json_atomic(payload, Path(str(snakemake.output.result)), overwrite=bool(section["overwrite"]), artifact_name="Benchmark result")  # type: ignore[name-defined]
