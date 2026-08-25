"""Snakemake bridge for exact-catalog TR benchmarking."""
from pathlib import Path
from hifivar.benchmark import TruthSet, compare_tr_vcfs
from hifivar.context import AnalysisContext
from hifivar.serialization import write_json_atomic

section = snakemake.config["benchmark"]  # type: ignore[name-defined]
track = section["tr"]
context = AnalysisContext.from_config(snakemake.config)  # type: ignore[name-defined]
truth = TruthSet(Path(track["truth_vcf"]), track["truth_version"], context.reference.build, track["truth_source"], track["catalog_id"])
result = compare_tr_vcfs(benchmark_id=section["benchmark_id"], sample_id=section["sample_id"], query_vcf=Path(track["query_vcf"]), truth_set=truth, query_catalog_id=track["catalog_id"], output_tsv=Path(str(snakemake.output.loci)), overwrite=bool(section["overwrite"]))  # type: ignore[name-defined]
write_json_atomic(result.to_dict(), Path(str(snakemake.output.result)), overwrite=bool(section["overwrite"]), artifact_name="Benchmark result")  # type: ignore[name-defined]
