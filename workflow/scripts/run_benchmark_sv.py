"""Snakemake bridge for read- or assembly-SV Truvari benchmarking."""
from pathlib import Path
from hifivar.benchmark import BenchmarkVariantClass, TruthSet
from hifivar.context import AnalysisContext
from hifivar.serialization import write_json_atomic
from hifivar.truvari import TruvariRequest, TruvariThresholds, TruvariWrapper, parse_truvari_summary, stratify_truvari_outputs

section = snakemake.config["benchmark"]  # type: ignore[name-defined]
track_name = getattr(snakemake.params, "track", "sv")  # type: ignore[name-defined]
track = section[track_name]
context = AnalysisContext.from_config(snakemake.config)  # type: ignore[name-defined]
thresholds = TruvariThresholds(**{key: track.get(key) for key in ("refdist", "pctseq", "pctsize", "pctovl", "sizemin", "sizemax", "bnddist")}, pass_only=track.get("pass_only", False))
request = TruvariRequest(section["sample_id"], context.reference, Path(track["truth_vcf"]), Path(track["query_vcf"]), Path(str(snakemake.output[0])), bool(section["overwrite"]), Path(track["confident_bed"]) if track.get("confident_bed") else None, thresholds)  # type: ignore[name-defined]
result = TruvariWrapper(executable=track["truvari_executable"]).run(request, stderr_path=Path(str(snakemake.log[0])))  # type: ignore[name-defined]
truth = TruthSet(Path(track["truth_vcf"]), track["truth_version"], context.reference.build, track["truth_source"])
variant_class = BenchmarkVariantClass.ASSEMBLY_SV if track_name == "assembly_sv" else BenchmarkVariantClass.SV
native_metrics = parse_truvari_summary(result.summary_path, variant_class=variant_class)
stratified, unsupported = stratify_truvari_outputs(result.request.output_directory, size_bins=tuple(track.get("size_bins", [])), variant_class=variant_class)
payload = result.to_dict(); payload["execution_status"] = payload.pop("status"); payload.update({"benchmark_id": section["benchmark_id"], "sample_id": section["sample_id"], "variant_class": track_name, "status": "PARTIAL" if unsupported else "PASS", "query_path": track["query_vcf"], "truth_set": truth.to_dict(), "tool": "truvari", "metrics": [m.to_dict() for m in native_metrics + stratified], "unsupported_length_stratification_types": list(unsupported), "region_class": track.get("region_class", "confident")})
write_json_atomic(payload, Path(str(snakemake.output.result)), overwrite=bool(section["overwrite"]), artifact_name="Benchmark result")  # type: ignore[name-defined]
