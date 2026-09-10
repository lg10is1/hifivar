"""Benchmark one Phase 15 DeepVariant call with explicit truth resources."""

from pathlib import Path

from hifivar.benchmark import BenchmarkRegion, TruthSet
from hifivar.happy import HappyRequest, HappyWrapper
from hifivar.reference import ReferenceGenome
from hifivar.serialization import write_json_atomic

section = snakemake.config["low_depth"]  # type: ignore[name-defined]
e2e = section["end_to_end"]
reference = ReferenceGenome.from_fasta(
    Path(str(snakemake.input.reference)),  # type: ignore[name-defined]
    build=snakemake.config["reference"]["build"],  # type: ignore[name-defined]
)
truth = TruthSet(
    Path(str(snakemake.input.truth)), str(e2e["truth_version"]),
    str(reference.build), str(e2e["truth_source"]),
)
region = BenchmarkRegion(
    "confident", Path(str(snakemake.input.confident)),  # type: ignore[name-defined]
    str(e2e["confident_bed_version"]),
)
request = HappyRequest(
    benchmark_id=str(section["study_id"]), sample_id=str(section["sample_id"]),
    reference=reference, query_vcf=Path(str(snakemake.input.vcf)),  # type: ignore[name-defined]
    truth_set=truth, confident_regions=region,
    output_prefix=Path(str(snakemake.params.prefix)),  # type: ignore[name-defined]
    threads=int(snakemake.threads),  # type: ignore[name-defined]
    overwrite=bool(section.get("overwrite", False)),
)
result = HappyWrapper(
    executable=str(e2e.get("happy_executable", "hap.py")),
    configured_version=e2e.get("happy_version"),
).run(request, stderr_path=Path(str(snakemake.log[0])))  # type: ignore[name-defined]
if result.metrics_artifact != Path(str(snakemake.output.metrics)):  # type: ignore[name-defined]
    raise RuntimeError("hap.py metrics artifact did not match the Phase 15 contract.")
write_json_atomic(
    result.to_dict(), Path(str(snakemake.output.result)),  # type: ignore[name-defined]
    overwrite=bool(section.get("overwrite", False)),
    artifact_name="Low-depth hap.py result",
)
