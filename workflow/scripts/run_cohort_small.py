"""Snakemake bridge for Phase 12 GLnexus joint genotyping."""

import json
from pathlib import Path

from hifivar.cohort import CohortDefinition, CohortTrack, read_cohort_input_manifest
from hifivar.context import AnalysisContext
from hifivar.glnexus import GLnexusRequest, GLnexusResources, GLnexusWrapper

config = snakemake.config  # type: ignore[name-defined]
section = config["cohort"]
track = section["small_variants"]
context = AnalysisContext.from_config(config)
cohort = CohortDefinition(str(section["cohort_id"]), context.sample_ids, context.reference)
inputs = read_cohort_input_manifest(Path(str(snakemake.input.manifest)), cohort, CohortTrack.SMALL_VARIANTS)  # type: ignore[name-defined]
request = GLnexusRequest(
    cohort, inputs,
    Path(str(config["paths"]["workdir"] or "work")) / "cohort" / cohort.cohort_id / "glnexus.DB",
    Path(str(snakemake.output.bcf)), Path(str(snakemake.output.vcf)),  # type: ignore[name-defined]
    str(track.get("preset", "DeepVariantWGS")),
    GLnexusResources(int(track.get("threads", 8)), int(track.get("memory_gb", 32))),
    bool(section.get("overwrite", False)),
)
result = GLnexusWrapper(executable=str(track.get("glnexus_executable", "glnexus_cli")), bcftools_executable=str(track.get("bcftools_executable", "bcftools"))).run(request, log_path=Path(str(snakemake.log[0])))  # type: ignore[name-defined]
Path(str(snakemake.output.result)).write_text(json.dumps(result.as_track_result().to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")  # type: ignore[name-defined]
