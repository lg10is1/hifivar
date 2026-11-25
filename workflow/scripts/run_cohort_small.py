"""Snakemake bridge for Phase 12 GLnexus joint genotyping."""

import json
from pathlib import Path

from hifivar.cohort import CohortDefinition, CohortTrack, read_cohort_input_manifest
from hifivar.attempt_receipts import attempt_step
from hifivar.context import AnalysisContext
from hifivar.glnexus import GLnexusRequest, GLnexusResources, GLnexusWrapper

config = snakemake.config  # type: ignore[name-defined]
section = config["cohort"]
track = section["small_variants"]
context = AnalysisContext.from_config(config)
cohort = CohortDefinition(str(section["cohort_id"]), context.sample_ids, context.reference)
manifest_path = Path(str(snakemake.input.manifest))  # type: ignore[name-defined]
inputs = read_cohort_input_manifest(manifest_path, cohort, CohortTrack.SMALL_VARIANTS)
request = GLnexusRequest(
    cohort, inputs,
    Path(str(config["paths"]["workdir"] or "work")) / "cohort" / cohort.cohort_id / "glnexus.DB",
    Path(str(snakemake.output.bcf)), Path(str(snakemake.output.vcf)),  # type: ignore[name-defined]
    str(track.get("preset", "DeepVariantWGS")),
    GLnexusResources(int(track.get("threads", 8)), int(track["memory_gb"])),
    bool(section.get("overwrite", False)),
)
receipt_outputs = {
    "bcf": request.output_bcf, "bcf_index": request.bcf_index,
    "vcf": request.output_vcf, "index": request.vcf_index,
    "provenance": Path(str(snakemake.output.result)),  # type: ignore[name-defined]
}
receipt_inputs = {item.sample_id: {"vcf": item.source_path, "index": item.index_path} for item in request.inputs}
with attempt_step("cohort.small_variants", tuple(context.sample_ids), outputs=receipt_outputs,
                  inputs=receipt_inputs, input_manifest=manifest_path) as receipt:
    result = GLnexusWrapper(executable=str(track.get("glnexus_executable", "glnexus_cli")), bcftools_executable=str(track.get("bcftools_executable", "bcftools"))).run(request, log_path=Path(str(snakemake.log[0])))  # type: ignore[name-defined]
    receipt_outputs["provenance"].write_text(json.dumps(result.as_track_result().to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    receipt.complete(
        outputs=receipt_outputs, tool_versions=result.versions,
        validation={"validated": True, "scope": "WRAPPER_INDEXED_GVCF_INPUT_HEADERS_NONEMPTY_BCF_CSI_VCF_TBI_AND_NAME_ASSOCIATED_OUTPUT_SAMPLE_QC"},
        metadata={"reference_build": context.reference.build, "reference_build_source": "EXPLICIT_REFERENCE_CONTEXT",
                  "declared_sample_order": list(context.sample_ids), "output_sample_order": result.qc["output_sample_order"],
                  "preset": request.preset, "memory_gb": request.resources.memory_gb,
                  "tool_version_sources": {tool: "COMMAND_PROBED" for tool in result.versions}},
    )
