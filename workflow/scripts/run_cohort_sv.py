"""Snakemake bridge for Phase 12 native SV cohort tables."""

from pathlib import Path

from hifivar import __version__
from hifivar.attempt_receipts import attempt_step
from hifivar.cohort import CohortDefinition, CohortTrack, read_cohort_input_manifest
from hifivar.cohort_tracks import build_sv_cohort_tables, write_track_result
from hifivar.context import AnalysisContext

config = snakemake.config  # type: ignore[name-defined]
section = config["cohort"]
context = AnalysisContext.from_config(config)
cohort = CohortDefinition(str(section["cohort_id"]), context.sample_ids, context.reference)
manifest_path = Path(str(snakemake.input.manifest))  # type: ignore[name-defined]
inputs = read_cohort_input_manifest(manifest_path, cohort, CohortTrack.SV)
receipt_outputs = {
    "sites_tsv": Path(str(snakemake.output.sites)),  # type: ignore[name-defined]
    "sample_matrix_tsv": Path(str(snakemake.output.matrix)),  # type: ignore[name-defined]
    "provenance": Path(str(snakemake.output.result)),  # type: ignore[name-defined]
}
receipt_inputs = {item.sample_id: {"vcf": item.source_path, "index": item.index_path} if item.callable else {} for item in inputs}
with attempt_step("cohort.sv", tuple(context.sample_ids), outputs=receipt_outputs,
                  inputs=receipt_inputs, input_manifest=manifest_path) as receipt:
    result = build_sv_cohort_tables(cohort, inputs, site_table=receipt_outputs["sites_tsv"], sample_matrix=receipt_outputs["sample_matrix_tsv"])
    write_track_result(receipt_outputs["provenance"], result)
    receipt.complete(
        outputs=receipt_outputs, tool_versions={"hifivar": __version__}, execution_kind="IN_PROCESS",
        validation={"validated": True, "scope": "NATIVE_REPRESENTATION_ROWS_INPUT_SAMPLE_CONTIG_AND_TABLE_BUILD_VALIDATION_NOT_CROSS_SAMPLE_CLUSTERING"},
        metadata={"reference_build": context.reference.build, "reference_build_source": "EXPLICIT_REFERENCE_CONTEXT",
                  "input_states": {item.sample_id: item.state.value for item in inputs},
                  "cross_sample_clustering": False, "tool_version_sources": {"hifivar": "RUNNING_PACKAGE_VERSION"}},
    )
