"""Snakemake bridge for Phase 12 TR cohort tables."""

from pathlib import Path

from hifivar import __version__
from hifivar.attempt_receipts import attempt_step
from hifivar.cohort import CohortDefinition, CohortTrack, read_cohort_input_manifest
from hifivar.cohort_tracks import build_tr_cohort_tables, write_track_result
from hifivar.context import AnalysisContext

config = snakemake.config  # type: ignore[name-defined]
section = config["cohort"]
context = AnalysisContext.from_config(config)
cohort = CohortDefinition(str(section["cohort_id"]), context.sample_ids, context.reference)
manifest_path = Path(str(snakemake.input.manifest))  # type: ignore[name-defined]
inputs = read_cohort_input_manifest(manifest_path, cohort, CohortTrack.TR)
scratch = Path(str(config["paths"]["workdir"] or "work")) / "cohort" / cohort.cohort_id / "tr.sqlite"
receipt_outputs = {
    "loci_tsv": Path(str(snakemake.output.loci)),  # type: ignore[name-defined]
    "sample_matrix_tsv": Path(str(snakemake.output.matrix)),  # type: ignore[name-defined]
    "provenance": Path(str(snakemake.output.result)),  # type: ignore[name-defined]
}
receipt_inputs = {item.sample_id: {"vcf": item.source_path, "index": item.index_path} if item.callable else {} for item in inputs}
with attempt_step("cohort.tr", tuple(context.sample_ids), outputs=receipt_outputs,
                  inputs=receipt_inputs, input_manifest=manifest_path) as receipt:
    result = build_tr_cohort_tables(cohort, inputs, locus_table=receipt_outputs["loci_tsv"], sample_matrix=receipt_outputs["sample_matrix_tsv"], scratch_database=scratch)
    write_track_result(receipt_outputs["provenance"], result)
    receipt.complete(
        outputs=receipt_outputs, tool_versions={"hifivar": __version__}, execution_kind="IN_PROCESS",
        validation={"validated": True, "scope": "SHARED_CATALOG_ID_INPUT_SAMPLE_CONTIG_AND_LOCUS_TABLE_BUILD_VALIDATION_NOT_TR_GENOTYPE_ACCURACY"},
        metadata={"reference_build": context.reference.build, "reference_build_source": "EXPLICIT_REFERENCE_CONTEXT",
                  "catalog_id": result.metrics["catalog_id"],
                  "catalog_id_source": "SHARED_EXPLICIT_INPUT_MANIFEST_NOT_A_NEW_CATALOG_FILE_HASH",
                  "input_states": {item.sample_id: item.state.value for item in inputs},
                  "tool_version_sources": {"hifivar": "RUNNING_PACKAGE_VERSION"}},
    )
