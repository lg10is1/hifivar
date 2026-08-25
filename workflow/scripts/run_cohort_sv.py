"""Snakemake bridge for Phase 12 native SV cohort tables."""

from pathlib import Path

from hifivar.cohort import CohortDefinition, CohortTrack, read_cohort_input_manifest
from hifivar.cohort_tracks import build_sv_cohort_tables, write_track_result
from hifivar.context import AnalysisContext

config = snakemake.config  # type: ignore[name-defined]
section = config["cohort"]
context = AnalysisContext.from_config(config)
cohort = CohortDefinition(str(section["cohort_id"]), context.sample_ids, context.reference)
inputs = read_cohort_input_manifest(Path(str(snakemake.input.manifest)), cohort, CohortTrack.SV)  # type: ignore[name-defined]
result = build_sv_cohort_tables(cohort, inputs, site_table=Path(str(snakemake.output.sites)), sample_matrix=Path(str(snakemake.output.matrix)))  # type: ignore[name-defined]
write_track_result(Path(str(snakemake.output.result)), result)  # type: ignore[name-defined]
