"""Snakemake bridge for Phase 12 TR cohort tables."""

from pathlib import Path

from hifivar.cohort import CohortDefinition, CohortTrack, read_cohort_input_manifest
from hifivar.cohort_tracks import build_tr_cohort_tables, write_track_result
from hifivar.context import AnalysisContext

config = snakemake.config  # type: ignore[name-defined]
section = config["cohort"]
context = AnalysisContext.from_config(config)
cohort = CohortDefinition(str(section["cohort_id"]), context.sample_ids, context.reference)
inputs = read_cohort_input_manifest(Path(str(snakemake.input.manifest)), cohort, CohortTrack.TR)  # type: ignore[name-defined]
scratch = Path(str(config["paths"]["workdir"] or "work")) / "cohort" / cohort.cohort_id / "tr.sqlite"
result = build_tr_cohort_tables(cohort, inputs, locus_table=Path(str(snakemake.output.loci)), sample_matrix=Path(str(snakemake.output.matrix)), scratch_database=scratch)  # type: ignore[name-defined]
write_track_result(Path(str(snakemake.output.result)), result)  # type: ignore[name-defined]
