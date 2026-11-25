"""Snakemake bridge for explicit review and its optional owned attempt receipt."""
from pathlib import Path

from hifivar.context import AnalysisContext
from hifivar.review import read_review_selection
from hifivar.review_handoff import run_review_step

context = AnalysisContext.from_config(snakemake.config)  # type: ignore[name-defined]
alignments = {record.sample.sample_id:record.sample.input.files[0] for record in context.samples}
section = snakemake.config['review']  # type: ignore[name-defined]
output_directory = Path(str(snakemake.output.manifest_json)).parent  # type: ignore[name-defined]
selection = Path(str(snakemake.input.selection))  # type: ignore[name-defined]
targets = read_review_selection(selection,alignments=alignments,
    reference_fasta=context.reference.fasta,output_directory=output_directory,
    default_flank_bp=int(section.get('flank_bp',500)),reference_build=context.reference.build)
run_review_step(snakemake.config,context,targets,output=output_directory,selection=selection)  # type: ignore[name-defined]
