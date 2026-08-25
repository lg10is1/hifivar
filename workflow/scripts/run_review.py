"""Snakemake bridge for the optional Phase 10 review branch."""

from pathlib import Path

from hifivar.context import AnalysisContext
from hifivar.igv import IgvWrapper
from hifivar.phase10 import run_phase10
from hifivar.review import read_review_selection


context = AnalysisContext.from_config(snakemake.config)  # type: ignore[name-defined]
alignments = {
    record.sample.sample_id: record.sample.input.files[0]
    for record in context.samples
}
section = snakemake.config["review"]  # type: ignore[name-defined]
output_directory = Path(str(snakemake.output.manifest_json)).parent  # type: ignore[name-defined]
targets = read_review_selection(
    Path(str(snakemake.input.selection)),  # type: ignore[name-defined]
    alignments=alignments,
    reference_fasta=context.reference.fasta,
    output_directory=output_directory,
    default_flank_bp=int(section.get("flank_bp", 500)),
)
run_phase10(
    targets,
    output_directory=output_directory,
    igv_wrapper=IgvWrapper(executable=str(section.get("igv_executable", "igv.sh"))),
    overwrite=bool(section.get("overwrite", False)),
)
