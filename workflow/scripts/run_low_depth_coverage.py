"""Measure post-alignment coverage for one Phase 15 replicate."""

from pathlib import Path

from hifivar.low_depth import SamtoolsCoverageWrapper
from hifivar.serialization import write_json_atomic

section = snakemake.config["low_depth"]  # type: ignore[name-defined]
result = SamtoolsCoverageWrapper(
    executable=str(section["end_to_end"].get("samtools_executable", "samtools")),
).run(
    sample_id=str(section["sample_id"]),
    alignment_bam=Path(str(snakemake.input.bam)),  # type: ignore[name-defined]
    alignment_index=Path(str(snakemake.input.index)),  # type: ignore[name-defined]
    stderr_path=Path(str(snakemake.log[0])),  # type: ignore[name-defined]
)
write_json_atomic(
    result.to_dict(), Path(str(snakemake.output.json)),  # type: ignore[name-defined]
    overwrite=bool(section.get("overwrite", False)),
    artifact_name="Low-depth observed coverage",
)
