"""Write the immutable Phase 15 experiment plan before materialization."""

from pathlib import Path


snakemake.params.plan.write(  # type: ignore[name-defined]
    Path(str(snakemake.output.plan)),  # type: ignore[name-defined]
    overwrite=bool(snakemake.config["low_depth"].get("overwrite", False)),  # type: ignore[name-defined]
)
