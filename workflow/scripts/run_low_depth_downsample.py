"""Materialize Phase 15 subsamples through SamtoolsDownsampleWrapper."""

from pathlib import Path

from hifivar.low_depth import SamtoolsDownsampleWrapper


wrapper = SamtoolsDownsampleWrapper(
    executable=str(snakemake.params.executable),  # type: ignore[name-defined]
)
source = Path(str(snakemake.input.source))  # type: ignore[name-defined]
wrapper.run(
    source,
    snakemake.params.run,  # type: ignore[name-defined]
    threads=int(snakemake.threads),  # type: ignore[name-defined]
    overwrite=bool(snakemake.params.overwrite),  # type: ignore[name-defined]
    stderr_path=Path(str(snakemake.log[0])),  # type: ignore[name-defined]
)
