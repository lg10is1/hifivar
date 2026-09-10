"""Materialize all nested raw FASTQ subsets in one streaming job."""

from hifivar.low_depth import FastqNestedDownsampler


FastqNestedDownsampler().run(
    snakemake.params.plan,  # type: ignore[name-defined]
    replicate=int(snakemake.params.replicate),  # type: ignore[name-defined]
    overwrite=bool(snakemake.params.overwrite),  # type: ignore[name-defined]
)
