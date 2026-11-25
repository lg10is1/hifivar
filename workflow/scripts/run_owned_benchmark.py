"""Owned launcher benchmark stage; standalone Phase13 rules stay separate."""
from hifivar.benchmark_handoff import execute
from hifivar.context import AnalysisContext

execute(snakemake.config, AnalysisContext.from_config(snakemake.config))  # type: ignore[name-defined]
