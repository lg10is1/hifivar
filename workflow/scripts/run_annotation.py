"""Snakemake bridge for one independent Phase 11 annotation input."""

from pathlib import Path

from hifivar.annotation import (
    RegionCategory,
    RegionDatabase,
    annotate_region_overlaps,
    read_annotation_inputs,
    read_selected_variant_loci,
)
from hifivar.annovar import AnnovarRequest, AnnovarWrapper
from hifivar.context import AnalysisContext
from hifivar.phase11 import AnnotationJob, run_phase11
from hifivar.vep import VepRequest, VepWrapper


config = snakemake.config  # type: ignore[name-defined]
section = config["annotation"]
context = AnalysisContext.from_config(config)
all_inputs = read_annotation_inputs(
    Path(str(snakemake.input.manifest)),  # type: ignore[name-defined]
    reference=context.reference,
)
sample = str(snakemake.wildcards.sample)  # type: ignore[name-defined]
category = str(snakemake.wildcards.category)  # type: ignore[name-defined]
selected = next(
    item for item in all_inputs
    if item.sample_id == sample and item.variant_category.value == category
)
output = Path(str(snakemake.output.result))  # type: ignore[name-defined]
overwrite = bool(section.get("overwrite", False))

annovar_request = None
if section.get("annovar_enabled") is True:
    annovar_request = AnnovarRequest(
        selected,
        Path(str(section["annovar_database_root"])),
        str(section["annovar_database_version"]),
        tuple(str(item) for item in section["annovar_protocols"]),
        tuple(str(item) for item in section["annovar_operations"]),
        output / f"{sample}.{category}.annovar",
        str(section["annovar_version"]),
        overwrite,
    )

vep_request = None
if section.get("vep_enabled") is True:
    vep_request = VepRequest(
        selected,
        Path(str(section["vep_cache_directory"])),
        str(section["vep_cache_version"]),
        str(section["vep_species"]),
        str(section["vep_assembly"]),
        output / f"{sample}.{category}.vep.tsv",
        int(section.get("threads", 4)),
        overwrite,
    )

overlap_results = []
if section.get("overlap_enabled") is True and category in {"sv", "tr"}:
    region_databases = []
    for key, region_category in (
        ("gene", RegionCategory.GENE), ("exon", RegionCategory.EXON),
        ("regulatory", RegionCategory.REGULATORY),
        ("repeat", RegionCategory.REPEAT), ("segdup", RegionCategory.SEGDUP),
    ):
        database_path = section.get(f"{key}_bed")
        if database_path:
            region_databases.append(
                RegionDatabase(
                    region_category, Path(str(database_path)),
                    str(section[f"{key}_version"]),
                    context.reference.build or "unknown",
                )
            )
    overlap_results.append(
        annotate_region_overlaps(
            read_selected_variant_loci(selected),
            region_databases,
            reference=context.reference,
            output_tsv=output / f"{sample}.{category}.region-overlap.tsv",
            overwrite=overwrite,
        )
    )

jobs = ()
if annovar_request is not None or vep_request is not None:
    jobs = (AnnotationJob(selected, annovar_request, vep_request),)
report = run_phase11(
    jobs,
    annovar_wrapper=AnnovarWrapper(executable=str(section["annovar_executable"])),
    vep_wrapper=VepWrapper(executable=str(section["vep_executable"])),
    region_overlaps=overlap_results,
    log_directory=output / "logs",
)
report.write_json(output / "phase11.provenance.json", overwrite=overwrite)
report.write_yaml(output / "phase11.provenance.yaml", overwrite=overwrite)
