"""Snakemake bridge for one independent Phase 11 annotation input."""

from pathlib import Path

from hifivar import __version__
from hifivar.attempt_receipts import attempt_step
from hifivar.annotation import (
    RegionCategory,
    RegionDatabase,
    AnnotationInput,
    annotate_region_overlaps,
    read_annotation_inputs,
    read_selected_variant_loci,
)
from hifivar.annovar import AnnovarRequest, AnnovarWrapper
from hifivar.context import AnalysisContext
from hifivar.phase11 import AnnotationJob, run_phase11
from hifivar.vep import VepRequest, VepWrapper


def receipt_sources(item: AnnotationInput) -> dict:
    """Declare the actual selected source and its existing conventional index."""
    roles = {"vcf": item.source_vcf,
            "index": Path(str(item.source_vcf) + ".tbi") if str(item.source_vcf).lower().endswith(".vcf.gz") else None}
    if item.bnd_bundle_manifest:
        roles["bnd_bundle_manifest"] = item.bnd_bundle_manifest
    return {item.sample_id: roles}


class ReceiptAnnovarWrapper(AnnovarWrapper):
    """Link one validated native backend execution without changing Phase 11."""

    def run(self, request, **kwargs):
        outputs = {"tsv": request.output_tsv, "vcf": request.output_vcf}
        with attempt_step("annotation.annovar", (request.input.sample_id,),
                          outputs=outputs, unit_id=request.input.variant_category.value,
                          inputs=receipt_sources(request.input), input_manifest=manifest_path) as receipt:
            result = super().run(request, **kwargs)
            receipt.complete(
                outputs=outputs, tool_versions={"annovar": result.tool_version},
                validation={"validated": True, "scope": "WRAPPER_SOURCE_IDENTITY_NONEMPTY_MULTIANNO_TSV_AND_VCF_BASIC_STRUCTURE_NOT_ANNOTATION_ACCURACY"},
                metadata={"reference_build": request.input.reference.build,
                          "reference_build_source": "EXPLICIT_REFERENCE_CONTEXT",
                          "variant_category": request.input.variant_category.value,
                          "tool_annotation_scope": "WHOLE_SOURCE_VCF",
                          "source_variant_ids_filter_applied": False,
                          "database_version": request.database_version,
                          "database_version_source": "EXPLICIT_REQUEST_NOT_INDEPENDENT_DATABASE_VERIFICATION",
                          "tool_version_sources": {"annovar": "USER_DECLARED_NOT_NATIVE_PROBED"}},
            )
            return result


class ReceiptVepWrapper(VepWrapper):
    """Keep VEP native execution and its command-probed version distinct."""

    def run(self, request, **kwargs):
        outputs = {"tsv": request.output_tsv}
        with attempt_step("annotation.vep", (request.input.sample_id,),
                          outputs=outputs, unit_id=request.input.variant_category.value,
                          inputs=receipt_sources(request.input), input_manifest=manifest_path) as receipt:
            result = super().run(request, **kwargs)
            receipt.complete(
                outputs=outputs, tool_versions={"vep": result.tool_version},
                validation={"validated": True, "scope": "WRAPPER_SOURCE_IDENTITY_AND_NONEMPTY_VEP_TSV_NOT_ANNOTATION_ACCURACY"},
                metadata={"reference_build": request.input.reference.build,
                          "reference_build_source": "EXPLICIT_REFERENCE_CONTEXT",
                          "variant_category": request.input.variant_category.value,
                          "tool_annotation_scope": "WHOLE_SOURCE_VCF",
                          "source_variant_ids_filter_applied": False,
                          "cache_version": request.cache_version,
                          "cache_version_source": "EXPLICIT_REQUEST_NOT_INDEPENDENT_CACHE_VERIFICATION",
                          "tool_version_sources": {"vep": "COMMAND_PROBED"}},
            )
            return result


config = snakemake.config  # type: ignore[name-defined]
section = config["annotation"]
context = AnalysisContext.from_config(config)
sample = str(snakemake.wildcards.sample)  # type: ignore[name-defined]
category = str(snakemake.wildcards.category)  # type: ignore[name-defined]
manifest_path = Path(str(snakemake.input.manifest))  # type: ignore[name-defined]
selected_inputs = read_annotation_inputs(
    manifest_path,
    reference=context.reference,
    sample_id=sample,
    variant_category=category,
)
selected = selected_inputs[0]
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
    overlap_output = output / f"{sample}.{category}.region-overlap.tsv"
    with attempt_step("annotation.overlap", (sample,), outputs={"tsv": overlap_output},
                      unit_id=category, inputs=receipt_sources(selected), input_manifest=manifest_path) as receipt:
        overlap = annotate_region_overlaps(
            read_selected_variant_loci(selected),
            region_databases,
            reference=context.reference,
            output_tsv=overlap_output,
            overwrite=overwrite,
        )
        overlap_results.append(overlap)
        receipt.complete(
            outputs={"tsv": overlap_output}, tool_versions={"hifivar": __version__},
            execution_kind="IN_PROCESS",
            validation={"validated": True, "scope": "SELECTED_LOCI_REFERENCE_COMPATIBILITY_REGION_DATABASE_AND_OVERLAP_TABLE_GENERATION_NOT_FUNCTIONAL_TRUTH"},
            metadata={"reference_build": context.reference.build,
                      "reference_build_source": "EXPLICIT_REFERENCE_CONTEXT", "variant_category": category,
                      "database_versions_source": "EXPLICIT_REGION_DATABASE_DECLARATIONS",
                      "databases": [{"category": item.category.value, "version": item.version,
                                     "reference_build": item.reference_build} for item in region_databases],
                      "tool_version_sources": {"hifivar": "RUNNING_PACKAGE_VERSION"}},
        )

jobs = ()
if annovar_request is not None or vep_request is not None:
    jobs = (AnnotationJob(selected, annovar_request, vep_request),)
report = run_phase11(
    jobs,
    annovar_wrapper=ReceiptAnnovarWrapper(executable=str(section["annovar_executable"])),
    vep_wrapper=ReceiptVepWrapper(executable=str(section["vep_executable"])),
    region_overlaps=overlap_results,
    log_directory=output / "logs",
)
report.write_json(output / "phase11.provenance.json", overwrite=overwrite)
report.write_yaml(output / "phase11.provenance.yaml", overwrite=overwrite)
