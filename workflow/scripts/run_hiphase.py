"""Snakemake bridge delegating Phase 6 exclusively to HiPhaseWrapper."""

from pathlib import Path

from hifivar.alignment import AlignmentOutputFormat
from hifivar.attempt_receipts import artifact_identity, attempt_step
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.alignment_postprocess import (
    AlignmentArtifact,
    AlignmentSortOrder,
    AlignmentSource,
)
from hifivar.hiphase import HiPhaseWrapper
from hifivar.phasing import PhasingRequest, PhasingResources
from hifivar.joint_phasing import PhasingTrackInput, phasing_options
from hifivar.reference import ReferenceGenome
from hifivar.small import SmallVariantArtifact


reference = ReferenceGenome.from_fasta(
    Path(str(snakemake.input.reference)),  # type: ignore[name-defined]
    build=snakemake.config["reference"].get("build"),  # type: ignore[name-defined]
)
alignment = AlignmentArtifact(
    sample_id=str(snakemake.wildcards.sample),  # type: ignore[name-defined]
    path=Path(str(snakemake.input.alignment)),  # type: ignore[name-defined]
    output_format=AlignmentOutputFormat.BAM,
    reference=reference,
    source=AlignmentSource.EXISTING,
    sort_order=AlignmentSortOrder.UNKNOWN,
    index_path=Path(str(snakemake.input.alignment_index)),  # type: ignore[name-defined]
)
small = SmallVariantArtifact(
    sample_id=str(snakemake.wildcards.sample),  # type: ignore[name-defined]
    reference_build=reference.build,
    vcf_path=Path(str(snakemake.input.vcf)),  # type: ignore[name-defined]
    gvcf_path=Path(str(snakemake.input.gvcf)),  # type: ignore[name-defined]
    vcf_index_path=Path(str(snakemake.input.vcf_index)),  # type: ignore[name-defined]
    gvcf_index_path=Path(str(snakemake.input.gvcf_index)),  # type: ignore[name-defined]
)
section = snakemake.config["phasing"]  # type: ignore[name-defined]
sv_caller, include_tr, global_mode = phasing_options(section)
selected = ([("sv", sv_caller)] if sv_caller else []) + ([("tr", "trgt")] if include_tr else [])
joint_inputs = list(getattr(snakemake.input, "joint_vcfs", []))  # type: ignore[name-defined]
joint_indexes = list(getattr(snakemake.input, "joint_indexes", []))  # type: ignore[name-defined]
if len(joint_inputs) != len(selected) or len(joint_indexes) != len(selected):
    raise InputValidationError("Joint phasing DAG sources do not match explicit category selection.")
additional = tuple(PhasingTrackInput(category, tool, str(snakemake.wildcards.sample),  # type: ignore[name-defined]
    reference.build, Path(str(source)), Path(str(index)))
    for (category, tool), source, index in zip(selected, joint_inputs, joint_indexes))
request = PhasingRequest(
    alignment=alignment,
    small_variants=small,
    output_vcf=Path(str(snakemake.output.vcf)),  # type: ignore[name-defined]
    resources=PhasingResources(
        threads=int(snakemake.threads),  # type: ignore[name-defined]
        memory_mb=int(snakemake.resources.mem_mb),  # type: ignore[name-defined]
        runtime_minutes=int(snakemake.resources.runtime_min),  # type: ignore[name-defined]
    ),
    overwrite=bool(section.get("overwrite", False)),
    additional_inputs=additional,
    global_realignment=global_mode,
)
wrapper = HiPhaseWrapper(
    executable=str(section.get("executable", "hiphase")),
    tabix_executable=str(section.get("tabix_executable", "tabix")),
)
receipt_outputs = {"vcf": request.output_vcf, "index": request.output_vcf_index}
if additional:
    actual_outputs = [Path(str(p)).resolve() for p in getattr(snakemake.output, "joint_vcfs", [])]  # type: ignore[name-defined]
    actual_indexes = [Path(str(p)).resolve() for p in getattr(snakemake.output, "joint_indexes", [])]  # type: ignore[name-defined]
    if actual_outputs != [p.resolve() for _, _, _, p, _, _ in request.vcf_pairs[1:]] or (
        actual_indexes != [p.resolve() for _, _, _, _, p, _ in request.vcf_pairs[1:]]
    ):
        raise InputValidationError("Joint phasing wrapper outputs differ from the declared DAG outputs.")
    for category, _, _, output, index, _ in request.vcf_pairs[1:]:
        receipt_outputs.update({category + "_vcf": output, category + "_index": index})
source_before = ({category: {"vcf": artifact_identity(source), "index": artifact_identity(index)}
    for category, source, index, _, _, _ in request.vcf_pairs} if additional else {})
with attempt_step("phasing", (request.sample_id,), outputs=receipt_outputs) as receipt:
    result = wrapper.run(request, stderr_path=Path(str(snakemake.log[0])))  # type: ignore[name-defined]
    if additional:
        source_after = {category: {"vcf": artifact_identity(source), "index": artifact_identity(index)}
            for category, source, index, _, _, _ in request.vcf_pairs}
        if source_after != source_before:
            raise OutputValidationError("Joint phasing source VCF/index changed during execution.")
    receipt.complete(
        outputs=receipt_outputs, tool_versions={"hiphase": result.hiphase_version},
        validation={"validated": True, "scope": "WRAPPER_PHASED_VCF_HEADER_PS_SAMPLE_CONTIGS_BGZF_AND_TBI_NOT_PHASE_ACCURACY",
                    "joint_records_identity_and_genotype_checked": bool(additional), "phase_accuracy": "NOT_ASSESSED"},
        metadata={"reference_build": reference.build, "reference_build_source": "EXPLICIT_REFERENCE_CONTEXT",
                  "joint_source_identities": source_before,
                  "global_realignment": request.use_global_realignment,
                  "category_order": [p[0] for p in request.vcf_pairs],
                  "tool_version_sources": {"hiphase": "COMMAND_PROBED"}},
    )
