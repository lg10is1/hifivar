"""Snakemake bridge delegating one rule exclusively to DeepVariantWrapper."""

from pathlib import Path

from hifivar.alignment import AlignmentOutputFormat
from hifivar.attempt_receipts import attempt_step
from hifivar.alignment_postprocess import (
    AlignmentArtifact,
    AlignmentSortOrder,
    AlignmentSource,
)
from hifivar.deepvariant import DeepVariantRuntime, DeepVariantWrapper
from hifivar.reference import ReferenceGenome
from hifivar.small import DeepVariantRequest, SmallVariantResources


alignment_path = Path(str(snakemake.input.alignment))  # type: ignore[name-defined]
output_format = (
    AlignmentOutputFormat.BAM
    if alignment_path.suffix.lower() == ".bam"
    else AlignmentOutputFormat.CRAM
)
reference = ReferenceGenome.from_fasta(
    Path(str(snakemake.input.reference)),  # type: ignore[name-defined]
    build=snakemake.config["reference"].get("build"),  # type: ignore[name-defined]
)
artifact = AlignmentArtifact(
    sample_id=str(snakemake.wildcards.sample),  # type: ignore[name-defined]
    path=alignment_path,
    output_format=output_format,
    reference=reference,
    source=AlignmentSource.EXISTING,
    sort_order=AlignmentSortOrder.UNKNOWN,
    index_path=Path(str(snakemake.input.alignment_index)),  # type: ignore[name-defined]
)
small_config = snakemake.config["small"]  # type: ignore[name-defined]
request = DeepVariantRequest(
    artifact=artifact,
    output_vcf=Path(str(snakemake.output.vcf)),  # type: ignore[name-defined]
    output_gvcf=Path(str(snakemake.output.gvcf)),  # type: ignore[name-defined]
    resources=SmallVariantResources(
        threads=int(snakemake.threads),  # type: ignore[name-defined]
        memory_mb=int(snakemake.resources.mem_mb),  # type: ignore[name-defined]
        runtime_minutes=int(snakemake.resources.runtime_min),  # type: ignore[name-defined]
    ),
    overwrite=bool(small_config.get("overwrite", False)),
    region=snakemake.params.region,  # type: ignore[name-defined]
    intermediate_directory=Path(str(snakemake.params.intermediate)),  # type: ignore[name-defined]
    logging_directory=Path(str(Path(str(snakemake.log[0])).parent)),  # type: ignore[name-defined]
    temporary_directory=Path(str(snakemake.params.tmpdir)),  # type: ignore[name-defined]
)
wrapper = DeepVariantWrapper(
    runtime=DeepVariantRuntime.from_config(snakemake.config),  # type: ignore[name-defined]
)
receipt_outputs = {
    "vcf": request.output_vcf, "index": request.output_vcf_index,
    "gvcf": request.output_gvcf, "gvcf_index": request.output_gvcf_index,
}
with attempt_step("small", (request.sample_id,), outputs=receipt_outputs) as receipt:
    result = wrapper.run(request, stderr_path=Path(str(snakemake.log[0])))  # type: ignore[name-defined]
    receipt.complete(
        outputs=receipt_outputs,
        tool_versions={"deepvariant": result.tool_version},
        validation={"validated": True, "scope": "WRAPPER_BGZF_VCF_GVCF_HEADERS_EXACT_SAMPLE_CONTIGS_AND_TBI_MAGIC_NOT_ALL_RECORDS"},
        metadata={"reference_build": reference.build, "reference_build_source": "EXPLICIT_REFERENCE_CONTEXT",
                  "region": request.region, "tool_version_sources": {"deepvariant": "COMMAND_PROBED"}},
    )
