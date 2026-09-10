"""Run DeepVariant for one independently aligned Phase 15 subset."""

from pathlib import Path

from hifivar.alignment import AlignmentOutputFormat, AlignmentTool
from hifivar.alignment_postprocess import AlignmentArtifact, AlignmentSortOrder, AlignmentSource
from hifivar.deepvariant import DeepVariantRuntime, DeepVariantWrapper
from hifivar.reference import ReferenceGenome
from hifivar.small import DeepVariantRequest, SmallVariantResources

section = snakemake.config["low_depth"]  # type: ignore[name-defined]
reference = ReferenceGenome.from_fasta(
    Path(str(snakemake.input.reference)),  # type: ignore[name-defined]
    build=snakemake.config["reference"]["build"],  # type: ignore[name-defined]
)
artifact = AlignmentArtifact(
    sample_id=str(section["sample_id"]),
    path=Path(str(snakemake.input.bam)),  # type: ignore[name-defined]
    output_format=AlignmentOutputFormat.BAM,
    reference=reference,
    source=AlignmentSource.GENERATED,
    sort_order=AlignmentSortOrder.COORDINATE,
    index_path=Path(str(snakemake.input.index)),  # type: ignore[name-defined]
    tool=AlignmentTool.PBMM2,
)
request = DeepVariantRequest(
    artifact=artifact,
    output_vcf=Path(str(snakemake.output.vcf)),  # type: ignore[name-defined]
    output_gvcf=Path(str(snakemake.output.gvcf)),  # type: ignore[name-defined]
    resources=SmallVariantResources(
        int(snakemake.threads),  # type: ignore[name-defined]
        int(snakemake.resources.mem_mb),  # type: ignore[name-defined]
        int(snakemake.resources.runtime_min),  # type: ignore[name-defined]
    ),
    overwrite=bool(section.get("overwrite", False)),
    intermediate_directory=Path(str(snakemake.params.intermediate)),  # type: ignore[name-defined]
    logging_directory=Path(str(Path(str(snakemake.log[0])).parent)),  # type: ignore[name-defined]
    temporary_directory=Path(str(snakemake.params.tmpdir)),  # type: ignore[name-defined]
)
DeepVariantWrapper(runtime=DeepVariantRuntime.from_config(snakemake.config)).run(
    request, stderr_path=Path(str(snakemake.log[0])),  # type: ignore[name-defined]
)
