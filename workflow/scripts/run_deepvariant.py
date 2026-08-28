"""Snakemake bridge delegating one rule exclusively to DeepVariantWrapper."""

from pathlib import Path

from hifivar.alignment import AlignmentOutputFormat
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
    intermediate_directory=Path(str(snakemake.params.intermediate)),  # type: ignore[name-defined]
    logging_directory=Path(str(Path(str(snakemake.log[0])).parent)),  # type: ignore[name-defined]
    temporary_directory=Path(str(snakemake.params.tmpdir)),  # type: ignore[name-defined]
)
wrapper = DeepVariantWrapper(
    runtime=DeepVariantRuntime.from_config(snakemake.config),  # type: ignore[name-defined]
)
wrapper.run(request, stderr_path=Path(str(snakemake.log[0])))  # type: ignore[name-defined]
