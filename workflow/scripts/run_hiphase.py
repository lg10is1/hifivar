"""Snakemake bridge delegating Phase 6 exclusively to HiPhaseWrapper."""

from pathlib import Path

from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import (
    AlignmentArtifact,
    AlignmentSortOrder,
    AlignmentSource,
)
from hifivar.hiphase import HiPhaseWrapper
from hifivar.phasing import PhasingRequest, PhasingResources
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
)
wrapper = HiPhaseWrapper(
    executable=str(section.get("executable", "hiphase")),
    tabix_executable=str(section.get("tabix_executable", "tabix")),
)
wrapper.run(request, stderr_path=Path(str(snakemake.log[0])))  # type: ignore[name-defined]
