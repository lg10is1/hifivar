from pathlib import Path

from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import AlignmentArtifact, AlignmentSortOrder, AlignmentSource
from hifivar.context import AnalysisContext
from hifivar.phase5 import resolve_trgt_karyotype
from hifivar.reference import ReferenceGenome
from hifivar.tr import TandemRepeatCatalog
from hifivar.trgt import TrgtPreset, TrgtRequest, TrgtResources, TrgtWrapper


sample = str(snakemake.wildcards.sample)  # type: ignore[name-defined]
alignment_path = Path(str(snakemake.input.alignment))  # type: ignore[name-defined]
reference = ReferenceGenome.from_fasta(
    Path(str(snakemake.input.reference)),  # type: ignore[name-defined]
    build=snakemake.config["reference"].get("build"),  # type: ignore[name-defined]
)
artifact = AlignmentArtifact(
    sample,
    alignment_path,
    AlignmentOutputFormat.BAM,
    reference,
    AlignmentSource.EXISTING,
    AlignmentSortOrder.UNKNOWN,
    Path(str(snakemake.input.alignment_index)),  # type: ignore[name-defined]
)
tr_config = snakemake.config["tr"]  # type: ignore[name-defined]
context = AnalysisContext.from_config(snakemake.config)  # type: ignore[name-defined]
record = next(item for item in context.samples if item.sample.sample_id == sample)
karyotype = resolve_trgt_karyotype(str(tr_config["karyotype"]), record.sex, sample)
catalog = TandemRepeatCatalog(
    Path(str(snakemake.input.catalog)),  # type: ignore[name-defined]
    tr_config.get("catalog_reference_build"),
)
request = TrgtRequest(
    artifact=artifact,
    catalog=catalog,
    raw_output_prefix=Path(str(snakemake.params.raw_prefix)),  # type: ignore[name-defined]
    final_vcf=Path(str(snakemake.output.vcf)),  # type: ignore[name-defined]
    final_spanning_bam=Path(str(snakemake.output.spanning_bam)),  # type: ignore[name-defined]
    karyotype=karyotype,
    resources=TrgtResources(
        int(snakemake.threads),  # type: ignore[name-defined]
        int(snakemake.resources.mem_mb),  # type: ignore[name-defined]
        int(snakemake.resources.runtime_min),  # type: ignore[name-defined]
    ),
    preset=TrgtPreset(str(tr_config["preset"]).lower()),
    overwrite=bool(tr_config.get("overwrite", False)),
)
TrgtWrapper(
    executable=str(tr_config["executable"]),
    bcftools_executable=str(tr_config["bcftools_executable"]),
    samtools_executable=str(tr_config["samtools_executable"]),
).run(request, stderr_path=Path(str(snakemake.log[0])))  # type: ignore[name-defined]
