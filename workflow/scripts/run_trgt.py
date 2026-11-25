import hashlib
from pathlib import Path

from hifivar.alignment import AlignmentOutputFormat
from hifivar.attempt_receipts import attempt_step
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
wrapper = TrgtWrapper(
    executable=str(tr_config["executable"]),
    bcftools_executable=str(tr_config["bcftools_executable"]),
    samtools_executable=str(tr_config["samtools_executable"]),
)
catalog_digest = hashlib.sha256()
with catalog.path.open("rb") as catalog_stream:
    for chunk in iter(lambda: catalog_stream.read(4 * 1024 * 1024), b""):
        catalog_digest.update(chunk)
receipt_outputs = {
    "vcf": request.final_vcf, "index": request.final_vcf_index,
    "spanning_bam": request.final_spanning_bam,
    "spanning_bam_index": request.final_spanning_bam_index,
}
with attempt_step("tr", (sample,), outputs=receipt_outputs) as receipt:
    result = wrapper.run(request, stderr_path=Path(str(snakemake.log[0])))  # type: ignore[name-defined]
    receipt.complete(
        outputs=receipt_outputs, tool_versions=result.tool_versions,
        validation={"validated": True, "scope": "WRAPPER_CATALOG_VCF_HEADER_TRGT_FIELDS_SAMPLE_CONTIGS_BGZF_TBI_AND_SPANNING_BAM_BAI_BASIC_VALIDATION"},
        metadata={"reference_build": reference.build, "reference_build_source": "EXPLICIT_REFERENCE_CONTEXT",
                  "catalog_id": "sha256:" + catalog_digest.hexdigest(),
                  "catalog_reference_build": catalog.reference_build,
                  "karyotype": karyotype,
                  "tool_version_sources": {tool: "COMMAND_PROBED" for tool in result.tool_versions}},
    )
