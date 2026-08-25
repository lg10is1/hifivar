"""Config-driven single-sample TRGT rules for Phase 5."""

from hifivar.alignment_postprocess import find_alignment_index
from hifivar.context import AnalysisContext
from hifivar.exceptions import WorkflowError
from hifivar.sample import InputType


TR_CONFIG = config.get("tr", {})
if not isinstance(TR_CONFIG, dict):
    raise WorkflowError("Effective config section tr must be a mapping.")
TR_ENABLED = TR_CONFIG.get("enabled", False)
if not isinstance(TR_ENABLED, bool):
    raise WorkflowError("tr.enabled must be boolean.")

TR_SAMPLE_IDS = []
TR_ALIGNMENTS = {}
TR_ALIGNMENT_INDEXES = {}
if TR_ENABLED:
    catalog = TR_CONFIG.get("catalog")
    if not isinstance(catalog, str) or not catalog.strip():
        raise WorkflowError("tr.enabled requires tr.catalog.")
    TR_CONTEXT = AnalysisContext.from_config(config)
    for record in TR_CONTEXT.samples:
        sample = record.sample
        if sample.input.input_type is not InputType.BAM:
            raise WorkflowError(
                f"TRGT Snakemake input for sample '{sample.sample_id}' must be aligned BAM."
            )
        alignment = sample.input.files[0]
        index = find_alignment_index(alignment)
        if index is None:
            raise WorkflowError(
                f"TRGT alignment index is missing for sample '{sample.sample_id}': '{alignment}'."
            )
        TR_SAMPLE_IDS.append(sample.sample_id)
        TR_ALIGNMENTS[sample.sample_id] = str(alignment)
        TR_ALIGNMENT_INDEXES[sample.sample_id] = str(index)


TR_VCF_PATTERN = (OUTPUT_ROOT / "tr" / "{sample}.tr.vcf.gz").as_posix()
TR_SPANNING_PATTERN = (OUTPUT_ROOT / "tr" / "{sample}.tr.spanning.bam").as_posix()
TR_VCFS = [TR_VCF_PATTERN.format(sample=sample) for sample in TR_SAMPLE_IDS]
TR_SPANNING_BAMS = [TR_SPANNING_PATTERN.format(sample=sample) for sample in TR_SAMPLE_IDS]
WORKFLOW_TARGETS += [
    *TR_VCFS,
    *(f"{path}.tbi" for path in TR_VCFS),
    *TR_SPANNING_BAMS,
    *(f"{path}.bai" for path in TR_SPANNING_BAMS),
]


rule tandem_repeat:
    input:
        reference=lambda wildcards: str(TR_CONTEXT.reference.fasta),
        reference_index=lambda wildcards: str(TR_CONTEXT.reference.fai),
        alignment=lambda wildcards: TR_ALIGNMENTS[wildcards.sample],
        alignment_index=lambda wildcards: TR_ALIGNMENT_INDEXES[wildcards.sample],
        catalog=lambda wildcards: str(TR_CONFIG["catalog"])
    output:
        vcf=TR_VCF_PATTERN,
        vcf_index=f"{TR_VCF_PATTERN}.tbi",
        spanning_bam=TR_SPANNING_PATTERN,
        spanning_bam_index=f"{TR_SPANNING_PATTERN}.bai"
    log:
        (LOG_ROOT / "tr" / "{sample}.trgt.log").as_posix()
    threads:
        int(TR_CONFIG.get("threads", 8))
    resources:
        mem_mb=int(TR_CONFIG.get("memory_mb", 16000)),
        runtime_min=int(TR_CONFIG.get("runtime_minutes", 720))
    params:
        raw_prefix=lambda wildcards: (
            WORK_ROOT / "tr" / wildcards.sample / f"{wildcards.sample}.trgt"
        ).as_posix()
    script:
        "../scripts/run_trgt.py"
