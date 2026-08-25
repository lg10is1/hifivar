"""Config-driven HiPhase rule for Phase 6."""

from hifivar.exceptions import WorkflowError


PHASING_CONFIG = config.get("phasing", {})
if not isinstance(PHASING_CONFIG, dict):
    raise WorkflowError("Effective config section phasing must be a mapping.")

PHASING_ENABLED = PHASING_CONFIG.get("enabled", False)
if not isinstance(PHASING_ENABLED, bool):
    raise WorkflowError("phasing.enabled must be boolean.")
if PHASING_ENABLED and not SMALL_ENABLED:
    raise WorkflowError("Phase 6 phasing requires small.enabled: true.")

PHASING_SAMPLE_IDS = list(SMALL_SAMPLE_IDS) if PHASING_ENABLED else []


def _phasing_output(sample, suffix):
    return (OUTPUT_ROOT / "phasing" / f"{sample}{suffix}").as_posix()


PHASED_VCFS = [
    _phasing_output(sample, ".phased.vcf.gz")
    for sample in PHASING_SAMPLE_IDS
]
PHASED_VCF_INDEXES = [f"{path}.tbi" for path in PHASED_VCFS]
PHASED_VCF_PATTERN = (
    OUTPUT_ROOT / "phasing" / "{sample}.phased.vcf.gz"
).as_posix()
WORKFLOW_TARGETS = [
    *WORKFLOW_TARGETS,
    *PHASED_VCFS,
    *PHASED_VCF_INDEXES,
]


rule hiphase_phasing:
    input:
        reference=lambda wildcards: str(SMALL_CONTEXT.reference.fasta),
        reference_index=lambda wildcards: str(SMALL_CONTEXT.reference.fai),
        alignment=lambda wildcards: SMALL_ALIGNMENTS[wildcards.sample],
        alignment_index=lambda wildcards: SMALL_ALIGNMENT_INDEXES[wildcards.sample],
        vcf=SMALL_VCF_PATTERN,
        vcf_index=f"{SMALL_VCF_PATTERN}.tbi",
        gvcf=SMALL_GVCF_PATTERN,
        gvcf_index=f"{SMALL_GVCF_PATTERN}.tbi"
    output:
        vcf=PHASED_VCF_PATTERN,
        vcf_index=f"{PHASED_VCF_PATTERN}.tbi"
    log:
        (LOG_ROOT / "phasing" / "{sample}.hiphase.log").as_posix()
    threads:
        PHASING_CONFIG.get("threads", 16)
    resources:
        mem_mb=PHASING_CONFIG.get("memory_mb", 32000),
        runtime_min=PHASING_CONFIG.get("runtime_minutes", 1440)
    script:
        "../scripts/run_hiphase.py"
