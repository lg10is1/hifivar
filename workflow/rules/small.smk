"""Config-driven single-sample DeepVariant rules for Phase 3."""

from hifivar.alignment_postprocess import find_alignment_index
from hifivar.context import AnalysisContext
from hifivar.exceptions import WorkflowError
from hifivar.sample import InputType


SMALL_CONFIG = config.get("small", {})
if not isinstance(SMALL_CONFIG, dict):
    raise WorkflowError("Effective config section small must be a mapping.")

DEEPVARIANT_MAX_CONCURRENT = int(
    SMALL_CONFIG.get("max_concurrent_samples", 1)
)
if DEEPVARIANT_MAX_CONCURRENT <= 0:
    raise WorkflowError("small.max_concurrent_samples must be positive.")
if "deepvariant_slots" not in workflow.global_resources:
    workflow.register_resource(
        "deepvariant_slots", DEEPVARIANT_MAX_CONCURRENT
    )

SMALL_ENABLED = SMALL_CONFIG.get("enabled", False)
if not isinstance(SMALL_ENABLED, bool):
    raise WorkflowError("small.enabled must be boolean.")

SMALL_SAMPLE_IDS = []
SMALL_ALIGNMENTS = {}
SMALL_ALIGNMENT_INDEXES = {}
if SMALL_ENABLED:
    SMALL_CONTEXT = AnalysisContext.from_config(config)
    for record in SMALL_CONTEXT.samples:
        sample = record.sample
        if sample.input.input_type not in (InputType.BAM, InputType.CRAM):
            raise WorkflowError(
                f"DeepVariant Snakemake input for sample '{sample.sample_id}' "
                "must be aligned BAM/CRAM, not raw FASTQ."
            )
        alignment = sample.input.files[0]
        index = find_alignment_index(alignment)
        if index is None:
            raise WorkflowError(
                f"DeepVariant alignment index is missing for sample "
                f"'{sample.sample_id}': '{alignment}'."
            )
        SMALL_SAMPLE_IDS.append(sample.sample_id)
        SMALL_ALIGNMENTS[sample.sample_id] = str(alignment)
        SMALL_ALIGNMENT_INDEXES[sample.sample_id] = str(index)


def _small_output(sample, suffix):
    return (OUTPUT_ROOT / "small" / f"{sample}{suffix}").as_posix()


def _small_tmpdir(sample):
    configured_root = RUNTIME_CONFIG.get("tmpdir")
    temporary_root = Path(configured_root) if configured_root else WORK_ROOT
    return (temporary_root / "deepvariant" / sample / "tmp").as_posix()


SMALL_VCFS = [_small_output(sample, ".small.vcf.gz") for sample in SMALL_SAMPLE_IDS]
SMALL_GVCFS = [_small_output(sample, ".g.vcf.gz") for sample in SMALL_SAMPLE_IDS]
SMALL_VCF_INDEXES = [f"{path}.tbi" for path in SMALL_VCFS]
SMALL_GVCF_INDEXES = [f"{path}.tbi" for path in SMALL_GVCFS]
SMALL_VCF_PATTERN = (OUTPUT_ROOT / "small" / "{sample}.small.vcf.gz").as_posix()
SMALL_GVCF_PATTERN = (OUTPUT_ROOT / "small" / "{sample}.g.vcf.gz").as_posix()
WORKFLOW_TARGETS = [
    PHASE0_SMOKE_MARKER,
    *SMALL_VCFS,
    *SMALL_GVCFS,
    *SMALL_VCF_INDEXES,
    *SMALL_GVCF_INDEXES,
]


rule deepvariant_small:
    input:
        reference=lambda wildcards: str(SMALL_CONTEXT.reference.fasta),
        reference_index=lambda wildcards: str(SMALL_CONTEXT.reference.fai),
        alignment=lambda wildcards: SMALL_ALIGNMENTS[wildcards.sample],
        alignment_index=lambda wildcards: SMALL_ALIGNMENT_INDEXES[wildcards.sample]
    output:
        vcf=SMALL_VCF_PATTERN,
        gvcf=SMALL_GVCF_PATTERN,
        vcf_index=f"{SMALL_VCF_PATTERN}.tbi",
        gvcf_index=f"{SMALL_GVCF_PATTERN}.tbi"
    log:
        (LOG_ROOT / "small" / "{sample}.deepvariant.log").as_posix()
    threads:
        SMALL_CONFIG.get("threads", 16)
    resources:
        mem_mb=SMALL_CONFIG.get("memory_mb", 64000),
        runtime_min=SMALL_CONFIG.get("runtime_minutes", 2880),
        deepvariant_slots=1
    params:
        intermediate=lambda wildcards: (
            WORK_ROOT / "deepvariant" / wildcards.sample
        ).as_posix(),
        tmpdir=lambda wildcards: _small_tmpdir(wildcards.sample)
    script:
        "../scripts/run_deepvariant.py"
