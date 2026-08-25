"""Config-driven independent read-based SV caller rules for Phase 4."""

from hifivar.alignment_postprocess import find_alignment_index
from hifivar.context import AnalysisContext
from hifivar.exceptions import WorkflowError
from hifivar.sample import InputType


SV_CONFIG = config.get("sv", {})
if not isinstance(SV_CONFIG, dict):
    raise WorkflowError("Effective config section sv must be a mapping.")
SV_ENABLED = SV_CONFIG.get("enabled", False)
if not isinstance(SV_ENABLED, bool):
    raise WorkflowError("sv.enabled must be boolean.")

SV_CALLER_ORDER = ("sawfish", "sniffles2", "pbsv", "cutesv")
SV_CALLERS = [
    caller for caller in SV_CALLER_ORDER
    if isinstance(SV_CONFIG.get(caller), dict)
    and SV_CONFIG[caller].get("enabled", True)
]
SV_SAMPLE_IDS = []
SV_ALIGNMENTS = {}
SV_ALIGNMENT_INDEXES = {}
if SV_ENABLED:
    if not SV_CALLERS:
        raise WorkflowError("sv.enabled requires at least one enabled caller.")
    SV_CONTEXT = AnalysisContext.from_config(config)
    for record in SV_CONTEXT.samples:
        sample = record.sample
        if sample.input.input_type not in (InputType.BAM, InputType.CRAM):
            raise WorkflowError(
                f"Read-based SV Snakemake input for sample '{sample.sample_id}' "
                "must be aligned BAM/CRAM, not raw FASTQ."
            )
        if sample.input.input_type is InputType.CRAM and any(
            caller in SV_CALLERS for caller in ("pbsv", "cutesv")
        ):
            raise WorkflowError(
                f"Sample '{sample.sample_id}' is CRAM, but enabled pbsv/cuteSV "
                "Phase 4 rules require BAM. Disable those callers or provide BAM."
            )
        alignment = sample.input.files[0]
        index = find_alignment_index(alignment)
        if index is None:
            raise WorkflowError(
                f"Read-based SV alignment index is missing for sample "
                f"'{sample.sample_id}': '{alignment}'."
            )
        SV_SAMPLE_IDS.append(sample.sample_id)
        SV_ALIGNMENTS[sample.sample_id] = str(alignment)
        SV_ALIGNMENT_INDEXES[sample.sample_id] = str(index)


def _sv_output(sample, caller):
    return (OUTPUT_ROOT / "sv" / f"{sample}.{caller}.sv.vcf.gz").as_posix()


def _sv_caller_config(wildcards):
    return SV_CONFIG[wildcards.caller]


SV_VCFS = [
    _sv_output(sample, caller)
    for sample in SV_SAMPLE_IDS
    for caller in SV_CALLERS
]
SV_VCF_INDEXES = [f"{path}.tbi" for path in SV_VCFS]
SV_VCF_PATTERN = (OUTPUT_ROOT / "sv" / "{sample}.{caller}.sv.vcf.gz").as_posix()
WORKFLOW_TARGETS += [*SV_VCFS, *SV_VCF_INDEXES]


rule read_based_sv:
    input:
        reference=lambda wildcards: str(SV_CONTEXT.reference.fasta),
        reference_index=lambda wildcards: str(SV_CONTEXT.reference.fai),
        alignment=lambda wildcards: SV_ALIGNMENTS[wildcards.sample],
        alignment_index=lambda wildcards: SV_ALIGNMENT_INDEXES[wildcards.sample]
    output:
        vcf=SV_VCF_PATTERN,
        index=f"{SV_VCF_PATTERN}.tbi"
    log:
        (LOG_ROOT / "sv" / "{sample}.{caller}.log").as_posix()
    threads:
        lambda wildcards: int(_sv_caller_config(wildcards).get("threads", 8))
    resources:
        mem_mb=lambda wildcards: int(_sv_caller_config(wildcards).get("memory_mb", 16000)),
        runtime_min=lambda wildcards: int(_sv_caller_config(wildcards).get("runtime_minutes", 1440))
    params:
        workdir=lambda wildcards: (
            WORK_ROOT / "sv" / wildcards.caller / wildcards.sample
        ).as_posix()
    wildcard_constraints:
        caller="sawfish|sniffles2|pbsv|cutesv"
    script:
        "../scripts/run_sv_caller.py"
