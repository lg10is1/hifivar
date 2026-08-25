"""Config-driven, reference-independent hifiasm branch for Phase 7."""

from hifivar.context import AnalysisContext
from hifivar.exceptions import WorkflowError
from hifivar.sample import InputType


ASSEMBLY_CONFIG = config.get("assembly", {})
if not isinstance(ASSEMBLY_CONFIG, dict):
    raise WorkflowError("Effective config section assembly must be a mapping.")

ASSEMBLY_ENABLED = ASSEMBLY_CONFIG.get("enabled", False)
if not isinstance(ASSEMBLY_ENABLED, bool):
    raise WorkflowError("assembly.enabled must be boolean.")

ASSEMBLY_SAMPLE_IDS = []
ASSEMBLY_READS = {}
if ASSEMBLY_ENABLED:
    ASSEMBLY_CONTEXT = AnalysisContext.from_config(config)
    for record in ASSEMBLY_CONTEXT.samples:
        sample = record.sample
        if sample.input.input_type is not InputType.FASTQ:
            raise WorkflowError(
                f"Phase 7 assembly sample '{sample.sample_id}' requires primary "
                "HiFi FASTQ; BAM/CRAM extraction is disabled."
            )
        ASSEMBLY_SAMPLE_IDS.append(sample.sample_id)
        ASSEMBLY_READS[sample.sample_id] = [
            str(path) for path in sample.input.files
        ]


def _assembly_fasta(sample, role):
    return (
        OUTPUT_ROOT / "assembly" / sample / f"{sample}.{role}.fa"
    ).as_posix()


ASSEMBLY_FASTAS = [
    _assembly_fasta(sample, role)
    for sample in ASSEMBLY_SAMPLE_IDS
    for role in ("primary", "hap1", "hap2")
]
ASSEMBLY_PREFIX_PATTERN = (
    WORK_ROOT / "hifiasm" / "{sample}" / "{sample}.asm"
).as_posix()
ASSEMBLY_DIRECTORY_PATTERN = (
    OUTPUT_ROOT / "assembly" / "{sample}"
).as_posix()
WORKFLOW_TARGETS = [*WORKFLOW_TARGETS, *ASSEMBLY_FASTAS]


rule hifiasm_assemble:
    input:
        reads=lambda wildcards: ASSEMBLY_READS[wildcards.sample]
    output:
        primary_gfa=f"{ASSEMBLY_PREFIX_PATTERN}.bp.p_ctg.gfa",
        hap1_gfa=f"{ASSEMBLY_PREFIX_PATTERN}.bp.hap1.p_ctg.gfa",
        hap2_gfa=f"{ASSEMBLY_PREFIX_PATTERN}.bp.hap2.p_ctg.gfa",
        primary_fasta=f"{ASSEMBLY_DIRECTORY_PATTERN}/{{sample}}.primary.fa",
        hap1_fasta=f"{ASSEMBLY_DIRECTORY_PATTERN}/{{sample}}.hap1.fa",
        hap2_fasta=f"{ASSEMBLY_DIRECTORY_PATTERN}/{{sample}}.hap2.fa"
    log:
        (LOG_ROOT / "assembly" / "{sample}.hifiasm.log").as_posix()
    threads:
        ASSEMBLY_CONFIG.get("threads", 32)
    resources:
        mem_mb=ASSEMBLY_CONFIG.get("memory_mb", 128000),
        runtime_min=ASSEMBLY_CONFIG.get("runtime_minutes", 4320)
    params:
        output_prefix=ASSEMBLY_PREFIX_PATTERN,
        assembly_directory=ASSEMBLY_DIRECTORY_PATTERN
    script:
        "../scripts/run_hifiasm.py"
