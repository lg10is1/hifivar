"""Phase 8 independent assembly-derived SV branches."""

from hifivar.exceptions import WorkflowError

ASSEMBLY_SV_CONFIG = config.get("assembly_sv", {})
if not isinstance(ASSEMBLY_SV_CONFIG, dict):
    raise WorkflowError("assembly_sv must be a mapping.")
ASSEMBLY_SV_ENABLED = ASSEMBLY_SV_CONFIG.get("enabled", False)
if not isinstance(ASSEMBLY_SV_ENABLED, bool):
    raise WorkflowError("assembly_sv.enabled must be boolean.")
PAV_CONFIG = ASSEMBLY_SV_CONFIG.get("pav", {})
SVIM_ASM_CONFIG = ASSEMBLY_SV_CONFIG.get("svim_asm", {})
if not isinstance(PAV_CONFIG, dict) or not isinstance(SVIM_ASM_CONFIG, dict):
    raise WorkflowError("assembly_sv caller sections must be mappings.")
if ASSEMBLY_SV_ENABLED and not ASSEMBLY_ENABLED:
    raise WorkflowError("Phase 8 assembly_sv requires assembly.enabled: true.")

PAV_ENABLED = ASSEMBLY_SV_ENABLED and PAV_CONFIG.get("enabled", True)
SVIM_ASM_ENABLED = ASSEMBLY_SV_ENABLED and SVIM_ASM_CONFIG.get("enabled", True)
ASSEMBLY_SV_TARGETS = []
if PAV_ENABLED:
    ASSEMBLY_SV_TARGETS.extend(
        (OUTPUT_ROOT / "assembly_sv" / sample / f"{sample}.pav.assembly.sv.vcf.gz").as_posix()
        for sample in ASSEMBLY_SAMPLE_IDS
    )
if SVIM_ASM_ENABLED:
    ASSEMBLY_SV_TARGETS.extend(
        (OUTPUT_ROOT / "assembly_sv" / sample / f"{sample}.svim_asm.assembly.sv.vcf.gz").as_posix()
        for sample in ASSEMBLY_SAMPLE_IDS
    )
WORKFLOW_TARGETS = [*WORKFLOW_TARGETS, *ASSEMBLY_SV_TARGETS]


rule pav_assembly_sv:
    input:
        reference=lambda wildcards: config["reference"]["fasta"],
        fai=lambda wildcards: f"{config['reference']['fasta']}.fai",
        hap1=f"{ASSEMBLY_DIRECTORY_PATTERN}/{{sample}}.hap1.fa",
        hap2=f"{ASSEMBLY_DIRECTORY_PATTERN}/{{sample}}.hap2.fa"
    output:
        vcf=(OUTPUT_ROOT / "assembly_sv" / "{sample}" / "{sample}.pav.assembly.sv.vcf.gz").as_posix(),
        index=(OUTPUT_ROOT / "assembly_sv" / "{sample}" / "{sample}.pav.assembly.sv.vcf.gz.tbi").as_posix()
    log:
        (LOG_ROOT / "assembly_sv" / "{sample}.pav.log").as_posix()
    threads:
        PAV_CONFIG.get("threads", 32)
    resources:
        mem_mb=PAV_CONFIG.get("memory_mb", 64000),
        runtime_min=PAV_CONFIG.get("runtime_minutes", 2880)
    params:
        work=(WORK_ROOT / "pav" / "{sample}").as_posix()
    script:
        "../scripts/run_pav.py"


rule svim_asm_call:
    input:
        reference=lambda wildcards: config["reference"]["fasta"],
        fai=lambda wildcards: f"{config['reference']['fasta']}.fai",
        hap1=f"{ASSEMBLY_DIRECTORY_PATTERN}/{{sample}}.hap1.fa",
        hap2=f"{ASSEMBLY_DIRECTORY_PATTERN}/{{sample}}.hap2.fa"
    output:
        vcf=(OUTPUT_ROOT / "assembly_sv" / "{sample}" / "{sample}.svim_asm.assembly.sv.vcf.gz").as_posix(),
        index=(OUTPUT_ROOT / "assembly_sv" / "{sample}" / "{sample}.svim_asm.assembly.sv.vcf.gz.tbi").as_posix()
    log:
        (LOG_ROOT / "assembly_sv" / "{sample}.svim_asm.log").as_posix()
    threads:
        SVIM_ASM_CONFIG.get("threads", 16)
    resources:
        mem_mb=SVIM_ASM_CONFIG.get("memory_mb", 32000),
        runtime_min=SVIM_ASM_CONFIG.get("runtime_minutes", 1440)
    params:
        work=(WORK_ROOT / "svim_asm" / "{sample}").as_posix()
    script:
        "../scripts/run_svim_asm.py"
