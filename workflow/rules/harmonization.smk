"""Phase 9 per-sample Jasmine harmonization and Truvari concordance."""

from hifivar.context import AnalysisContext
from hifivar.exceptions import WorkflowError

HARMONIZATION_CONFIG = SV_CONFIG.get("harmonization", {})
if not isinstance(HARMONIZATION_CONFIG, dict):
    raise WorkflowError("sv.harmonization must be a mapping.")
HARMONIZATION_ENABLED = HARMONIZATION_CONFIG.get("enabled", False)
if not isinstance(HARMONIZATION_ENABLED, bool):
    raise WorkflowError("sv.harmonization.enabled must be boolean.")

HARMONIZATION_SAMPLES = []
EXTERNAL_HARMONIZATION_INPUTS = HARMONIZATION_CONFIG.get("input_vcfs", {})
if not isinstance(EXTERNAL_HARMONIZATION_INPUTS, dict):
    raise WorkflowError("sv.harmonization.input_vcfs must be a mapping.")

HARMONIZATION_INPUTS = {}
if HARMONIZATION_ENABLED:
    HARMONIZATION_CONTEXT = AnalysisContext.from_config(config)
    read_callers = [
        name for name in ("sawfish", "sniffles2", "pbsv", "cutesv")
        if SV_ENABLED and isinstance(SV_CONFIG.get(name), dict) and SV_CONFIG[name].get("enabled", True)
    ]
    assembly_callers = [
        name for name in ("pav", "svim_asm")
        if ASSEMBLY_SV_ENABLED and isinstance(ASSEMBLY_SV_CONFIG.get(name), dict)
        and ASSEMBLY_SV_CONFIG[name].get("enabled", True)
    ]
    if not read_callers and not assembly_callers and not EXTERNAL_HARMONIZATION_INPUTS:
        raise WorkflowError("Phase 9 requires at least one explicitly enabled caller.")
    for record in HARMONIZATION_CONTEXT.samples:
        sample = record.sample.sample_id
        HARMONIZATION_SAMPLES.append(sample)
        paths = [
            (OUTPUT_ROOT / "sv" / f"{sample}.{caller}.sv.vcf.gz").as_posix()
            for caller in read_callers
        ]
        paths.extend(
            (OUTPUT_ROOT / "assembly_sv" / sample / f"{sample}.{caller}.assembly.sv.vcf.gz").as_posix()
            for caller in assembly_callers
        )
        external = EXTERNAL_HARMONIZATION_INPUTS.get(sample, {})
        if not isinstance(external, dict):
            raise WorkflowError(f"External Phase 9 inputs for '{sample}' must be a mapping.")
        for caller, configured_path in external.items():
            if caller in read_callers or caller in assembly_callers:
                raise WorkflowError(
                    f"Phase 9 caller '{caller}' for sample '{sample}' is both generated and external."
                )
            path = Path(configured_path)
            suffix = ".assembly.sv.vcf.gz" if caller in {"pav", "svim_asm"} else ".sv.vcf.gz"
            expected = f"{sample}.{caller}{suffix}"
            if path.name != expected:
                raise WorkflowError(
                    f"External Phase 9 input for '{caller}' must be named '{expected}'."
                )
            paths.append(path.as_posix())
        if not paths:
            raise WorkflowError(f"Phase 9 has no enabled or external evidence for sample '{sample}'.")
        HARMONIZATION_INPUTS[sample] = paths


def _harmonization_inputs(wildcards):
    vcfs = HARMONIZATION_INPUTS[wildcards.sample]
    return [item for path in vcfs for item in (path, f"{path}.tbi")]


HARMONIZATION_TARGETS = [
    (OUTPUT_ROOT / "sv_harmonized" / sample / f"{sample}.harmonized.sv.vcf.gz").as_posix()
    for sample in HARMONIZATION_SAMPLES
]
WORKFLOW_TARGETS = [*WORKFLOW_TARGETS, *HARMONIZATION_TARGETS]


rule harmonize_sv:
    input:
        sources=_harmonization_inputs,
        reference=lambda wildcards: config["reference"]["fasta"],
        fai=lambda wildcards: f"{config['reference']['fasta']}.fai"
    output:
        vcf=(OUTPUT_ROOT / "sv_harmonized" / "{sample}" / "{sample}.harmonized.sv.vcf.gz").as_posix(),
        index=(OUTPUT_ROOT / "sv_harmonized" / "{sample}" / "{sample}.harmonized.sv.vcf.gz.tbi").as_posix(),
        evidence=(OUTPUT_ROOT / "sv_harmonized" / "{sample}" / "{sample}.sv.evidence.tsv").as_posix(),
        truvari=directory((OUTPUT_ROOT / "sv_harmonized" / "{sample}" / "truvari").as_posix()),
        provenance=(OUTPUT_ROOT / "sv_harmonized" / "{sample}" / "{sample}.phase9.provenance.json").as_posix()
    log:
        (LOG_ROOT / "sv_harmonization" / "{sample}.jasmine.log").as_posix()
    threads:
        HARMONIZATION_CONFIG.get("threads", 8)
    resources:
        mem_mb=HARMONIZATION_CONFIG.get("memory_mb", 32000),
        runtime_min=HARMONIZATION_CONFIG.get("runtime_minutes", 1440)
    params:
        work_root=(WORK_ROOT / "sv_harmonization").as_posix()
    script:
        "../scripts/run_harmonization.py"
