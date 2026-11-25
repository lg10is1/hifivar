"""Phase 9 per-sample Jasmine harmonization and Truvari concordance."""

from hifivar.context import AnalysisContext
from hifivar.breakends import BndPolicy
from hifivar.phase9 import phase9_bnd_paths
from hifivar.sv_support_selection import support_selection_paths
from hifivar.exceptions import WorkflowError
from hifivar.sv_input_identity import SourceIdentityPolicy

HARMONIZATION_CONFIG = SV_CONFIG.get("harmonization", {})
if not isinstance(HARMONIZATION_CONFIG, dict):
    raise WorkflowError("sv.harmonization must be a mapping.")
HARMONIZATION_ENABLED = HARMONIZATION_CONFIG.get("enabled", False)
HARMONIZATION_IDENTITY_POLICY = SourceIdentityPolicy.from_config(
    HARMONIZATION_CONFIG.get("input_identity", {"enabled": False})
)
if HARMONIZATION_IDENTITY_POLICY is not None and not HARMONIZATION_ENABLED:
    raise WorkflowError("SV source identity requires sv.harmonization.enabled: true.")
if not isinstance(HARMONIZATION_ENABLED, bool):
    raise WorkflowError("sv.harmonization.enabled must be boolean.")
HARMONIZATION_BND_POLICY = BndPolicy.from_config(
    HARMONIZATION_CONFIG.get("bnd_representation", {"enabled": False})
)
if HARMONIZATION_BND_POLICY is not None and not HARMONIZATION_ENABLED:
    raise WorkflowError("BND representation requires sv.harmonization.enabled: true.")

HARMONIZATION_MINIMUM_CALLER_SUPPORT = HARMONIZATION_CONFIG.get("minimum_caller_support", 1)
if type(HARMONIZATION_MINIMUM_CALLER_SUPPORT) is not int or not 1 <= HARMONIZATION_MINIMUM_CALLER_SUPPORT <= 6:
    raise WorkflowError("sv.harmonization.minimum_caller_support must be an integer between 1 and 6.")
if HARMONIZATION_MINIMUM_CALLER_SUPPORT > 1 and not HARMONIZATION_ENABLED:
    raise WorkflowError("Caller-support selection requires sv.harmonization.enabled: true.")

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
    if set(EXTERNAL_HARMONIZATION_INPUTS) - set(HARMONIZATION_CONTEXT.sample_ids):
        raise WorkflowError("External Phase 9 sample is not declared in the sample sheet.")
    for record in HARMONIZATION_CONTEXT.samples:
        sample = record.sample.sample_id
        if HARMONIZATION_IDENTITY_POLICY is not None:
            identity_bundle = OUTPUT_ROOT / "sv_harmonized" / sample / "source-identity-v1"
            if identity_bundle.exists() or identity_bundle.is_symlink():
                raise WorkflowError("SV source identity requires a fresh per-sample bundle before scheduling; overwrite is unsupported.")
        if HARMONIZATION_BND_POLICY is not None:
            # Snakemake may delete declared directory outputs before executing
            # a job. Reject an existing bundle while constructing the DAG,
            # before that cleanup can bypass the wrapper's no-overwrite guard.
            bundle = phase9_bnd_paths(OUTPUT_ROOT / "sv_harmonized", sample)["directory"]
            if bundle.exists() or bundle.is_symlink():
                raise WorkflowError("BND representation requires a fresh per-sample bundle before scheduling; overwrite is unsupported.")
        if HARMONIZATION_MINIMUM_CALLER_SUPPORT > 1:
            selected_bundle = support_selection_paths(OUTPUT_ROOT / "sv_harmonized", sample, HARMONIZATION_MINIMUM_CALLER_SUPPORT)["directory"]
            if selected_bundle.exists() or selected_bundle.is_symlink():
                raise WorkflowError("Caller-support selection requires a fresh per-sample bundle before scheduling; overwrite is unsupported.")
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
        from hifivar.harmonization_inputs import CALLER_ORDER
        if set(external) - set(CALLER_ORDER):
            raise WorkflowError("Unsupported external Phase 9 caller.")
        for caller in (name for name in CALLER_ORDER if name in external):
            configured_path = external[caller]
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
        if HARMONIZATION_MINIMUM_CALLER_SUPPORT > len(paths):
            raise WorkflowError(f"minimum_caller_support exceeds the actual input callers for sample '{sample}'.")
        HARMONIZATION_INPUTS[sample] = paths


def _harmonization_inputs(wildcards):
    vcfs = HARMONIZATION_INPUTS[wildcards.sample]
    return [item for path in vcfs for item in (path, f"{path}.tbi")]


HARMONIZATION_TARGETS = [
    (OUTPUT_ROOT / "sv_harmonized" / sample / f"{sample}.harmonized.sv.vcf.gz").as_posix()
    for sample in HARMONIZATION_SAMPLES
]
if HARMONIZATION_MINIMUM_CALLER_SUPPORT > 1:
    HARMONIZATION_TARGETS.extend(
        support_selection_paths(OUTPUT_ROOT / "sv_harmonized", sample, HARMONIZATION_MINIMUM_CALLER_SUPPORT)["directory"].as_posix()
        for sample in HARMONIZATION_SAMPLES
    )
WORKFLOW_TARGETS = [*WORKFLOW_TARGETS, *HARMONIZATION_TARGETS]


if HARMONIZATION_ENABLED:
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
            provenance=(OUTPUT_ROOT / "sv_harmonized" / "{sample}" / "{sample}.phase9.provenance.json").as_posix(),
            support_bundle=directory((OUTPUT_ROOT / "sv_harmonized" / "{sample}" / f"support-ge-{HARMONIZATION_MINIMUM_CALLER_SUPPORT}").as_posix()) if HARMONIZATION_MINIMUM_CALLER_SUPPORT > 1 else [],
            bnd_bundle=directory((OUTPUT_ROOT / "sv_harmonized" / "{sample}" / "bndrepr-v1").as_posix()) if HARMONIZATION_BND_POLICY is not None else [],
            input_identity_bundle=directory((OUTPUT_ROOT / "sv_harmonized" / "{sample}" / "source-identity-v1").as_posix()) if HARMONIZATION_IDENTITY_POLICY is not None else []
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
