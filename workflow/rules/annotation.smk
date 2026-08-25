"""Optional Phase 11 annotation branch from an explicit input manifest."""

from pathlib import Path

from hifivar.annotation import VariantCategory, read_annotation_input_rows
from hifivar.context import AnalysisContext
from hifivar.exceptions import WorkflowError


ANNOTATION_CONFIG = config.get("annotation", {})
if not isinstance(ANNOTATION_CONFIG, dict):
    raise WorkflowError("Effective config section annotation must be a mapping.")
ANNOTATION_ENABLED = ANNOTATION_CONFIG.get("enabled", False)
if not isinstance(ANNOTATION_ENABLED, bool):
    raise WorkflowError("annotation.enabled must be boolean.")

ANNOTATION_MANIFEST = None
ANNOTATION_ROWS = ()
ANNOTATION_INPUTS = {}
ANNOTATION_OUTPUTS = []
if ANNOTATION_ENABLED:
    configured_manifest = ANNOTATION_CONFIG.get("input_manifest")
    if not isinstance(configured_manifest, str) or not configured_manifest.strip():
        raise WorkflowError("annotation.enabled requires annotation.input_manifest.")
    if ANNOTATION_CONFIG.get("functional_enabled") is True:
        raise WorkflowError(
            "AlphaGenome is currently an interface-only boundary; run it with an "
            "explicit injected backend after Linux/cloud verification."
        )
    ANNOTATION_MANIFEST = Path(configured_manifest).expanduser()
    ANNOTATION_ROWS = read_annotation_input_rows(ANNOTATION_MANIFEST)
    ANNOTATION_CONTEXT = AnalysisContext.from_config(config)
    known_samples = set(ANNOTATION_CONTEXT.sample_ids)
    manifest_root = ANNOTATION_MANIFEST.parent
    for row in ANNOTATION_ROWS:
        sample = row["sample"].strip()
        category = VariantCategory(row["variant_category"].strip().lower()).value
        if sample not in known_samples:
            raise WorkflowError(f"Annotation manifest references unknown sample '{sample}'.")
        source = Path(row["source_vcf"].strip()).expanduser()
        if not source.is_absolute():
            source = manifest_root / source
        key = (sample, category)
        inputs = [str(source)]
        if str(source).lower().endswith(".vcf.gz"):
            inputs.append(f"{source}.tbi")
        ANNOTATION_INPUTS[key] = inputs
        ANNOTATION_OUTPUTS.append(
            (OUTPUT_ROOT / "annotation" / sample / category).as_posix()
        )

WORKFLOW_TARGETS += ANNOTATION_OUTPUTS


def _annotation_inputs(wildcards):
    return ANNOTATION_INPUTS[(wildcards.sample, wildcards.category)]


rule annotate_variants:
    input:
        manifest=lambda wildcards: str(ANNOTATION_MANIFEST),
        reference=lambda wildcards: str(ANNOTATION_CONTEXT.reference.fasta),
        fai=lambda wildcards: str(ANNOTATION_CONTEXT.reference.fai),
        source=_annotation_inputs
    output:
        result=directory((OUTPUT_ROOT / "annotation" / "{sample}" / "{category}").as_posix())
    log:
        (LOG_ROOT / "annotation" / "{sample}.{category}.log").as_posix()
    threads:
        int(ANNOTATION_CONFIG.get("threads", 4))
    resources:
        mem_mb=int(ANNOTATION_CONFIG.get("memory_mb", 16000)),
        runtime_min=int(ANNOTATION_CONFIG.get("runtime_minutes", 480))
    script:
        "../scripts/run_annotation.py"
