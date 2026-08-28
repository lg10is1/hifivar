"""Optional, independent Phase 12 cohort tracks."""

from pathlib import Path

from hifivar.cohort import CohortDefinition, CohortTrack, read_cohort_input_manifest
from hifivar.context import AnalysisContext
from hifivar.exceptions import WorkflowError


COHORT_CONFIG = config.get("cohort", {})
if not isinstance(COHORT_CONFIG, dict):
    raise WorkflowError("Effective config section cohort must be a mapping.")
COHORT_ENABLED = COHORT_CONFIG.get("enabled", False)
COHORT_TRACKS = {}
COHORT_INPUTS = {}
COHORT_ID = "_disabled"
COHORT_MANIFEST = None
COHORT_CONTEXT = None
COHORT_ROOT = OUTPUT_ROOT / "cohort"
COHORT_RESULT_FILES = []
COHORT_RUN_MANIFEST = (COHORT_ROOT / COHORT_ID / "cohort_manifest.json").as_posix()
COHORT_RUN_MANIFEST_YAML = (COHORT_ROOT / COHORT_ID / "cohort_manifest.yaml").as_posix()

if COHORT_ENABLED:
    COHORT_ID = COHORT_CONFIG.get("cohort_id")
    configured_manifest = COHORT_CONFIG.get("input_manifest")
    if not isinstance(COHORT_ID, str) or not COHORT_ID.strip():
        raise WorkflowError("cohort.enabled requires cohort.cohort_id.")
    if not isinstance(configured_manifest, str) or not configured_manifest.strip():
        raise WorkflowError("cohort.enabled requires cohort.input_manifest.")
    COHORT_MANIFEST = Path(configured_manifest).expanduser()
    COHORT_CONTEXT = AnalysisContext.from_config(config)
    COHORT_DEFINITION = CohortDefinition(COHORT_ID, COHORT_CONTEXT.sample_ids, COHORT_CONTEXT.reference)
    for track in CohortTrack:
        subsection = COHORT_CONFIG.get(track.value, {})
        if not isinstance(subsection, dict):
            raise WorkflowError(f"cohort.{track.value} must be a mapping.")
        COHORT_TRACKS[track] = subsection.get("enabled", False)
        if COHORT_TRACKS[track]:
            items = read_cohort_input_manifest(COHORT_MANIFEST, COHORT_DEFINITION, track)
            COHORT_INPUTS[track] = [str(value) for item in items for value in (item.source_path, item.index_path) if value is not None]


def _cohort_track_inputs(track):
    return lambda wildcards: COHORT_INPUTS.get(track, [])


SMALL_ROOT = COHORT_ROOT / str(COHORT_ID) / "small"
SV_ROOT = COHORT_ROOT / str(COHORT_ID) / "sv"
TR_ROOT = COHORT_ROOT / str(COHORT_ID) / "tr"


def _cohort_small_memory_gb():
    if not COHORT_TRACKS.get(CohortTrack.SMALL_VARIANTS):
        return 1
    value = COHORT_CONFIG.get("small_variants", {}).get("memory_gb")
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise WorkflowError(
            "Enabled cohort small-variant calling requires an explicit positive "
            "cohort.small_variants.memory_gb."
        )
    return value


COHORT_SMALL_MEMORY_GB = _cohort_small_memory_gb()

if COHORT_TRACKS.get(CohortTrack.SMALL_VARIANTS):
    SMALL_RESULTS = [
        (SMALL_ROOT / f"{COHORT_ID}.small.bcf").as_posix(),
        (SMALL_ROOT / f"{COHORT_ID}.small.bcf.csi").as_posix(),
        (SMALL_ROOT / f"{COHORT_ID}.small.vcf.gz").as_posix(),
        (SMALL_ROOT / f"{COHORT_ID}.small.vcf.gz.tbi").as_posix(),
        (SMALL_ROOT / "track_result.json").as_posix(),
    ]
    WORKFLOW_TARGETS += SMALL_RESULTS
    COHORT_RESULT_FILES.append(SMALL_RESULTS[-1])

if COHORT_TRACKS.get(CohortTrack.SV):
    SV_RESULTS = [
        (SV_ROOT / f"{COHORT_ID}.sv.sites.tsv").as_posix(),
        (SV_ROOT / f"{COHORT_ID}.sv.sample-matrix.tsv").as_posix(),
        (SV_ROOT / "track_result.json").as_posix(),
    ]
    WORKFLOW_TARGETS += SV_RESULTS
    COHORT_RESULT_FILES.append(SV_RESULTS[-1])

if COHORT_TRACKS.get(CohortTrack.TR):
    TR_RESULTS = [
        (TR_ROOT / f"{COHORT_ID}.tr.loci.tsv").as_posix(),
        (TR_ROOT / f"{COHORT_ID}.tr.sample-matrix.tsv").as_posix(),
        (TR_ROOT / "track_result.json").as_posix(),
    ]
    WORKFLOW_TARGETS += TR_RESULTS
    COHORT_RESULT_FILES.append(TR_RESULTS[-1])

if COHORT_RESULT_FILES:
    COHORT_RUN_MANIFEST = (COHORT_ROOT / str(COHORT_ID) / "cohort_manifest.json").as_posix()
    COHORT_RUN_MANIFEST_YAML = (COHORT_ROOT / str(COHORT_ID) / "cohort_manifest.yaml").as_posix()
    WORKFLOW_TARGETS += [COHORT_RUN_MANIFEST, COHORT_RUN_MANIFEST_YAML]


rule cohort_small_variants:
    input:
        manifest=lambda wildcards: str(COHORT_MANIFEST),
        reference=lambda wildcards: str(COHORT_CONTEXT.reference.fasta),
        fai=lambda wildcards: str(COHORT_CONTEXT.reference.fai),
        artifacts=_cohort_track_inputs(CohortTrack.SMALL_VARIANTS)
    output:
        bcf=(SMALL_ROOT / f"{COHORT_ID}.small.bcf").as_posix(),
        bcf_index=(SMALL_ROOT / f"{COHORT_ID}.small.bcf.csi").as_posix(),
        vcf=(SMALL_ROOT / f"{COHORT_ID}.small.vcf.gz").as_posix(),
        vcf_index=(SMALL_ROOT / f"{COHORT_ID}.small.vcf.gz.tbi").as_posix(),
        result=(SMALL_ROOT / "track_result.json").as_posix()
    log:
        (LOG_ROOT / "cohort" / str(COHORT_ID) / "glnexus.log").as_posix()
    threads:
        int(COHORT_CONFIG.get("small_variants", {}).get("threads", 8))
    resources:
        mem_mb=COHORT_SMALL_MEMORY_GB * 1024,
        runtime_min=int(COHORT_CONFIG.get("small_variants", {}).get("runtime_minutes", 1440))
    conda:
        "../envs/glnexus.yaml"
    script:
        "../scripts/run_cohort_small.py"


rule cohort_sv:
    input:
        manifest=lambda wildcards: str(COHORT_MANIFEST),
        reference=lambda wildcards: str(COHORT_CONTEXT.reference.fasta),
        fai=lambda wildcards: str(COHORT_CONTEXT.reference.fai),
        artifacts=_cohort_track_inputs(CohortTrack.SV)
    output:
        sites=(SV_ROOT / f"{COHORT_ID}.sv.sites.tsv").as_posix(),
        matrix=(SV_ROOT / f"{COHORT_ID}.sv.sample-matrix.tsv").as_posix(),
        result=(SV_ROOT / "track_result.json").as_posix()
    threads: 1
    resources:
        mem_mb=int(COHORT_CONFIG.get("sv", {}).get("memory_mb", 8000)),
        runtime_min=int(COHORT_CONFIG.get("sv", {}).get("runtime_minutes", 240))
    script:
        "../scripts/run_cohort_sv.py"


rule cohort_tr:
    input:
        manifest=lambda wildcards: str(COHORT_MANIFEST),
        reference=lambda wildcards: str(COHORT_CONTEXT.reference.fasta),
        fai=lambda wildcards: str(COHORT_CONTEXT.reference.fai),
        artifacts=_cohort_track_inputs(CohortTrack.TR)
    output:
        loci=(TR_ROOT / f"{COHORT_ID}.tr.loci.tsv").as_posix(),
        matrix=(TR_ROOT / f"{COHORT_ID}.tr.sample-matrix.tsv").as_posix(),
        result=(TR_ROOT / "track_result.json").as_posix()
    threads: 1
    resources:
        mem_mb=int(COHORT_CONFIG.get("tr", {}).get("memory_mb", 8000)),
        runtime_min=int(COHORT_CONFIG.get("tr", {}).get("runtime_minutes", 240))
    script:
        "../scripts/run_cohort_tr.py"


rule cohort_manifest:
    input:
        lambda wildcards: COHORT_RESULT_FILES
    output:
        json=COHORT_RUN_MANIFEST,
        yaml=COHORT_RUN_MANIFEST_YAML
    script:
        "../scripts/run_cohort_manifest.py"
