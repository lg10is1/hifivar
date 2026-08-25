"""Optional Phase 10 IGV/manual-review branch from explicit selections."""

from pathlib import Path

from hifivar.alignment_postprocess import find_alignment_index
from hifivar.context import AnalysisContext
from hifivar.exceptions import WorkflowError
from hifivar.review import read_review_selection_rows
from hifivar.sample import InputType


REVIEW_CONFIG = config.get("review", {})
if not isinstance(REVIEW_CONFIG, dict):
    raise WorkflowError("Effective config section review must be a mapping.")
REVIEW_ENABLED = REVIEW_CONFIG.get("enabled", False)
if not isinstance(REVIEW_ENABLED, bool):
    raise WorkflowError("review.enabled must be boolean.")

REVIEW_SELECTION = None
REVIEW_ROWS = ()
REVIEW_ALIGNMENTS = {}
REVIEW_ALIGNMENT_INDEXES = {}
REVIEW_SOURCE_VCFS = []
REVIEW_SOURCE_INDEXES = []
REVIEW_TRGT_VISUALIZATIONS = []
if REVIEW_ENABLED:
    configured_selection = REVIEW_CONFIG.get("selection_file")
    if not isinstance(configured_selection, str) or not configured_selection.strip():
        raise WorkflowError("review.enabled requires review.selection_file.")
    REVIEW_SELECTION = Path(configured_selection).expanduser()
    REVIEW_ROWS = read_review_selection_rows(REVIEW_SELECTION)
    REVIEW_CONTEXT = AnalysisContext.from_config(config)
    sample_map = {record.sample.sample_id: record.sample for record in REVIEW_CONTEXT.samples}
    selected_samples = sorted({row["sample"].strip() for row in REVIEW_ROWS})
    for sample_id in selected_samples:
        sample = sample_map.get(sample_id)
        if sample is None:
            raise WorkflowError(f"Review selection references unknown sample '{sample_id}'.")
        if sample.input.input_type not in (InputType.BAM, InputType.CRAM):
            raise WorkflowError(
                f"Review sample '{sample_id}' requires aligned BAM/CRAM, not raw FASTQ."
            )
        alignment = sample.input.files[0]
        index = find_alignment_index(alignment)
        if index is None:
            raise WorkflowError(f"Review alignment index is missing for '{sample_id}'.")
        REVIEW_ALIGNMENTS[sample_id] = str(alignment)
        REVIEW_ALIGNMENT_INDEXES[sample_id] = str(index)
    selection_root = REVIEW_SELECTION.parent
    for row in REVIEW_ROWS:
        source = Path(row["source_vcf"].strip()).expanduser()
        if not source.is_absolute():
            source = selection_root / source
        REVIEW_SOURCE_VCFS.append(str(source))
        if str(source).lower().endswith(".vcf.gz"):
            REVIEW_SOURCE_INDEXES.append(f"{source}.tbi")
        visualization = (row.get("trgt_visualization") or "").strip()
        if visualization:
            path = Path(visualization).expanduser()
            REVIEW_TRGT_VISUALIZATIONS.append(str(path if path.is_absolute() else selection_root / path))


REVIEW_OUTPUT = OUTPUT_ROOT / "review"
REVIEW_MANIFEST_JSON = (REVIEW_OUTPUT / "review_manifest.json").as_posix()
REVIEW_MANIFEST_YAML = (REVIEW_OUTPUT / "review_manifest.yaml").as_posix()
REVIEW_MANIFEST_TSV = (REVIEW_OUTPUT / "review_manifest.tsv").as_posix()
REVIEW_BATCH = (REVIEW_OUTPUT / "review.igv.batch").as_posix()
REVIEW_SCREENSHOTS = (REVIEW_OUTPUT / "screenshots").as_posix()
if REVIEW_ENABLED:
    WORKFLOW_TARGETS += [
        REVIEW_MANIFEST_JSON,
        REVIEW_MANIFEST_YAML,
        REVIEW_MANIFEST_TSV,
        REVIEW_BATCH,
        REVIEW_SCREENSHOTS,
    ]


rule manual_review:
    input:
        selection=lambda wildcards: str(REVIEW_SELECTION),
        reference=lambda wildcards: str(REVIEW_CONTEXT.reference.fasta),
        reference_index=lambda wildcards: str(REVIEW_CONTEXT.reference.fai),
        alignments=lambda wildcards: list(REVIEW_ALIGNMENTS.values()),
        alignment_indexes=lambda wildcards: list(REVIEW_ALIGNMENT_INDEXES.values()),
        sources=lambda wildcards: REVIEW_SOURCE_VCFS,
        source_indexes=lambda wildcards: REVIEW_SOURCE_INDEXES,
        trgt_visualizations=lambda wildcards: REVIEW_TRGT_VISUALIZATIONS
    output:
        manifest_json=REVIEW_MANIFEST_JSON,
        manifest_yaml=REVIEW_MANIFEST_YAML,
        manifest_tsv=REVIEW_MANIFEST_TSV,
        batch=REVIEW_BATCH,
        screenshots=directory(REVIEW_SCREENSHOTS)
    log:
        (LOG_ROOT / "review" / "igv.log").as_posix()
    threads:
        int(REVIEW_CONFIG.get("threads", 1))
    resources:
        mem_mb=int(REVIEW_CONFIG.get("memory_mb", 8000)),
        runtime_min=int(REVIEW_CONFIG.get("runtime_minutes", 240))
    script:
        "../scripts/run_review.py"
