"""Optional Phase 15 deterministic aligned-BAM or raw-FASTQ branch."""

from pathlib import Path

from hifivar.exceptions import WorkflowError
from hifivar.low_depth import (
    build_low_depth_benchmark_plan,
    build_low_depth_plan_from_config,
    build_raw_fastq_subsampling_plan_from_config,
)


LOW_DEPTH_CONFIG = config.get("low_depth", {})
if not isinstance(LOW_DEPTH_CONFIG, dict):
    raise WorkflowError("Effective config section low_depth must be a mapping.")
LOW_DEPTH_ENABLED = LOW_DEPTH_CONFIG.get("enabled", False)
if not isinstance(LOW_DEPTH_ENABLED, bool):
    raise WorkflowError("low_depth.enabled must be boolean.")

LOW_DEPTH_PLAN = None
LOW_DEPTH_PLAN_JSON = None
LOW_DEPTH_RUN_BY_OUTPUT = {}
LOW_DEPTH_BAMS = []
LOW_DEPTH_INDEXES = []
LOW_DEPTH_BAM_PATTERN = None
LOW_DEPTH_INDEX_PATTERN = None
LOW_DEPTH_INPUT_MODE = str(LOW_DEPTH_CONFIG.get("input_mode", "aligned_bam"))
LOW_DEPTH_RAW_FASTQS = []
LOW_DEPTH_RAW_PATTERN = []
LOW_DEPTH_E2E_ENABLED = False
LOW_DEPTH_E2E_CONFIG = {}
LOW_DEPTH_BENCHMARK_PLAN = None
LOW_DEPTH_BENCHMARK_RUNS = {}
if LOW_DEPTH_ENABLED:
    study_id = LOW_DEPTH_CONFIG.get("study_id")
    if not isinstance(study_id, str) or not study_id.strip():
        raise WorkflowError("low_depth.enabled requires low_depth.study_id.")
    LOW_DEPTH_ROOT = OUTPUT_ROOT / "low_depth" / study_id
    if LOW_DEPTH_INPUT_MODE == "aligned_bam":
        LOW_DEPTH_PLAN = build_low_depth_plan_from_config(
            LOW_DEPTH_CONFIG,
            LOW_DEPTH_ROOT / "subsamples",
        )
    elif LOW_DEPTH_INPUT_MODE == "raw_fastq":
        LOW_DEPTH_PLAN = build_raw_fastq_subsampling_plan_from_config(
            LOW_DEPTH_CONFIG,
            LOW_DEPTH_ROOT / "raw_subsamples",
        )
    else:
        raise WorkflowError(f"Unsupported low_depth.input_mode: {LOW_DEPTH_INPUT_MODE!r}.")
    LOW_DEPTH_PLAN_JSON = (LOW_DEPTH_ROOT / "low_depth_plan.json").as_posix()
    if LOW_DEPTH_INPUT_MODE == "aligned_bam":
        LOW_DEPTH_BAM_PATTERN = (
            LOW_DEPTH_ROOT / "subsamples" / "{coverage}" /
            "replicate-{replicate}.bam"
        ).as_posix()
        LOW_DEPTH_INDEX_PATTERN = f"{LOW_DEPTH_BAM_PATTERN}.bai"
        LOW_DEPTH_BAMS = list(dict.fromkeys(run.output_bam.as_posix() for run in LOW_DEPTH_PLAN.runs))
        LOW_DEPTH_INDEXES = list(dict.fromkeys(run.output_index.as_posix() for run in LOW_DEPTH_PLAN.runs))
        LOW_DEPTH_RUN_BY_ID = {
            (run.output_bam.parent.name, str(run.replicate)): run
            for run in LOW_DEPTH_PLAN.runs
        }
        WORKFLOW_TARGETS += [LOW_DEPTH_PLAN_JSON, *LOW_DEPTH_BAMS, *LOW_DEPTH_INDEXES]
    else:
        LOW_DEPTH_RAW_FASTQS = [item.output_fastq.as_posix() for item in LOW_DEPTH_PLAN.subsets]
        LOW_DEPTH_RAW_PATTERN = expand(
            (
                LOW_DEPTH_ROOT / "raw_subsamples" / "{coverage}" /
                "replicate-{replicate}.fastq.gz"
            ).as_posix(),
            coverage=list(dict.fromkeys(item.output_fastq.parent.name for item in LOW_DEPTH_PLAN.subsets)),
            allow_missing=True,
        )
        WORKFLOW_TARGETS += [LOW_DEPTH_PLAN_JSON, *LOW_DEPTH_RAW_FASTQS]
        LOW_DEPTH_E2E_CONFIG = LOW_DEPTH_CONFIG.get("end_to_end", {})
        if not isinstance(LOW_DEPTH_E2E_CONFIG, dict):
            raise WorkflowError("low_depth.end_to_end must be a mapping.")
        LOW_DEPTH_E2E_ENABLED = bool(LOW_DEPTH_E2E_CONFIG.get("enabled", False))
        if LOW_DEPTH_E2E_ENABLED:
            LOW_DEPTH_BENCHMARK_PLAN = build_low_depth_benchmark_plan(
                LOW_DEPTH_PLAN,
                LOW_DEPTH_ROOT,
                metrics_compression=str(
                    LOW_DEPTH_E2E_CONFIG.get("metrics_compression", "gzip")
                ),
            )
            LOW_DEPTH_BENCHMARK_RUNS = {
                (run.alignment_bam.parent.name, str(run.replicate)): run
                for run in LOW_DEPTH_BENCHMARK_PLAN.runs
            }
            WORKFLOW_TARGETS += [
                run.observation_json.as_posix()
                for run in LOW_DEPTH_BENCHMARK_PLAN.runs
            ]


if LOW_DEPTH_ENABLED:
    rule low_depth_plan:
        output:
            plan=LOW_DEPTH_PLAN_JSON
        params:
            plan=LOW_DEPTH_PLAN
        script:
            "../scripts/write_low_depth_plan.py"


if LOW_DEPTH_ENABLED and LOW_DEPTH_INPUT_MODE == "aligned_bam":
    rule low_depth_downsample:
        input:
            source=str(LOW_DEPTH_PLAN.source_bam),
            plan=LOW_DEPTH_PLAN_JSON
        output:
            bam=LOW_DEPTH_BAM_PATTERN,
            index=LOW_DEPTH_INDEX_PATTERN
        log:
            (LOG_ROOT / "low_depth" / "{coverage}.r{replicate}.log").as_posix()
        threads:
            int(LOW_DEPTH_CONFIG.get("threads", 1))
        resources:
            mem_mb=int(LOW_DEPTH_CONFIG.get("memory_mb", 8000)),
            runtime_min=int(LOW_DEPTH_CONFIG.get("runtime_minutes", 120))
        params:
            executable=str(LOW_DEPTH_CONFIG.get("samtools_executable", "samtools")),
            overwrite=bool(LOW_DEPTH_CONFIG.get("overwrite", False)),
            run=lambda wildcards: LOW_DEPTH_RUN_BY_ID[
                (wildcards.coverage, str(wildcards.replicate))
            ]
        script:
            "../scripts/run_low_depth_downsample.py"


if LOW_DEPTH_ENABLED and LOW_DEPTH_INPUT_MODE == "raw_fastq":
    rule low_depth_raw_fastq_nested_subsamples:
        input:
            sources=[path.as_posix() for path in LOW_DEPTH_PLAN.source_fastqs],
            plan=LOW_DEPTH_PLAN_JSON
        output:
            LOW_DEPTH_RAW_PATTERN
        log:
            (LOG_ROOT / "low_depth" / "raw_fastq.r{replicate}.log").as_posix()
        resources:
            mem_mb=int(LOW_DEPTH_CONFIG.get("memory_mb", 8000)),
            runtime_min=int(LOW_DEPTH_CONFIG.get("runtime_minutes", 120))
        params:
            plan=LOW_DEPTH_PLAN,
            replicate=lambda wildcards: int(wildcards.replicate),
            overwrite=bool(LOW_DEPTH_CONFIG.get("overwrite", False))
        script:
            "../scripts/run_low_depth_raw_fastq.py"


if LOW_DEPTH_ENABLED and LOW_DEPTH_INPUT_MODE == "raw_fastq" and LOW_DEPTH_E2E_ENABLED:
    _ld_alignment = (LOW_DEPTH_ROOT / "alignment" / "{coverage}" / "replicate-{replicate}.bam").as_posix()
    _ld_small_vcf = (LOW_DEPTH_ROOT / "small" / "{coverage}" / "replicate-{replicate}.small.vcf.gz").as_posix()
    _ld_small_gvcf = (LOW_DEPTH_ROOT / "small" / "{coverage}" / "replicate-{replicate}.g.vcf.gz").as_posix()
    _ld_coverage = (LOW_DEPTH_ROOT / "coverage" / "{coverage}" / "replicate-{replicate}.coverage.json").as_posix()
    _ld_happy_prefix = (LOW_DEPTH_ROOT / "benchmark" / "{coverage}" / "replicate-{replicate}.happy").as_posix()
    _ld_happy_metrics = _ld_happy_prefix + (
        ".metrics.json.gz"
        if LOW_DEPTH_BENCHMARK_PLAN.metrics_compression == "gzip"
        else ".metrics.json"
    )
    _ld_observation = (LOW_DEPTH_ROOT / "observations" / "{coverage}" / "replicate-{replicate}.json").as_posix()

    rule low_depth_align:
        input:
            fastq=lambda wc: LOW_DEPTH_BENCHMARK_RUNS[(wc.coverage, str(wc.replicate))].input_fastq,
            reference=config["reference"]["fasta"],
            fai=lambda wc: str(config["reference"]["fasta"]) + ".fai"
        output:
            bam=_ld_alignment,
            index=_ld_alignment + ".bai",
            provenance=(LOW_DEPTH_ROOT / "provenance" / "alignment" / "{coverage}" / "replicate-{replicate}.json").as_posix()
        log:
            (LOG_ROOT / "low_depth" / "alignment.{coverage}.r{replicate}.log").as_posix()
        threads:
            int(LOW_DEPTH_E2E_CONFIG.get("alignment_threads", config["alignment"]["threads"]))
        resources:
            mem_mb=int(LOW_DEPTH_E2E_CONFIG.get("alignment_memory_mb", config["alignment"]["memory_mb"])),
            runtime_min=int(LOW_DEPTH_E2E_CONFIG.get("alignment_runtime_minutes", config["alignment"]["runtime_minutes"]))
        script:
            "../scripts/run_low_depth_alignment.py"

    rule low_depth_measure_coverage:
        input:
            bam=_ld_alignment,
            index=_ld_alignment + ".bai"
        output:
            json=_ld_coverage
        log:
            (LOG_ROOT / "low_depth" / "coverage.{coverage}.r{replicate}.log").as_posix()
        threads: 1
        resources:
            mem_mb=2000,
            runtime_min=int(LOW_DEPTH_E2E_CONFIG.get("coverage_runtime_minutes", 120))
        script:
            "../scripts/run_low_depth_coverage.py"

    rule low_depth_deepvariant:
        input:
            reference=config["reference"]["fasta"],
            fai=lambda wc: str(config["reference"]["fasta"]) + ".fai",
            bam=_ld_alignment,
            index=_ld_alignment + ".bai"
        output:
            vcf=_ld_small_vcf,
            gvcf=_ld_small_gvcf,
            vcf_index=_ld_small_vcf + ".tbi",
            gvcf_index=_ld_small_gvcf + ".tbi"
        log:
            (LOG_ROOT / "low_depth" / "deepvariant.{coverage}.r{replicate}.log").as_posix()
        threads:
            int(LOW_DEPTH_E2E_CONFIG.get("deepvariant_threads", config["small"]["threads"]))
        resources:
            mem_mb=int(LOW_DEPTH_E2E_CONFIG.get("deepvariant_memory_mb", config["small"]["memory_mb"])),
            runtime_min=int(LOW_DEPTH_E2E_CONFIG.get("deepvariant_runtime_minutes", config["small"]["runtime_minutes"])),
            deepvariant_slots=1
        params:
            intermediate=lambda wc: (WORK_ROOT / "low_depth" / LOW_DEPTH_CONFIG["study_id"] / "deepvariant" / wc.coverage / f"replicate-{wc.replicate}").as_posix(),
            tmpdir=lambda wc: (WORK_ROOT / "low_depth" / LOW_DEPTH_CONFIG["study_id"] / "deepvariant" / wc.coverage / f"replicate-{wc.replicate}" / "tmp").as_posix()
        script:
            "../scripts/run_low_depth_deepvariant.py"

    rule low_depth_happy:
        input:
            vcf=_ld_small_vcf,
            vcf_index=_ld_small_vcf + ".tbi",
            truth=LOW_DEPTH_E2E_CONFIG["truth_vcf"],
            truth_index=lambda wc: str(LOW_DEPTH_E2E_CONFIG["truth_vcf"]) + ".tbi",
            confident=LOW_DEPTH_E2E_CONFIG["confident_bed"],
            reference=config["reference"]["fasta"],
            fai=lambda wc: str(config["reference"]["fasta"]) + ".fai"
        output:
            summary=_ld_happy_prefix + ".summary.csv",
            metrics=_ld_happy_metrics,
            result=(LOW_DEPTH_ROOT / "benchmark" / "{coverage}" / "replicate-{replicate}.result.json").as_posix()
        log:
            (LOG_ROOT / "low_depth" / "happy.{coverage}.r{replicate}.log").as_posix()
        threads:
            int(LOW_DEPTH_E2E_CONFIG.get("happy_threads", 2))
        resources:
            mem_mb=int(LOW_DEPTH_E2E_CONFIG.get("happy_memory_mb", 16000)),
            runtime_min=int(LOW_DEPTH_E2E_CONFIG.get("happy_runtime_minutes", 480))
        params:
            prefix=_ld_happy_prefix
        script:
            "../scripts/run_low_depth_happy.py"

    rule low_depth_observation:
        input:
            coverage=_ld_coverage,
            summary=_ld_happy_prefix + ".summary.csv",
            benchmark=(LOW_DEPTH_ROOT / "benchmark" / "{coverage}" / "replicate-{replicate}.result.json").as_posix()
        output:
            json=_ld_observation
        params:
            run_id=lambda wc: LOW_DEPTH_BENCHMARK_RUNS[(wc.coverage, str(wc.replicate))].run_id,
            target_coverage=lambda wc: LOW_DEPTH_BENCHMARK_RUNS[(wc.coverage, str(wc.replicate))].target_coverage
        script:
            "../scripts/write_low_depth_observation.py"
