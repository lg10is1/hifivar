"""Cross-platform rules proving the Phase 0 workflow DAG."""


rule phase0_prepare:
    output:
        PHASE0_PREPARE_MARKER
    log:
        PHASE0_PREPARE_LOG
    threads:
        RUNTIME_THREADS
    resources:
        mem_mb=DEFAULT_MEM_MB,
        runtime_min=DEFAULT_RUNTIME_MIN
    params:
        marker_text="HiFiVar effective config accepted\n",
        rule_name="phase0_prepare"
    script:
        "../scripts/write_phase0_marker.py"


rule phase0_smoke:
    input:
        PHASE0_PREPARE_MARKER
    output:
        PHASE0_SMOKE_MARKER
    log:
        PHASE0_SMOKE_LOG
    threads:
        RUNTIME_THREADS
    resources:
        mem_mb=DEFAULT_MEM_MB,
        runtime_min=DEFAULT_RUNTIME_MIN
    params:
        marker_text=(
            "HiFiVar Snakemake infrastructure OK\n"
            f"project: {PROJECT_NAME}\n"
            f"preset: {WORKFLOW_PRESET}\n"
        ),
        rule_name="phase0_smoke"
    script:
        "../scripts/write_phase0_marker.py"
