"""Optional Phase 13 truth-set benchmark tracks; never upstream dependencies."""

from pathlib import Path
from hifivar.exceptions import WorkflowError

BENCHMARK_CONFIG = config.get("benchmark", {})
if not isinstance(BENCHMARK_CONFIG, dict):
    raise WorkflowError("Effective config section benchmark must be a mapping.")
BENCHMARK_ENABLED = BENCHMARK_CONFIG.get("enabled", False)
BENCHMARK_ID = BENCHMARK_CONFIG.get("benchmark_id") or "_disabled"
BENCHMARK_ROOT = OUTPUT_ROOT / "benchmark" / str(BENCHMARK_ID)
BENCHMARK_RESULTS = []

if BENCHMARK_ENABLED:
    for _track in ("small_variants", "sv", "assembly_sv", "tr"):
        if BENCHMARK_CONFIG.get(_track, {}).get("enabled", False):
            _result = (BENCHMARK_ROOT / _track / "benchmark_result.json").as_posix()
            BENCHMARK_RESULTS.append(_result)
            WORKFLOW_TARGETS.append(_result)
    if BENCHMARK_RESULTS:
        BENCHMARK_MANIFEST = (BENCHMARK_ROOT / "benchmark_manifest.json").as_posix()
        BENCHMARK_METRICS = (BENCHMARK_ROOT / "benchmark_metrics.tsv").as_posix()
        WORKFLOW_TARGETS += [BENCHMARK_MANIFEST, BENCHMARK_METRICS]

if BENCHMARK_ENABLED and BENCHMARK_CONFIG.get("small_variants", {}).get("enabled", False):
    _small = BENCHMARK_CONFIG["small_variants"]
    _happy_metrics_suffix = ".metrics.json.gz" if _small.get("metrics_compression", "gzip") == "gzip" else ".metrics.json"
    rule benchmark_small_variants:
        input:
            query=_small["query_vcf"], query_index=lambda wc: str(_small["query_vcf"]) + ".tbi",
            truth=_small["truth_vcf"], truth_index=lambda wc: str(_small["truth_vcf"]) + ".tbi",
            confident=_small["confident_bed"], reference=config["reference"]["fasta"],
            fai=lambda wc: str(config["reference"]["fasta"]) + ".fai",
            stratifications=lambda wc: [item["path"] for item in _small.get("stratifications", [])]
        output:
            summary=(BENCHMARK_ROOT / "small_variants" / "happy.summary.csv").as_posix(),
            metrics=(BENCHMARK_ROOT / "small_variants" / f"happy{_happy_metrics_suffix}").as_posix(),
            result=(BENCHMARK_ROOT / "small_variants" / "benchmark_result.json").as_posix()
        log: (LOG_ROOT / "benchmark" / str(BENCHMARK_ID) / "happy.log").as_posix()
        threads: int(_small.get("threads", 8))
        resources:
            mem_mb=int(_small.get("memory_mb", 16000)), runtime_min=int(_small.get("runtime_minutes", 480))
        script: "../scripts/run_benchmark_small.py"

if BENCHMARK_ENABLED and BENCHMARK_CONFIG.get("sv", {}).get("enabled", False):
    _svb = BENCHMARK_CONFIG["sv"]
    rule benchmark_sv:
        input:
            query=_svb["query_vcf"], query_index=lambda wc: str(_svb["query_vcf"]) + ".tbi",
            truth=_svb["truth_vcf"], truth_index=lambda wc: str(_svb["truth_vcf"]) + ".tbi",
            reference=config["reference"]["fasta"], fai=lambda wc: str(config["reference"]["fasta"]) + ".fai",
            confident=lambda wc: [_svb["confident_bed"]] if _svb.get("confident_bed") else []
        output:
            directory((BENCHMARK_ROOT / "sv" / "truvari").as_posix()),
            result=(BENCHMARK_ROOT / "sv" / "benchmark_result.json").as_posix()
        log: (LOG_ROOT / "benchmark" / str(BENCHMARK_ID) / "truvari-sv.log").as_posix()
        threads: int(_svb.get("threads", 1))
        resources:
            mem_mb=int(_svb.get("memory_mb", 16000)), runtime_min=int(_svb.get("runtime_minutes", 480))
        script: "../scripts/run_benchmark_sv.py"

if BENCHMARK_ENABLED and BENCHMARK_CONFIG.get("assembly_sv", {}).get("enabled", False):
    _asvb = BENCHMARK_CONFIG["assembly_sv"]
    rule benchmark_assembly_sv:
        input:
            query=_asvb["query_vcf"], query_index=lambda wc: str(_asvb["query_vcf"]) + ".tbi",
            truth=_asvb["truth_vcf"], truth_index=lambda wc: str(_asvb["truth_vcf"]) + ".tbi",
            reference=config["reference"]["fasta"], fai=lambda wc: str(config["reference"]["fasta"]) + ".fai",
            confident=lambda wc: [_asvb["confident_bed"]] if _asvb.get("confident_bed") else []
        output:
            directory((BENCHMARK_ROOT / "assembly_sv" / "truvari").as_posix()),
            result=(BENCHMARK_ROOT / "assembly_sv" / "benchmark_result.json").as_posix()
        log: (LOG_ROOT / "benchmark" / str(BENCHMARK_ID) / "truvari-assembly-sv.log").as_posix()
        threads: int(_asvb.get("threads", 1))
        resources:
            mem_mb=int(_asvb.get("memory_mb", 16000)), runtime_min=int(_asvb.get("runtime_minutes", 480))
        params: track="assembly_sv"
        script: "../scripts/run_benchmark_sv.py"

if BENCHMARK_ENABLED and BENCHMARK_CONFIG.get("tr", {}).get("enabled", False):
    _trb = BENCHMARK_CONFIG["tr"]
    rule benchmark_tr:
        input: query=_trb["query_vcf"], truth=_trb["truth_vcf"]
        output:
            loci=(BENCHMARK_ROOT / "tr" / "locus_comparison.tsv").as_posix(),
            result=(BENCHMARK_ROOT / "tr" / "benchmark_result.json").as_posix()
        threads: 1
        resources:
            mem_mb=int(_trb.get("memory_mb", 4000)), runtime_min=int(_trb.get("runtime_minutes", 120))
        script: "../scripts/run_benchmark_tr.py"

if BENCHMARK_ENABLED and BENCHMARK_RESULTS:
    rule benchmark_manifest:
        input: BENCHMARK_RESULTS
        output: json=BENCHMARK_MANIFEST, tsv=BENCHMARK_METRICS
        script: "../scripts/run_benchmark_manifest.py"
