"""Phase 13 benchmark contracts and streaming truth-set comparisons.

Benchmark evidence is descriptive.  It never changes a caller VCF and never
turns functional or clinical interpretation into a calling-confidence claim.
"""

from __future__ import annotations

import csv
import gzip
import sqlite3
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Iterator, TextIO

from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.serialization import utc_now_iso8601, write_json_atomic, write_yaml_atomic


class BenchmarkStatus(str, Enum):
    PASS = "PASS"
    PARTIAL = "PARTIAL"
    NOT_RUN = "NOT_RUN"
    UNSUPPORTED = "UNSUPPORTED"


class BenchmarkVariantClass(str, Enum):
    SMALL_VARIANT = "small_variant"
    SV = "sv"
    TR = "tr"
    ASSEMBLY_SV = "assembly_sv"


@dataclass(frozen=True, slots=True)
class TruthSet:
    path: Path
    version: str
    reference_build: str
    source: str
    catalog_id: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "path", Path(self.path))
        for name in ("version", "reference_build", "source"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise InputValidationError(f"Truth-set {name} must be explicit and non-empty.")
        if self.version.strip().casefold() == "latest":
            raise InputValidationError("Truth-set version must be pinned; 'latest' is not reproducible.")

    def to_dict(self) -> dict[str, object]:
        return {"path": str(self.path), "version": self.version,
                "reference_build": self.reference_build, "source": self.source,
                "catalog_id": self.catalog_id}


@dataclass(frozen=True, slots=True)
class BenchmarkRegion:
    name: str
    path: Path
    version: str
    region_class: str = "confident"

    def __post_init__(self) -> None:
        object.__setattr__(self, "path", Path(self.path))
        for name in ("name", "version", "region_class"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise InputValidationError(f"Benchmark region {name} must be non-empty.")
        if self.version.strip().casefold() == "latest":
            raise InputValidationError("Benchmark region version must be pinned, not 'latest'.")

    def to_dict(self) -> dict[str, str]:
        return {"name": self.name, "path": str(self.path), "version": self.version,
                "region_class": self.region_class}


@dataclass(frozen=True, slots=True)
class BenchmarkMetric:
    name: str
    value: float | int | None
    variant_class: BenchmarkVariantClass
    stratum: str = "ALL"
    source_field: str | None = None

    def to_dict(self) -> dict[str, object]:
        return {"name": self.name, "value": self.value,
                "variant_class": self.variant_class.value, "stratum": self.stratum,
                "source_field": self.source_field}


@dataclass(frozen=True, slots=True)
class BenchmarkResult:
    benchmark_id: str
    sample_id: str
    variant_class: BenchmarkVariantClass
    status: BenchmarkStatus
    query_path: Path
    truth_set: TruthSet
    tool: str
    tool_version: str | None
    outputs: tuple[Path, ...] = ()
    metrics: tuple[BenchmarkMetric, ...] = ()
    command: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()
    regions: tuple[BenchmarkRegion, ...] = ()
    effective_config: dict[str, object] = field(default_factory=dict)

    def to_dict(self) -> dict[str, object]:
        return {
            "benchmark_id": self.benchmark_id, "sample_id": self.sample_id,
            "variant_class": self.variant_class.value, "status": self.status.value,
            "query_path": str(self.query_path), "truth_set": self.truth_set.to_dict(),
            "tool": self.tool, "tool_version": self.tool_version,
            "outputs": [str(path) for path in self.outputs],
            "metrics": [metric.to_dict() for metric in self.metrics],
            "command": list(self.command), "notes": list(self.notes),
            "regions": [region.to_dict() for region in self.regions],
            "effective_config": self.effective_config,
            "scientific_semantics": {
                "truth_set_is_context_bound": True,
                "benchmark_performance_is_not_pathogenicity": True,
                "functional_impact_is_not_call_confidence": True,
                "f1_is_not_clinical_utility": True,
                "caller_count_is_not_benchmark_truth": True,
                "benchmark_concordance_is_not_absolute_truth": True,
                "raw_query_modified": False,
            },
        }


@dataclass(frozen=True, slots=True)
class BenchmarkManifest:
    benchmark_id: str
    reference_build: str
    results: tuple[BenchmarkResult, ...]
    config: dict[str, object] = field(default_factory=dict)
    created_at: str = field(default_factory=utc_now_iso8601)
    git_commit: str | None = None

    def to_dict(self) -> dict[str, object]:
        return {"benchmark_id": self.benchmark_id, "reference_build": self.reference_build,
                "created_at": self.created_at, "git_commit": self.git_commit,
                "config": self.config, "results": [result.to_dict() for result in self.results]}

    def write(self, root: Path, *, overwrite: bool = False, markdown: bool = True) -> tuple[Path, ...]:
        root = Path(root)
        json_path, yaml_path, tsv_path = root / "benchmark_manifest.json", root / "benchmark_manifest.yaml", root / "benchmark_metrics.tsv"
        write_json_atomic(self.to_dict(), json_path, overwrite=overwrite, artifact_name="Benchmark manifest")
        write_yaml_atomic(self.to_dict(), yaml_path, overwrite=overwrite, artifact_name="Benchmark manifest")
        _write_metrics_tsv(tsv_path, self.results, overwrite=overwrite)
        outputs: list[Path] = [json_path, yaml_path, tsv_path]
        if markdown:
            md_path = root / "benchmark_summary.md"
            _write_markdown(md_path, self, overwrite=overwrite)
            outputs.append(md_path)
        return tuple(outputs)


def compare_tr_vcfs(
    *, benchmark_id: str, sample_id: str, query_vcf: Path, truth_set: TruthSet,
    query_catalog_id: str, output_tsv: Path, overwrite: bool = False,
) -> BenchmarkResult:
    """Compare TRGT-like single-sample VCFs by explicit TRID without rematching loci."""
    if not truth_set.catalog_id or truth_set.catalog_id != query_catalog_id:
        raise InputValidationError("TR benchmark requires identical explicit truth/query catalog IDs.")
    query_vcf, output_tsv = Path(query_vcf), Path(output_tsv)
    for path in (query_vcf, truth_set.path):
        if not path.is_file():
            raise InputValidationError(f"TR benchmark input does not exist: '{path}'.")
    if output_tsv.exists() and not overwrite:
        raise OutputValidationError(f"TR benchmark output exists: '{output_tsv}'.")
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    scratch = output_tsv.with_name(f".{output_tsv.name}.benchmark.sqlite")
    if scratch.exists():
        raise OutputValidationError(f"TR benchmark scratch path exists: '{scratch}'.")
    connection=sqlite3.connect(scratch)
    compared = exact_allele = genotype_equal = no_call = total = 0
    try:
        connection.execute("CREATE TABLE calls (source TEXT, trid TEXT, alleles TEXT, genotype TEXT, PRIMARY KEY(source,trid))")
        for source,path in (("truth",truth_set.path),("query",query_vcf)):
            for trid,alleles,genotype in _tr_records(path,sample_id):
                try: connection.execute("INSERT INTO calls VALUES (?,?,?,?)",(source,trid,alleles,genotype))
                except sqlite3.IntegrityError as error: raise OutputValidationError(f"Duplicate TRID '{trid}' in '{path}'.") from error
        connection.commit()
        mode = "w" if overwrite else "x"
        with output_tsv.open(mode, encoding="utf-8", newline="\n") as handle:
            handle.write("trid\ttruth_alleles\tquery_alleles\ttruth_gt\tquery_gt\texact_allele_agreement\tgenotype_agreement\tquery_no_call\n")
            rows=connection.execute("SELECT t.trid,t.alleles,t.genotype,q.alleles,q.genotype FROM calls t LEFT JOIN calls q ON q.source='query' AND q.trid=t.trid WHERE t.source='truth' ORDER BY t.trid")
            for trid,truth_al,truth_gt,query_al,query_gt in rows:
                total += 1; query_al=query_al or ""; query_gt=query_gt or ""
                missing = query_gt in {"", ".", "./.", ".|."}
                no_call += int(missing)
                if not missing:
                    compared += 1; exact_allele += int(bool(truth_al) and truth_al == query_al); genotype_equal += int(bool(truth_gt) and truth_gt == query_gt)
                handle.write("\t".join((trid, truth_al, query_al, truth_gt, query_gt,
                    str(bool(truth_al) and truth_al == query_al).lower(), str(bool(truth_gt) and truth_gt == query_gt).lower(), str(missing).lower())) + "\n")
    finally:
        connection.close()
        scratch.unlink(missing_ok=True)
    metrics = (
        BenchmarkMetric("truth_locus_count", total, BenchmarkVariantClass.TR),
        BenchmarkMetric("compared_locus_count", compared, BenchmarkVariantClass.TR),
        BenchmarkMetric("exact_allele_agreement", exact_allele / compared if compared else None, BenchmarkVariantClass.TR),
        BenchmarkMetric("genotype_agreement", genotype_equal / compared if compared else None, BenchmarkVariantClass.TR),
        BenchmarkMetric("query_no_call_rate", no_call / total if total else None, BenchmarkVariantClass.TR),
    )
    status = BenchmarkStatus.PASS if total else BenchmarkStatus.UNSUPPORTED
    return BenchmarkResult(benchmark_id, sample_id, BenchmarkVariantClass.TR, status,
                           query_vcf, truth_set, "hifivar-tr-exact", "1", (output_tsv,), metrics,
                           notes=("Exact catalog/TRID comparison; no clinical interpretation.",))


def _tr_records(path: Path, sample_id: str) -> Iterator[tuple[str, str, str]]:
    sample_column: int | None = None
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("#CHROM\t"):
                samples = line.rstrip().split("\t")[9:]
                if sample_id not in samples:
                    raise InputValidationError(f"TR VCF '{path}' lacks sample '{sample_id}'.")
                sample_column = 9 + samples.index(sample_id)
            elif not line.startswith("#") and line.strip():
                if sample_column is None:
                    raise InputValidationError(f"TR VCF '{path}' has no sample header.")
                fields = line.rstrip().split("\t")
                info = {item.partition("=")[0]: item.partition("=")[2] for item in fields[7].split(";")}
                trid = info.get("TRID") or fields[2]
                if not trid or trid == ".": raise InputValidationError(f"TR VCF '{path}' record lacks TRID/ID.")
                keys, values = fields[8].split(":"), fields[sample_column].split(":")
                mapping = dict(zip(keys, values))
                alleles = mapping.get("AL", mapping.get("MC", ""))
                yield trid, alleles, mapping.get("GT", "")


def _open_text(path: Path) -> TextIO:
    return gzip.open(path, "rt", encoding="utf-8", newline="") if str(path).endswith(".gz") else path.open("r", encoding="utf-8", newline="")


def _write_metrics_tsv(path: Path, results: tuple[BenchmarkResult, ...], *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise OutputValidationError(f"Benchmark metrics output exists: '{path}'.")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w" if overwrite else "x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("benchmark_id", "sample", "variant_class", "status", "tool", "tool_version", "stratum", "metric", "value"))
        for result in results:
            if not result.metrics:
                writer.writerow((result.benchmark_id, result.sample_id, result.variant_class.value, result.status.value, result.tool, result.tool_version or "", "", "", ""))
            for metric in result.metrics:
                writer.writerow((result.benchmark_id, result.sample_id, result.variant_class.value, result.status.value, result.tool, result.tool_version or "", metric.stratum, metric.name, "" if metric.value is None else metric.value))


def _write_markdown(path: Path, manifest: BenchmarkManifest, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise OutputValidationError(f"Benchmark summary output exists: '{path}'.")
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"# Benchmark {manifest.benchmark_id}", "", "Benchmark metrics are truth-set/context specific and are not pathogenicity claims.", "", "| Track | Sample | Status | Tool |", "|---|---|---|---|"]
    lines.extend(f"| {r.variant_class.value} | {r.sample_id} | {r.status.value} | {r.tool} {r.tool_version or ''} |" for r in manifest.results)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


__all__ = ["BenchmarkManifest", "BenchmarkMetric", "BenchmarkRegion", "BenchmarkResult",
           "BenchmarkStatus", "BenchmarkVariantClass", "TruthSet", "compare_tr_vcfs"]
