"""Explicit, sequential benchmarks of existing immutable SV callsets.

This entry runs no caller and performs no chromosome or allele repair. Initial
bench and refined-variant metrics remain distinct, context-bound artifacts.
"""

from __future__ import annotations

from dataclasses import dataclass, fields, replace
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import time
from typing import Any

import yaml

from hifivar import __version__
from hifivar.command import CommandRunner, format_command
from hifivar.exceptions import (
    CommandExecutionError, HiFiVarError, InputValidationError,
    OutputValidationError, ToolVersionError,
)
from hifivar.reference import ReferenceGenome, normalize_reference_build
from hifivar.serialization import utc_now_iso8601, write_json_atomic
from hifivar.truvari import (
    TruvariRequest, TruvariThresholds, TruvariWrapper,
    parse_truvari_refined_summary, parse_truvari_summary,
)
from hifivar.validation import compute_sha256, validate_vcf


_IDENTIFIER = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")
_VERSION = re.compile(r"\d+(?:\.\d+)+(?:[A-Za-z0-9.+-]*)")
_TOP_KEYS = {"schema_version", "benchmark_id", "sample_id", "reference", "truth", "tools", "thresholds", "queries"}


class SourceImmutabilityError(HiFiVarError):
    """An input changed; no further query may be benchmarked in this run."""


class _UniqueLoader(yaml.SafeLoader):
    """Reject ambiguous YAML duplicate keys rather than taking the last one."""


def _unique_mapping(loader: yaml.SafeLoader, node: yaml.MappingNode) -> dict:
    result: dict = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=True)
        if not isinstance(key, str) or key in result:
            raise InputValidationError("Manifest mapping keys must be unique strings.")
        result[key] = loader.construct_object(value_node, deep=True)
    return result


_UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _unique_mapping)


@dataclass(frozen=True, slots=True)
class SVBenchmarkQuery:
    """One named existing callset; caller labels are provenance, not confidence."""

    id: str
    vcf: Path
    caller: str


@dataclass(frozen=True, slots=True)
class SVBenchmarkJob:
    """Validated schema-1 job whose relative paths resolve beside its manifest."""

    manifest: Path
    benchmark_id: str
    sample_id: str
    reference: ReferenceGenome
    truth_vcf: Path
    confident_bed: Path
    truth_version: str
    truth_source: str
    truvari_executable: str
    bcftools_executable: str
    truvari_version: str
    thresholds: TruvariThresholds
    queries: tuple[SVBenchmarkQuery, ...]
    manifest_sha256: str

    @property
    def source_paths(self) -> tuple[Path, ...]:
        """All data/index/config inputs, without loading large data in memory."""
        paths = [self.manifest, self.reference.fasta, self.reference.fai,
                 self.truth_vcf, Path(f"{self.truth_vcf}.tbi"), self.confident_bed]
        for query in self.queries:
            paths.extend((query.vcf, Path(f"{query.vcf}.tbi")))
        return tuple(dict.fromkeys(paths))


def _mapping(value: object, keys: set[str], label: str) -> dict:
    if not isinstance(value, dict) or set(value) != keys:
        actual = set(value) if isinstance(value, dict) else set()
        raise InputValidationError(
            f"{label} requires exactly {sorted(keys)}; missing={sorted(keys - actual)}, "
            f"unknown={sorted(actual - keys)}."
        )
    return value


def _text(value: object, label: str, *, pinned: bool = False) -> str:
    if not isinstance(value, str) or not value.strip() or value != value.strip():
        raise InputValidationError(f"{label} must be a non-empty, trimmed string.")
    if "\x00" in value or "\n" in value or "\r" in value:
        raise InputValidationError(f"{label} contains an invalid control character.")
    if pinned and value.casefold() in {"latest", "unknown", "tbd", "unversioned"}:
        raise InputValidationError(f"{label} must be a pinned version, not {value!r}.")
    return value


def _identifier(value: object, label: str) -> str:
    text = _text(value, label)
    if not _IDENTIFIER.fullmatch(text) or text in {".", ".."}:
        raise InputValidationError(f"{label} is not a safe single-component identifier.")
    return text


def _path(value: object, parent: Path, label: str, *, existing: bool = True) -> Path:
    text = _text(value, label)
    if ".." in PurePosixPath(text.replace("\\", "/")).parts or text.startswith("~"):
        raise InputValidationError(f"{label} must not contain traversal or implicit home expansion.")
    path = Path(text)
    path = (path if path.is_absolute() else parent / path).resolve()
    if existing and (not path.is_file() or path.stat().st_size == 0):
        raise InputValidationError(f"{label} is missing, empty, or not a file: '{path}'.")
    return path


def _tool(value: object, parent: Path, label: str) -> str:
    text = _text(value, label)
    if "/" in text or "\\" in text:
        return str(_path(text, parent, label, existing=False))
    if text.startswith("-") or any(c.isspace() for c in text):
        raise InputValidationError(f"{label} must be one executable, not a shell command.")
    return text


def load_sv_benchmark_job(manifest: str | Path) -> SVBenchmarkJob:
    """Load a strict schema-1 manifest; never create outputs or run a tool."""
    source = _path(str(manifest), Path.cwd(), "manifest")
    try:
        before = _stat(source)
        raw_manifest = source.read_bytes()
        if _stat(source) != before:
            raise SourceImmutabilityError("Manifest changed while being read.")
        payload = yaml.load(raw_manifest.decode("utf-8-sig"), Loader=_UniqueLoader)
    except (OSError, UnicodeError, yaml.YAMLError) as error:
        raise InputValidationError(f"Unable to read SV benchmark manifest: {error}") from error
    data = _mapping(payload, _TOP_KEYS, "SV benchmark manifest")
    if type(data["schema_version"]) is not int or data["schema_version"] != 1:
        raise InputValidationError("SV benchmark schema_version must be integer 1.")
    ref = _mapping(data["reference"], {"fasta", "build"}, "reference")
    truth = _mapping(data["truth"], {"vcf", "confident_bed", "version", "source", "reference_build"}, "truth")
    tools = _mapping(data["tools"], {"truvari_executable", "bcftools_executable", "truvari_version"}, "tools")
    reference = ReferenceGenome.from_fasta(
        _path(ref["fasta"], source.parent, "reference.fasta"),
        build=_text(ref["build"], "reference.build", pinned=True),
    )
    truth_build = normalize_reference_build(_text(truth["reference_build"], "truth.reference_build", pinned=True))
    if truth_build != reference.build:
        raise InputValidationError("Truth/reference build mismatch; chromosome renaming cannot change genome build.")
    thresholds = data["thresholds"]
    known = {field.name for field in fields(TruvariThresholds)}
    if not isinstance(thresholds, dict) or set(thresholds) - known:
        raise InputValidationError("thresholds must contain only TruvariThresholds fields.")
    policy = TruvariThresholds(**thresholds)
    if any(isinstance(v, float) and not math.isfinite(v) for v in thresholds.values()):
        raise InputValidationError("Truvari thresholds must be finite.")
    pinned_version = _text(tools["truvari_version"], "tools.truvari_version", pinned=True)
    if not _VERSION.fullmatch(pinned_version):
        raise InputValidationError("tools.truvari_version must be a specific numeric release.")
    truth_vcf = _path(truth["vcf"], source.parent, "truth.vcf")
    queries: list[SVBenchmarkQuery] = []
    ids: set[str] = set()
    identities = {_file_identity(truth_vcf)}
    if not isinstance(data["queries"], list) or not data["queries"]:
        raise InputValidationError("queries must be a non-empty list.")
    for raw in data["queries"]:
        item = _mapping(raw, {"id", "vcf", "caller"}, "query")
        query_id = _identifier(item["id"], "query.id")
        if query_id.casefold() in ids:
            raise InputValidationError(f"Duplicate query id: '{query_id}'.")
        ids.add(query_id.casefold())
        vcf = _path(item["vcf"], source.parent, f"query {query_id}.vcf")
        identity = _file_identity(vcf)
        if identity in identities:
            raise InputValidationError("Truth/query and query/query source aliases are not allowed.")
        identities.add(identity)
        queries.append(SVBenchmarkQuery(query_id, vcf, _text(item["caller"], "query.caller")))
    return SVBenchmarkJob(
        source, _identifier(data["benchmark_id"], "benchmark_id"),
        _identifier(data["sample_id"], "sample_id"), reference, truth_vcf,
        _path(truth["confident_bed"], source.parent, "truth.confident_bed"),
        _text(truth["version"], "truth.version", pinned=True),
        _text(truth["source"], "truth.source"),
        _tool(tools["truvari_executable"], source.parent, "tools.truvari_executable"),
        _tool(tools["bcftools_executable"], source.parent, "tools.bcftools_executable"),
        pinned_version, policy, tuple(queries), hashlib.sha256(raw_manifest).hexdigest(),
    )


def _file_identity(path: Path) -> tuple[int, int]:
    stat = path.stat()
    return stat.st_dev, stat.st_ino


def _stat(path: Path) -> dict[str, int]:
    try:
        stat = path.stat()
    except OSError as error:
        raise SourceImmutabilityError(f"Source disappeared/unreadable: '{path}': {error}") from error
    return {"device": stat.st_dev, "inode": stat.st_ino, "bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns, "ctime_ns": stat.st_ctime_ns}


def _snapshot(paths: tuple[Path, ...]) -> dict[str, dict[str, Any]]:
    snapshot: dict[str, dict[str, Any]] = {}
    for path in paths:
        before = _stat(path)
        digest = compute_sha256(path)
        if _stat(path) != before:
            raise SourceImmutabilityError(f"Source changed while hashing: '{path}'.")
        snapshot[str(path)] = {"stat": before, "sha256": digest}
    return snapshot


def _assert_unchanged(snapshot: dict[str, dict[str, Any]]) -> None:
    for name, item in snapshot.items():
        if _stat(Path(name)) != item["stat"]:
            raise SourceImmutabilityError(f"Source mutation detected; remaining queries stopped: '{name}'.")


def _output_root(output_dir: str | Path, paths: tuple[Path, ...]) -> Path:
    root = _path(str(output_dir), Path.cwd(), "output directory", existing=False)
    if root.exists() or Path(output_dir).is_symlink():
        raise OutputValidationError(f"SV benchmark requires a fresh output directory: '{root}'.")
    if root == Path(root.anchor):
        raise OutputValidationError("A filesystem root cannot be a benchmark output directory.")
    for path in paths:
        resolved = path.resolve()
        if root == resolved or root in resolved.parents:
            raise OutputValidationError("Benchmark output would alias or contain an immutable input.")
    return root


def _ref_command(job: SVBenchmarkJob, vcf: Path) -> tuple[str, ...]:
    return (job.bcftools_executable, "norm", "-N", "-c", "e", "-f",
            str(job.reference.fasta), "-Ou", "-o", os.devnull, str(vcf))


class _RecordingRunner:
    """Add owned logs and structured exit evidence around shared CommandRunner."""

    def __init__(self, runner: CommandRunner, root: Path, records: list[dict], label: str):
        self.runner, self.root, self.records, self.label = runner, root, records, label

    def require_executable(self, executable: str) -> Path:
        return self.runner.require_executable(executable)

    def run(self, command, **kwargs):
        prefix = self.root / "logs" / f"{self.label}.{len(self.records):03d}"
        kwargs.setdefault("cwd", self.root)
        if kwargs.get("stdout_path") is None:
            kwargs["stdout_path"] = Path(f"{prefix}.stdout.log")
        if kwargs.get("stderr_path") is None:
            kwargs["stderr_path"] = Path(f"{prefix}.stderr.log")
        if tuple(map(str, command[1:])) in {("version",), ("--version",)}:
            kwargs.setdefault("timeout", 30)
        check = kwargs.pop("check", True)
        record = {"command": list(map(str, command)), "display_command": format_command(command),
                  "stdout_path": str(kwargs["stdout_path"]), "stderr_path": str(kwargs["stderr_path"]),
                  "cwd": str(kwargs["cwd"]), "return_code": None, "executed": False}
        self.records.append(record)
        started = time.monotonic()
        try:
            result = self.runner.run(command, check=False, **kwargs)
        except HiFiVarError as error:
            record.update(error_type=type(error).__name__, error=str(error), duration_seconds=time.monotonic() - started)
            raise
        record.update(return_code=result.returncode, executed=result.executed, duration_seconds=result.duration_seconds)
        if check and result.returncode not in {0, None}:
            raise CommandExecutionError(
                f"Benchmark command failed with return code {result.returncode}; stderr: {kwargs['stderr_path']}"
            )
        if tuple(map(str, command[1:])) in {("version",), ("--version",)}:
            # Existing tool adapters parse captured text. Keep bounded version
            # readback while retaining the complete version logs on disk.
            return replace(result, stdout=_logged_text(result), stderr=None)
        return result


def _logged_text(result) -> str:
    parts = [item for item in (result.stdout, result.stderr) if item]
    for path in (result.stdout_path, result.stderr_path):
        if path is not None and Path(path).is_file():
            with Path(path).open("r", encoding="utf-8", errors="replace") as stream:
                parts.append(stream.read(8192))
    return "\n".join(parts)


def _metrics(path: Path, *, refined: bool = False) -> list[dict]:
    if not refined:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as error:
            raise OutputValidationError(f"Unreadable native Truvari summary: {error}") from error
        required = {"TP-base", "TP-comp", "FP", "FN", "base cnt", "comp cnt", "precision", "recall", "f1"}
        if not isinstance(payload, dict) or not required.issubset(payload):
            raise OutputValidationError("New batch executions require native Truvari summary counts and denominators.")
    parsed = parse_truvari_refined_summary(path) if refined else parse_truvari_summary(path)
    result = [metric.to_dict() for metric in parsed]
    for item in result:
        value = item["value"]
        if item["name"] in {"tp_base", "tp_call", "fp", "fn"}:
            if type(value) is not int or value < 0:
                raise OutputValidationError("Truvari counts must be non-negative integers.")
        elif value is not None and (isinstance(value, bool) or not math.isfinite(value) or not 0 <= value <= 1):
            raise OutputValidationError("Truvari rates must be finite values in [0, 1] or null.")
    return result


def run_sv_benchmark(
    manifest: str | Path, output_dir: str | Path, *, dry_run: bool = False,
    runner: CommandRunner | None = None,
) -> dict[str, Any]:
    """Evaluate existing queries sequentially, persisting failures without retry.

    Dry-run reads metadata/headers only and creates no files, directories or
    executable probes. It never marks sequence-level REF compatibility PASS.
    Any source mutation is a global failure, not an ordinary per-query failure.
    """
    from hifivar.benchmark_inputs import audit_benchmark_regions, audit_sv_benchmark_input

    if not isinstance(dry_run, bool):
        raise InputValidationError("dry_run must be boolean.")
    job = load_sv_benchmark_job(manifest)
    root = _output_root(output_dir, job.source_paths)
    truth_audit = audit_sv_benchmark_input(job.truth_vcf, job.sample_id, job.reference)
    regions_audit = audit_benchmark_regions(job.confident_bed, job.reference)
    records: list[dict] = []
    report: dict[str, Any] = {
        "schema_version": 1, "benchmark_id": job.benchmark_id, "sample_id": job.sample_id,
        "hifivar_version": __version__, "created_at": utc_now_iso8601(),
        "manifest": str(job.manifest), "output_directory": str(root),
        "status": "PLANNED" if dry_run else "RUNNING", "dry_run": dry_run,
        "reference": job.reference.to_dict(),
        "truth": {"vcf": str(job.truth_vcf), "version": job.truth_version,
                  "source": job.truth_source, "reference_build": job.reference.build,
                  "confident_bed": str(job.confident_bed)},
        "thresholds": job.thresholds.to_dict(), "truth_audit": truth_audit,
        "regions_audit": regions_audit, "commands": records, "queries": [],
        "tools": {"truvari_executable": job.truvari_executable,
                  "bcftools_executable": job.bcftools_executable,
                  "expected_truvari_version": job.truvari_version},
        "ref_sequence_check": "NOT_RUN_DRY_RUN" if dry_run else "NOT_RUN",
        "source_immutability": "NOT_CHECKED_DRY_RUN" if dry_run else "PENDING",
        "scientific_scope": "Declared truth/BED only; no caller rerun, confidence, pathogenicity or cross-sample inference.",
    }
    for query in job.queries:
        request = TruvariRequest(job.sample_id, job.reference, job.truth_vcf, query.vcf,
                                 root / "queries" / query.id, confident_regions=job.confident_bed,
                                 thresholds=job.thresholds)
        report["queries"].append({"id": query.id, "caller": query.caller, "vcf": str(query.vcf),
                                  "status": "NOT_RUN", "command": list(TruvariWrapper(executable=job.truvari_executable).plan_command(request)),
                                  "ref_sequence_check": "NOT_RUN_DRY_RUN" if dry_run else "NOT_RUN",
                                  "bench_metrics": [], "refined_metrics": [],
                                  "refined_stratification_status": "NOT_RUN_NEEDS_REFINED_ASSIGNMENTS" if job.thresholds.refine else "NOT_REQUESTED"})
    if dry_run:
        for query, entry in zip(job.queries, report["queries"]):
            try:
                entry["input_audit"] = audit_sv_benchmark_input(query.vcf, job.sample_id, job.reference)
                entry["status"] = "PLANNED"
                entry["ref_check_command"] = list(_ref_command(job, query.vcf))
            except HiFiVarError as error:
                entry.update(status="FAILED", error_type=type(error).__name__, error=str(error))
        if any(entry["status"] == "FAILED" for entry in report["queries"]):
            report["status"] = "FAILED"
        return report

    underlying = runner or CommandRunner()
    tool_paths = tuple(underlying.require_executable(exe).resolve() for exe in (job.truvari_executable, job.bcftools_executable))
    _output_root(root, job.source_paths + tool_paths)
    for executable in tool_paths:
        installation = executable.parent.parent if executable.parent.name.casefold() in {"bin", "scripts"} else executable.parent
        if installation == root or installation in root.parents:
            raise OutputValidationError("Benchmark outputs must not be created inside a tool installation.")
    snapshot = _snapshot(tuple(dict.fromkeys(job.source_paths + tool_paths)))
    if snapshot[str(job.manifest)]["sha256"] != job.manifest_sha256:
        raise SourceImmutabilityError("Manifest changed between parsing and input hashing.")
    report["input_checksums"] = snapshot
    report["reference"]["sha256"] = snapshot[str(job.reference.fasta)]["sha256"]
    root.mkdir(parents=True, exist_ok=False)
    matrix_path = root / "benchmark_matrix.json"

    def persist() -> None:
        write_json_atomic(report, matrix_path, overwrite=True, artifact_name="SV benchmark matrix")

    persist()
    try:
        shared = _RecordingRunner(underlying, root, records, "shared")
        # Pin matching version before starting comparisons. Keep all probe logs.
        version_result = shared.run([job.truvari_executable, "version"])
        match = re.search(r"\d+(?:\.\d+)+", _logged_text(version_result))
        measured = match.group(0) if match else None
        if measured != job.truvari_version:
            raise ToolVersionError(f"Truvari version mismatch: expected {job.truvari_version}, observed {measured!r}.")
        report["tools"]["observed_truvari_version"] = measured
        bcftools_version = shared.run([job.bcftools_executable, "--version"])
        text = _logged_text(bcftools_version)
        match = re.search(r"bcftools\s+(\S+)", text)
        if match is None:
            raise ToolVersionError("Unable to identify bcftools version.")
        report["tools"]["observed_bcftools_version"] = match.group(1)
        _assert_unchanged(snapshot)
        shared.run(_ref_command(job, job.truth_vcf))
        _assert_unchanged(snapshot)
        report["ref_sequence_check"] = "PASS_BCFTOOLS_NORM_CHECK_ONLY"
        for query, entry in zip(job.queries, report["queries"]):
            _assert_unchanged(snapshot)
            active = _RecordingRunner(underlying, root, records, query.id)
            entry["status"] = "RUNNING"
            try:
                entry["input_audit"] = audit_sv_benchmark_input(query.vcf, job.sample_id, job.reference)
                active.run(_ref_command(job, query.vcf))
                _assert_unchanged(snapshot)
                entry["ref_sequence_check"] = "PASS_BCFTOOLS_NORM_CHECK_ONLY"
                request = TruvariRequest(job.sample_id, job.reference, job.truth_vcf, query.vcf,
                                         root / "queries" / query.id, confident_regions=job.confident_bed,
                                         thresholds=job.thresholds)
                result = TruvariWrapper(executable=job.truvari_executable, runner=active).run(request)
                if result.version != job.truvari_version:
                    raise ToolVersionError("Per-query Truvari version differs from the pinned job version.")
                entry["result"] = result.to_dict()
                for field in ("refinement_requested", "refinement_status", "refinement_proof",
                              "effective_summary_kind", "effective_summary_path", "effective_summary_sha256",
                              "refined_stratification_status"):
                    entry[field] = entry["result"][field]
                for name in ("tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz"):
                    validate_vcf(request.output_directory / name, require_index=True)
                entry["bench_metrics"] = _metrics(result.summary_path)
                if result.refined_summary_path is not None:
                    entry["refined_metrics"] = _metrics(result.refined_summary_path, refined=True)
                entry["status"] = "COMPLETE"
            except SourceImmutabilityError:
                raise
            except (HiFiVarError, OSError) as error:
                entry.update(status="FAILED", error_type=type(error).__name__, error=str(error))
            # Check even after a failed tool: mutation invalidates the entire matrix.
            _assert_unchanged(snapshot)
            persist()
        report["status"] = "COMPLETE" if all(e["status"] == "COMPLETE" for e in report["queries"]) else "PARTIAL"
        _assert_unchanged(snapshot)
        final_snapshot = _snapshot(tuple(Path(name) for name in snapshot))
        report["input_checksums_after"] = final_snapshot
        if final_snapshot != snapshot:
            raise SourceImmutabilityError("Final complete source SHA256/stat snapshot differs from baseline.")
        report["source_immutability"] = "PASS_FULL_SHA256_BEFORE_AFTER_AND_UNCHANGED_STAT"
    except SourceImmutabilityError as error:
        report.update(status="FAILED", source_immutability="FAIL", error_type=type(error).__name__, error=str(error))
        for entry in report["queries"]:
            if entry["status"] == "NOT_RUN":
                entry["status"] = "NOT_RUN_SOURCE_MUTATION"
            elif entry["status"] == "RUNNING":
                entry["status"] = "FAILED_SOURCE_MUTATION"
            elif entry["status"] == "COMPLETE":
                entry["status"] = "INVALIDATED_SOURCE_MUTATION"
        report["metrics_valid_for_claims"] = False
    except (HiFiVarError, OSError) as error:
        report.update(status="FAILED", error_type=type(error).__name__, error=str(error))
        try:
            _assert_unchanged(snapshot)
            report["source_immutability"] = "PASS_FULL_SHA256_BASELINE_AND_UNCHANGED_STAT"
        except SourceImmutabilityError as mutation:
            report.update(source_immutability="FAIL", source_error=str(mutation), metrics_valid_for_claims=False)
    report["completed_at"] = utc_now_iso8601()
    report["output_checksums"] = {
        str(path.relative_to(root)): {"bytes": path.stat().st_size, "sha256": compute_sha256(path)}
        for path in sorted(root.rglob("*")) if path.is_file() and path != matrix_path
    }
    try:
        _assert_unchanged(snapshot)
    except SourceImmutabilityError as error:
        report.update(status="FAILED", source_immutability="FAIL", error_type=type(error).__name__,
                      error=str(error), metrics_valid_for_claims=False)
        for entry in report["queries"]:
            if entry["status"] == "COMPLETE":
                entry["status"] = "INVALIDATED_SOURCE_MUTATION"
    persist()
    return report


__all__ = ["SVBenchmarkJob", "SVBenchmarkQuery", "SourceImmutabilityError", "load_sv_benchmark_job", "run_sv_benchmark"]
