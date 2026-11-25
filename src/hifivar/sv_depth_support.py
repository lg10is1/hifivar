"""Explicit four-read-caller support-response analysis of sealed SV-depth outputs.

No caller, sampling or alignment is run. Native raw calls and the Jasmine union
remain canonical inputs; >=1/2/3/4 views are separate research comparisons.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Mapping, Sequence

import yaml

from hifivar import __version__
from hifivar.assembly_sv import SVEvidenceSource
from hifivar.command import CommandRunner
from hifivar.exceptions import HiFiVarError, InputValidationError, OutputValidationError
from hifivar.harmonization import EvidenceRunStatus, SVEvidenceSourceArtifact, SVHarmonizationRequest
from hifivar.jasmine import JasmineWrapper
from hifivar.serialization import utc_now_iso8601, write_json_atomic
from hifivar.sv import BgzipTabixWrapper, SvCaller, VcfFinalizeRequest
from hifivar.sv_benchmark import SourceImmutabilityError, _metrics, run_sv_benchmark
from hifivar.sv_benchmark_selection import derive_sv_benchmark_query
from hifivar.sv_depth_study import (
    SVDepthStudy, _Loader, _commands, _request, _selection, _stat,
    _TimedRunner, _verify_receipt, load_sv_depth_study,
)
from hifivar.sv_input_identity import SourceIdentityPolicy
from hifivar.sv_support_selection import SupportSelectionRequest, SupportSelectionWrapper, support_selection_paths
from hifivar.validation import compute_sha256

CALLERS = ("sawfish", "sniffles2", "pbsv", "cutesv")
_THRESHOLDS = (1, 2, 3, 4)
_HEX = re.compile(r"[0-9a-f]{64}\Z")
_METRICS = ("tp_base", "tp_call", "fp", "fn", "precision", "recall", "f1")


class SupportContractError(InputValidationError):
    """A static contract reason safe for the user-facing analysis report."""
    def __init__(self, reason):
        self.reason = reason
        super().__init__(reason)


def _require(condition, message):
    if not condition:
        raise SupportContractError(message)


def _failure(error):
    return {"error_type": type(error).__name__,
            "error_code": "SUPPORT_CONTRACT_REJECTED" if isinstance(error, SupportContractError) else "SOURCE_IMMUTABILITY_FAIL" if isinstance(error, SourceImmutabilityError) else "EXTERNAL_VALIDATION_FAILURE",
            "error_reason": error.reason if isinstance(error, SupportContractError) else "Protected source integrity failed." if isinstance(error, SourceImmutabilityError) else "External validation failed; inspect retained private evidence."}



def _plain_path(value, parent, *, exists=True):
    _require(isinstance(value, (str, Path)) and str(value).strip(), "Expected an explicit nonempty path.")
    raw = Path(value).expanduser()
    _require(".." not in raw.parts, "Path traversal is unsupported.")
    raw = raw if raw.is_absolute() else parent / raw
    raw = raw.absolute()
    for part in (raw, *raw.parents):
        _require(not part.is_symlink() and not getattr(part, "is_junction", lambda: False)(), "Symlink/junction paths are unsupported.")
    if exists:
        _require(raw.exists(), "Required source path is missing.")
    return raw


def _mapping(value, keys, label):
    _require(isinstance(value, dict) and set(value) == set(keys), label + " requires exactly its documented keys.")
    return value


def _json_unique(pairs):
    result = {}
    for key, value in pairs:
        _require(key not in result, "Duplicate JSON control key.")
        result[key] = value
    return result


class _ReadSet:
    """Track parsed control bytes, sealed outputs, and explicitly weaker BAM stat guards."""
    def __init__(self):
        self.files = {}
        self.inherited = {}
        self.shared_paths = set()

    def add(self, path, *, fresh=True, sealed_sha=None):
        p = _plain_path(path, Path.cwd())
        _require(p.is_file(), "Source is not a regular file.")
        before = list(_stat(p))
        old = self.files.get(str(p))
        if old is not None and old["stat"] != before:
            raise SourceImmutabilityError("Source metadata changed while reading the analysis.")
        digest = compute_sha256(p) if fresh else sealed_sha
        if list(_stat(p)) != before:
            raise SourceImmutabilityError("Source changed while hashing.")
        if old and old.get("sha256") and digest and old["sha256"] != digest:
            raise SourceImmutabilityError("Source content identity changed.")
        self.files[str(p)] = {"stat": before, "sha256": digest or (old or {}).get("sha256")}
        return p

    def json(self, path):
        p = _plain_path(path, Path.cwd())
        _require(p.is_file() and p.stat().st_size <= 16 * 1024 * 1024, "Control JSON is missing or too large.")
        before = list(_stat(p)); raw = p.read_bytes()
        if list(_stat(p)) != before:
            raise SourceImmutabilityError("Control changed during read.")
        digest = hashlib.sha256(raw).hexdigest()
        old = self.files.get(str(p))
        if old and (old["stat"] != before or old.get("sha256") not in (None, digest)):
            raise SourceImmutabilityError("Control changed between reads.")
        self.files[str(p)] = {"stat": before, "sha256": digest}
        try:
            return json.loads(raw, object_pairs_hook=_json_unique)
        except (ValueError, UnicodeError) as error:
            raise InputValidationError("Invalid source JSON.") from error

    def history(self, path, saved, *, fresh=False):
        p = self.add(path, fresh=fresh)
        _require(isinstance(saved, dict) and isinstance(saved.get("sha256"), str)
                 and _HEX.fullmatch(saved["sha256"]) and saved.get("stat") == list(_stat(p)),
                 "Historical input checksum/stat does not match the declared current input.")
        current = self.files[str(p)]
        if current["sha256"] is not None and current["sha256"] != saved["sha256"]:
            raise SourceImmutabilityError("Historical source SHA differs from freshly verified content.")
        previous = self.inherited.get(str(p))
        if previous is not None and previous != saved:
            raise InputValidationError("Explicit attempts disagree on the same historical source identity.")
        self.inherited[str(p)] = dict(saved)

    def check(self, *, full=False, include_shared=False):
        for raw, item in self.files.items():
            p = _plain_path(raw, Path.cwd())
            if list(_stat(p)) != item["stat"] or full and item["sha256"] is not None and (include_shared or raw not in self.shared_paths) and compute_sha256(p) != item["sha256"]:
                raise SourceImmutabilityError("Source changed during support-response analysis.")


@dataclass(frozen=True)
class SVDepthSupportPlan:
    manifest: Path
    manifest_sha256: str
    job: SVDepthStudy
    source_roots: Mapping[str, Mapping[str, Path]]
    harmonization: Mapping[str, Any]
    identity_policy: SourceIdentityPolicy


def load_sv_depth_support_plan(manifest: str | Path) -> SVDepthSupportPlan:
    """Load only the explicit plan; no output, native probe or caller is executed."""
    p = _plain_path(manifest, Path.cwd())
    _require(p.is_file() and p.stat().st_size <= 16 * 1024 * 1024, "Support plan must be a small file.")
    before = _stat(p); raw = p.read_bytes()
    if _stat(p) != before:
        raise SourceImmutabilityError("Support plan changed during parsing.")
    try:
        data = yaml.load(raw.decode("utf-8-sig"), Loader=_Loader)
    except (UnicodeError, yaml.YAMLError) as error:
        raise InputValidationError("Invalid SV-depth support plan.") from error
    _mapping(data, {"schema_version", "study_manifest", "source_roots", "harmonization", "support_thresholds"}, "Support plan")
    _require(type(data["schema_version"]) is int and data["schema_version"] == 1, "Support schema_version must be integer1.")
    values = data["support_thresholds"]
    _require(isinstance(values, list) and all(type(k) is int for k in values) and values == list(_THRESHOLDS),
             "support_thresholds must be exactly [1,2,3,4].")
    job = load_sv_depth_study(_plain_path(data["study_manifest"], p.parent))
    _require(len({u.coverage_reference_bases for u in job.units}) == 1, "Depth response requires one consistent recorded coverage denominator.")
    _require(job.caller_order == CALLERS, "The original study must declare all four callers in canonical order.")
    roots = data["source_roots"]
    _require(isinstance(roots, dict) and roots and set(roots) <= {u.id for u in job.units}, "source_roots must name declared study units.")
    parsed = {}
    for unit, mapping in roots.items():
        _mapping(mapping, CALLERS, "Each unit source_roots")
        parsed[unit] = {c: _plain_path(mapping[c], p.parent) for c in CALLERS}
        _require(all(v.is_dir() for v in parsed[unit].values()), "Each explicit source root must be an existing study directory.")
    h = _mapping(data["harmonization"], {"jasmine_executable", "bgzip_executable", "tabix_executable", "max_dist", "distance_type", "input_identity_policy"}, "harmonization").copy()
    for key in ("jasmine_executable", "bgzip_executable", "tabix_executable"):
        value = h[key]
        _require(isinstance(value, str) and value.strip() and "\0" not in value, "Harmonization executable must be explicit.")
        if Path(value).is_absolute() or "/" in value or "\\" in value:
            h[key] = str(_plain_path(value, p.parent))
    _require(type(h["max_dist"]) is int and h["max_dist"] > 0 and h["distance_type"] in {"linear", "nonlinear"}, "Invalid frozen Jasmine distance policy.")
    policy = SourceIdentityPolicy.from_config(h["input_identity_policy"])
    _require(policy is not None, "Auditable source identity must be explicitly enabled.")
    return SVDepthSupportPlan(p, hashlib.sha256(raw).hexdigest(), job, parsed, h, policy)


def _source_master(plan, root, reads, cache):
    if root in cache:
        return cache[root]
    _require(not (root / "study.lock").exists(), "Source study is active or has an unresolved lock.")
    job = plan.job
    owner = reads.json(root / "study_owner.json")
    _require(owner == {"schema_version": 1, "study_id": job.study_id, "manifest_sha256": job.manifest_sha256}
             and type(owner.get("schema_version")) is int, "Source study owner identity differs.")
    report = reads.json(root / "study_report.json")
    _require(isinstance(report, dict), "Source master must be a JSON object.")
    _require(report.get("schema_version") == 1 and report.get("study_id") == job.study_id
             and report.get("sample_id") == job.sample_id and report.get("manifest_sha256") == job.manifest_sha256
             and report.get("caller_order") == list(CALLERS), "Source master identity differs from the original study.")
    _require(report.get("status") in {"COMPLETE", "INCOMPLETE", "PARTIAL"} and report.get("completed_at")
             and report.get("metrics_valid_for_claims", True) is True and report.get("output_validation", "PASS") == "PASS"
             and report.get("source_immutability") == "PASS_FULL_SHA256_SELECTED_INPUTS_BEFORE_AFTER_STAT_BETWEEN_CALLERS", "Source master has no valid terminal immutable attempt.")
    rows = report.get("runs")
    expected = {u.id + "/" + c for u in job.units for c in CALLERS}
    _require(isinstance(rows, list) and all(isinstance(r, dict) for r in rows) and len(rows) == len(expected)
             and {r.get("id") for r in rows} == expected, "Source master matrix membership is incomplete or duplicated.")
    for u in job.units:
        for c in CALLERS:
            row = next(r for r in rows if r["id"] == u.id + "/" + c)
            expect = {"unit_id": u.id, "caller": c, "nominal_depth": u.nominal_depth,
                      "replicate": u.replicate, "sampling_fraction": u.sampling_fraction,
                      "replicate_kind": u.replicate_kind, "observed_coverage": u.observed_coverage,
                      "coverage_reference_bases": u.coverage_reference_bases, "coverage_sha256": u.coverage_sha256}
            _require(all(row.get(k) == v for k, v in expect.items()), "Source row depth/replicate/coverage identity differs.")
    shared = (job.manifest, job.reference.fasta, job.reference.fai, job.truth["vcf"],
              Path(str(job.truth["vcf"]) + ".tbi"), job.truth["confident_bed"])
    for p in shared:
        reads.history(p, report.get("shared_input_checksums", {}).get(str(p)), fresh=False)
    cache[root] = report
    return report


def _admit(plan, unit, reads, cache):
    job = plan.job; sources = []; lineage = []
    raw_files = set()
    for caller in CALLERS:
        source_root = plan.source_roots[unit.id][caller]
        report = _source_master(plan, source_root, reads, cache)
        row = next(r for r in report["runs"] if r["id"] == unit.id + "/" + caller)
        _require(row.get("status") == "COMPLETE", "All four explicit caller receipts must be COMPLETE.")
        cell = source_root / "units" / unit.id / caller
        reads.add(cell / "unit_receipt.json")
        receipt = _verify_receipt(cell, row, job.manifest_sha256, expected_job=job)
        _require(reads.json(cell / "unit_receipt.json") == receipt, "Receipt changed while admitting the source.")
        for rel, identity in receipt["output_checksums"].items():
            reads.add(cell / rel, fresh=False, sealed_sha=identity["sha256"])
        # Preserve the producer's full BAM SHA, but do not read a huge BAM again.
        for p in (unit.bam, unit.index, unit.coverage_json):
            reads.history(p, report.get("unit_input_checksums", {}).get(unit.id, {}).get(str(p)), fresh=p == unit.coverage_json)
        _, wrapper, request = _request(job, unit, caller, cell, None)
        result = receipt.get("caller_provenance", {}).get("caller", {})
        expected_commands = _commands(wrapper, request)
        actual_commands = result.get("commands") if caller in {"sawfish", "pbsv"} else [result.get("command")]
        _require(result.get("status") == "completed" and result.get("executed") is True
                 and result.get("request") == request.to_dict() and isinstance(result.get("tool_version"), str)
                 and len(result["tool_version"]) <= 64 and re.fullmatch(r"[0-9]+(?:\.[0-9]+)+[A-Za-z0-9._+-]*", result["tool_version"])
                 and isinstance(actual_commands, list) and len(actual_commands) == len(expected_commands)
                 and all(isinstance(x, dict) for x in actual_commands)
                 and [x.get("args") for x in actual_commands] == expected_commands
                 and all(x.get("shell") is False for x in actual_commands if isinstance(x, dict)),
                 "Caller provenance does not bind the original BAM, parameters and executed commands.")
        raw = _plain_path(cell / "raw" / f"{job.sample_id}.{caller}.sv.vcf.gz", Path.cwd())
        index = _plain_path(str(raw) + ".tbi", Path.cwd())
        sealed = {Path(k).as_posix(): v for k, v in receipt["output_checksums"].items()}
        _require(len(sealed) == len(receipt["output_checksums"]), "Ambiguous canonical output paths.")
        for p in (raw, index):
            key = p.relative_to(cell).as_posix()
            _require(key in sealed and p.stat().st_size == sealed[key]["bytes"], "Canonical raw caller VCF/index is not sealed.")
            reads.add(p, fresh=False, sealed_sha=sealed[key]["sha256"])
        inode = (_stat(raw)[0], _stat(raw)[1])
        _require(inode not in raw_files, "Different callers alias one physical raw VCF.")
        raw_files.add(inode)
        derivation = receipt.get("derivation", {})
        _require(derivation.get("status") == "COMPLETE" and derivation.get("source_vcf") == str(raw)
                 and derivation.get("output_paths", {}).get("vcf") == str(cell / "query/query.ins_del.vcf.gz"),
                 "Benchmark derivation is not bound to the canonical raw caller VCF.")
        for p in (raw, index):
            item = derivation.get("input_checksums", {}).get(str(p), {})
            key = p.relative_to(cell).as_posix()
            _require(item.get("sha256") == sealed[key]["sha256"] and item.get("stat") == list(_stat(p)),
                     "Derived-query source identity differs from sealed raw VCF/index.")
        selection = cell / "query/selection_manifest.json"
        _require(selection.relative_to(cell).as_posix() in sealed and reads.json(selection) == derivation,
                 "Sealed query selection manifest differs from the recorded derivation.")
        if caller in {"pbsv", "cutesv"}:
            final = receipt["caller_provenance"].get("finalization", {})
            artifact = final.get("artifact", {})
            finalize = VcfFinalizeRequest(SvCaller(caller), job.sample_id, job.reference, request.raw_vcf, raw,
                                          result["tool_version"], tuple(tuple(x) for x in expected_commands))
            fw = BgzipTabixWrapper(bgzip_executable=job.settings.bgzip_executable, tabix_executable=job.settings.tabix_executable)
            _require(final.get("status") == "completed"
                     and artifact.get("vcf_path") == str(raw) and artifact.get("index_path") == str(index)
                     and artifact.get("sample_id") == job.sample_id and artifact.get("reference_fasta") == str(job.reference.fasta)
                     and artifact.get("reference_build") == job.reference.build and artifact.get("caller") == caller
                     and [x.get("args") for x in final.get("commands", [])] == [list(x) for x in fw.build_commands(finalize)],
                     "Caller finalization is not bound to the expected native/raw outputs.")
        sources.append(SVEvidenceSourceArtifact(job.sample_id, SVEvidenceSource.READ, caller, raw, index, EvidenceRunStatus.COMPLETED))
        lineage.append({"caller": caller, "source_root": str(source_root), "receipt": str(cell / "unit_receipt.json"),
                        "receipt_sha256": reads.files[str(cell / "unit_receipt.json")]["sha256"],
                        "study_report_sha256": reads.files[str(source_root / "study_report.json")]["sha256"],
                        "producer_hifivar_version": report.get("hifivar_version"), "caller_version": result["tool_version"],
                        "raw_vcf": str(raw), "raw_sha256": sealed[raw.relative_to(cell).as_posix()]["sha256"]})
    return tuple(sources), lineage


def _output_root(output, plan):
    root = _plain_path(output, Path.cwd(), exists=False)
    _require(not root.exists() and root != Path(root.anchor), "Support analysis requires a fresh output root.")
    job = plan.job
    files = (plan.manifest, job.manifest, job.reference.fasta, job.reference.fai, job.truth["vcf"],
             Path(str(job.truth["vcf"]) + ".tbi"), job.truth["confident_bed"])
    protected = {Path(sys.prefix).resolve(), Path(__file__).resolve().parent}
    protected.update(v for m in plan.source_roots.values() for v in m.values())
    for path in files:
        _require(root != path and root not in path.parents, "Output aliases or contains immutable input.")
    for directory in protected:
        _require(root != directory and directory not in root.parents and root not in directory.parents,
                 "Output overlaps a source data, installation or study directory.")
    return root


def _summary(root, report):
    columns = ("receipt_path", "receipt_sha256", "unit_id", "support_ge", "status", "nominal_depth", "observed_coverage", "replicate",
               "sampling_fraction", "replicate_kind", "coverage_reference_bases", *_METRICS,
               "summary_kind", "refinement_status", "summary_path", "summary_sha256", "benchmark_report")
    path = root / "support_study_summary.tsv"
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for unit in report["units"]:
            for row in unit["thresholds"]:
                values = {k: unit.get(k, "") for k in columns}
                values.update({k: row.get(k, "") for k in columns if k in row})
                if row["status"] == "COMPLETE":
                    values.update(row["metrics"])
                writer.writerow({k: "" if values.get(k) is None else values.get(k, "") for k in columns})


def run_sv_depth_support(manifest: str | Path, output_dir: str | Path, *, dry_run: bool = False,
                         unit_ids: Sequence[str] | None = None, runner: CommandRunner | None = None) -> dict[str, Any]:
    """Combine explicit receipt-complete units; never call, retry or repair inputs."""
    _require(type(dry_run) is bool, "dry_run must be boolean.")
    plan = load_sv_depth_support_plan(manifest); job = plan.job
    ids = _selection(unit_ids, list(plan.source_roots), "unit_ids")
    root = _output_root(output_dir, plan)
    reads = _ReadSet(); reads.add(plan.manifest); reads.add(job.manifest); reads.add(Path(__file__).resolve()); cache = {}
    _require(reads.files[str(plan.manifest)]["sha256"] == plan.manifest_sha256
             and reads.files[str(job.manifest)]["sha256"] == job.manifest_sha256, "Manifest changed after parsing.")
    shared = (job.reference.fasta, job.reference.fai, job.truth["vcf"],
              Path(str(job.truth["vcf"]) + ".tbi"), job.truth["confident_bed"])
    reads.shared_paths.update(str(p) for p in shared)
    for p in shared:
        reads.add(p, fresh=not dry_run)
    report = {"schema_version": 1, "hifivar_version": __version__, "status": "PLANNED" if dry_run else "RUNNING",
              "selection_status": "PLANNED" if dry_run else "RUNNING", "manifest_sha256": plan.manifest_sha256,
              "study_manifest_sha256": job.manifest_sha256, "study_id": job.study_id,
              "validator_source_sha256": reads.files[str(Path(__file__).resolve())]["sha256"], "sample_id": job.sample_id,
              "parent_unit_ids": [u.id for u in job.units], "declared_unit_ids": list(plan.source_roots),
              "selected_unit_ids": list(ids), "caller_order": list(CALLERS), "support_thresholds": list(_THRESHOLDS),
              "truth_prefilter": "NONE_ORIGINAL_MIXED_TRUTH_PRESERVED", "fixed_thresholds": job.thresholds.to_dict(),
              "harmonization": {**dict(plan.harmonization), "input_identity_policy": plan.identity_policy.to_dict()},
              "metrics_valid_for_claims": False, "selection_metrics_valid": False,
              "per_command_timeout_policy": "Inherited maximum caller.runtime_minutes*60; not whole-job/memory enforcement",
              "biological_sample_count": 1, "study_design": dict(job.design), "caller_support_is_confidence": False,
              "coverage_source": "INHERITED_MEASURED_JSON_NOT_REMEASURED",
              "coverage_denominator_bases": job.units[0].coverage_reference_bases,
              "coverage_scope": "DECLARED_REFERENCE_DENOMINATOR_NOT_REMEASURED_BENCHMARK_BED_DEPTH", "source_immutability": "PENDING",
              "bam_integrity": "INHERITED_FULL_SHA_WITH_CURRENT_STAT_ONLY_NOT_FRESH_CONTENT_HASH",
              "minimum_depth_claim": "NOT_ESTABLISHED", "cost_reduction_claim": "NOT_ESTABLISHED",
              "created_at": utc_now_iso8601(), "units": []}
    ready = {}
    for unit in job.units:
        selected = unit.id in ids
        state = "PENDING" if selected else "NOT_SELECTED" if unit.id in plan.source_roots else "NOT_DECLARED"
        item = {"unit_id": unit.id, "status": state, "nominal_depth": unit.nominal_depth,
                "observed_coverage": unit.observed_coverage, "replicate": unit.replicate,
                "sampling_fraction": unit.sampling_fraction, "replicate_kind": unit.replicate_kind,
                "coverage_reference_bases": unit.coverage_reference_bases,
                "thresholds": [{"support_ge": k, "status": state, "metrics": None} for k in _THRESHOLDS]}
        report["units"].append(item)
        if selected:
            try:
                ready[unit.id], item["sources"] = _admit(plan, unit, reads, cache)
                item["status"] = "READY"
                for r in item["thresholds"]: r["status"] = "PLANNED" if dry_run else "PENDING"
            except SourceImmutabilityError:
                raise
            except (HiFiVarError, OSError, ValueError, TypeError, KeyError, AttributeError) as error:
                item.update(status="BLOCKED_REQUIRED_CALLERS", **_failure(error), failure_policy="PRESERVED_NO_AUTOMATIC_RETRY")
                for r in item["thresholds"]: r.update(status="BLOCKED_REQUIRED_CALLERS", metrics=None)
    versions = {c: sorted({source["caller_version"] for item in report["units"] if item.get("status") == "READY"
                           for source in item["sources"] if source["caller"] == c}) for c in CALLERS}
    report["caller_version_observations"] = versions
    report["caller_versions"] = {c: values[0] if len(values) == 1 else None for c, values in versions.items()}
    if any(len(values) > 1 for values in versions.values()):
        ready.clear()
        for item in report["units"]:
            if item["status"] == "READY":
                item.update(status="BLOCKED_CALLER_VERSION_MISMATCH", error_code="CALLER_VERSION_MISMATCH",
                            error_reason="The same caller has different recorded versions across selected units.")
                for row in item["thresholds"]:
                    row.update(status="BLOCKED_CALLER_VERSION_MISMATCH", metrics=None)
    reads.check(full=bool(dry_run))
    if dry_run:
        report["selection_status"] = "PLANNED" if all(u["status"] == "READY" for u in report["units"] if u["unit_id"] in ids) else "BLOCKED"
        report["status"] = "PLANNED" if report["selection_status"] == "PLANNED" and len(ids) == len(job.units) else "INCOMPLETE"
        report["source_immutability"] = "READ_ONLY_CONTROLS_AND_SEALED_OUTPUTS_CHECKED_BAM_AND_SHARED_LARGE_INPUTS_STAT_ONLY"
        report["will_execute"] = ["Jasmine once per ready unit", "separate support>=1/2/3/4 views",
                                  "native INS_DEL queries", "one fixed mixed-truth multiquery Truvari job"]
        return report
    root.mkdir(parents=True, exist_ok=False)
    def persist():
        write_json_atomic(report, root / "support_study_report.json", overwrite=True)
        _summary(root, report)
    persist(); queries = []; pending = {}
    timeout = max(getattr(job.settings, c + "_resources").runtime_minutes for c in CALLERS) * 60
    timed = _TimedRunner(runner or CommandRunner(), timeout)
    try:
        for item in report["units"]:
            uid = item["unit_id"]
            if uid not in ready:
                continue
            target = root / "units" / uid
            try:
                reads.check()
                request = SVHarmonizationRequest(job.sample_id, job.reference, ready[uid], target / "work",
                     target / "union" / f"{job.sample_id}.harmonized.sv.vcf.gz",
                     target / "union" / f"{job.sample_id}.sv.evidence.tsv",
                     plan.harmonization["max_dist"], plan.harmonization["distance_type"], False, plan.identity_policy)
                result = JasmineWrapper(executable=plan.harmonization["jasmine_executable"],
                     bgzip_executable=plan.harmonization["bgzip_executable"], tabix_executable=plan.harmonization["tabix_executable"], runner=timed).run(request, stderr_path=target / "logs/jasmine.log")
                _require(result.artifact is not None, "Jasmine did not produce a validated artifact.")
                item["jasmine"] = result.to_dict()
                for p in (result.artifact.vcf_path, result.artifact.index_path, result.artifact.evidence_table):
                    reads.add(p)
                selector = SupportSelectionWrapper(bgzip_executable=plan.harmonization["bgzip_executable"],
                              tabix_executable=plan.harmonization["tabix_executable"], runner=timed)
                for row in item["thresholds"]:
                    k = row["support_ge"]
                    try:
                        selected = selector.run(SupportSelectionRequest(result.artifact, job.reference, k,
                            support_selection_paths(target / "selected", job.sample_id, k)["directory"]),
                            stderr_path=target / "logs" / f"support-ge-{k}.log")
                        row["selection"] = selected.to_dict()
                        for p in (selected.vcf_path, selected.index_path, selected.evidence_table, selected.manifest_path):
                            reads.add(p)
                        view = derive_sv_benchmark_query(selected.vcf_path, target / "query" / f"support-ge-{k}", sample_id=job.sample_id)
                        row["derivation"] = view
                        for p in (Path(view["output_paths"]["vcf"]), Path(view["output_paths"]["vcf"] + ".tbi"), Path(view["output_paths"]["manifest"])):
                            reads.add(p)
                        qid = uid + ".support-ge-" + str(k)
                        queries.append({"id": qid, "vcf": view["output_paths"]["vcf"], "caller": "jasmine_support_ge_" + str(k)})
                        pending[qid] = row; row.update(status="AWAITING_BENCHMARK", query_id=qid)
                    except SourceImmutabilityError:
                        raise
                    except (HiFiVarError, OSError, ValueError, TypeError, KeyError, AttributeError) as error:
                        row.update(status="FAILED_SELECTION_OR_QUERY", metrics=None, **_failure(error), failure_policy="PRESERVED_NO_AUTOMATIC_RETRY")
                    reads.check(); persist()
                item["status"] = "AWAITING_BENCHMARK"
            except SourceImmutabilityError:
                raise
            except (HiFiVarError, OSError, ValueError, TypeError, KeyError, AttributeError) as error:
                item.update(status="FAILED_HARMONIZATION", **_failure(error), failure_policy="PRESERVED_NO_AUTOMATIC_RETRY")
                for row in item["thresholds"]:
                    row.update(status="BLOCKED_HARMONIZATION", metrics=None)
            reads.check(); persist()
        if queries:
            bm = {"schema_version": 1, "benchmark_id": job.study_id + ".support-response",
                  "sample_id": job.sample_id, "reference": {"fasta": str(job.reference.fasta), "build": job.reference.build},
                  "truth": {k: str(v) if isinstance(v, Path) else v for k, v in job.truth.items()},
                  "tools": {k: v for k, v in job.tools.items() if k != "samtools_executable"},
                  "thresholds": job.thresholds.to_dict(), "queries": queries}
            bmfile = root / "benchmark_manifest.json"; write_json_atomic(bm, bmfile)
            result = run_sv_benchmark(bmfile, root / "benchmark", runner=timed)
            report["benchmark_report"] = str(root / "benchmark/benchmark_matrix.json")
            if result.get("source_immutability") == "FAIL":
                raise SourceImmutabilityError("Native benchmark recorded immutable-source failure.")
            if result.get("status") == "FAILED" or result.get("metrics_valid_for_claims", True) is not True:
                raise OutputValidationError("Native benchmark failed its global integrity contract.")
            _require(reads.json(root / "benchmark/benchmark_matrix.json") == result, "Native benchmark report changed while reading.")
            native_seals = result.get("output_checksums")
            _require(isinstance(native_seals, dict) and native_seals, "Native benchmark has no sealed output identities.")
            normalized_seals = {Path(k).as_posix(): v for k, v in native_seals.items()}
            _require(len(normalized_seals) == len(native_seals), "Native benchmark output paths are ambiguous.")
            for rel, identity in normalized_seals.items():
                path = _plain_path(root / "benchmark" / rel, root)
                _require(path.is_relative_to(root / "benchmark") and isinstance(identity, dict)
                         and type(identity.get("bytes")) is int and identity["bytes"] >= 0
                         and isinstance(identity.get("sha256"), str) and _HEX.fullmatch(identity["sha256"]),
                         "Native output checksum contract is malformed.")
                reads.add(path)
                _require(path.stat().st_size == identity["bytes"] and reads.files[str(path)]["sha256"] == identity["sha256"],
                         "Native output differs from its sealed checksum.")
            native_rows = result.get("queries", [])
            _require(len(native_rows) == len(pending) and {r.get("id") for r in native_rows} == set(pending), "Native benchmark query membership changed.")
            if result.get("source_immutability") == "FAIL":
                raise SourceImmutabilityError("Native benchmark recorded immutable-source failure.")
            if result.get("status") == "FAILED" or result.get("metrics_valid_for_claims", True) is not True:
                raise OutputValidationError("Native benchmark failed its global integrity contract.")
            for native in native_rows:
                row = pending[native["id"]]
                if native.get("status") != "COMPLETE":
                    row.update(status="FAILED_BENCHMARK", metrics=None, error_type=native.get("error_type"), failure_policy="NATIVE_FAILURE_PRESERVED")
                    continue
                kind = native.get("effective_summary_kind", "refined" if job.thresholds.refine else "bench")
                if job.thresholds.refine and kind == "bench":
                    from hifivar.truvari_refinement import verify_noop_result
                    proof = verify_noop_result(native.get("result", {}))
                    _require(native.get("refinement_status") == "NOT_NEEDED" and native.get("refined_metrics") == []
                             and native.get("refinement_proof") == proof, "Bench metrics cannot impersonate requested refinement.")
                elif job.thresholds.refine:
                    _require(kind == "refined" and native.get("refinement_status") == "COMPLETE", "Refined output did not complete.")
                else:
                    _require(kind == "bench" and native.get("refinement_status") == "NOT_REQUESTED", "Unexpected refinement branch.")
                summary = root / "benchmark/queries" / native["id"] / ("refine.variant_summary.json" if kind == "refined" else "summary.json")
                _require(kind in {"bench", "refined"} and summary.relative_to(root / "benchmark").as_posix() in normalized_seals, "Unknown or unsealed native metric source.")
                metrics = {v["name"]: v["value"] for v in _metrics(summary, refined=kind == "refined")}
                _require(set(_METRICS) <= set(metrics) and metrics == {v["name"]: v["value"] for v in native[kind + "_metrics"]}, "Native metric row differs from its actual summary.")
                for name in ("tp_base", "tp_call", "fp", "fn"):
                    _require(type(metrics[name]) is int and metrics[name] >= 0, "Invalid native assignment count.")
                for name in ("precision", "recall", "f1"):
                    _require(metrics[name] is None or type(metrics[name]) in {int, float} and math.isfinite(metrics[name]) and 0 <= metrics[name] <= 1, "Invalid native rate.")
                for name, a, b in (("precision", "tp_call", "fp"), ("recall", "tp_base", "fn")):
                    if metrics[a] + metrics[b]:
                        _require(metrics[name] is not None and math.isclose(metrics[name], metrics[a] / (metrics[a] + metrics[b]), rel_tol=1e-5, abs_tol=1e-8), "Native TP denominator does not reconcile.")
                if metrics["precision"] is not None and metrics["recall"] is not None and metrics["precision"] + metrics["recall"]:
                    expected_f1 = 2 * metrics["precision"] * metrics["recall"] / (metrics["precision"] + metrics["recall"])
                    _require(metrics["f1"] is not None and math.isclose(metrics["f1"], expected_f1, rel_tol=1e-5, abs_tol=1e-8), "Native F1 does not reconcile.")
                row.update(status="COMPLETE", metrics=metrics, summary_kind=kind,
                           refinement_status=native.get("refinement_status"), summary_path=str(summary),
                           summary_sha256=compute_sha256(summary), benchmark_report=report["benchmark_report"],
                           precision_formula="TP_COMP/(TP_COMP+FP)", recall_formula="TP_BASE/(TP_BASE+FN)")
            if result.get("source_immutability") == "FAIL":
                raise SourceImmutabilityError("Native benchmark source mutation invalidates all comparisons.")
        reads.check()
        for item in report["units"]:
            if item["unit_id"] in ids and item["status"] == "AWAITING_BENCHMARK":
                item["status"] = "COMPLETE" if all(r["status"] == "COMPLETE" for r in item["thresholds"]) else "PARTIAL"
        report["selection_status"] = "COMPLETE" if all(u["status"] == "COMPLETE" for u in report["units"] if u["unit_id"] in ids) else "FAILED"
        report["status"] = "COMPLETE" if all(u["status"] == "COMPLETE" for u in report["units"]) else "PARTIAL" if report["selection_status"] == "FAILED" else "INCOMPLETE"
        report["source_immutability"] = "PASS_FRESH_CONTROL_AND_OUTPUT_SHA256_WITH_INHERITED_BAM_SHA_CURRENT_STAT_ONLY"
    except SourceImmutabilityError as error:
        report.update(status="FAILED", selection_status="FAILED", source_immutability="FAIL", **_failure(error), failure_policy="PRESERVED_NO_AUTOMATIC_RETRY")
        for item in report["units"]:
            if item["unit_id"] in ids:
                item["status"] = "INVALIDATED_SOURCE_MUTATION"
            for row in item["thresholds"]:
                row.update(status="INVALIDATED_SOURCE_MUTATION", metrics=None)
    except (HiFiVarError, OSError, ValueError, TypeError, KeyError, AttributeError) as error:
        report.update(status="FAILED", selection_status="FAILED", output_validation="FAIL", **_failure(error), failure_policy="PRESERVED_NO_AUTOMATIC_RETRY")
        for item in report["units"]:
            if item["unit_id"] in ids:
                item["status"] = "INVALIDATED_PIPELINE_FAILURE"
            for row in item["thresholds"]:
                row.update(status="INVALIDATED_PIPELINE_FAILURE", metrics=None)
    # Every branch, including a tool/serialization error, ends with immutable-input readback.
    source_failure_recorded = report.get("source_immutability") == "FAIL"
    try:
        reads.check(full=True, include_shared=True)
        for mapping in plan.source_roots.values():
            for source in mapping.values():
                if (source / "study.lock").exists():
                    raise SourceImmutabilityError("Source writer appeared during support analysis.")
        if not source_failure_recorded:
            report["source_immutability"] = "PASS_FRESH_CONTROL_AND_OUTPUT_SHA256_WITH_INHERITED_BAM_SHA_CURRENT_STAT_ONLY"
    except (HiFiVarError, OSError, ValueError) as error:
        report.update(status="FAILED", selection_status="FAILED", source_immutability="FAIL", **_failure(error))
        for item in report["units"]:
            if item["unit_id"] in ids:
                item["status"] = "INVALIDATED_SOURCE_MUTATION"
            for row in item["thresholds"]:
                row.update(status="INVALIDATED_SOURCE_MUTATION", metrics=None)
    report["metrics_valid_for_claims"] = report["status"] == "COMPLETE"
    report["selection_metrics_valid"] = report["selection_status"] == "COMPLETE" and report["source_immutability"].startswith("PASS_")
    try:
        if report["source_immutability"].startswith("PASS_"):
            for item in report["units"]:
                for row in item["thresholds"]:
                    if row["status"] != "COMPLETE":
                        continue
                    k = row["support_ge"]; uid = item["unit_id"]
                    paths = [Path(row["summary_path"])]
                    selection = row["selection"]
                    for key in ("vcf", "index", "evidence"):
                        paths.append(Path(selection["request"]["paths"][key]))
                    paths.append(Path(selection["manifest_path"]))
                    for key in ("vcf", "index", "manifest"):
                        paths.append(Path(row["derivation"]["output_paths"][key]))
                    native = root / "benchmark/queries" / row["query_id"]
                    paths.extend(p for p in native.rglob("*") if p.is_file())
                    checks = {}
                    for path in dict.fromkeys(paths):
                        path = _plain_path(path, root)
                        _require(path.is_relative_to(root), "Threshold receipt output escaped the owned root.")
                        digest = compute_sha256(path)
                        if str(path) in reads.files and reads.files[str(path)]["sha256"] is not None:
                            _require(digest == reads.files[str(path)]["sha256"] and list(_stat(path)) == reads.files[str(path)]["stat"],
                                     "Threshold output changed after metric extraction.")
                        checks[path.relative_to(root).as_posix()] = {"sha256": digest, "bytes": path.stat().st_size}
                    receipt = {"schema_version": 1, "status": "COMPLETE", "unit_id": uid, "support_ge": k,
                               "plan_sha256": plan.manifest_sha256, "study_manifest_sha256": job.manifest_sha256,
                               "caller_order": list(CALLERS), "source_lineage": item["sources"],
                               "comparison": dict(row), "output_checksums": checks,
                               "parent_matrix_complete": report["status"] == "COMPLETE"}
                    target = root / "units" / uid / "threshold_receipts" / ("support-ge-" + str(k) + ".json")
                    write_json_atomic(receipt, target)
                    row.update(receipt_path=str(target), receipt_sha256=compute_sha256(target))
        reads.check(full=True)
    except (HiFiVarError, OSError, ValueError, TypeError, KeyError) as error:
        report.update(status="FAILED", selection_status="FAILED", output_validation="FAIL",
                      metrics_valid_for_claims=False, selection_metrics_valid=False, **_failure(error))
        for item in report["units"]:
            if item["unit_id"] in ids:
                item["status"] = "INVALIDATED_OUTPUT_SEAL"
            for row in item["thresholds"]:
                row.update(status="INVALIDATED_OUTPUT_SEAL", metrics=None)
        try:
            reads.check()
        except (HiFiVarError, OSError):
            report["source_immutability"] = "FAIL"
    report["inherited_input_checksums"] = reads.inherited
    report["current_validator_inputs"] = reads.files
    report["completed_at"] = utc_now_iso8601(); persist()
    return report


__all__ = ["SVDepthSupportPlan", "load_sv_depth_support_plan", "run_sv_depth_support"]
