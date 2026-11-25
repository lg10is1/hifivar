"""Depth/replicate read-SV calling and matched-truth benchmarking.

This separate manifest entry reuses measured, indexed BAM units. It neither
resamples reads nor silently repairs BAMs. Caller failures are independent;
resume verifies completed immutable bundles and preserves failed/incomplete
bundles without retrying them. Support selection remains a separate analysis.
"""
from __future__ import annotations

import csv
from contextlib import contextmanager
from dataclasses import dataclass, fields
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import re
from uuid import uuid4
from typing import Any, Mapping, Sequence

import yaml

from hifivar import __version__
from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import AlignmentArtifact, AlignmentSortOrder, AlignmentSource
from hifivar.command import CommandRunner
from hifivar.benchmark_inputs import audit_sv_benchmark_input, audit_benchmark_regions
from hifivar.cutesv import CuteSvRequest, CuteSvWrapper
from hifivar.exceptions import HiFiVarError, InputValidationError, OutputValidationError
from hifivar.launch_inputs import IndexedBamPreflight, unique_bam_index
from hifivar.pbsv import PbsvRequest, PbsvWrapper
from hifivar.phase4 import Phase4Settings
from hifivar.reference import ReferenceGenome, normalize_reference_build
from hifivar.sawfish import SawfishRequest, SawfishWrapper
from hifivar.serialization import utc_now_iso8601, write_json_atomic
from hifivar.sniffles2 import Sniffles2Request, Sniffles2Wrapper
from hifivar.sv import BgzipTabixWrapper, SvCaller, VcfFinalizeRequest, create_structural_variant_artifact
from hifivar.sv_benchmark import SourceImmutabilityError, run_sv_benchmark, _metrics
from hifivar.sv_benchmark_selection import derive_sv_benchmark_query
from hifivar.truvari import TruvariThresholds
from hifivar.validation import compute_sha256

_SCHEMA = {"schema_version", "study_id", "sample_id", "reference", "truth",
           "tools", "thresholds", "sv", "design", "units"}
_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
_CALLERS = ("sawfish", "sniffles2", "pbsv", "cutesv")


class _Loader(yaml.SafeLoader):
    pass


def _unique(loader, node):
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=True)
        if not isinstance(key, str) or key in result:
            raise InputValidationError("SV-depth YAML mapping keys must be unique strings.")
        result[key] = loader.construct_object(value_node, deep=True)
    return result


_Loader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _unique)


def _map(value, required, label):
    if not isinstance(value, dict) or set(value) != set(required):
        raise InputValidationError(f"{label} requires exactly {sorted(required)}.")
    return value


def _text(value, label):
    if not isinstance(value, str) or not value.strip() or value != value.strip() or any(c in value for c in "\0\r\n"):
        raise InputValidationError(f"{label} must be an explicit trimmed string.")
    return value


def _id(value, label):
    value = _text(value, label)
    if not _ID.fullmatch(value):
        raise InputValidationError(f"{label} must be one safe identifier.")
    return value


def _number(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise InputValidationError(f"{label} must be finite and positive.")
    return float(value)


def _path(value, parent, label):
    raw = _text(value, label)
    if ".." in PurePosixPath(raw.replace("\\", "/")).parts or raw.startswith("~"):
        raise InputValidationError(f"{label} must not use traversal or implicit home.")
    p = Path(raw)
    p = (p if p.is_absolute() else parent / p).resolve()
    if not p.is_file() or p.stat().st_size == 0:
        raise InputValidationError(f"{label} requires a nonempty file: {p}")
    return p


def _exe(value, parent, label):
    value = _text(value, label)
    if value.startswith("-") or not ("/" in value or "\\" in value) and any(c.isspace() for c in value):
        raise InputValidationError(f"{label} must be one executable, not a shell command.")
    if "/" in value or "\\" in value:
        p = Path(value)
        return str((p if p.is_absolute() else parent / p).resolve())
    return value


def _strict_sv(sv):
    """Reject silent Phase4 int/bool coercion in this standalone manifest."""
    _map(sv, {"enabled", "overwrite", "sawfish", "sniffles2", "pbsv", "cutesv", "finalization"}, "sv")
    base = {"enabled", "executable", "threads", "memory_mb", "runtime_minutes"}
    extras = {"sawfish": {"disable_cnv"},
              "sniffles2": {"minimum_support", "minimum_sv_length"}, "pbsv": set(),
              "cutesv": {"minimum_support", "minimum_sv_size", "max_cluster_bias_ins",
                         "diff_ratio_merging_ins", "max_cluster_bias_del",
                         "diff_ratio_merging_del", "genotype"}}
    _map(sv["finalization"], {"bgzip_executable", "tabix_executable"}, "sv.finalization")
    for name in _CALLERS:
        block = _map(sv[name], base | extras[name], "sv." + name)
        for key in ("enabled", "disable_cnv", "genotype"):
            if key in block and type(block[key]) is not bool:
                raise InputValidationError("sv." + name + "." + key + " must be boolean.")
        for key in ("threads", "memory_mb", "runtime_minutes", "minimum_support", "minimum_sv_length",
                    "minimum_sv_size", "max_cluster_bias_ins", "max_cluster_bias_del"):
            if key not in block or key == "minimum_support" and name == "sniffles2" and block[key] is None:
                continue
            if type(block[key]) is not int or block[key] <= 0:
                raise InputValidationError("sv." + name + "." + key + " must be a positive integer.")
        for key in ("diff_ratio_merging_ins", "diff_ratio_merging_del"):
            if key in block:
                value = block[key]
                if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 <= value <= 1:
                    raise InputValidationError("sv." + name + "." + key + " must be a finite rate.")


@dataclass(frozen=True)
class SVDepthUnit:
    id: str
    nominal_depth: float
    replicate: int
    sampling_fraction: float
    bam: Path
    index: Path
    coverage_json: Path
    observed_coverage: float
    coverage_reference_bases: int
    coverage_sha256: str

    @property
    def replicate_kind(self) -> str:
        return "FULL_INPUT_TECHNICAL_REPEAT" if self.sampling_fraction == 1 else "NESTED_TECHNICAL_SUBSAMPLE"


@dataclass(frozen=True)
class SVDepthStudy:
    manifest: Path
    manifest_sha256: str
    study_id: str
    sample_id: str
    reference: ReferenceGenome
    reference_fai_sha256: str
    truth: Mapping[str, Any]
    tools: Mapping[str, str]
    thresholds: TruvariThresholds
    settings: Phase4Settings
    sv_config: Mapping[str, Any]
    design: Mapping[str, Any]
    units: tuple[SVDepthUnit, ...]

    @property
    def caller_order(self) -> tuple[str, ...]:
        return tuple(c.value for c in self.settings.enabled_callers)


def _coverage(path: Path, bam: Path, sample_id: str) -> tuple[float, int, str]:
    try:
        before = _stat(path)
        raw = path.read_bytes()
        if _stat(path) != before:
            raise SourceImmutabilityError("Coverage metadata changed while reading.")
        data = json.loads(raw.decode("utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise InputValidationError("Unreadable observed-coverage JSON.") from error
    if not isinstance(data, dict) or data.get("sample_id") != sample_id:
        raise InputValidationError("Coverage sample identity must match the study.")
    binding = Path(_text(data.get("alignment_bam"), "coverage.alignment_bam"))
    binding = (binding if binding.is_absolute() else path.parent / binding).resolve()
    if binding != bam:
        raise InputValidationError("Coverage measurement is bound to a different BAM.")
    if type(data.get("reference_bases")) is not int or data["reference_bases"] <= 0:
        raise InputValidationError("Coverage must retain its positive reference-base denominator.")
    return _number(data.get("observed_coverage"), "observed_coverage"), data["reference_bases"], hashlib.sha256(raw).hexdigest()


def load_sv_depth_study(manifest: str | Path) -> SVDepthStudy:
    """Read a strict study plan without native probes, hashes of BAMs or outputs."""
    source = _path(str(manifest), Path.cwd(), "manifest")
    stamp = _stat(source)
    raw = source.read_bytes()
    if _stat(source) != stamp:
        raise SourceImmutabilityError("Study manifest changed while reading.")
    try:
        data = yaml.load(raw.decode("utf-8-sig"), Loader=_Loader)
    except (UnicodeError, yaml.YAMLError) as error:
        raise InputValidationError("Invalid SV-depth YAML.") from error
    data = _map(data, _SCHEMA, "SV-depth manifest")
    if type(data["schema_version"]) is not int or data["schema_version"] != 1:
        raise InputValidationError("SV-depth schema_version must be integer 1.")
    sample = _id(data["sample_id"], "sample_id")
    ref = _map(data["reference"], {"fasta", "build"}, "reference")
    fasta = _path(ref["fasta"], source.parent, "reference.fasta")
    fai = _path(str(fasta) + ".fai", source.parent, "reference.fai")
    fai_stamp = _stat(fai)
    reference = ReferenceGenome.from_fasta(fasta, build=_text(ref["build"], "reference.build"))
    fai_digest = compute_sha256(fai)
    if _stat(fai) != fai_stamp:
        raise SourceImmutabilityError("FAI metadata changed during reference construction.")
    truth = _map(data["truth"], {"vcf", "confident_bed", "version", "source", "reference_build"}, "truth").copy()
    if normalize_reference_build(_text(truth["reference_build"], "truth.reference_build")) != reference.build:
        raise InputValidationError("Truth/reference build mismatch; contig renaming does not change build.")
    for name in ("version", "source"):
        truth[name] = _text(truth[name], "truth." + name)
    if truth["version"].casefold() in {"latest", "unknown", "tbd"}:
        raise InputValidationError("Truth version must be pinned.")
    truth["vcf"] = _path(truth["vcf"], source.parent, "truth.vcf")
    _path(str(truth["vcf"]) + ".tbi", source.parent, "truth.index")
    truth["confident_bed"] = _path(truth["confident_bed"], source.parent, "truth.confident_bed")
    tools = _map(data["tools"], {"samtools_executable", "truvari_executable", "bcftools_executable", "truvari_version"}, "tools").copy()
    if tools["truvari_version"] != "5.4.0":
        raise InputValidationError("SV-depth query selection requires pinned Truvari 5.4.0.")
    for name in ("samtools_executable", "truvari_executable", "bcftools_executable"):
        tools[name] = _exe(tools[name], source.parent, "tools." + name)
    policy = data["thresholds"]
    known = {f.name for f in fields(TruvariThresholds)}
    required = {"sizemin", "sizemax", "refdist", "pctseq", "pctsize", "pass_only", "pick", "refine"}
    if not isinstance(policy, dict) or set(policy) - known or not required <= set(policy):
        raise InputValidationError("SV-depth thresholds require explicit size/distance/sequence/filter/pick/refine policy.")
    thresholds = TruvariThresholds(**policy)
    if thresholds.sizemin is None or thresholds.sizemin < 50:
        raise InputValidationError("Primary SV-depth benchmark requires sizemin >= 50.")
    design = _map(data["design"], {"sampling_unit", "nested", "biological_samples", "full_input_nominal_depth"}, "design").copy()
    if design["sampling_unit"] != "pacbio_molecule" or design["nested"] is not True or type(design["biological_samples"]) is not int or design["biological_samples"] != 1:
        raise InputValidationError("This entry evaluates one matched-truth biological sample with declared nested molecule subsampling.")
    _number(design["full_input_nominal_depth"], "design.full_input_nominal_depth")
    if not isinstance(data["sv"], dict):
        raise InputValidationError("sv must contain the existing complete Phase4 configuration.")
    sv = data["sv"]
    _strict_sv(sv)
    if sv.get("enabled") is not True or sv.get("overwrite") is not False:
        raise InputValidationError("SV-depth requires sv.enabled=true and overwrite=false.")
    for section in ("sawfish", "sniffles2", "pbsv", "cutesv"):
        if not isinstance(sv.get(section), dict) or type(sv[section].get("enabled")) is not bool:
            raise InputValidationError("Every caller requires explicit boolean enabled.")
    for name in _CALLERS:
        sv[name]["executable"] = _exe(sv[name].get("executable"), source.parent, "sv." + name + ".executable")
    for name in ("bgzip_executable", "tabix_executable"):
        if not isinstance(sv.get("finalization"), dict):
            raise InputValidationError("sv.finalization must be explicit.")
        sv["finalization"][name] = _exe(sv["finalization"].get(name), source.parent, "sv.finalization." + name)
    settings = Phase4Settings.from_config({"sv": sv})
    if not settings.enabled_callers:
        raise InputValidationError("At least one read-SV caller must be enabled.")
    units = []
    seen_ids, seen_pairs, bam_ids = set(), set(), set()
    if not isinstance(data["units"], list) or not data["units"]:
        raise InputValidationError("units must be a nonempty list.")
    for entry in data["units"]:
        unit = _map(entry, {"id", "nominal_depth", "replicate", "sampling_fraction", "bam", "coverage_json"}, "unit")
        identifier = _id(unit["id"], "unit.id")
        nominal = _number(unit["nominal_depth"], "unit.nominal_depth")
        fraction = _number(unit["sampling_fraction"], "unit.sampling_fraction")
        replicate = unit["replicate"]
        if type(replicate) is not int or replicate <= 0 or fraction > 1:
            raise InputValidationError("Replicate must be positive integer and fraction <= 1.")
        if not math.isclose(fraction, nominal / float(design["full_input_nominal_depth"]), rel_tol=1e-9, abs_tol=1e-12):
            raise InputValidationError("sampling_fraction must equal nominal_depth/full_input_nominal_depth; observed depth is separate.")
        if nominal > float(design["full_input_nominal_depth"]) or (fraction == 1) != (nominal == float(design["full_input_nominal_depth"])):
            raise InputValidationError("Full-input units must match declared full-input nominal depth.")
        bam = _path(unit["bam"], source.parent, "unit.bam")
        if bam.suffix.casefold() != ".bam":
            raise InputValidationError("SV-depth input is indexed BAM only; CRAM is deferred.")
        index = unique_bam_index(bam)
        identity = (bam.stat().st_dev, bam.stat().st_ino)
        if identifier.casefold() in seen_ids or (nominal, replicate) in seen_pairs or identity in bam_ids:
            raise InputValidationError("Unit IDs, depth/replicate pairs and BAM file identities must be unique.")
        seen_ids.add(identifier.casefold());seen_pairs.add((nominal, replicate));bam_ids.add(identity)
        coverage = _path(unit["coverage_json"], source.parent, "unit.coverage_json")
        observed, denominator, coverage_digest = _coverage(coverage, bam, sample)
        units.append(SVDepthUnit(identifier, nominal, replicate, fraction, bam, index, coverage, observed, denominator, coverage_digest))
    audit_sv_benchmark_input(truth["vcf"], sample, reference)
    audit_benchmark_regions(truth["confident_bed"], reference)
    for replicate in {u.replicate for u in units}:
        ordered = sorted((u for u in units if u.replicate == replicate), key=lambda u: u.nominal_depth)
        if any(a.sampling_fraction >= b.sampling_fraction for a, b in zip(ordered, ordered[1:])):
            raise InputValidationError("Declared nested fractions must increase with nominal depth.")
    return SVDepthStudy(source, hashlib.sha256(raw).hexdigest(), _id(data["study_id"], "study_id"),
                        sample, reference, fai_digest, truth, tools, thresholds, settings, sv, design, tuple(units))


def _stat(path: Path) -> tuple[int, ...]:
    x = path.stat()
    return x.st_dev, x.st_ino, x.st_size, x.st_mtime_ns, x.st_ctime_ns


def _snapshot(paths):
    result = {}
    for path in dict.fromkeys(paths):
        before = _stat(path)
        digest = compute_sha256(path)
        if _stat(path) != before:
            raise SourceImmutabilityError("Study input changed while hashing.")
        result[str(path)] = {"sha256": digest, "stat": list(before)}
    return result


def _unchanged(snapshot, *, full=False):
    for raw, value in snapshot.items():
        p = Path(raw)
        if list(_stat(p)) != value["stat"] or full and compute_sha256(p) != value["sha256"]:
            raise SourceImmutabilityError(f"SV-depth source changed: {p}")


def _root(path, inputs, resume):
    raw = Path(path)
    if ".." in raw.parts or raw.is_symlink():
        raise OutputValidationError("Output must not use traversal or a symlink.")
    root = raw.resolve()
    if root == Path(root.anchor) or any(root == p or root in p.parents for p in inputs):
        raise OutputValidationError("Study output cannot alias/contain immutable input files.")
    if root.exists() and (not resume or not root.is_dir()):
        raise OutputValidationError("A fresh output directory is required unless resume is explicit.")
    return root


def _selection(values, known, label):
    if values is None:
        return tuple(known)
    if isinstance(values, (str, bytes)) or not values or len(set(values)) != len(values) or set(values) - set(known):
        raise InputValidationError(f"{label} must contain unique declared identifiers.")
    return tuple(x for x in known if x in values)


def _request(job, unit, caller, root, runner):
    a = AlignmentArtifact(job.sample_id, unit.bam, AlignmentOutputFormat.BAM, job.reference,
                          AlignmentSource.EXISTING, AlignmentSortOrder.COORDINATE, unit.index)
    s = job.settings
    outputs, work = root / "raw", root / "work"
    if caller == "sawfish":
        w = SawfishWrapper(executable=s.sawfish_executable, runner=runner)
        r = SawfishRequest.create(a, outputs, work, resources=s.sawfish_resources, disable_cnv=s.sawfish_disable_cnv)
    elif caller == "sniffles2":
        w = Sniffles2Wrapper(executable=s.sniffles2_executable, runner=runner)
        r = Sniffles2Request.create(a, outputs, resources=s.sniffles2_resources, minimum_support=s.sniffles2_minimum_support, minimum_sv_length=s.sniffles2_minimum_sv_length)
    elif caller == "pbsv":
        w = PbsvWrapper(executable=s.pbsv_executable, runner=runner)
        r = PbsvRequest.create(a, work / "native", work / "signatures", resources=s.pbsv_resources)
    else:
        w = CuteSvWrapper(executable=s.cutesv_executable, runner=runner)
        r = CuteSvRequest.create(a, work / "native", work / "tmp", resources=s.cutesv_resources,
             minimum_support=s.cutesv_minimum_support, minimum_sv_size=s.cutesv_minimum_sv_size,
             max_cluster_bias_ins=s.cutesv_max_cluster_bias_ins, diff_ratio_merging_ins=s.cutesv_diff_ratio_merging_ins,
             max_cluster_bias_del=s.cutesv_max_cluster_bias_del, diff_ratio_merging_del=s.cutesv_diff_ratio_merging_del,
             genotype=s.cutesv_genotype)
    return a, w, r


def _commands(wrapper, request):
    if hasattr(wrapper, "build_commands"):
        return [list(c) for c in wrapper.build_commands(request)]
    return [list(wrapper.build_command(request))]


def _run_caller(job, unit, caller, root, runner):
    a, wrapper, request = _request(job, unit, caller, root, runner)
    result = wrapper.run(request, stderr_path=root / "logs" / (caller + ".stderr.log"))
    commands = tuple(c.args for c in result.commands) if hasattr(result, "commands") else (result.command.args,)
    if caller in {"pbsv", "cutesv"}:
        final = BgzipTabixWrapper(runner=runner, bgzip_executable=job.settings.bgzip_executable,
                                  tabix_executable=job.settings.tabix_executable).run(
            VcfFinalizeRequest(SvCaller(caller), job.sample_id, job.reference, request.raw_vcf,
                 root / "raw" / f"{job.sample_id}.{caller}.sv.vcf.gz", result.tool_version, commands),
            stderr_path=root / "logs" / (caller + ".finalize.stderr.log"))
        if final.artifact is None:
            raise OutputValidationError("Caller finalization produced no validated SV artifact.")
        artifact = final.artifact
        caller_provenance = {"caller": result.to_dict(), "finalization": final.to_dict()}
    else:
        artifact = create_structural_variant_artifact(caller=SvCaller(caller), sample_id=job.sample_id,
                    reference=job.reference, vcf_path=request.output_vcf,
                    caller_version=result.tool_version, commands=commands)
        caller_provenance = {"caller": result.to_dict()}
    return artifact, caller_provenance


def _benchmark_manifest(job, unit, caller, query):
    return {"schema_version": 1, "benchmark_id": job.study_id + "." + unit.id + "." + caller,
            "sample_id": job.sample_id,
            "reference": {"fasta": str(job.reference.fasta), "build": job.reference.build},
            "truth": {k: str(v) if isinstance(v, Path) else v for k, v in job.truth.items()},
            "tools": {k: v for k, v in job.tools.items() if k != "samtools_executable"},
            "thresholds": job.thresholds.to_dict(),
            "queries": [{"id": caller, "vcf": str(query), "caller": caller}]}


def _output_checksums(root):
    # Seal final/raw VCFs, native summaries and logs. Heavy scratch/signature
    # binaries are retained but are not falsely represented as sealed outputs.
    result = {}
    for p in sorted(root.rglob("*")):
        if not p.is_file() or p.name == "unit_receipt.json":
            continue
        rel = p.relative_to(root)
        if rel.parts[0] == "work" and not p.name.endswith((".vcf", ".vcf.gz", ".tbi", ".json", ".log")):
            continue
        result[str(rel)] = {"sha256": compute_sha256(p), "bytes": p.stat().st_size}
    return result



def _verify_receipt(root, entry, manifest_sha, *, expected_job=None):
    receipt = root / "unit_receipt.json"
    data = json.loads(receipt.read_text(encoding="utf-8"))
    if data.get("status") != "COMPLETE" or data.get("manifest_sha256") != manifest_sha or data.get("id") != entry["id"]:
        raise OutputValidationError("Completed receipt identity does not match resume.")
    if (not isinstance(data.get("output_checksums"), dict) or not data["output_checksums"]
            or "benchmark/benchmark_matrix.json" not in {Path(k).as_posix() for k in data["output_checksums"]}):
        raise OutputValidationError("Completed receipt lacks sealed benchmark artifacts.")
    for rel, value in data["output_checksums"].items():
        p = root / rel
        if not p.resolve().is_relative_to(root) or p.is_symlink() or not p.is_file() or compute_sha256(p) != value["sha256"]:
            raise OutputValidationError("Completed study output changed; resume refuses reuse.")
    for name in ("unit_id", "caller", "nominal_depth", "observed_coverage", "replicate",
                 "sampling_fraction", "replicate_kind", "coverage_reference_bases", "coverage_sha256", "metrics", "summary_kind",
                 "benchmark_report", "caller_provenance", "derivation",
                  "precision_formula", "recall_formula", "refinement_requested", "refinement_status",
                  "refinement_proof", "effective_summary_kind", "effective_summary_path",
                  "effective_summary_sha256"):
        if data.get(name) != entry.get(name):
            raise OutputValidationError("Master report differs from sealed unit receipt: " + name)
    kind = data.get("summary_kind")
    caller = data.get("caller")
    if kind not in {"bench", "refined"} or caller not in _CALLERS:
        raise OutputValidationError("Sealed unit metric protocol is invalid.")
    control = root / "benchmark_manifest.json"
    sealed = {Path(k).as_posix() for k in data["output_checksums"]}
    if control.relative_to(root).as_posix() not in sealed:
        raise OutputValidationError("Completed receipt lacks its sealed benchmark protocol.")
    protocol = json.loads(control.read_text(encoding="utf-8"))
    requested = protocol.get("thresholds", {}).get("refine")
    if type(requested) is not bool:
        raise OutputValidationError("Sealed benchmark protocol lacks explicit refinement policy.")
    if expected_job is not None:
        unit = next((u for u in expected_job.units if u.id == data["unit_id"]), None)
        expected_protocol = _benchmark_manifest(expected_job, unit, caller, root / "query/query.ins_del.vcf.gz") if unit else None
        if protocol != expected_protocol:
            raise OutputValidationError("Sealed benchmark protocol differs from the declared study.")
    if data.get("refinement_requested", requested) != requested:
        raise OutputValidationError("Receipt refinement flag differs from the sealed protocol.")
    allowed_status = "NOT_NEEDED" if requested and kind == "bench" else "COMPLETE" if requested and kind == "refined" else "NOT_REQUESTED" if not requested and kind == "bench" else None
    # Historical receipts may lack the new field; explicit contradictory values fail.
    if allowed_status is None or data.get("refinement_status", allowed_status) != allowed_status:
        raise OutputValidationError("Receipt refinement state conflicts with its requested metric protocol.")
    if requested and kind == "bench":
        from hifivar.truvari_refinement import verify_noop_result
        benchmark = json.loads((root / "benchmark/benchmark_matrix.json").read_text(encoding="utf-8"))
        query_rows = benchmark.get("queries", [])
        if benchmark.get("status") != "COMPLETE" or len(query_rows) != 1:
            raise OutputValidationError("Native no-op benchmark evidence is incomplete.")
        row = query_rows[0]
        if row.get("status") != "COMPLETE":
            raise OutputValidationError("Native no-op query row did not complete.")
        native_result = row.get("result", {})
        command = native_result.get("command")
        records = benchmark.get("commands")
        matching = [r for r in records if isinstance(r, dict) and r.get("command") == command] if isinstance(records, list) else []
        if len(matching) != 1 or row.get("command") != command:
            raise OutputValidationError("Native no-op invocation lacks one matching command record.")
        record = matching[0]
        execution = (native_result.get("refinement_proof") or {}).get("execution", {})
        if record.get("executed") is not True or type(record.get("return_code")) is not int or record["return_code"] != 0 or record.get("stderr_path") != execution.get("stderr_path"):
            raise OutputValidationError("Native no-op command record conflicts with its current execution proof.")
        expected_request = {"sample_id": protocol["sample_id"], "reference": protocol["reference"],
                            "base_vcf": protocol["truth"]["vcf"], "comparison_vcf": protocol["queries"][0]["vcf"],
                            "output_directory": str(root / "benchmark/queries" / caller),
                            "confident_regions": protocol["truth"]["confident_bed"], "thresholds": protocol["thresholds"]}
        proof = verify_noop_result(row.get("result", {}), expected_request=expected_request)
        if data.get("refinement_status") != "NOT_NEEDED" or data.get("refinement_proof") != proof or row.get("refinement_proof") != proof:
            raise OutputValidationError("Receipt no-op status/proof differs from native evidence.")
        if row.get("refined_metrics") != [] or row.get("refinement_status") != "NOT_NEEDED":
            raise OutputValidationError("Native no-op falsely declares refined metrics.")
    elif data.get("refinement_status") == "NOT_NEEDED":
        raise OutputValidationError("No-op status cannot replace a different metric protocol.")
    summary_name = "refine.variant_summary.json" if kind == "refined" else "summary.json"
    summary = root / "benchmark" / "queries" / caller / summary_name
    sealed = {Path(k).as_posix() for k in data["output_checksums"]}
    if summary.relative_to(root).as_posix() not in sealed:
        raise OutputValidationError("Native metric summary is not sealed.")
    native_metrics = {m["name"]: m["value"] for m in _metrics(summary, refined=kind == "refined")}
    if data.get("metrics") != native_metrics:
        raise OutputValidationError("Receipt metrics differ from canonical native Truvari summary.")
    if data.get("effective_summary_kind", kind) != kind:
        raise OutputValidationError("Receipt effective summary kind differs.")
    if data.get("effective_summary_path", str(summary)) != str(summary) or data.get("effective_summary_sha256", compute_sha256(summary)) != compute_sha256(summary):
        raise OutputValidationError("Receipt effective summary source changed.")
    return data


class _TimedRunner:
    """Apply an explicit caller resource timeout through CommandRunner only."""
    def __init__(self, runner, seconds):
        self.runner, self.seconds = runner, float(seconds)

    def require_executable(self, executable, **kwargs):
        return self.runner.require_executable(executable, **kwargs)

    def run(self, command, **kwargs):
        probe = len(command) > 1 and command[1] in {"--version", "version", "-v"}
        limit = min(self.seconds, 30.0) if probe else self.seconds
        requested = kwargs.get("timeout")
        kwargs["timeout"] = min(float(requested), limit) if requested is not None else limit
        return self.runner.run(command, **kwargs)


def _write_metrics(root, report):
    columns = ("unit_id", "caller", "status", "nominal_depth", "observed_coverage", "replicate",
               "sampling_fraction", "replicate_kind", "coverage_reference_bases", "tp_base", "tp_call", "fp", "fn",
               "precision", "recall", "f1", "summary_kind", "benchmark_report")
    path = root / "sv_depth_metrics.tsv"
    partial = root / ("sv_depth_metrics." + uuid4().hex + ".partial")
    with partial.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for entry in report["runs"]:
            row = {k: entry.get(k, "") for k in columns}
            if entry["status"] == "COMPLETE":
                row.update(entry.get("metrics", {}))
            writer.writerow(row)
    partial.replace(path)


@contextmanager
def _study_session(root, job, resume):
    """Hold exclusive ownership before reading any mutable study state."""
    created = False
    owner = root / "study_owner.json"
    identity = {"schema_version": 1, "study_id": job.study_id, "manifest_sha256": job.manifest_sha256}
    if root.exists():
        if not resume:
            raise OutputValidationError("Output was created by another writer; explicit resume is required.")
        try:
            current = json.loads(owner.read_text(encoding="utf-8"))
        except (OSError, ValueError) as error:
            raise OutputValidationError("Resume requires an immutable study ownership marker.") from error
        if current != identity:
            raise OutputValidationError("Study ownership identity differs from this manifest.")
    else:
        if compute_sha256(job.manifest) != job.manifest_sha256:
            raise SourceImmutabilityError("Manifest changed after parsing.")
        root.mkdir(parents=True, exist_ok=False)
        created = True
    lock = root / "study.lock"
    token = uuid4().hex
    try:
        with lock.open("x", encoding="utf-8") as handle:
            json.dump({"token": token, "started_at": utc_now_iso8601()}, handle)
    except FileExistsError as error:
        raise OutputValidationError("Study lock already exists; no concurrent writer or automatic stale-lock removal is allowed.") from error
    try:
        if created:
            write_json_atomic(identity, owner)
        yield created
    finally:
        try:
            if json.loads(lock.read_text(encoding="utf-8")).get("token") != token:
                raise OutputValidationError("Study lock ownership changed; lock preserved.")
            lock.unlink()
        except (OSError, ValueError) as error:
            raise OutputValidationError("Unable to release the owned study lock safely.") from error


def run_sv_depth_study(manifest: str | Path, output_dir: str | Path, *, dry_run: bool = False,
                       resume: bool = False, unit_ids: Sequence[str] | None = None,
                       callers: Sequence[str] | None = None, runner: CommandRunner | None = None) -> dict[str, Any]:
    """Run selected depth/caller units; no failed/incomplete unit is retried."""
    if type(dry_run) is not bool or type(resume) is not bool:
        raise InputValidationError("dry_run and resume must be booleans.")
    job = load_sv_depth_study(manifest)
    ids = _selection(unit_ids, [u.id for u in job.units], "unit_ids")
    chosen_callers = _selection(callers, job.caller_order, "callers")
    selected = [u for u in job.units if u.id in ids]
    shared = (job.manifest, job.reference.fasta, job.reference.fai, job.truth["vcf"],
              Path(str(job.truth["vcf"]) + ".tbi"), job.truth["confident_bed"])
    inputs = shared + tuple(p for u in selected for p in (u.bam, u.index, u.coverage_json))
    root = _root(output_dir, inputs, resume)
    if root.exists() and (root / "study.lock").exists():
        raise OutputValidationError("Study lock exists; concurrent or stale attempts require manual review, not automatic removal.")
    entries = []
    for u in job.units:
        for caller in job.caller_order:
            _, wrapper, request = _request(job, u, caller, root / "units" / u.id / caller, runner)
            entries.append({"id": u.id + "/" + caller, "unit_id": u.id, "caller": caller,
                "nominal_depth": u.nominal_depth, "observed_coverage": u.observed_coverage, "replicate": u.replicate,
                "coverage_reference_bases": u.coverage_reference_bases, "coverage_sha256": u.coverage_sha256,
                "sampling_fraction": u.sampling_fraction, "replicate_kind": u.replicate_kind,
                "status": "PLANNED" if dry_run else "PENDING", "command_plan": _commands(wrapper, request)})
    planned = {"schema_version": 1, "study_id": job.study_id, "sample_id": job.sample_id,
        "hifivar_version": __version__, "manifest_sha256": job.manifest_sha256,
        "status": "PLANNED" if dry_run else "RUNNING", "created_at": utc_now_iso8601(),
        "caller_order": list(job.caller_order), "protocol": _benchmark_manifest(job, selected[0], chosen_callers[0], Path("NOT_GENERATED")),
        "design": dict(job.design), "runs": entries, "unit_input_checksums": {},
        "coverage_source": "PREEXISTING_MEASURED_JSON_BAM_PATH_BOUND_NOT_REMEASURED",
        "coverage_scope": "RECORDED_REFERENCE_DENOMINATOR_NOT_BENCHMARK_BED_LOCAL_DEPTH",
        "nested_membership": "DECLARED_FROM_EXISTING_STUDY_NOT_REESTABLISHED",
        "accuracy_scope": "NATIVE_TRUVARI_5_4_DECLARED_TRUTH_RULES_WITH_INS_DEL_QUERY",
        "truth_prefilter": "NONE_ORIGINAL_MIXED_TRUTH_PRESERVED",
        "per_command_timeout_policy": "caller.runtime_minutes*60; version probes <=30s; not whole-study or memory enforcement",
        "support_threshold_comparison": "NOT_RUN_SEPARATE_EXPLICIT_ANALYSIS",
        "minimum_depth_claim": "NOT_ESTABLISHED", "sequencing_cost_reduction_claim": "NOT_ESTABLISHED",
        "biological_sample_count": 1, "intermediate_signature_checksums": "NOT_SEALED_RETAINED_ONLY", "source_immutability": "NOT_CHECKED_DRY_RUN" if dry_run else "PENDING"}
    if dry_run:
        planned["selected_run_ids"] = [r["id"] for r in entries if r["unit_id"] in ids and r["caller"] in chosen_callers]
        return planned
    # Resolve the pinned native classifier before launching any costly caller.
    from hifivar.sv_benchmark_selection import _load_native_runtime
    _load_native_runtime()
    with _study_session(root, job, resume) as created:
        if not created:
            try:
                report = json.loads((root / "study_report.json").read_text(encoding="utf-8"))
            except (OSError, ValueError) as error:
                raise OutputValidationError("Resume requires an owned valid study_report.json.") from error
            if report.get("study_id") != job.study_id or report.get("manifest_sha256") != job.manifest_sha256:
                raise OutputValidationError("Resume manifest/study identity changed.")
            if (not isinstance(report.get("runs"), list)
                    or {r.get("id") for r in report["runs"]} != {r["id"] for r in entries}
                    or len(report["runs"]) != len(entries)):
                raise OutputValidationError("Resume study matrix membership changed.")
            planned_by_id = {r["id"]: r for r in entries}
            for row in report["runs"]:
                expected = planned_by_id[row["id"]]
                if any(row.get(k) != expected[k] for k in ("unit_id", "caller", "nominal_depth", "observed_coverage", "replicate", "sampling_fraction", "replicate_kind", "coverage_reference_bases", "coverage_sha256")):
                    raise OutputValidationError("Resume study unit metadata changed.")
            _unchanged(report["shared_input_checksums"], full=True)
            for row in report["runs"]:
                if row["status"] == "COMPLETE":
                    _verify_receipt(root / "units" / row["unit_id"] / row["caller"], row, job.manifest_sha256, expected_job=job)
            for saved in report["unit_input_checksums"].values():
                _unchanged(saved)
        else:
            report = planned
            snapshot = _snapshot(shared)
            if snapshot[str(job.manifest)]["sha256"] != job.manifest_sha256:
                raise SourceImmutabilityError("Manifest changed after parsing.")
            if snapshot[str(job.reference.fai)]["sha256"] != job.reference_fai_sha256:
                raise SourceImmutabilityError("FAI changed between parsing and input hashing.")
            report["shared_input_checksums"] = snapshot
        report["selected_run_ids"] = [row["id"] for row in report["runs"]
                                      if row["unit_id"] in ids and row["caller"] in chosen_callers]
        def persist():
            write_json_atomic(report, root / "study_report.json", overwrite=True, artifact_name="SV depth study")
            _write_metrics(root, report)
        persist()
        selected_unit_snapshots = {}
        actual_runner = runner or CommandRunner()
        try:
            for unit in selected:
                unit_paths = (unit.bam, unit.index, unit.coverage_json)
                old = report["unit_input_checksums"].get(unit.id)
                if old is None:
                    old = _snapshot(unit_paths)
                    report["unit_input_checksums"][unit.id] = old
                else:
                    _unchanged(old, full=True)
                if old[str(unit.coverage_json)]["sha256"] != unit.coverage_sha256:
                    raise SourceImmutabilityError("Coverage changed between parsing and input hashing.")
                selected_unit_snapshots[unit.id] = old
                preflight = None
                for entry in report["runs"]:
                    if entry["unit_id"] != unit.id or entry["caller"] not in chosen_callers:
                        continue
                    _unchanged(report["shared_input_checksums"]);_unchanged(old)
                    caller = entry["caller"]
                    target = root / "units" / unit.id / caller
                    if entry["status"] == "COMPLETE":
                        _verify_receipt(target, entry, job.manifest_sha256, expected_job=job)
                        continue
                    if entry["status"] != "PENDING" or target.exists():
                        entry["status"] = "FAILED_PRESERVED_NO_AUTOMATIC_RETRY"
                        persist()
                        continue
                    target.mkdir(parents=True, exist_ok=False)
                    timeout_seconds = getattr(job.settings, caller + "_resources").runtime_minutes * 60
                    entry.update(status="RUNNING", per_command_timeout_seconds=timeout_seconds,
                                 version_probe_timeout_seconds=min(30, timeout_seconds))
                    timed_runner = _TimedRunner(actual_runner, timeout_seconds)
                    persist()
                    try:
                        if preflight is None:
                            preflight = IndexedBamPreflight(timed_runner, job.tools["samtools_executable"]).run(
                                unit.bam, unit.index, job.sample_id, job.reference, target / "logs" / "preflight")
                        artifact, caller_result = _run_caller(job, unit, caller, target, timed_runner)
                        _unchanged(old);_unchanged(report["shared_input_checksums"])
                        derivation = derive_sv_benchmark_query(artifact.vcf_path, target / "query", sample_id=job.sample_id)
                        query = Path(derivation["output_paths"]["vcf"])
                        bench_manifest = target / "benchmark_manifest.json"
                        write_json_atomic(_benchmark_manifest(job, unit, caller, query), bench_manifest)
                        bench = run_sv_benchmark(bench_manifest, target / "benchmark", runner=timed_runner)
                        if bench["status"] != "COMPLETE" or len(bench["queries"]) != 1 or bench["queries"][0]["status"] != "COMPLETE":
                            raise OutputValidationError("Caller benchmark did not complete; raw outputs and benchmark failure retained.")
                        row = bench["queries"][0]
                        kind = "refined" if job.thresholds.refine else "bench"
                        if job.thresholds.refine and row.get("refinement_status") == "NOT_NEEDED":
                            from hifivar.truvari_refinement import verify_noop_result
                            if verify_noop_result(row.get("result", {})) != row.get("refinement_proof"):
                                raise OutputValidationError("Benchmark no-op proof is inconsistent.")
                            kind = "bench"
                        metrics = {m["name"]: m["value"] for m in row[kind + "_metrics"]}
                        if not {"tp_base", "tp_call", "fp", "fn", "precision", "recall", "f1"} <= set(metrics):
                            raise OutputValidationError("Benchmark lacks separately named truth/query TP counts and rates.")
                        _unchanged(old);_unchanged(report["shared_input_checksums"])
                        entry.update(status="COMPLETE", metrics=metrics, summary_kind=kind,
                            benchmark_report=str(target / "benchmark" / "benchmark_matrix.json"),
                            caller_provenance=caller_result, derivation=derivation, bam_preflight=preflight,
                            precision_formula="TP_COMP/(TP_COMP+FP)", recall_formula="TP_BASE/(TP_BASE+FN)",
                            refinement_requested=job.thresholds.refine,
                            refinement_status=row.get("refinement_status", "COMPLETE" if kind == "refined" else "NOT_REQUESTED"),
                            refinement_proof=row.get("refinement_proof"),
                            effective_summary_kind=kind,
                            effective_summary_path=str(target / "benchmark/queries" / caller / ("refine.variant_summary.json" if kind == "refined" else "summary.json")),
                            effective_summary_sha256=compute_sha256(target / "benchmark/queries" / caller / ("refine.variant_summary.json" if kind == "refined" else "summary.json")))
                        receipt = {**entry, "manifest_sha256": job.manifest_sha256,
                            "output_checksums": _output_checksums(target), "completed_at": utc_now_iso8601()}
                        write_json_atomic(receipt, target / "unit_receipt.json")
                    except SourceImmutabilityError:
                        raise
                    except (HiFiVarError, OSError, ValueError, TypeError, KeyError) as error:
                        entry.update(status="FAILED", error_type=type(error).__name__, error=str(error),
                                     failure_policy="PRESERVE_NO_AUTOMATIC_RETRY")
                        write_json_atomic({**entry, "manifest_sha256": job.manifest_sha256}, target / "failure.json")
                    persist()
            _unchanged(report["shared_input_checksums"], full=True)
            for snapshot in selected_unit_snapshots.values():_unchanged(snapshot, full=True)
            for row in report["runs"]:
                if row["status"] == "COMPLETE":
                    _verify_receipt(root / "units" / row["unit_id"] / row["caller"], row, job.manifest_sha256, expected_job=job)
            for saved in report["unit_input_checksums"].values():
                _unchanged(saved)
            selected_rows = [r for r in report["runs"] if r["unit_id"] in ids and r["caller"] in chosen_callers]
            report["selection_status"] = "COMPLETE" if all(r["status"] == "COMPLETE" for r in selected_rows) else "FAILED"
            report["status"] = ("COMPLETE" if all(r["status"] == "COMPLETE" for r in report["runs"])
                                else "PARTIAL" if any(r["status"].startswith("FAILED") for r in report["runs"]) else "INCOMPLETE")
            report["source_immutability"] = "PASS_FULL_SHA256_SELECTED_INPUTS_BEFORE_AFTER_STAT_BETWEEN_CALLERS"
        except (SourceImmutabilityError, OSError) as error:
            report.update(status="FAILED", selection_status="FAILED", source_immutability="FAIL",
                          metrics_valid_for_claims=False, error_type=type(error).__name__, error=str(error))
            for entry in report["runs"]:
                if entry["status"] == "COMPLETE":entry["status"] = "INVALIDATED_SOURCE_MUTATION"
        except OutputValidationError as error:
            report.update(status="FAILED", selection_status="FAILED", output_validation="FAIL",
                          metrics_valid_for_claims=False, error_type=type(error).__name__, error=str(error))
            for entry in report["runs"]:
                if entry["status"] == "COMPLETE":entry["status"] = "INVALIDATED_OUTPUT_INTEGRITY"
        report["completed_at"] = utc_now_iso8601()
        persist()
        return report



__all__ = ["SVDepthUnit", "SVDepthStudy", "load_sv_depth_study", "run_sv_depth_study"]
