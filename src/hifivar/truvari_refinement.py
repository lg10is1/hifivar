"""Strict evidence for the Truvari 5.4.0 perfect-bench zero-region refine exit.

This is not a general missing-output fallback. Large input identity remains the
benchmark runner's responsibility; this module seals native output evidence.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any, Mapping, Sequence

from hifivar.benchmark import BenchmarkVariantClass
from hifivar.exceptions import OutputValidationError, ToolNotFoundError
from hifivar.validation import compute_sha256

_POLICY = "TRUVARI_5_4_0_PERFECT_BENCH_ZERO_REGION_NOOP_V1"
_ROLES = {"tp-base": "TP-base", "tp-comp": "TP-comp", "fp": "FP", "fn": "FN"}
_REFINED = ("refine.variant_summary.json", "refine.region_summary.json",
            "refine.base.vcf.gz", "refine.base.vcf.gz.tbi",
            "refine.comp.vcf.gz", "refine.comp.vcf.gz.tbi")


def _require(ok: bool, message: str) -> None:
    if not ok:
        raise OutputValidationError("Unverified Truvari refinement no-op: " + message)


def _unique(pairs: list[tuple[str, Any]]) -> dict:
    result: dict = {}
    for key, value in pairs:
        _require(key not in result, "duplicate JSON key.")
        result[key] = value
    return result


def _object(path: Path) -> dict:
    try:
        _require(path.stat().st_size <= 8 * 1024 * 1024, "JSON evidence exceeds supported size.")
        value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_unique)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise OutputValidationError("Unreadable native no-op JSON evidence.") from error
    _require(isinstance(value, dict), "native JSON evidence is not an object.")
    return value


def _absolute(value: Any) -> str:
    _require(isinstance(value, (str, Path)) and bool(str(value)), "missing path binding.")
    return str(Path(value).absolute())


def _context(request: Mapping[str, Any]) -> dict:
    from hifivar.truvari import TruvariThresholds
    try:
        thresholds = TruvariThresholds(**request["thresholds"]).to_dict()
        ref = request["reference"]["fasta"]
        sample = request["sample_id"]
        context = {"base": _absolute(request["base_vcf"]),
                   "comp": _absolute(request["comparison_vcf"]),
                   "reference": _absolute(ref),
                   "output": _absolute(request["output_directory"]),
                   "includebed": (_absolute(request["confident_regions"])
                                  if request.get("confident_regions") is not None else None),
                   "sample_id": sample, "thresholds": thresholds}
    except (KeyError, TypeError, ValueError) as error:
        raise OutputValidationError("Malformed no-op request binding.") from error
    _require(isinstance(sample, str) and bool(sample), "sample identity is absent.")
    _require(thresholds["refine"] is True, "refinement was not requested.")
    return context


def _validate_binding(context: dict, params: dict, command: Sequence[str]) -> None:
    _require(isinstance(command, (tuple, list)) and len(command) >= 2
             and all(isinstance(x, str) for x in command) and command[1] == "bench",
             "not a native bench command.")
    values: dict[str, str] = {}
    flags: set[str] = set()
    iterator = iter(command[2:])
    for option in iterator:
        _require(option not in values and option not in flags, "duplicate command option.")
        if option in {"--passonly", "--refine"}:
            flags.add(option)
        else:
            try:
                values[option] = next(iterator)
            except StopIteration as error:
                raise OutputValidationError("Incomplete no-op command.") from error
    wanted = {"-b": context["base"], "-c": context["comp"],
              "-f": context["reference"], "-o": context["output"]}
    if context["includebed"] is not None:
        wanted["--includebed"] = context["includebed"]
    threshold_flags = {"refdist": "--refdist", "pctseq": "--pctseq",
                       "pctsize": "--pctsize", "pctovl": "--pctovl",
                       "sizemin": "--sizemin", "sizemax": "--sizemax",
                       "bnddist": "--bnddist", "pick": "--pick", "chunksize": "--chunksize"}
    thresholds = context["thresholds"]
    for name, option in threshold_flags.items():
        if thresholds[name] is not None:
            wanted[option] = str(thresholds[name])
    wanted_flags = {"--refine"} | ({"--passonly"} if thresholds["pass_only"] else set())
    _require(values == wanted and flags == wanted_flags, "command differs from requested inputs or thresholds.")
    for name in ("base", "comp", "reference", "output", "includebed"):
        _require(params.get(name) == context[name], "native params path differs: " + name)
    _require(params.get("bSample") == context["sample_id"] == params.get("cSample"),
             "native base/comparison sample differs.")
    _require(params.get("skip_gt") is False, "native sample-aware refinement was not used.")
    _require(type(params.get("passonly")) is bool
             and params["passonly"] == thresholds["pass_only"], "native filter policy differs.")
    for name in threshold_flags:
        expected = thresholds[name]
        if expected is None:
            continue
        actual = params.get(name)
        if isinstance(expected, int):
            valid = type(actual) is int and actual == expected
        elif isinstance(expected, float):
            valid = (not isinstance(actual, bool) and isinstance(actual, (int, float))
                     and math.isfinite(actual) and actual == expected)
        else:
            valid = type(actual) is type(expected) and actual == expected
        _require(valid, "native threshold differs: " + name)


def _identity(path: Path) -> dict:
    _require(path.is_file() and not path.is_symlink(), "proof file is absent or symlinked: " + path.name)
    return {"bytes": path.stat().st_size, "sha256": compute_sha256(path)}


def _partition_counts(vcf: Path, index: Path, sample_id: str) -> dict:
    """Read full gzip/VCF streams and require an equivalent native Tabix stream."""
    try:
        import pysam
    except ImportError as error:
        raise ToolNotFoundError("Strict Truvari no-op validation requires the optional pysam dependency.") from error
    count, header, raw_hash = 0, None, hashlib.sha256()
    try:
        with gzip.open(vcf, "rt", encoding="utf-8", newline="") as stream:
            for line in stream:
                if line.startswith("#CHROM"):
                    _require(header is None and count == 0, "duplicate/late VCF header.")
                    header = line.rstrip("\r\n").split("\t")[9:]
                elif line.startswith("#"):
                    continue
                elif line.strip():
                    _require(header is not None and len(line.rstrip("\r\n").split("\t")) >= 10,
                             "malformed partition VCF record.")
                    count += 1
                    raw_hash.update((line.rstrip("\r\n") + "\n").encode("utf-8"))
        _require(header == [sample_id], "partition sample header differs.")
        with pysam.VariantFile(str(vcf)) as variants:
            _require(list(variants.header.samples) == [sample_id], "native VCF sample differs.")
            native_count = sum(1 for _ in variants)
        indexed_count, indexed_hash = 0, hashlib.sha256()
        with pysam.TabixFile(str(vcf), index=str(index)) as tabix:
            for contig in tabix.contigs:
                for line in tabix.fetch(contig):
                    indexed_count += 1
                    indexed_hash.update((line.rstrip("\r\n") + "\n").encode("utf-8"))
        _require(count == native_count == indexed_count
                 and raw_hash.digest() == indexed_hash.digest(), "partition index/stream counts or records differ.")
    except (OSError, ValueError, EOFError, UnicodeError) as error:
        raise OutputValidationError("Invalid native Truvari partition VCF/BGZF/TBI.") from error
    return {"records": count, "indexed_records": indexed_count,
            "record_stream_sha256": raw_hash.hexdigest()}


def _log_evidence(path: Path, params: dict, context: dict) -> None:
    _require(path.stat().st_size <= 8 * 1024 * 1024, "execution log exceeds supported no-op size.")
    try:
        text = path.read_text(encoding="utf-8")
        markers = list(re.finditer(r"Params:\s*\n", text))
        _require(len(markers) == 2, "current native bench/refine parameter logs are absent or ambiguous.")
        decoder = json.JSONDecoder(object_pairs_hook=_unique)
        bench_params, _ = decoder.raw_decode(text[markers[0].end():].lstrip())
        refine_params, _ = decoder.raw_decode(text[markers[1].end():].lstrip())
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise OutputValidationError("Unreadable current native no-op execution log.") from error
    _require(bench_params == params, "execution log params differ from native params.json.")
    _require(isinstance(refine_params, dict)
             and refine_params.get("benchdir") == context["output"]
             and refine_params.get("reference") == context["reference"]
             and refine_params.get("regions") == str(Path(context["output"]) / "candidate.refine.bed")
             and refine_params.get("use_original_vcfs") is False
             and refine_params.get("coords") == "R", "refine hook context differs.")
    first = re.search(r"\[INFO\]\s+0 regions to be refined\s*(?:\n|$)", text)
    last = re.search(r"\[CRITICAL\]\s+No regions to be refined\. Exiting\s*(?:\n|$)", text)
    _require(first is not None and last is not None
             and markers[1].start() < first.start() < last.start(),
             "both ordered native zero-region markers are required.")
    _require("Truvari v5.4.0" in text and "Finished bench" in text[:markers[1].start()]
             and "Traceback" not in text and "[ERROR]" not in text, "current native execution log is not clean/bound.")


def verify_noop_refinement(
    request: Mapping[str, Any], command: Sequence[str], *, version: str | None,
    returncode: Any, executed: Any, stderr_path: str | Path | None,
    output_absent_before: bool,
) -> dict:
    """Accept only a fully evidenced native 5.4.0 perfect-bench no-op."""
    from hifivar.truvari import _parse_summary_metrics
    _require(version == "5.4.0" and type(returncode) is int and returncode == 0
             and executed is True and output_absent_before is True,
             "requires this fresh native 5.4.0 --refine invocation with exit zero.")
    context = _context(request)
    root = Path(context["output"])
    _require(not any((root / name).exists() or (root / name).is_symlink() for name in _REFINED),
             "partial refinement outputs are not a no-op.")
    _require(not any(p.name.startswith("phab")
                     or (p.name.startswith("refine.") and p.name != "refine.log.txt")
                     for p in root.rglob("*")), "partial phab/refinement artifacts are present.")
    params_path, summary_path = root / "params.json", root / "summary.json"
    _require(stderr_path is not None, "this invocation has no persisted stderr log.")
    log = Path(stderr_path).absolute()
    files = [params_path, summary_path, root / "candidate.refine.bed", log]
    if (root / "refine.log.txt").exists():
        files.append(root / "refine.log.txt")
    before = {str(p.absolute()): _identity(p) for p in files}
    params, summary = _object(params_path), _object(summary_path)
    _validate_binding(context, params, command)
    _parse_summary_metrics(summary_path, variant_class=BenchmarkVariantClass.SV, require_totals=True)
    _require(summary["FP"] == summary["FN"] == 0
             and summary["TP-base"] > 0 and summary["TP-comp"] > 0,
             "only positive-TP perfect-bench FP=FN=0 no-ops are supported; empty/nonperfect cases remain unsupported.")
    _log_evidence(log, params, context)
    partitions = {}
    for role, field in _ROLES.items():
        vcf = root / (role + ".vcf.gz")
        index = Path(str(vcf) + ".tbi")
        before[str(vcf.absolute())], before[str(index.absolute())] = _identity(vcf), _identity(index)
        counts = _partition_counts(vcf, index, context["sample_id"])
        _require(counts["records"] == summary[field], "partition count differs from native summary: " + role)
        partitions[role] = counts
    for path, identity in before.items():
        _require(_identity(Path(path)) == identity, "native no-op evidence changed while validating.")
    return {"schema_version": 1, "policy": _POLICY, "refinement_status": "NOT_NEEDED",
            "request_context": context,
            "execution": {"command": list(command), "version": version, "returncode": 0,
                          "executed": True, "stderr_path": str(log), "output_absent_before": True},
            "proof_files": before, "partition_counts": partitions,
            "effective_summary": {"kind": "bench", "path": str(summary_path),
                                  "sha256": before[str(summary_path)]["sha256"]},
            "limitations": "PERFECT_BENCH_ONLY_NOT_ALL_NATIVE_ZERO_REGION_CASES"}


def verify_noop_result(result: Mapping[str, Any], *, expected_request: Mapping[str, Any] | None = None) -> dict:
    """Recompute the evidence on reuse; a stored status/flag never grants success."""
    _require(result.get("status") == "completed", "stored invocation did not complete.")
    _require(result.get("refinement_status") == "NOT_NEEDED"
             and result.get("refinement_requested") is True
             and result.get("effective_summary_kind") == "bench", "stored no-op status is inconsistent.")
    _require(all(result.get(key) is None for key in
                 ("refined_summary_path", "refined_region_summary_path",
                  "refined_base_vcf_path", "refined_comparison_vcf_path")),
             "stored no-op falsely claims refined outputs.")
    proof = result.get("refinement_proof")
    _require(isinstance(proof, dict) and isinstance(proof.get("execution"), dict), "sealed no-op proof is absent.")
    request = result.get("request")
    _require(isinstance(request, dict), "stored request is absent.")
    if expected_request is not None:
        _require(_context(request) == _context(expected_request), "stored request differs from the sealed study protocol.")
    execution = proof["execution"]
    actual = verify_noop_refinement(
        request, result.get("command"), version=result.get("version"),
        returncode=execution.get("returncode"), executed=execution.get("executed"),
        stderr_path=execution.get("stderr_path"),
        output_absent_before=execution.get("output_absent_before"),
    )
    _require(actual == proof, "sealed no-op proof changed or was forged.")
    source = actual["effective_summary"]
    _require(result.get("effective_summary_path") == source["path"]
             and result.get("effective_summary_sha256") == source["sha256"],
             "stored effective summary identity differs from native evidence.")
    return actual
