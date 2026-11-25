"""Count native refined BD assignments, never initial bench-assigned VCFs."""

from __future__ import annotations

from collections import Counter, defaultdict
import json
from pathlib import Path
from typing import Any

from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.serialization import utc_now_iso8601, write_json_atomic
from hifivar.truvari import parse_truvari_refined_summary
from hifivar.sv_benchmark_selection import (
    _checked_fields, _check_native_header, _check_index, _fresh_output, _header,
    _input, _load_native_runtime, _native_type, _open_raw, _raw_rows, _record,
    _runtime_paths, _runtime_snapshot, _sample_id, _unchanged,
)


DEFAULT_SIZE_BOUNDARIES = (30, 50, 100, 500, 1000, 10000, 50001)
_ROLES = ("TP-base", "TP-comp", "FN", "FP")


def _rates(values: dict[str, int]) -> dict[str, Any]:
    """Match Truvari5.4 performance_metrics, including its null convention."""
    counts = {role: values.get(role, 0) for role in _ROLES}
    base = counts["TP-base"] + counts["FN"]
    comp = counts["TP-comp"] + counts["FP"]
    precision = recall = f1 = None
    if base and comp:
        precision, recall = counts["TP-comp"] / comp, counts["TP-base"] / base
        if precision + recall:
            f1 = 2 * precision * recall / (precision + recall)
    return {"counts": {**counts, "base cnt": base, "comp cnt": comp},
            "precision": precision, "recall": recall, "f1": f1}


def _bin(size: int, boundaries: tuple[int, ...]) -> str:
    lower = 0
    for upper in boundaries:
        if size < upper:
            return f"{lower}-{upper - 1}bp"
        lower = upper
    return f">={lower}bp"


def _json(path: Path) -> dict[str, Any]:
    try:
        result = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise InputValidationError(f"Invalid native protocol JSON: '{path}': {error}") from error
    if not isinstance(result, dict):
        raise InputValidationError(f"Native protocol JSON must be an object: '{path}'.")
    return result


def _assignment_header(header: Any, sample: str) -> None:
    _check_native_header(header, sample)
    for name in ("BD", "BK"):
        if name not in header.formats or header.formats[name].number != 1 or header.formats[name].type != "String":
            raise InputValidationError(f"Refined VCF requires FORMAT/{name},Number=1,Type=String.")
    if "IsRefined" not in header.info or header.info["IsRefined"].number != 0 or header.info["IsRefined"].type != "Flag":
        raise InputValidationError("Refined VCF requires native INFO/IsRefined provenance declaration.")


def stratify_refined_sv_benchmark(
    benchmark_dir: str | Path, output_dir: str | Path, *, sample_id: str,
    size_boundaries: tuple[int, ...] = DEFAULT_SIZE_BOUNDARIES,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Read both native refined assignment files and reconcile before publishing.

    Sample names are read explicitly from params.bSample/cSample and checked
    against the requested biological sample and each single-sample header.
    Strata classify each truth/query side independently; they are not matched
    record pairs. A 30--49 bp query TP therefore remains visible.
    """
    sample_id = _sample_id(sample_id)
    if not isinstance(dry_run, bool):
        raise InputValidationError("dry_run must be boolean.")
    if not isinstance(size_boundaries, (tuple, list)):
        raise InputValidationError("Size boundaries must be an explicit sequence of integers.")
    boundaries = tuple(size_boundaries)
    if (any(type(x) is not int or x <= 0 for x in boundaries)
            or boundaries != tuple(sorted(set(boundaries))) or 30 not in boundaries or 50 not in boundaries):
        raise InputValidationError("Size boundaries must be positive, unique, ascending integers including30 and50.")
    directory = Path(benchmark_dir).resolve()
    base = _input(directory / "refine.base.vcf.gz")
    comp = _input(directory / "refine.comp.vcf.gz")
    params_path = _input(directory / "params.json")
    summary_path = _input(directory / "refine.variant_summary.json")
    indexes = (_input(Path(f"{base}.tbi")), _input(Path(f"{comp}.tbi")))
    params = _json(params_path)
    sample_map = {"base": params.get("bSample"), "comparison": params.get("cSample")}
    if any(value != sample_id for value in sample_map.values()):
        raise InputValidationError("Native params.bSample/cSample must explicitly match the requested sample; no inferred alias mapping.")
    if params.get("dup_to_ins", False) is not False:
        raise InputValidationError("Refined stratification does not silently adopt DUP-to-INS conversion.")
    parse_truvari_refined_summary(summary_path)
    summary = _json(summary_path)
    paths = (base, comp, *indexes, params_path, summary_path)
    root = _fresh_output(output_dir, paths)
    for path in (base, comp):
        with _open_raw(path) as raw:
            _header(raw, sample_id)
    report: dict[str, Any] = {
        "schema_version": 1, "status": "PLANNED" if dry_run else "RUNNING",
        "created_at": utc_now_iso8601(), "sample_id": sample_id,
        "benchmark_directory": str(directory), "sample_mapping": sample_map,
        "source_kind": "NATIVE_REFINED_BASE_COMP_BD_ASSIGNMENTS",
        "size_boundaries": list(boundaries), "native_record_pairing_assumed": False,
        "output_paths": {"report": str(root / "refined_stratification.json")},
        "strata_are_post_assignment_descriptive": True,
        "scientific_scope": "Declared native refined assignments only, not all SV types or clinical validation.",
        "global_reconciliation": "NOT_RUN", "artifact_finalization": "NOT_FINALIZED",
        "public_release_modified": False,
    }
    if dry_run:
        return report
    runtime = _load_native_runtime()
    paths += _runtime_paths(runtime)
    _fresh_output(root, paths)
    before = _runtime_snapshot(runtime, paths)
    report.update(classifier=runtime.identity, input_checksums=before)
    root.mkdir(parents=True, exist_ok=False)
    output = root / "refined_stratification.json"
    write_json_atomic(report, output)
    global_counts: Counter[str] = Counter()
    strata: dict[str, Counter[str]] = defaultdict(Counter)
    provenance: Counter[str] = Counter()
    try:
        for side, path in (("base", base), ("comparison", comp)):
            mapped_sample = sample_map[side]
            _check_index(runtime, path, mapped_sample)
            with runtime.pysam.VariantFile(str(path)) as native, _open_raw(path) as raw:
                _assignment_header(native.header, mapped_sample)
                _header(raw, mapped_sample)
                records = iter(native)
                for line in _raw_rows(raw):
                    fields = _checked_fields(line, native.header)
                    parsed = next(records, None)
                    if parsed is None:
                        raise InputValidationError("Native/raw refined streams disagree.")
                    record = _record(runtime, parsed, fields)
                    formats, values = fields[8].split(":"), fields[9].split(":")
                    values += ["."] * (len(formats) - len(values))
                    if "BD" not in formats or "BK" not in formats:
                        raise InputValidationError("Every refined record needs explicit BD/BK sample values.")
                    bd, bk = values[formats.index("BD")], values[formats.index("BK")]
                    allowed_bd = {"TP", "FN"} if side == "base" else {"TP", "FP"}
                    if bd not in allowed_bd or (bd == "TP" and bk not in {"gm", "am"}) or (bd != "TP" and bk not in {"lm", "."}):
                        raise InputValidationError(f"Unsupported native refined {side} BD/BK: {bd}/{bk}.")
                    role = ("TP-base" if side == "base" else "TP-comp") if bd == "TP" else bd
                    # Read the explicitly mapped sample, never a guessed column0.
                    if parsed.samples[mapped_sample].get("BD") != bd:
                        raise InputValidationError("Native/raw refined BD values disagree.")
                    if (parsed.samples[mapped_sample].get("BK") or ".") != bk:
                        raise InputValidationError("Native/raw refined BK values disagree.")
                    kind = _native_type(record)
                    global_counts[role] += 1
                    strata[f"SVTYPE:{kind}"][role] += 1
                    if kind in {"INS", "DEL"}:
                        size = record.var_size()
                        if type(size) is not int or size < 0:
                            raise InputValidationError("Native INS/DEL size is not a nonnegative integer.")
                        label = _bin(size, boundaries)
                        strata[f"SIZE:{label}"][role] += 1
                        strata[f"SVTYPE_SIZE:{kind}:{label}"][role] += 1
                    else:
                        strata["SIZE:UNCLASSIFIED_NON_INS_DEL"][role] += 1
                    flag = bool(parsed.info.get("IsRefined", False))
                    provenance[f"{side}:IsRefined={str(flag).lower()}"] += 1
                    provenance[f"{side}:BK={bk}"] += 1
                if next(records, None) is not None:
                    raise InputValidationError("Native refined parser returned extra records.")
        for role in _ROLES:
            if global_counts[role] != summary[role]:
                raise OutputValidationError(f"REFINED_COUNT_RECONCILIATION_FAIL: {role} assignments={global_counts[role]}, summary={summary[role]}.")
        totals = _rates(global_counts)
        if totals["counts"]["base cnt"] != summary["base cnt"] or totals["counts"]["comp cnt"] != summary["comp cnt"]:
            raise OutputValidationError("Refined truth/query totals do not reconcile.")
        after = _unchanged(before)
        report.update(status="COMPLETE", global_reconciliation="PASS_EXACT_NATIVE_COUNTS",
                      global_metrics=totals, native_summary_rates={name: summary[name] for name in ("precision", "recall", "f1")},
                      strata=[{"stratum": name, **_rates(values)} for name, values in sorted(strata.items())],
                      assignment_provenance_counts=dict(sorted(provenance.items())),
                      input_checksums_after=after, source_immutability="PASS_FULL_SHA256_AND_STAT_BEFORE_AFTER",
                      artifact_finalization="FINALIZED_NEW_DESCRIPTIVE_STRATIFICATION",
                      semantics={"BK_gm": "Alternate-allele dosage match, not phase/clinical confidence.",
                                 "BK_am": "Allele match, not a confidence score.",
                                 "IsRefined": "Decision-path provenance flag, not confidence.",
                                 "size_bins": "Side-specific native sizes; matched30-49bp queryTP retained; no additional size filtering.",
                                 "rates": "Native NULL when either denominator is zero; F1NULL whenP+R=0.",
                                 "regions": "No new region filtering; assignments inherit native benchmark BED."})
    except (OSError, ValueError, TypeError, KeyError, IndexError, EOFError, AttributeError,
            RuntimeError, InputValidationError, OutputValidationError) as error:
        report.update(status="FAILED", global_reconciliation="FAIL_OR_NOT_ESTABLISHED",
                      artifact_finalization="NOT_FINALIZED_COMPLETE", error_type=type(error).__name__, error=str(error))
        try:
            report["input_checksums_after"] = _unchanged(before)
            report["source_immutability"] = "PASS_FULL_SHA256_AND_STAT_BEFORE_AFTER"
        except (OSError, InputValidationError, OutputValidationError) as mutation:
            report.update(source_immutability="FAIL", source_error=str(mutation))
        write_json_atomic(report, output, overwrite=True)
        raise OutputValidationError(f"Refined stratification failed; evidence: '{output}': {error}") from error
    report["completed_at"] = utc_now_iso8601()
    write_json_atomic(report, output, overwrite=True)
    return report


__all__ = ["DEFAULT_SIZE_BOUNDARIES", "stratify_refined_sv_benchmark"]
