"""Auditable INS/DEL query derivation using optional, pinned native Truvari.

This is a scope restriction, not normalization or a calling-confidence filter.
Canonical VCFs and the mixed benchmark truth are never rewritten. No external
process is launched: optional pysam APIs supply BGZF and tabix operations.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import gzip
import hashlib
import importlib
from pathlib import Path
import re
from typing import Any, BinaryIO, Iterator

from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar.serialization import utc_now_iso8601, write_json_atomic
from hifivar.validation import compute_sha256


POLICY_ID = "NATIVE_TRUVARI_5_4_INS_DEL_ONLY_V1"
NATIVE_VERSION = "5.4.0"
_SAMPLE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")


@dataclass(frozen=True)
class _NativeRuntime:
    truvari: Any
    pysam: Any
    identity: dict[str, Any]


def _load_native_runtime() -> _NativeRuntime:
    """Import optional dependencies only when an actual operation needs them."""
    try:
        native = importlib.import_module("truvari")
        pysam = importlib.import_module("pysam")
    except ImportError as error:
        raise ToolVersionError("SV protocol requires optional Truvari 5.4.0 and compatible pysam in this Python environment.") from error
    if getattr(native, "__version__", None) != NATIVE_VERSION:
        raise ToolVersionError("Native SV classifier must be exactly Truvari 5.4.0.")
    sources: dict[str, dict[str, str]] = {}
    for name in ("truvari", "truvari.variant_record", "truvari.variant_params",
                 "truvari.variant_file", "truvari.vcf2df", "truvari.utils"):
        try:
            module = importlib.import_module(name)
            path = Path(module.__file__).resolve()
        except (ImportError, AttributeError, TypeError) as error:
            raise ToolVersionError(f"Native classifier source identity is unavailable: {name}.") from error
        if path.suffix != ".py" or not path.is_file():
            raise ToolVersionError(f"Readable native Python source is required: {name}.")
        sources[name] = {"path": str(path), "sha256": compute_sha256(path)}
    return _NativeRuntime(native, pysam, {
        "truvari_version": NATIVE_VERSION, "pysam_version": str(getattr(pysam, "__version__", "UNKNOWN")),
        "sources": sources, "dup_to_ins": False, "classification": "native VariantRecord.var_type",
    })


def _sample_id(value: str) -> str:
    if not isinstance(value, str) or not _SAMPLE.fullmatch(value) or value.casefold() in {"unknown", "unknownsamplename"}:
        raise InputValidationError("An explicit non-placeholder sample_id is required.")
    return value


def _input(path: str | Path) -> Path:
    if not isinstance(path, (str, Path)) or not str(path).strip():
        raise InputValidationError("A nonempty existing input path is required.")
    result = Path(path).resolve()
    if not result.is_file() or result.stat().st_size == 0:
        raise InputValidationError(f"Input is missing or empty: '{result}'.")
    return result


def _fresh_output(path: str | Path, inputs: tuple[Path, ...]) -> Path:
    if not isinstance(path, (str, Path)) or not str(path).strip() or ".." in Path(path).parts:
        raise OutputValidationError("Output must be an explicit fresh directory without traversal.")
    supplied = Path(path)
    result = supplied.resolve()
    if supplied.is_symlink() or result.exists() or result == Path(result.anchor):
        raise OutputValidationError(f"A fresh output directory is required: '{result}'.")
    if any(result == p.resolve() or result in p.resolve().parents for p in inputs):
        raise OutputValidationError("Output cannot alias or contain an immutable input.")
    return result


def _stat(path: Path) -> tuple[int, int, int, int, int]:
    item = path.stat()
    return item.st_dev, item.st_ino, item.st_size, item.st_mtime_ns, item.st_ctime_ns


def _snapshot(paths: tuple[Path, ...]) -> dict[str, dict[str, Any]]:
    result = {}
    for path in dict.fromkeys(paths):
        before = _stat(path)
        digest = compute_sha256(path)
        if _stat(path) != before:
            raise OutputValidationError(f"Input changed while hashing: '{path}'.")
        result[str(path)] = {"stat": list(before), "sha256": digest}
    return result


def _unchanged(before: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    after = _snapshot(tuple(Path(p) for p in before))
    if after != before:
        raise OutputValidationError("SOURCE_IMMUTABILITY_FAIL: protocol input or classifier source changed.")
    return after


def _runtime_paths(runtime: _NativeRuntime) -> tuple[Path, ...]:
    return tuple(Path(value["path"]) for value in runtime.identity["sources"].values())


def _runtime_snapshot(runtime: _NativeRuntime, paths: tuple[Path, ...]) -> dict[str, dict[str, Any]]:
    before = _snapshot(paths)
    for item in runtime.identity["sources"].values():
        if before[str(Path(item["path"]))]["sha256"] != item["sha256"]:
            raise ToolVersionError("Native classifier source changed after identity capture.")
    return before


def _open_raw(path: Path) -> BinaryIO:
    return gzip.open(path, "rb") if path.name.endswith((".gz", ".bgz")) else path.open("rb")


def _raw_rows(handle: BinaryIO) -> Iterator[bytes]:
    for line in handle:
        if line.startswith(b"#"):
            raise InputValidationError("VCF header line occurs after #CHROM.")
        if not line.strip():
            raise InputValidationError("Blank VCF records are not accepted by the SV protocol.")
        yield line


def _header(handle: BinaryIO, sample_id: str) -> list[bytes]:
    header: list[bytes] = []
    for raw in handle:
        if not raw.startswith(b"#"):
            break
        header.append(raw)
        if raw.startswith(b"#CHROM\t"):
            try:
                fields = raw.rstrip(b"\r\n").decode("utf-8").split("\t")
            except UnicodeError as error:
                raise InputValidationError("VCF header must be UTF-8.") from error
            if len(fields) != 10 or fields[-1] != sample_id:
                raise InputValidationError("VCF must have exactly the explicitly mapped sample column.")
            if not header[0].startswith(b"##fileformat=VCFv"):
                raise InputValidationError("VCF fileformat header is missing.")
            return header
    raise InputValidationError("VCF lacks a valid single-sample #CHROM header.")


def _checked_fields(raw: bytes, header: Any) -> list[str]:
    try:
        fields = raw.rstrip(b"\r\n").decode("utf-8").split("\t")
    except UnicodeError as error:
        raise InputValidationError("VCF records must be UTF-8.") from error
    if len(fields) != 10 or any("\x00" in value for value in fields):
        raise InputValidationError("Malformed VCF record: exactly ten non-NUL fields required.")
    if not fields[1].isascii() or not fields[1].isdecimal():
        raise InputValidationError("Malformed VCF POS.")
    # VCF 4.3/4.4 fixed fields permit virtual telomeric positions 0 and N+1.
    # Preserve these coordinates verbatim; format eligibility is not call truth.
    if fields[0] not in header.contigs or int(fields[1]) > header.contigs[fields[0]].length + 1:
        raise InputValidationError("VCF record contig/POS is outside declared header lengths.")
    if not fields[3] or fields[3] == "." or not fields[4] or fields[4] == ".":
        raise InputValidationError("VCF REF/ALT must be present.")
    if "," in fields[4]:
        raise InputValidationError("MULTIALLELIC_NOT_SUPPORTED: explicitly split separately, never silently drop an ALT.")
    info_keys = [part.partition("=")[0] for part in fields[7].split(";")] if fields[7] != "." else []
    if len(info_keys) != len(set(info_keys)) or any(key not in header.info for key in info_keys):
        raise InputValidationError("VCF INFO is duplicate or undeclared.")
    for part in fields[7].split(";"):
        key, separator, _ = part.partition("=")
        if key in header.info and header.info[key].type == "Flag" and separator:
            raise InputValidationError("A VCF INFO Flag must not carry an assigned value.")
    formats, values = fields[8].split(":"), fields[9].split(":")
    if (len(formats) != len(set(formats)) or len(values) > len(formats)
            or "GT" not in formats or any(key not in header.formats for key in formats)):
        raise InputValidationError("VCF FORMAT/sample values are malformed or undeclared.")
    values += ["."] * (len(formats) - len(values))
    gt = values[formats.index("GT")]
    if not gt or any(a not in {"0", "1", "."} for a in re.split(r"[/|]", gt)):
        raise InputValidationError("VCF GT is invalid for one ALT.")
    if any(value not in {".", "PASS"} and value not in header.filters for value in fields[6].split(";")):
        raise InputValidationError("VCF FILTER is undeclared.")
    return fields


def _record(runtime: _NativeRuntime, parsed: Any, fields: list[str]) -> Any:
    if (parsed.contig, parsed.pos, parsed.id or ".", parsed.ref, tuple(parsed.alts or ())) != (
            fields[0], int(fields[1]), fields[2], fields[3], (fields[4],)):
        raise InputValidationError("Native and raw VCF record streams do not agree.")
    params = runtime.truvari.VariantParams(dup_to_ins=False)
    return runtime.truvari.VariantRecord(parsed, params)


def _native_type(record: Any) -> str:
    result = record.var_type()
    value = getattr(result, "name", None)
    if not isinstance(value, str) or not value:
        raise InputValidationError("Native Truvari did not return an explicit SV type.")
    return value


def _check_native_header(header: Any, sample_id: str) -> None:
    if list(header.samples) != [sample_id]:
        raise InputValidationError("Native header sample differs from explicit sample mapping.")
    if not header.contigs or any(type(value.length) is not int or value.length <= 0 for value in header.contigs.values()):
        raise InputValidationError("VCF requires declared positive contig lengths.")


def _check_index(runtime: _NativeRuntime, path: Path, sample_id: str) -> None:
    _input(Path(f"{path}.tbi"))
    with runtime.pysam.VariantFile(str(path)) as handle:
        _check_native_header(handle.header, sample_id)
    with runtime.pysam.TabixFile(str(path), index=str(Path(f"{path}.tbi"))) as index:
        # A header-only derived query legitimately has an empty indexed-contig
        # set. Do not mistake that for an unreadable index.
        if index.contigs:
            next(iter(index.fetch(index.contigs[0], 0, 1)), None)


def derive_sv_benchmark_query(
    source_vcf: str | Path, output_dir: str | Path, *, sample_id: str,
    policy_id: str = POLICY_ID, dry_run: bool = False,
) -> dict[str, Any]:
    """Write a new indexed INS/DEL-only query without changing retained rows.

    All FILTER states and short INS/DEL records are retained. Size, genotype,
    BED and pass-only decisions belong to the later native benchmark. The
    function does not accept truth inputs or produce benchmark accuracy.
    """
    sample_id = _sample_id(sample_id)
    if policy_id != POLICY_ID or not isinstance(dry_run, bool):
        raise InputValidationError("Only the frozen native INS/DEL policy and boolean dry_run are supported.")
    source = _input(source_vcf)
    index = _input(Path(f"{source}.tbi"))
    root = _fresh_output(output_dir, (source, index))
    with _open_raw(source) as stream:
        _header(stream, sample_id)
    report: dict[str, Any] = {
        "schema_version": 1, "status": "PLANNED" if dry_run else "RUNNING",
        "policy_id": policy_id, "sample_id": sample_id, "source_vcf": str(source),
        "created_at": utc_now_iso8601(), "canonical_query_modified": False,
        "truth_modified": False, "extra_query_size_filter": "NONE",
        "output_paths": {"vcf": str(root / "query.ins_del.vcf.gz"),
                         "index": str(root / "query.ins_del.vcf.gz.tbi"),
                         "manifest": str(root / "selection_manifest.json")},
        "artifact_finalization": "NOT_FINALIZED", "public_release_modified": False,
        "scientific_accuracy": "NOT_ASSESSED",
    }
    if dry_run:
        report["classification"] = "NOT_RUN_DRY_RUN_NO_DEPENDENCY_PROBE"
        return report
    runtime = _load_native_runtime()
    paths = (source, index) + _runtime_paths(runtime)
    _fresh_output(root, paths)
    before = _runtime_snapshot(runtime, paths)
    report.update(classifier=runtime.identity, input_checksums=before)
    root.mkdir(parents=True, exist_ok=False)
    manifest = root / "selection_manifest.json"
    write_json_atomic(report, manifest)
    partial = root / "query.ins_del.partial.vcf.gz"
    counts: Counter[str] = Counter()
    excluded: Counter[str] = Counter()
    raw_digest = hashlib.sha256()
    try:
        _check_index(runtime, source, sample_id)
        with runtime.pysam.VariantFile(str(source)) as native, _open_raw(source) as raw, runtime.pysam.BGZFile(str(partial), "wb") as target:
            _check_native_header(native.header, sample_id)
            header = _header(raw, sample_id)
            for line in header[:-1]:
                target.write(line)
            target.write(f'##HiFiVarBenchmarkQuerySelection=<Policy={policy_id},Canonical=false,Truvari=5.4.0,Types="INS|DEL",ExtraSizeFilter=NONE>\n'.encode())
            target.write(header[-1])
            records = iter(native)
            for ordinal, line in enumerate(_raw_rows(raw), start=1):
                fields = _checked_fields(line, native.header)
                parsed = next(records, None)
                if parsed is None:
                    raise InputValidationError("Native parser stopped before the raw VCF stream.")
                record = _record(runtime, parsed, fields)
                kind = _native_type(record)
                counts["source_records"] += 1
                if kind in {"INS", "DEL"}:
                    target.write(line)
                    raw_digest.update(line)
                    counts[f"retained_{kind}"] += 1
                else:
                    excluded[kind] += 1
            if next(records, None) is not None:
                raise InputValidationError("Native parser returned extra VCF records.")
        runtime.pysam.tabix_index(str(partial), preset="vcf", force=False)
        _check_index(runtime, partial, sample_id)
        readback = hashlib.sha256()
        retained_count = 0
        with _open_raw(partial) as stream:
            _header(stream, sample_id)
            for line in _raw_rows(stream):
                readback.update(line)
                retained_count += 1
        if readback.digest() != raw_digest.digest() or retained_count != counts["retained_INS"] + counts["retained_DEL"]:
            raise OutputValidationError("Retained VCF records changed during BGZF/index finalization.")
        after = _unchanged(before)
        if counts["source_records"] != counts["retained_INS"] + counts["retained_DEL"] + sum(excluded.values()):
            raise OutputValidationError("Selection count reconciliation failed.")
        output = Path(report["output_paths"]["vcf"])
        partial.rename(output)
        Path(f"{partial}.tbi").rename(Path(f"{output}.tbi"))
        report.update(status="COMPLETE", artifact_finalization="FINALIZED_NEW_DERIVED_COPY_ONLY",
                      counts={name: counts[name] for name in ("source_records", "retained_INS", "retained_DEL")},
                      excluded_by_native_type=dict(sorted(excluded.items())), retained_raw_record_sha256=raw_digest.hexdigest(),
                      input_checksums_after=after, source_immutability="PASS_FULL_SHA256_AND_STAT_BEFORE_AFTER",
                      derived_sha256=compute_sha256(output), index_sha256=compute_sha256(Path(f"{output}.tbi")),
                      execution_api="pysam.BGZFile/tabix_index; no external commands")
    except (OSError, ValueError, TypeError, IndexError, KeyError, EOFError, AttributeError,
            RuntimeError, InputValidationError, OutputValidationError) as error:
        report.update(status="FAILED", error_type=type(error).__name__, error=str(error),
                      counts=dict(counts), excluded_by_native_type=dict(excluded), artifact_finalization="NOT_FINALIZED_COMPLETE")
        try:
            report["input_checksums_after"] = _unchanged(before)
            report["source_immutability"] = "PASS_FULL_SHA256_AND_STAT_BEFORE_AFTER"
        except (OSError, InputValidationError, OutputValidationError) as mutation:
            report.update(source_immutability="FAIL", source_error=str(mutation))
        write_json_atomic(report, manifest, overwrite=True)
        raise OutputValidationError(f"SV query derivation failed; evidence: '{manifest}': {error}") from error
    report["completed_at"] = utc_now_iso8601()
    write_json_atomic(report, manifest, overwrite=True)
    return report


__all__ = ["POLICY_ID", "derive_sv_benchmark_query"]
