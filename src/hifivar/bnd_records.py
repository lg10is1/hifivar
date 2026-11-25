"""Streaming, reversible BND INFO derivation with explicit junction QC.

This is a narrow representation contract, not a general VCF validator. Native
records, alleles, genotypes, membership and breakpoints are never normalized.
"""
from __future__ import annotations

import base64
from collections import Counter
from contextlib import contextmanager
import csv
import gzip
import hashlib
from itertools import zip_longest
import json
from pathlib import Path
import re
import sqlite3
from typing import BinaryIO, Iterator, Mapping

from hifivar.breakends import BndPolicy, BreakendJunction, parse_breakend
from hifivar.exceptions import InputValidationError, OutputValidationError


_KEY = re.compile(r"[A-Za-z_][A-Za-z0-9_.]*\Z")
_INTEGER = re.compile(r"-?[0-9]+\Z")
_RESERVED = {"HIFIVAR_BND_REPR", "HIFIVAR_BND_LEDGER"}
_JUNCTION_COLUMNS = ("ledger_key", "record_ordinal", "record_id", "contig", "position",
    "local_side", "local_virtual", "mate_contig", "mate_position", "mate_side",
    "mate_virtual", "shape", "mate_qc_json")


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


@contextmanager
def _open_binary(path: Path) -> Iterator[BinaryIO]:
    with path.open("rb") as raw:
        magic = raw.read(2)
        raw.seek(0)
        if magic == b"\x1f\x8b":
            with gzip.GzipFile(fileobj=raw) as handle:
                yield handle
        else:
            yield raw


def _line(raw: bytes) -> tuple[list[str], bytes]:
    if raw.endswith(b"\r\n"):
        body, ending = raw[:-2], b"\r\n"
    elif raw.endswith(b"\n"):
        body, ending = raw[:-1], b"\n"
    else:
        raise InputValidationError("VCF records require an LF or CRLF terminator.")
    try:
        text = body.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise InputValidationError("VCF is not valid UTF-8.") from error
    if any(ord(char) < 32 and char != "\t" for char in text) or "\x7f" in text:
        raise InputValidationError("Unexpected control character in VCF.")
    return text.split("\t"), ending


def _attributes(text: str, prefix: str) -> dict[str, str]:
    if not text.startswith(prefix) or not text.endswith(">"):
        raise InputValidationError("Malformed structured VCF header.")
    tokens, buffer, quoted, escaped = [], [], False, False
    for char in text[len(prefix):-1]:
        if char == "," and not quoted:
            tokens.append("".join(buffer))
            buffer = []
        else:
            buffer.append(char)
            if char == '"' and not escaped:
                quoted = not quoted
        if char == "\\" and quoted and not escaped:
            escaped = True
        else:
            escaped = False
    if quoted or escaped:
        raise InputValidationError("Unclosed quote in VCF header.")
    tokens.append("".join(buffer))
    result: dict[str, str] = {}
    for token in tokens:
        key, separator, value = token.partition("=")
        if not separator or not key or not value or key in result:
            raise InputValidationError("Malformed or duplicate VCF header attribute.")
        result[key] = value
    return result


def _info(text: str) -> tuple[list[str], dict[str, str | None]]:
    if text == ".":
        return [], {}
    tokens, values = text.split(";"), {}
    for token in tokens:
        key, separator, value = token.partition("=")
        if _KEY.fullmatch(key) is None or key in values or (separator and not value):
            raise InputValidationError("Malformed or duplicate VCF INFO token.")
        if any(char.isspace() for char in token) or (separator and "=" in value):
            raise InputValidationError("Unescaped reserved character in VCF INFO token.")
        values[key] = value if separator else None
    if _RESERVED.intersection(values):
        raise InputValidationError("VCF already contains reserved BND derivation fields.")
    return tokens, values


def _numeric_info(values: Mapping[str, str | None]) -> None:
    for key in ("END", "SVLEN"):
        if key in values and values[key] != ".":
            if values[key] is None or _INTEGER.fullmatch(values[key]) is None:
                raise InputValidationError(f"Invalid scalar {key} value under the single-ALT policy.")
    if "CIPOS" in values:
        value = values["CIPOS"]
        fields = value.split(",") if isinstance(value, str) else []
        if value != "." and (len(fields) != 2 or any(_INTEGER.fullmatch(v) is None for v in fields)):
            raise InputValidationError("CIPOS must be missing or two explicit integers.")
        if len(fields) == 2 and int(fields[0]) > int(fields[1]):
            raise InputValidationError("CIPOS lower bound exceeds upper bound.")


def _header(handle: BinaryIO, writer: BinaryIO, lengths: Mapping[str, int],
            sample: str, policy: BndPolicy, native_sha256: str) -> tuple[int, set[str], set[str]]:
    definitions, contigs, formats = {}, {}, []
    for raw in handle:
        fields, ending = _line(raw)
        text = "\t".join(fields)
        if text.startswith("##fileformat="):
            formats.append(text)
        elif text.startswith("##INFO="):
            attributes = _attributes(text, "##INFO=<")
            name = attributes.get("ID")
            if not name or name in definitions or name in _RESERVED:
                raise InputValidationError("Duplicate or reserved INFO header definition.")
            if "Number" not in attributes or "Type" not in attributes:
                raise InputValidationError("INFO header lacks Number/Type.")
            definitions[name] = attributes
        elif text.startswith("##contig="):
            attributes = _attributes(text, "##contig=<")
            name, length = attributes.get("ID"), attributes.get("length")
            if (name not in lengths or name in contigs or length is None
                    or not length.isascii() or not length.isdigit() or int(length) != lengths[name]):
                raise InputValidationError("REFERENCE_CONTIG_MISMATCH: exact VCF/FAI name and length required.")
            contigs[name] = int(length)
        elif text.startswith("##HIFIVAR_BND_"):
            raise InputValidationError("Native VCF already declares BND derivation provenance.")
        if text.startswith("#CHROM"):
            if formats != ["##fileformat=VCFv4.4"] or not contigs:
                raise InputValidationError("BND policy requires native VCFv4.4 and declared contigs.")
            if fields != ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample]:
                raise InputValidationError("BND policy requires the exact single-sample native header.")
            for key, number in (("END", "1"), ("SVLEN", "A")):
                definition = definitions.get(key)
                if definition is not None and (definition["Number"] != number or definition["Type"] != "Integer"):
                    raise InputValidationError(f"Reserved {key} header has incompatible Number/Type.")
            if "CIPOS" in definitions and (definitions["CIPOS"]["Type"] != "Integer" or definitions["CIPOS"]["Number"] not in {"2", "."}):
                raise InputValidationError("CIPOS header must describe two integer bounds or a variable integer list.")
            added = [f"##HIFIVAR_BND_POLICY={policy.policy_id}:{policy.policy_version}:{policy.policy_sha256}",
                     f"##HIFIVAR_BND_NATIVE_SHA256={native_sha256}",
                     '##INFO=<ID=HIFIVAR_BND_REPR,Number=1,Type=String,Description="Explicit derived BND INFO policy; not call confidence">',
                     '##INFO=<ID=HIFIVAR_BND_LEDGER,Number=1,Type=String,Description="Mandatory native-record ledger key">']
            for item in added:
                writer.write(item.encode("ascii") + ending)
            writer.write(raw)
            return len(fields), set(definitions), set(contigs)
        if not text.startswith("##"):
            raise InputValidationError("Data or malformed header before #CHROM.")
        writer.write(raw)
    raise InputValidationError("VCF lacks #CHROM.")


def _mate_qc(db: sqlite3.Connection, row: dict[str, object]) -> dict[str, object]:
    junction = row["junction"]
    mates = row["mate_ids"]
    if not mates:
        return {"id_resolution": "NO_MATEID", "links": []}
    links = []
    for identifier in mates:
        number = db.execute("SELECT count(*) FROM bnd WHERE record_id=?", (identifier,)).fetchone()[0]
        link: dict[str, object] = {"mate_id": identifier, "matching_record_count": number}
        if number != 1:
            link["id_resolution"] = "UNRESOLVED" if not number else "AMBIGUOUS_DUPLICATE_ID"
        else:
            ordinal, payload = db.execute("SELECT ordinal, payload FROM bnd WHERE record_id=?", (identifier,)).fetchone()
            target = json.loads(payload)
            other = target["junction"]
            exact = (junction["mate_contig"], junction["mate_position"]) == (other["contig"], other["position"])
            cipos = target["cipos"]
            uncertainty = "EXACT" if exact else "NOT_ESTABLISHED"
            if (not exact and cipos is not None and junction["mate_contig"] == other["contig"]
                    and junction["mate_position"] is not None):
                uncertainty = ("WITHIN_EXPLICIT_MATE_CIPOS" if other["position"] + cipos[0] <= junction["mate_position"] <= other["position"] + cipos[1]
                               else "OUTSIDE_EXPLICIT_MATE_CIPOS")
            link.update(id_resolution="RESOLVED", target_ordinal=ordinal,
                reciprocal_id_link=row["record_id"] in target["mate_ids"],
                self_link=ordinal == row["record_ordinal"], exact_mate_coordinate=exact,
                orientation_consistent=(junction["mate_side"] is not None and other["mate_side"] is not None
                    and junction["local_side"] == other["mate_side"] and junction["mate_side"] == other["local_side"]),
                coordinate_consistency=uncertainty)
        links.append(link)
    return {"id_resolution": "PER_LINK_QC_ONLY", "links": links, "not_event_validity_or_confidence": True}


def derive_record_stream(native: Path, derived: Path, ledger: Path, junctions: Path, *,
                         reference_lengths: Mapping[str, int], sample_id: str,
                         native_sha256: str, policy: BndPolicy, work_directory: Path) -> dict[str, object]:
    """Write to exclusive owned paths; preserve failed intermediates for diagnosis.

    The caller owns the fresh directory. No existing artifact is replaced or
    removed. Mate lookup uses disk-backed SQLite rather than retaining a VCF.
    """
    spool = work_directory / "ledger.spool.jsonl"
    database = work_directory / "mates.sqlite3"
    if any(path.exists() or path.is_symlink() for path in (derived, ledger, junctions, spool, database)):
        raise OutputValidationError("BND derivation refuses existing output or intermediate paths.")
    native_digest = hashlib.sha256()
    with native.open("rb") as reader:
        for chunk in iter(lambda: reader.read(1024 * 1024), b""):
            native_digest.update(chunk)
    if native_digest.hexdigest() != native_sha256:
        raise InputValidationError("Native VCF bytes differ from the explicit derivation identity.")
    counts: Counter[str] = Counter(records=0, bnd_records=0, unchanged_non_bnd_records=0, unresolved_mate_links=0)
    with sqlite3.connect(database) as db:
        db.execute("CREATE TABLE bnd (ordinal INTEGER PRIMARY KEY, record_id TEXT, payload TEXT)")
        db.execute("CREATE INDEX bnd_id ON bnd(record_id)")
        with _open_binary(native) as reader, derived.open("xb") as writer, spool.open("x", encoding="utf-8", newline="\n") as receipt:
            columns, definitions, contigs = _header(reader, writer, reference_lengths, sample_id, policy, native_sha256)
            for ordinal, raw in enumerate(reader, 1):
                fields, ending = _line(raw)
                if len(fields) != columns or any(not value for value in fields) or fields[0].startswith("#"):
                    raise InputValidationError(f"Malformed VCF data record at ordinal {ordinal}.")
                if not re.fullmatch(r"[0-9]+", fields[1]):
                    raise InputValidationError("Invalid VCF POS token.")
                if fields[0] not in contigs:
                    raise InputValidationError("Record contig lacks a compatible native header declaration.")
                tokens, values = _info(fields[7])
                if "SVTYPE" in values and values["SVTYPE"] is None:
                    raise InputValidationError("SVTYPE cannot be a flag.")
                junction = parse_breakend(fields[0], int(fields[1]), fields[3], fields[4], reference_lengths,
                                          svtype=values.get("SVTYPE"))
                key = "bnd1_" + _digest(f"{native_sha256}:{ordinal}:{_digest(raw)}".encode("ascii"))
                row: dict[str, object] = {"record_ordinal": ordinal, "record_id": fields[2], "ledger_key": key,
                    "native_record_sha256": _digest(raw), "action": "BND_INFO_DERIVED" if junction else "UNCHANGED_NON_BND",
                    "eol": "CRLF" if ending == b"\r\n" else "LF", "junction": junction.to_dict() if junction else None}
                if junction:
                    _numeric_info(values)
                    if "MATEID" in values and values["MATEID"] is None:
                        raise InputValidationError("MATEID cannot be a flag.")
                    if any(k in values and k not in definitions for k in ("END", "SVLEN", "CIPOS")):
                        raise InputValidationError("BND numeric INFO is used without its native header definition.")
                    removed = [{"token_index": i, "raw_token": token} for i, token in enumerate(tokens)
                               if token.partition("=")[0] in {"END", "SVLEN"}]
                    kept = [token for token in tokens if token.partition("=")[0] not in {"END", "SVLEN"}]
                    row["native_info_b64"] = base64.b64encode(fields[7].encode("utf-8")).decode("ascii")
                    row["removed_info_tokens"] = removed
                    row["mate_ids"] = [] if values.get("MATEID") in {None, "."} else values["MATEID"].split(",")
                    if any(not token for token in row["mate_ids"]):
                        raise InputValidationError("Empty MATEID token.")
                    row["cipos"] = [int(v) for v in values["CIPOS"].split(",")] if values.get("CIPOS") not in {None, "."} else None
                    fields[7] = ";".join([*kept, f"HIFIVAR_BND_REPR={policy.policy_id}", f"HIFIVAR_BND_LEDGER={key}"])
                    changed = "\t".join(fields).encode("utf-8") + ending
                    db.execute("INSERT INTO bnd VALUES (?, ?, ?)", (ordinal, fields[2], json.dumps(row, sort_keys=True)))
                    counts["bnd_records"] += 1
                else:
                    changed = raw
                    counts["unchanged_non_bnd_records"] += 1
                row["derived_record_sha256"] = _digest(changed)
                writer.write(changed)
                receipt.write(json.dumps(row, sort_keys=True) + "\n")
                counts["records"] += 1
        db.commit()
        with spool.open(encoding="utf-8") as reader, ledger.open("x", encoding="utf-8", newline="\n") as writer, junctions.open("x", encoding="utf-8", newline="") as table:
            tab = csv.DictWriter(table, fieldnames=_JUNCTION_COLUMNS, delimiter="\t", lineterminator="\n")
            tab.writeheader()
            for line in reader:
                row = json.loads(line)
                if row["junction"]:
                    qc = _mate_qc(db, row)
                    row["mate_qc"] = qc
                    entry = {k: row[k] for k in ("ledger_key", "record_ordinal", "record_id")}
                    entry.update({k: v for k, v in row["junction"].items() if k in _JUNCTION_COLUMNS})
                    entry["mate_qc_json"] = json.dumps(qc, sort_keys=True, separators=(",", ":"))
                    tab.writerow(entry)
                    counts["unresolved_mate_links"] += sum(v["id_resolution"] == "UNRESOLVED" for v in qc["links"])
                writer.write(json.dumps(row, sort_keys=True) + "\n")
    return {**dict(counts), "round_trip": verify_record_round_trip(native, derived, ledger),
            "format_scope": "NARROW_BND_POLICY_NOT_FULL_VCF_CONFORMANCE",
            "mate_QC_is_biological_validation": False, "index_semantics": "LOCAL_CHROM_ONLY_NOT_ALT_MATE"}


def _data_lines(path: Path) -> Iterator[bytes]:
    with _open_binary(path) as reader:
        for raw in reader:
            if not raw.startswith(b"#"):
                yield raw


def verify_record_round_trip(native: Path, derived: Path, ledger: Path) -> dict[str, object]:
    """Independently reconstruct and compare every original data-record byte."""
    count = 0
    native_stream, reconstructed_stream = hashlib.sha256(), hashlib.sha256()
    with ledger.open(encoding="utf-8") as reader:
        for ordinal, (original, changed, line) in enumerate(zip_longest(_data_lines(native), _data_lines(derived), reader), 1):
            if original is None or changed is None or line is None:
                raise OutputValidationError("BND round-trip record/ledger counts differ.")
            row = json.loads(line)
            if row["record_ordinal"] != ordinal or _digest(changed) != row["derived_record_sha256"]:
                raise OutputValidationError("BND ledger ordinal or derived record digest differs.")
            fields, ending = _line(changed)
            if row["action"] == "BND_INFO_DERIVED":
                fields[7] = base64.b64decode(row["native_info_b64"], validate=True).decode("utf-8")
                reconstructed = "\t".join(fields).encode("utf-8") + ending
            elif row["action"] == "UNCHANGED_NON_BND":
                reconstructed = changed
            else:
                raise OutputValidationError("Unknown BND ledger action.")
            if reconstructed != original or _digest(original) != row["native_record_sha256"]:
                raise OutputValidationError("BND round-trip changed protected native record bytes.")
            native_stream.update(original)
            reconstructed_stream.update(reconstructed)
            count += 1
    return {"status": "PASS", "records": count, "native_records_sha256": native_stream.hexdigest(),
            "reconstructed_records_sha256": reconstructed_stream.hexdigest(), "includes_original_line_terminators": True}


__all__ = ["derive_record_stream", "verify_record_round_trip"]
