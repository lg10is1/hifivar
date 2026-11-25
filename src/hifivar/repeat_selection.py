"""Explicit TRGT INFO/TRID selection without rewriting raw VCF IDs or alleles."""
from __future__ import annotations

from dataclasses import dataclass
import gzip
import hashlib
from pathlib import Path
from typing import Sequence

from hifivar.exceptions import InputValidationError
from hifivar.reference import ReferenceGenome


PREFIX = "TRID:"


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _lines(path: Path):
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    with opener(path, "rb") as stream:
        yield from stream


@dataclass(frozen=True, slots=True)
class SelectedRepeatRecord:
    """One selected raw record; repeat coordinates exclude the VCF left anchor."""

    source_vcf: Path
    source_sha256: str
    ordinal: int
    record_sha256: str
    selector: str
    raw_vcf_id: str
    trid: str
    contig: str
    anchor_position: int
    start: int
    end: int

    def to_dict(self) -> dict[str, object]:
        return {"source_vcf": str(self.source_vcf), "source_sha256": self.source_sha256,
                "source_record_key": f"sha256:{self.source_sha256}:record:{self.ordinal}",
                "record_ordinal": self.ordinal, "selected_record_sha256": self.record_sha256,
                "selection_namespace": "TRID", "selector": self.selector,
                "raw_vcf_id": self.raw_vcf_id, "trid": self.trid,
                "contig": self.contig, "anchor_position": self.anchor_position,
                "start": self.start, "end": self.end,
                "coordinate_semantics": "TRGT_REPEAT_1_BASED_INCLUSIVE_EXCLUDING_LEFT_ANCHOR",
                "raw_vcf_modified": False, "annotation_is_call_confidence": False}


def read_selected_repeat_records(source_vcf: Path, selectors: Sequence[str], *,
                                 reference: ReferenceGenome, sample_id: str) -> tuple[SelectedRepeatRecord, ...]:
    """Stream unique explicit TRID selectors, retaining only named records.

    TRGT's repeat occupies (POS, END] in one-based reference coordinates;
    POS itself is the preceding VCF anchor. No ordinary-ID fallback occurs.
    """
    from hifivar.annotation import audit_annotation_vcf_identity

    if isinstance(selectors, (str, bytes)) or not isinstance(selectors, Sequence):
        raise InputValidationError("TRGT selectors must be an explicit sequence, not a string.")
    selected = tuple(selectors)
    if (not selected or len(selected) > 4096
            or any(not isinstance(x, str) or not x.startswith(PREFIX)
                   or not x[len(PREFIX):] or x[len(PREFIX):] == "."
                   or any(c.isspace() or c in ";,=" for c in x[len(PREFIX):]) for x in selected)
            or len(set(selected)) != len(selected)):
        raise InputValidationError("TRGT selection requires 1-4096 unique explicit TRID:<locus> selectors.")
    source = Path(source_vcf).expanduser().absolute()
    selected_set = set(selected)
    audit_annotation_vcf_identity(source, sample_id, reference)
    lengths = {c.name: c.length for c in reference.contigs}
    try:
        before = _sha(source)
        found = {}
        ordinal = 0
        for raw in _lines(source):
            if raw.startswith(b"#") or not raw.strip():
                continue
            ordinal += 1
            fields = raw.decode("utf-8").rstrip("\r\n").split("\t")
            if len(fields) < 8:
                raise InputValidationError("Malformed TRGT VCF record.")
            entries = [x.split("=", 1) for x in fields[7].split(";") if "=" in x]
            trids = [v for k, v in entries if k == "TRID"]
            if not any(PREFIX + v in selected_set for v in trids):
                continue
            if len(trids) != 1 or len({k for k, _ in entries}) != len(entries):
                raise InputValidationError("Selected TRGT INFO keys/TRID are ambiguous.")
            info = dict(entries)
            selector = PREFIX + info["TRID"]
            if selector in found:
                raise InputValidationError("Selected TRID is ambiguous: duplicate source records.")
            if not info.get("MOTIFS") or not info.get("STRUC"):
                raise InputValidationError("Selected TRID requires explicit TRGT MOTIFS and STRUC.")
            anchor, end = int(fields[1]), int(info["END"])
            if fields[0] not in lengths or not 1 <= anchor < end <= lengths[fields[0]]:
                raise InputValidationError("Selected TRGT repeat interval/reference contig is invalid.")
            if len(fields[3]) != end - anchor + 1 or any(x not in "ACGTNacgtn" for x in fields[3]):
                raise InputValidationError("Selected TRGT REF must include the declared repeat and one left anchor.")
            found[selector] = SelectedRepeatRecord(source, before, ordinal, hashlib.sha256(raw).hexdigest(),
                selector, fields[2], info["TRID"], fields[0], anchor, anchor + 1, end)
        if set(found) != set(selected):
            raise InputValidationError("Explicit selected TRID is absent from VCF; ordinary ID fallback is disabled.")
        if _sha(source) != before:
            raise InputValidationError("TRGT source changed during selection.")
        return tuple(found[x] for x in selected)
    except (OSError, EOFError, UnicodeError, KeyError, ValueError) as error:
        raise InputValidationError(f"Cannot inspect selected TRGT source: {error}") from error


def revalidate_repeat_record(record: SelectedRepeatRecord) -> None:
    """Reject changed sources before publishing overlap results."""
    try:
        if _sha(record.source_vcf) != record.source_sha256:
            raise InputValidationError("Selected TRGT source bytes changed after selection.")
    except OSError as error:
        raise InputValidationError(f"Cannot revalidate selected TRGT source: {error}") from error


def revalidate_repeat_records(records: Sequence[SelectedRepeatRecord]) -> None:
    """Hash each distinct source identity once per validation boundary."""
    seen: set[tuple[Path, str]] = set()
    for record in records:
        key = (record.source_vcf, record.source_sha256)
        if key not in seen:
            revalidate_repeat_record(record)
            seen.add(key)


__all__ = ["SelectedRepeatRecord", "read_selected_repeat_records", "revalidate_repeat_record", "revalidate_repeat_records"]
