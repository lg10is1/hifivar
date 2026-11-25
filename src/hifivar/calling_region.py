"""Explicit single-interval calling scope; never infer or rename contigs."""
from __future__ import annotations

from dataclasses import dataclass
import re
from typing import TYPE_CHECKING

from hifivar.exceptions import InputValidationError, ReferenceError

if TYPE_CHECKING:
    from hifivar.reference import ReferenceGenome


@dataclass(frozen=True, slots=True)
class CallingRegion:
    """One 1-based, inclusive reference interval parsed without normalization."""

    contig: str
    start: int
    end: int


def parse_calling_region(value: object) -> CallingRegion | None:
    """Accept null or exactly contig:start-end, not BED paths or interval lists."""
    if value is None:
        return None
    match = re.fullmatch(r"([A-Za-z0-9_][A-Za-z0-9_.-]*):([1-9][0-9]*)-([1-9][0-9]*)", value) if isinstance(value, str) else None
    if match is None:
        raise InputValidationError(
            "small.region must be null or one contig:start-end interval "
            "(1-based inclusive); BED files, lists and implicit normalization are unsupported."
        )
    region = CallingRegion(match[1], int(match[2]), int(match[3]))
    if region.start > region.end:
        raise InputValidationError("small.region start must not exceed end.")
    return region


def validate_calling_region(value: object, reference: ReferenceGenome) -> CallingRegion | None:
    """Check exact contig identity and bounds against the configured FASTA index."""
    region = parse_calling_region(value)
    if region is None:
        return None
    length = next((c.length for c in reference.contigs if c.name == region.contig), None)
    if length is None:
        raise ReferenceError(f"REFERENCE_CONTIG_MISMATCH: small.region contig '{region.contig}' is absent from the configured reference.")
    if region.end > length:
        raise ReferenceError(f"small.region end {region.end} exceeds {region.contig} length {length}.")
    return region
