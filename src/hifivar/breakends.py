"""Explicit BND grammar and an immutable, opt-in representation policy.

Coordinates are VCF anchors, not rearrangement spans. This module neither
normalizes coordinates nor infers biological validity or missing mate links.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import hashlib
from importlib.resources import files
import json
import re

from hifivar.exceptions import ConfigurationError, InputValidationError


POLICY = json.loads(files("hifivar").joinpath("resources/policies/bnd_repr_v1.json").read_text(encoding="utf-8"))
POLICY_SHA256 = hashlib.sha256(json.dumps(POLICY, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")).hexdigest()
_DNA = re.compile(r"[ACGTNacgtn]+\Z")
_PAIRED = re.compile(r"([^\[\]]*)([\[\]])([^\[\]]+)\2([^\[\]]*)\Z")


@dataclass(frozen=True, slots=True)
class BndPolicy:
    policy_id: str
    policy_version: int
    policy_sha256: str
    mode: str
    timeout_seconds: int = 300

    def __post_init__(self) -> None:
        if (self.policy_id != POLICY["policy_id"] or type(self.policy_version) is not int
                or self.policy_version != POLICY["policy_version"]
                or self.policy_sha256 != POLICY_SHA256 or self.mode != POLICY["mode"]):
            raise ConfigurationError("BND representation requires the exact supported policy identity and mode.")
        if type(self.timeout_seconds) is not int or self.timeout_seconds <= 0:
            raise ConfigurationError("BND representation timeout_seconds must be positive.")

    @classmethod
    def from_config(cls, section: object) -> BndPolicy | None:
        if not isinstance(section, Mapping):
            raise ConfigurationError("sv.harmonization.bnd_representation must be a mapping.")
        expected = {"enabled", "policy_id", "policy_version", "policy_sha256", "mode", "timeout_seconds"}
        if set(section) - expected or any(not isinstance(k, str) for k in section):
            raise ConfigurationError("Unknown BND representation configuration key.")
        enabled = section.get("enabled", False)
        if type(enabled) is not bool:
            raise ConfigurationError("BND representation enabled must be boolean.")
        timeout = section.get("timeout_seconds", 300)
        if type(timeout) is not int or timeout <= 0:
            raise ConfigurationError("BND representation timeout_seconds must be positive.")
        if not enabled:
            if any(section.get(k) is not None for k in ["policy_id", "policy_version", "policy_sha256", "mode"]):
                raise ConfigurationError("Disabled BND representation must not declare an active policy.")
            return None
        try:
            return cls(section["policy_id"], section["policy_version"], section["policy_sha256"],
                       section["mode"], timeout)
        except KeyError as error:
            raise ConfigurationError("Enabled BND representation requires an explicit policy id/version/hash/mode.") from error

    def to_dict(self) -> dict[str, object]:
        return {"enabled": True, "policy_id": self.policy_id, "policy_version": self.policy_version,
                "policy_sha256": self.policy_sha256, "mode": self.mode, "timeout_seconds": self.timeout_seconds}


@dataclass(frozen=True, slots=True)
class BreakendJunction:
    contig: str
    position: int
    local_side: str
    local_virtual: bool
    mate_contig: str | None
    mate_position: int | None
    mate_side: str | None
    mate_virtual: bool | None
    shape: str

    def to_dict(self) -> dict[str, object]:
        return {"contig": self.contig, "position": self.position, "local_side": self.local_side,
                "local_virtual": self.local_virtual, "mate_contig": self.mate_contig,
                "mate_position": self.mate_position, "mate_side": self.mate_side,
                "mate_virtual": self.mate_virtual, "shape": self.shape,
                "coordinate_semantics": "VCF_ANCHORS_NOT_EVENT_SPAN"}


def _bounds(contig: str, position: int, lengths: Mapping[str, int]) -> bool:
    if contig not in lengths:
        raise InputValidationError(f"REFERENCE_CONTIG_MISMATCH: BND contig '{contig}'.")
    if type(position) is not int or not 0 <= position <= lengths[contig] + 1:
        raise InputValidationError("BND coordinate exceeds the exact reference's 0..N+1 bounds.")
    return position in {0, lengths[contig] + 1}


def parse_breakend(contig: str, position: int, ref: str, alt: str,
                   lengths: Mapping[str, int], *, svtype: str | None = None) -> BreakendJunction | None:
    """Parse four paired orientations or two single breakends, without repair.

    A paired virtual local anchor accepts exactly t='.' and REF='N'. A mate
    coordinate at 0/N+1 is retained; it is not an ordinary reference interval.
    Unsupported representations fail rather than being dropped or rewritten.
    """
    virtual = _bounds(contig, position, lengths)
    if not isinstance(ref, str) or _DNA.fullmatch(ref) is None:
        raise InputValidationError("Invalid BND REF sequence.")
    if not virtual and position + len(ref) - 1 > lengths[contig]:
        raise InputValidationError("Local REF anchor exceeds contig length.")
    if not isinstance(alt, str) or not alt or "," in alt:
        raise InputValidationError("BND representation v1 requires exactly one ALT; no silent splitting.")
    bracketed = "[" in alt or "]" in alt
    single = alt != "." and (alt.startswith(".") or alt.endswith("."))
    if not bracketed and not single:
        if svtype == "BND":
            raise InputValidationError("SVTYPE=BND lacks a supported breakend ALT.")
        return None
    if svtype not in {None, "BND"}:
        raise InputValidationError("SVTYPE and breakend ALT contradict each other.")
    if bracketed:
        match = _PAIRED.fullmatch(alt)
        if match is None:
            raise InputValidationError("Malformed or multiple BND mate tokens.")
        before, bracket, locus, after = match.groups()
        if bool(before) == bool(after):
            raise InputValidationError("BND replacement must occur on exactly one side.")
        replacement = before or after
        if virtual:
            if ref.upper() != "N" or replacement != ".":
                raise InputValidationError("Unsupported virtual-local BND replacement; no coordinate clamp.")
        elif _DNA.fullmatch(replacement) is None:
            raise InputValidationError("Unsupported ordinary BND replacement sequence.")
        mate_contig, separator, mate_raw = locus.rpartition(":")
        if not separator or not mate_contig or not re.fullmatch(r"[0-9]+", mate_raw):
            raise InputValidationError("Malformed BND mate coordinate.")
        mate_position = int(mate_raw)
        mate_virtual = _bounds(mate_contig, mate_position, lengths)
        return BreakendJunction(contig, position, "LEFT" if before else "RIGHT", virtual,
                mate_contig, mate_position, "RIGHT" if bracket == "[" else "LEFT", mate_virtual,
                ("t" + bracket + "p" + bracket) if before else (bracket + "p" + bracket + "t"))
    replacement = alt[1:] if alt.startswith(".") else alt[:-1]
    if virtual or _DNA.fullmatch(replacement) is None:
        raise InputValidationError("Unsupported single-breakend ALT or virtual local anchor.")
    return BreakendJunction(contig, position, "RIGHT" if alt.startswith(".") else "LEFT", False,
                           None, None, None, None, ".t" if alt.startswith(".") else "t.")


__all__ = ["BndPolicy", "BreakendJunction", "POLICY", "POLICY_SHA256", "parse_breakend"]
