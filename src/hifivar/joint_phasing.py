"""Explicit optional SV/TR inputs for HiPhase; no merged variant representation."""
from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass
import gzip
from itertools import zip_longest
from pathlib import Path
import re
from typing import TYPE_CHECKING

from hifivar.exceptions import ConfigurationError, InputValidationError, OutputValidationError

if TYPE_CHECKING:
    from hifivar.phasing import PhasingRequest


@dataclass(frozen=True, slots=True)
class PhasingTrackInput:
    """One explicitly declared, supported caller output, separate from small VCF."""

    category: str
    source_tool: str
    sample_id: str
    reference_build: str | None
    vcf_path: Path
    index_path: Path

    def __post_init__(self) -> None:
        supported = {"sv": {"pbsv", "sawfish"}, "tr": {"trgt"}}
        if (not isinstance(self.category, str) or not isinstance(self.source_tool, str)
                or self.category not in supported or self.source_tool not in supported[self.category]):
            raise InputValidationError("Joint HiPhase requires declared pbsv/Sawfish SV or TRGT TR input.")
        if not isinstance(self.sample_id, str) or not self.sample_id.strip():
            raise InputValidationError("Joint phasing sample_id must be non-empty.")
        vcf, index = Path(self.vcf_path).expanduser(), Path(self.index_path).expanduser()
        if not vcf.name.endswith(".vcf.gz") or index != Path(f"{vcf}.tbi"):
            raise InputValidationError("Joint phasing requires BGZF VCF and explicit matching .tbi.")
        object.__setattr__(self, "vcf_path", vcf)
        object.__setattr__(self, "index_path", index)

    def to_dict(self) -> dict[str, object]:
        return {"category": self.category, "source_tool": self.source_tool,
                "sample_id": self.sample_id, "reference_build": self.reference_build,
                "vcf_path": str(self.vcf_path), "index_path": str(self.index_path)}


@dataclass(frozen=True, slots=True)
class PhasedTrackArtifact:
    """A separate phased category output; not an independently called variant set."""

    category: str
    source_tool: str
    source_vcf: Path
    vcf_path: Path
    index_path: Path
    records_identity_and_genotype_checked: int

    def to_dict(self) -> dict[str, object]:
        return {"category": self.category, "source_tool": self.source_tool,
                "source_vcf": str(self.source_vcf), "vcf_path": str(self.vcf_path),
                "index_path": str(self.index_path),
                "records_identity_and_genotype_checked": self.records_identity_and_genotype_checked,
                "phase_accuracy": "NOT_ASSESSED"}


def phasing_options(section: Mapping[str, object]) -> tuple[str | None, bool, str]:
    """Validate options without silently enabling a producer or choosing a caller."""
    caller = section.get("sv_caller")
    include_tr = section.get("include_tr", False)
    global_mode = section.get("global_realignment", "auto")
    if caller is not None and caller not in ("pbsv", "sawfish"):
        raise ConfigurationError("phasing.sv_caller must be null, pbsv, or sawfish.")
    if not isinstance(include_tr, bool):
        raise ConfigurationError("phasing.include_tr must be boolean.")
    if global_mode not in ("auto", "enabled", "disabled"):
        raise ConfigurationError("phasing.global_realignment must be auto, enabled, or disabled.")
    return caller, include_tr, global_mode  # type: ignore[return-value]


def validate_phasing_producers(config: Mapping[str, object]) -> None:
    """Launcher/DAG dependency checks; independent API can consume ready artifacts."""
    section = config.get("phasing", {})
    if not isinstance(section, Mapping):
        raise ConfigurationError("phasing must be a mapping.")
    caller, include_tr, _ = phasing_options(section)
    if not section.get("enabled", False):
        return
    small, sv, tr = (config.get(key, {}) for key in ("small", "sv", "tr"))
    if not isinstance(small, Mapping) or not small.get("enabled", False):
        raise ConfigurationError("phasing requires small.enabled for this workflow.")
    if caller is not None and (
        not isinstance(sv, Mapping) or not sv.get("enabled", False)
        or not isinstance(sv.get(caller), Mapping) or not sv[caller].get("enabled", False)
    ):
        raise ConfigurationError(f"phasing.sv_caller={caller} requires that enabled SV producer.")
    if include_tr and (not isinstance(tr, Mapping) or not tr.get("enabled", False)):
        raise ConfigurationError("phasing.include_tr requires tr.enabled.")


def validate_joint_source_headers(request: PhasingRequest) -> None:
    """Check exact identity and category headers, BGZF/TBI; not sequence accuracy."""
    from hifivar.annotation import audit_annotation_vcf_identity
    from hifivar.small import _validate_bgzf, _validate_tabix_index

    for category, source, index, _, _, _ in request.vcf_pairs:
        _validate_bgzf(source, f"joint phasing {category} input")
        _validate_tabix_index(index)
        audit_annotation_vcf_identity(source, request.sample_id, request.alignment.reference)
        if category == "small":
            continue
        required = {"SVTYPE"} if category == "sv" else {"TRID", "MOTIFS", "STRUC"}
        observed: set[str] = set()
        try:
            with gzip.open(source, "rt", encoding="utf-8") as handle:
                for line in handle:
                    match = re.match(r"##INFO=<ID=([^,>]+)", line)
                    if match:
                        observed.add(match.group(1))
                    if line.startswith("#CHROM\t"):
                        break
        except (OSError, EOFError, UnicodeError) as error:
            raise InputValidationError(f"Unable to inspect joint phasing source: {source}") from error
        if not required <= observed:
            raise InputValidationError(f"Joint {category} input lacks category INFO headers: {required - observed}.")


def _records(path: Path) -> Iterator[list[str]]:
    try:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if not line.startswith("#"):
                    fields = line.rstrip("\r\n").split("\t")
                    if len(fields) != 10:
                        raise OutputValidationError(f"Joint phasing requires one-sample VCF records: {path}")
                    yield fields
    except (OSError, EOFError, UnicodeError) as error:
        raise OutputValidationError(f"Unable to stream joint phasing VCF: {path}") from error


def _alleles(fields: list[str]) -> tuple[str, ...]:
    keys, values = fields[8].split(":"), fields[9].split(":")
    if len(keys) != len(set(keys)) or "GT" not in keys:
        raise OutputValidationError("Joint phasing record lacks unique FORMAT/GT.")
    if keys.index("GT") >= len(values):
        raise OutputValidationError("Joint phasing record has missing GT value.")
    alleles = re.split(r"[/|]", values[keys.index("GT")])
    if not alleles or any(a != "." and not a.isdecimal() for a in alleles):
        raise OutputValidationError("Malformed joint phasing GT.")
    return tuple(sorted(alleles))


def validate_record_preservation(source: Path, output: Path) -> int:
    """Stream ordered CHROM/POS/ID/REF/ALT and unordered GT alleles; no truth claim.

    Phase orientation and PS/PF may change by design. This check does not claim
    byte-identical non-phase annotations or correctness of assigned phase.
    """
    count = 0
    for before, after in zip_longest(_records(source), _records(output)):
        if before is None or after is None:
            raise OutputValidationError("HiPhase changed the number of source records.")
        if before[:5] != after[:5] or _alleles(before) != _alleles(after):
            raise OutputValidationError(f"HiPhase changed variant identity or genotype alleles at record {count + 1}.")
        count += 1
    return count
