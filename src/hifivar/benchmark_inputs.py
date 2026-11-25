"""Read-only, fail-closed reference gates for structural-variant benchmarks.

Header compatibility is not reference-sequence identity. Callers must separately
check VCF REF alleles against the intended FASTA before scientific benchmarking.
No function in this module renames contigs, indexes or rewrites an input.
"""

from __future__ import annotations

import gzip
import hashlib
import re
from pathlib import Path
from typing import TextIO

from hifivar.exceptions import InputValidationError, OutputValidationError, ReferenceError
from hifivar.reference import ReferenceGenome
from hifivar.sv import validate_sv_vcf


_EXPLICIT_BUILD = re.compile(
    r"(?<![A-Za-z0-9])(?:GRCh(?P<grch>37|38)|hg(?P<hg>19|38))(?![0-9])",
    re.IGNORECASE,
)
_UNIDENTIFIED_SAMPLES = frozenset({"unknown", "unknownsamplename"})


def audit_sv_benchmark_input(
    vcf: Path,
    sample_id: str,
    reference: ReferenceGenome,
    *,
    require_contig_lengths: bool = True,
) -> dict[str, object]:
    """Audit one indexed single-sample VCF without reading variant records.

    The default benchmark gate requires positive lengths for every declared
    contig. An explicitly non-strict diagnostic call may return ``UNVERIFIED``
    for absent lengths, but never upgrades them to a passed length check.
    Exact-name contig subsets are allowed; a truth VCF need not list reference
    decoys. A missing or unrecognized ``##reference`` is recorded, not guessed
    from contig names. Explicit GRCh37/hg19 versus GRCh38/hg38 conflicts fail.

    This reuses the SV BGZF/TBI/header/sample checks. It does not establish TBI
    queryability, scan record coordinates or verify REF alleles. The returned
    ``sequence_check`` always requires a separate external sequence audit.
    """
    path = _input_path(vcf, "SV benchmark VCF")
    _check_reference(reference)
    if (
        not isinstance(sample_id, str)
        or not sample_id
        or any(character.isspace() for character in sample_id)
        or sample_id.casefold() in _UNIDENTIFIED_SAMPLES
    ):
        raise InputValidationError(
            "SAMPLE_IDENTITY_UNVERIFIED: benchmark requires an explicit, "
            "non-placeholder sample_id."
        )
    if not isinstance(require_contig_lengths, bool):
        raise InputValidationError("require_contig_lengths must be a boolean.")
    index = Path(f"{path}.tbi")
    try:
        validate_sv_vcf(path, index, sample_id=sample_id, reference=reference)
    except OutputValidationError as error:
        raise InputValidationError(f"SV_BENCHMARK_INPUT_INVALID: {error}") from error

    lengths: dict[str, int | None] = {}
    reference_headers: list[str] = []
    header_digest = hashlib.sha256()
    try:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            for raw_line in handle:
                header_digest.update(raw_line.encode("utf-8"))
                line = raw_line.rstrip("\r\n")
                if line.startswith("##contig="):
                    attributes = _contig_attributes(line, path)
                    name = attributes.get("ID", "")
                    if not name or any(character.isspace() for character in name):
                        raise InputValidationError(f"Invalid VCF contig ID in '{path}'.")
                    if name in lengths:
                        raise InputValidationError(
                            f"Duplicate VCF contig declaration '{name}' in '{path}'."
                        )
                    length_text = attributes.get("length")
                    length = None
                    if length_text is not None:
                        if not length_text.isascii() or not length_text.isdecimal():
                            raise InputValidationError(
                                f"Invalid VCF contig length for '{name}' in '{path}'."
                            )
                        length = int(length_text)
                        if length <= 0:
                            raise InputValidationError(
                                f"Invalid VCF contig length for '{name}' in '{path}'."
                            )
                    lengths[name] = length
                elif line.startswith("##reference="):
                    reference_headers.append(line.partition("=")[2])
                elif line.startswith("#CHROM\t"):
                    break
    except (OSError, EOFError, UnicodeError) as error:
        raise InputValidationError(f"Unable to read VCF header '{path}': {error}") from error

    if not lengths:
        raise InputValidationError(f"VCF header has no contig declarations: '{path}'.")
    explicit_builds = _explicit_builds(reference_headers)
    declared_builds = _explicit_builds([reference.build] if reference.build else [])
    if len(explicit_builds) > 1 or (
        explicit_builds and declared_builds and explicit_builds != declared_builds
    ):
        raise ReferenceError(
            "REFERENCE_BUILD_MISMATCH: VCF ##reference identifies "
            f"{sorted(explicit_builds)!r}, declared reference is {reference.build!r}."
        )
    reference_lengths = {contig.name: contig.length for contig in reference.contigs}
    for name, length in lengths.items():
        if name not in reference_lengths:
            raise ReferenceError(
                f"REFERENCE_CONTIG_MISMATCH: '{name}' is absent from the reference."
            )
        if length is not None and length != reference_lengths[name]:
            raise ReferenceError(
                f"CONTIG_LENGTH_MISMATCH: '{name}' VCF length {length} != "
                f"reference length {reference_lengths[name]}."
            )
    missing = [name for name, length in lengths.items() if length is None]
    if missing and require_contig_lengths:
        raise InputValidationError(
            "CONTIG_LENGTH_UNVERIFIED: benchmark requires VCF header lengths; "
            f"missing for {missing!r}."
        )
    build_check = "UNVERIFIED_MISSING_OR_UNRECOGNIZED_REFERENCE_METADATA"
    if explicit_builds and declared_builds:
        build_check = "PASS"
    elif explicit_builds:
        build_check = "UNVERIFIED_DECLARED_REFERENCE_BUILD"
    return {
        "audit_status": "UNVERIFIED" if missing else "PASS",
        "scope": "HEADER_AND_BGZF_TBI_SIGNATURE_ONLY",
        "vcf": str(path),
        "index": str(index),
        "sample_id": sample_id,
        "sample_check": "PASS",
        "contig_name_check": "PASS_SUBSET",
        "contig_length_check": "UNVERIFIED_MISSING_LENGTHS" if missing else "PASS",
        "declared_contig_count": len(lengths),
        "reference_contig_count": len(reference_lengths),
        "contig_lengths": lengths,
        "missing_contig_lengths": missing,
        "reference_headers": reference_headers,
        "explicit_vcf_builds": sorted(explicit_builds),
        "declared_reference_build": reference.build,
        "reference_build_check": build_check,
        "header_sha256": header_digest.hexdigest(),
        "sequence_check": "NOT_RUN_REQUIRES_BCFTOOLS",
        "record_validation": "NOT_RUN_HEADER_ONLY",
        "input_modified": False,
    }


def audit_benchmark_regions(bed: Path, reference: ReferenceGenome) -> dict[str, object]:
    """Stream a nonempty BED and count its zero-based half-open interval union.

    Contigs must be grouped (in any order), with nondecreasing starts within a
    contig. This explicit ordering contract permits an exact union with constant
    interval memory; unsorted inputs fail rather than being silently sorted.
    Overlapping and duplicate intervals are retained and reported, not rewritten.
    Blank, comment, track and browser lines are ignored. Extra BED columns are
    preserved on disk and have no effect on coordinate validation.
    """
    path = _input_path(bed, "Benchmark BED")
    _check_reference(reference)
    reference_lengths = {contig.name: contig.length for contig in reference.contigs}
    counts: dict[str, dict[str, int]] = {}
    current: str | None = None
    previous_start = -1
    union_end = -1
    interval_count = input_bases = union_bases = overlap_count = adjacent_count = 0
    try:
        with _open_bed(path) as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()
                if not line or line.startswith("#") or line.split(maxsplit=1)[0] in {"track", "browser"}:
                    continue
                fields = line.split()
                if len(fields) < 3:
                    raise InputValidationError(f"BED line {line_number} requires at least three columns.")
                contig, start_text, end_text = fields[:3]
                if contig not in reference_lengths:
                    raise ReferenceError(
                        f"REFERENCE_CONTIG_MISMATCH: BED line {line_number} contig '{contig}' "
                        "is absent from the reference."
                    )
                if any(not value.isascii() or not value.isdecimal() for value in (start_text, end_text)):
                    raise InputValidationError(f"Invalid BED integer coordinates at line {line_number}.")
                start, end = int(start_text), int(end_text)
                if start >= end or end > reference_lengths[contig]:
                    raise InputValidationError(
                        f"BED_COORDINATE_OUT_OF_BOUNDS: line {line_number} has "
                        f"{contig}:{start}-{end}; require 0 <= start < end <= "
                        f"{reference_lengths[contig]}."
                    )
                if contig != current:
                    if contig in counts:
                        raise InputValidationError(
                            f"BED_NOT_SORTED: contig '{contig}' reappears at line {line_number}."
                        )
                    current, previous_start, union_end = contig, -1, -1
                    counts[contig] = {"interval_count": 0, "input_bases": 0, "union_bases": 0}
                if start < previous_start:
                    raise InputValidationError(
                        f"BED_NOT_SORTED: decreasing start at line {line_number}."
                    )
                span = end - start
                added = max(0, end - max(start, union_end))
                overlap_count += int(start < union_end)
                adjacent_count += int(start == union_end)
                interval_count += 1
                input_bases += span
                union_bases += added
                counts[contig]["interval_count"] += 1
                counts[contig]["input_bases"] += span
                counts[contig]["union_bases"] += added
                previous_start, union_end = start, max(end, union_end)
    except (OSError, EOFError, UnicodeError) as error:
        raise InputValidationError(f"Unable to stream benchmark BED '{path}': {error}") from error
    if interval_count == 0:
        raise InputValidationError(f"Benchmark BED contains no intervals: '{path}'.")
    return {
        "audit_status": "PASS",
        "bed": str(path),
        "coordinate_system": "ZERO_BASED_HALF_OPEN",
        "ordering": "CONTIG_GROUPED_NONDECREASING_START",
        "contig_name_check": "PASS_SUBSET",
        "coordinate_bounds_check": "PASS",
        "interval_count": interval_count,
        "contig_count": len(counts),
        "input_bases": input_bases,
        "union_bases": union_bases,
        "overlap_bases_counted_once": input_bases - union_bases,
        "overlapping_interval_count": overlap_count,
        "adjacent_interval_count": adjacent_count,
        "per_contig": counts,
        "intervals_rewritten": False,
        "input_modified": False,
    }


def _input_path(value: Path, label: str) -> Path:
    if not isinstance(value, (str, Path)) or (isinstance(value, str) and not value.strip()):
        raise InputValidationError(f"{label} must be a nonempty path.")
    path = Path(value)
    if not path.is_file() or path.stat().st_size == 0:
        raise InputValidationError(f"{label} is missing, empty or not a file: '{path}'.")
    return path


def _check_reference(reference: ReferenceGenome) -> None:
    if not isinstance(reference, ReferenceGenome):
        raise InputValidationError("Benchmark audit requires ReferenceGenome.")


def _explicit_builds(values: list[str]) -> set[str]:
    builds: set[str] = set()
    for value in values:
        for match in _EXPLICIT_BUILD.finditer(value):
            builds.add("GRCh37" if (match.group("grch") or match.group("hg")) in {"37", "19"} else "GRCh38")
    return builds


def _contig_attributes(line: str, path: Path) -> dict[str, str]:
    if not line.startswith("##contig=<") or not line.endswith(">"):
        raise InputValidationError(f"Malformed VCF contig declaration in '{path}'.")
    body = line[len("##contig=<") : -1]
    fields: list[str] = []
    quoted = escaped = False
    start = 0
    for index, character in enumerate(body):
        if escaped:
            escaped = False
        elif character == "\\" and quoted:
            escaped = True
        elif character == '"':
            quoted = not quoted
        elif character == "," and not quoted:
            fields.append(body[start:index])
            start = index + 1
    if quoted or escaped:
        raise InputValidationError(f"Malformed quoted VCF contig declaration in '{path}'.")
    fields.append(body[start:])
    attributes: dict[str, str] = {}
    for field in fields:
        key, separator, value = field.partition("=")
        if not separator or not key or not value or key in attributes:
            raise InputValidationError(f"Malformed VCF contig attributes in '{path}'.")
        attributes[key] = value.strip('"')
    return attributes


def _open_bed(path: Path) -> TextIO:
    if path.suffix.lower() in {".gz", ".bgz"}:
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


__all__ = ["audit_benchmark_regions", "audit_sv_benchmark_input"]
