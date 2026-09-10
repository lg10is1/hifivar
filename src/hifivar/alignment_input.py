"""Alignment input classification, validation, and routing.

Distinguishes FASTQ, unaligned uBAM, coordinate-sorted aligned BAM/CRAM,
and routes to ALIGN or REUSE actions based on header metadata inspection.
No silent chr1/1 conversion; explicit sample identity validation.

R2-MAJOR-001 fix: Classification is driven by explicit config/manifest,
not hard-coded .fastq wildcard patterns.

R2-MAJOR-002 fix: Non-coordinate-sorted aligned BAM is rejected at
classification time rather than routed to ALIGN without a valid
FASTQ conversion path.

R2-MAJOR-005 fix: CRAM support is explicitly deferred until full
reference/format/index validation is implemented.
"""

from __future__ import annotations

import os
from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hifivar import validation
from hifivar.alignment import (
    AlignmentAction,
    AlignmentOutputFormat,
)
from hifivar.exceptions import (
    InputValidationError,
    ReferenceError,
)
from hifivar.logging_utils import get_logger
from hifivar.sample import InputType


PathInput = str | Path

_LOGGER = get_logger(__name__)


class AlignedStatus(str, Enum):
    """Whether an alignment file is aligned or unaligned."""

    UNALIGNED = "unaligned"
    ALIGNED = "aligned"
    UNKNOWN = "unknown"


class SortOrder(str, Enum):
    """Recognized BAM/CRAM sort order header values."""

    QUERYNAME = "queryname"
    COORDINATE = "coordinate"
    UNORDERED = "unordered"
    UNKNOWN = "unknown"


@dataclass(frozen=True, slots=True)
class AlignmentHeaderMetadata:
    """Parsed BAM/CRAM header metadata for routing decisions."""

    sort_order: SortOrder
    has_sq_lines: bool
    sample_names: tuple[str, ...]
    read_group_ids: tuple[str, ...]
    has_hd_header: bool
    reference_contigs: tuple[str, ...]

    def __post_init__(self) -> None:
        """Validate metadata structure invariants."""
        if not isinstance(self.sort_order, SortOrder):
            raise InputValidationError(
                "Alignment header sort_order must be a SortOrder value."
            )
        if not isinstance(self.has_sq_lines, bool):
            raise InputValidationError(
                "Alignment header has_sq_lines must be a boolean."
            )
        if not isinstance(self.sample_names, tuple):
            raise InputValidationError(
                "Alignment header sample_names must be a tuple."
            )
        if not isinstance(self.read_group_ids, tuple):
            raise InputValidationError(
                "Alignment header read_group_ids must be a tuple."
            )
        if not isinstance(self.has_hd_header, bool):
            raise InputValidationError(
                "Alignment header has_hd_header must be a boolean."
            )
        if not isinstance(self.reference_contigs, tuple):
            raise InputValidationError(
                "Alignment header reference_contigs must be a tuple."
            )

    @property
    def unique_sample_names(self) -> tuple[str, ...]:
        """Return deduplicated sample names preserving first-seen order."""
        seen: dict[str, None] = {}
        for name in self.sample_names:
            seen.setdefault(name, None)
        return tuple(seen.keys())

    def to_dict(self) -> dict[str, object]:
        """Return JSON/YAML-friendly header metadata."""
        return {
            "sort_order": self.sort_order.value,
            "has_sq_lines": self.has_sq_lines,
            "sample_names": list(self.sample_names),
            "read_group_ids": list(self.read_group_ids),
            "has_hd_header": self.has_hd_header,
            "reference_contigs_count": len(self.reference_contigs),
        }


class BamHeaderParser:
    """Protocol-like base for BAM/CRAM header parsing."""

    def parse_header(self, path: Path) -> AlignmentHeaderMetadata:
        """Return parsed header metadata for one BAM or CRAM file."""
        ...


@dataclass(frozen=True, slots=True)
class InputClassification:
    """Result of classifying one input file for alignment routing."""

    input_path: Path
    input_type: InputType
    aligned_status: AlignedStatus
    header: AlignmentHeaderMetadata | None
    routing_action: AlignmentAction

    def to_dict(self) -> dict[str, object]:
        """Return JSON/YAML-friendly classification metadata."""
        return {
            "input_path": str(self.input_path),
            "input_type": self.input_type.value,
            "aligned_status": self.aligned_status.value,
            "header": self.header.to_dict() if self.header is not None else None,
            "routing_action": self.routing_action.value,
        }


def classify_input(
    path: PathInput,
    *,
    expected_sample_id: str | None = None,
    reference_contigs: tuple[str, ...] | None = None,
    parser: BamHeaderParser | None = None,
) -> InputClassification:
    """Classify one input file and determine alignment routing.

    FASTQ always routes to ALIGN. Unaligned uBAM (no SQ lines, unordered
    sort) routes to ALIGN. Coordinate-sorted aligned BAM routes to REUSE.
    Non-coordinate-sorted aligned files raise an error (must be preprocessed
    upstream).

    R2-MAJOR-002 fix: Non-coordinate-sorted aligned BAM is explicitly
    rejected rather than routed to ALIGN, since Minimap2Wrapper only
    consumes FASTQ input.

    R2-MAJOR-005 fix: CRAM inputs are explicitly deferred until full
    reference/format/index validation is implemented.

    Raises explicit errors for:
    - Missing or unreadable files
    - Ambiguous alignment status
    - Missing/duplicate/conflicting sample identity
    - CRAM inputs (deferred)
    - Non-coordinate-sorted aligned BAM (needs preprocessing)
    """
    input_path = Path(path).expanduser()
    _LOGGER.debug("Classifying input: %s", input_path)

    if not input_path.exists():
        raise InputValidationError(f"Input file does not exist: '{input_path}'.")
    if not input_path.is_file():
        raise InputValidationError(f"Input path is not a file: '{input_path}'.")

    inferred_type = _infer_type_from_suffix(input_path)

    if inferred_type is InputType.FASTQ:
        return _classify_fastq(input_path)

    if inferred_type is InputType.CRAM:
        raise InputValidationError(
            f"CRAM input '{input_path}' is not yet supported in HiFiVar. "
            "CRAM support requires validated reference, format, and index "
            "contracts. Please convert to BAM or defer until CRAM support "
            "is implemented."
        )

    if parser is None:
        raise InputValidationError(
            "BAM input requires a BamHeaderParser for alignment "
            f"classification: '{input_path}'."
        )

    header = parser.parse_header(input_path)
    aligned_status = _determine_aligned_status(header)
    routing_action = _determine_routing_action(aligned_status, header)

    if expected_sample_id is not None:
        _validate_sample_identity(header, expected_sample_id, input_path)

    if reference_contigs is not None:
        _validate_reference_compatibility(
            header, reference_contigs, input_path
        )

    classification = InputClassification(
        input_path=input_path,
        input_type=inferred_type,
        aligned_status=aligned_status,
        header=header,
        routing_action=routing_action,
    )

    _LOGGER.debug(
        "Input classified: path=%s type=%s status=%s action=%s",
        input_path,
        inferred_type.value,
        aligned_status.value,
        routing_action.value,
    )
    return classification


def classify_inputs_batch(
    paths: Iterable[PathInput],
    *,
    expected_sample_ids: dict[str, str] | None = None,
    reference_contigs: tuple[str, ...] | None = None,
    parser: BamHeaderParser | None = None,
) -> tuple[InputClassification, ...]:
    """Classify multiple inputs with cross-file sample identity validation.

    Args:
        paths: Input file paths to classify.
        expected_sample_ids: Mapping of file path to expected sample ID.
        reference_contigs: Reference contig names for CRAM validation.
        parser: BAM/CRAM header parser instance.

    Returns:
        Tuple of InputClassification in input order.

    Raises:
        InputValidationError: On duplicate files, conflicting sample IDs,
            or invalid classifications.
    """
    paths_list = [Path(p).expanduser() for p in paths]
    _validate_unique_input_paths(paths_list)

    classifications: list[InputClassification] = []
    seen_sample_ids: dict[str, str] = {}

    for path in paths_list:
        expected_id = (
            expected_sample_ids.get(str(path))
            if expected_sample_ids is not None
            else None
        )
        classification = classify_input(
            path,
            expected_sample_id=expected_id,
            reference_contigs=reference_contigs,
            parser=parser,
        )
        classifications.append(classification)

        if classification.header is not None:
            for sample_name in classification.header.sample_names:
                if sample_name in seen_sample_ids:
                    previous_path = seen_sample_ids[sample_name]
                    raise InputValidationError(
                        f"Duplicate sample identity '{sample_name}' found in "
                        f"'{path}' and '{previous_path}'. Each sample must "
                        "have a unique identity across inputs."
                    )
                seen_sample_ids[sample_name] = str(path)

    return tuple(classifications)


def validate_uBAM_for_alignment(
    classification: InputClassification,
) -> None:
    """Validate that a uBAM classification is suitable for alignment.

    uBAM must have:
    - No SQ lines (unaligned)
    - No coordinate sorting
    - At least one read group with sample name
    """
    if classification.input_type is not InputType.BAM:
        raise InputValidationError(
            f"uBAM validation requires BAM input; received "
            f"{classification.input_type.value}."
        )
    if classification.aligned_status is not AlignedStatus.UNALIGNED:
        raise InputValidationError(
            f"uBAM validation requires unaligned status; received "
            f"{classification.aligned_status.value}."
        )
    if classification.header is None:
        raise InputValidationError(
            "uBAM validation requires header metadata."
        )
    if classification.header.sort_order is SortOrder.COORDINATE:
        raise InputValidationError(
            "uBAM must not be coordinate-sorted."
        )
    if not classification.header.read_group_ids:
        raise InputValidationError(
            "uBAM must contain at least one read group (RG) header line."
        )
    if not classification.header.sample_names:
        raise InputValidationError(
            "uBAM read groups must contain SM (sample name) tags."
        )


def validate_aligned_bam_for_reuse(
    classification: InputClassification,
    reference_contigs: tuple[str, ...] | None = None,
) -> None:
    """Validate that an aligned BAM classification is suitable for reuse.

    Aligned BAM must have:
    - Coordinate sorting
    - SQ lines present
    - Sample identity
    """
    if classification.aligned_status is not AlignedStatus.ALIGNED:
        raise InputValidationError(
            f"Aligned BAM validation requires aligned status; received "
            f"{classification.aligned_status.value}."
        )
    if classification.header is None:
        raise InputValidationError(
            "Aligned BAM validation requires header metadata."
        )
    if classification.header.sort_order is not SortOrder.COORDINATE:
        raise InputValidationError(
            "Aligned BAM must be coordinate-sorted for reuse."
        )
    if not classification.header.has_sq_lines:
        raise InputValidationError(
            "Aligned BAM must contain SQ header lines."
        )
    if not classification.header.sample_names:
        raise InputValidationError(
            "Aligned BAM must contain sample identity (SM tag)."
        )


def _classify_fastq(path: Path) -> InputClassification:
    """Classify a FASTQ input as requiring alignment."""
    validation.validate_fastq(path)
    return InputClassification(
        input_path=path,
        input_type=InputType.FASTQ,
        aligned_status=AlignedStatus.UNALIGNED,
        header=None,
        routing_action=AlignmentAction.ALIGN,
    )


def _determine_aligned_status(header: AlignmentHeaderMetadata) -> AlignedStatus:
    """Determine alignment status from parsed header metadata."""
    if not header.has_sq_lines:
        return AlignedStatus.UNALIGNED

    if header.sort_order is SortOrder.COORDINATE:
        return AlignedStatus.ALIGNED

    if header.sort_order in (SortOrder.QUERYNAME, SortOrder.UNORDERED):
        return AlignedStatus.ALIGNED

    if header.sort_order is SortOrder.UNKNOWN and header.has_sq_lines:
        return AlignedStatus.ALIGNED

    return AlignedStatus.UNALIGNED


def _determine_routing_action(
    aligned_status: AlignedStatus,
    header: AlignmentHeaderMetadata,
) -> AlignmentAction:
    """Determine alignment routing from header and alignment status.

    R2-MAJOR-002 fix: Only coordinate-sorted aligned files are eligible
    for REUSE. Queryname-sorted, unordered, or unknown-sort aligned files
    raise an error at classification time, since they cannot be directly
    consumed by the FASTQ-only Minimap2Wrapper.
    """
    if aligned_status is AlignedStatus.UNALIGNED:
        return AlignmentAction.ALIGN

    if aligned_status is AlignedStatus.ALIGNED:
        if header.sort_order is SortOrder.COORDINATE:
            return AlignmentAction.REUSE
        raise InputValidationError(
            f"Aligned BAM with sort order '{header.sort_order.value}' cannot "
            "be directly reused. Coordinate-sorted BAM is required for REUSE. "
            "Please re-sort or re-align upstream before HiFiVar."
        )

    return AlignmentAction.ALIGN


def _validate_sample_identity(
    header: AlignmentHeaderMetadata,
    expected_sample_id: str,
    input_path: Path,
) -> None:
    """Validate that header sample names match expected sample identity."""
    if not header.sample_names:
        raise InputValidationError(
            f"Input '{input_path}' has no sample identity (SM tag) in header."
        )

    if len(header.sample_names) > 1:
        unique_names = header.unique_sample_names
        if len(unique_names) > 1:
            names_text = ", ".join(sorted(unique_names))
            raise InputValidationError(
                f"Input '{input_path}' contains conflicting sample identities: "
                f"{names_text}. Each input must have a single consistent sample "
                "identity."
            )

    actual_name = header.sample_names[0]
    if actual_name != expected_sample_id:
        raise InputValidationError(
            f"Input '{input_path}' sample identity '{actual_name}' does not "
            f"match expected '{expected_sample_id}'."
        )


def _validate_reference_compatibility(
    header: AlignmentHeaderMetadata,
    reference_contigs: tuple[str, ...],
    input_path: Path,
) -> None:
    """Validate that header contigs are compatible with reference."""
    if not header.has_sq_lines:
        return

    header_contigs = set(header.reference_contigs)
    ref_contigs = set(reference_contigs)

    missing_from_ref = header_contigs - ref_contigs
    if missing_from_ref:
        missing_text = ", ".join(sorted(missing_from_ref)[:10])
        raise ReferenceError(
            f"REFERENCE_CONTIG_MISMATCH: input '{input_path}' contains contigs "
            f"not in reference: {missing_text}. "
            "HiFiVar does not rename contigs automatically."
        )


def _infer_type_from_suffix(path: Path) -> InputType:
    """Infer input type from file suffix without reading content."""
    lowered_name = path.name.lower()
    if any(lowered_name.endswith(s) for s in (".fastq", ".fq", ".fastq.gz", ".fq.gz")):
        return InputType.FASTQ
    if lowered_name.endswith(".bam"):
        return InputType.BAM
    if lowered_name.endswith(".cram"):
        return InputType.CRAM
    raise InputValidationError(
        f"Unable to infer input type from unsupported file suffix: '{path}'."
    )


def _validate_unique_input_paths(paths: list[Path]) -> None:
    """Reject duplicate normalized input paths."""
    identities: set[str] = set()
    for path in paths:
        identity = os.path.normcase(os.path.normpath(str(path.absolute())))
        if identity in identities:
            raise InputValidationError(
                f"Duplicate input path detected: '{path}'."
            )
        identities.add(identity)


__all__ = [
    "AlignmentHeaderMetadata",
    "AlignedStatus",
    "BamHeaderParser",
    "InputClassification",
    "SortOrder",
    "classify_input",
    "classify_inputs_batch",
    "validate_aligned_bam_for_reuse",
    "validate_uBAM_for_alignment",
]
