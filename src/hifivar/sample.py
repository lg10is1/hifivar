"""Sample identity and primary HiFi input metadata."""

from __future__ import annotations

import os
import re
from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hifivar import validation
from hifivar.exceptions import InputValidationError
from hifivar.logging_utils import get_logger


PathInput = str | Path

_FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
_SAMPLE_ID_PATTERN = re.compile(r"(?!\.)[A-Za-z0-9_.-]+\Z")
_LOGGER = get_logger(__name__)


class InputType(str, Enum):
    """Supported primary sequencing-input types."""

    FASTQ = "fastq"
    BAM = "bam"
    CRAM = "cram"


@dataclass(frozen=True, slots=True)
class InputDataset:
    """Immutable metadata for one sample's primary sequencing input.

    Use :meth:`from_files` to validate files on disk. Direct construction is
    intended for already validated metadata and still enforces type, suffix,
    uniqueness, and file-count invariants without accessing file contents.
    """

    input_type: InputType
    files: tuple[Path, ...]

    def __post_init__(self) -> None:
        """Normalize immutable paths and enforce structural invariants."""
        if not isinstance(self.input_type, InputType):
            raise InputValidationError("input_type must be an InputType value.")

        normalized_files = _coerce_files(self.files)
        _validate_unique_files(normalized_files)
        inferred_type = _infer_normalized_input_type(normalized_files)
        if inferred_type is not self.input_type:
            raise InputValidationError(
                f"Explicit input type '{self.input_type.value}' does not match "
                f"file suffix type '{inferred_type.value}'."
            )
        _validate_file_count(self.input_type, normalized_files)
        object.__setattr__(self, "files", normalized_files)

    @classmethod
    def from_files(
        cls,
        files: Iterable[PathInput],
        *,
        input_type: InputType | None = None,
    ) -> InputDataset:
        """Create and validate one primary input dataset.

        FASTQ datasets contain one or more files. BAM and CRAM datasets contain
        exactly one file. Validation is lightweight: FASTQ checks only its first
        record, while BAM/CRAM checks path, suffix, and non-empty status.
        """
        normalized_files = _coerce_files(files)
        _validate_unique_files(normalized_files)
        inferred_type = _infer_normalized_input_type(normalized_files)

        if input_type is not None and not isinstance(input_type, InputType):
            raise InputValidationError(
                "Explicit input_type must be an InputType value."
            )
        selected_type = inferred_type if input_type is None else input_type
        if selected_type is not inferred_type:
            raise InputValidationError(
                f"Explicit input type '{selected_type.value}' does not match "
                f"file suffix type '{inferred_type.value}'."
            )
        _validate_file_count(selected_type, normalized_files)

        _LOGGER.debug(
            "Creating input dataset: input_type=%s, file_count=%d",
            selected_type.value,
            len(normalized_files),
        )
        if selected_type is InputType.FASTQ:
            for path in normalized_files:
                validation.validate_fastq(path)
        else:
            validation.validate_alignment_file(
                normalized_files[0],
                require_index=False,
            )

        return cls(input_type=selected_type, files=normalized_files)

    def validate_index(self) -> None:
        """Require an existing BAM/CRAM index without creating one."""
        if self.input_type is InputType.FASTQ:
            raise InputValidationError(
                "Index validation is not applicable to FASTQ input."
            )
        validation.validate_alignment_file(self.files[0], require_index=True)

    def to_dict(self) -> dict[str, object]:
        """Return JSON/YAML-friendly primary-input metadata."""
        return {
            "type": self.input_type.value,
            "files": [str(path) for path in self.files],
        }


@dataclass(frozen=True, slots=True)
class Sample:
    """Immutable machine sample identity and its one primary input dataset."""

    sample_id: str
    input: InputDataset

    def __post_init__(self) -> None:
        """Validate the machine-safe sample identifier and input model."""
        validate_sample_id(self.sample_id)
        if not isinstance(self.input, InputDataset):
            raise InputValidationError(
                "Sample input must be an InputDataset instance."
            )
        _LOGGER.debug(
            "Creating sample: sample_id=%s, input_type=%s, file_count=%d",
            self.sample_id,
            self.input.input_type.value,
            len(self.input.files),
        )

    def to_dict(self) -> dict[str, object]:
        """Return JSON/YAML-friendly sample and primary-input metadata."""
        return {
            "sample_id": self.sample_id,
            "input": self.input.to_dict(),
        }


def infer_input_type(files: Iterable[PathInput]) -> InputType:
    """Infer one unambiguous input type from supported filename suffixes."""
    return _infer_normalized_input_type(_coerce_files(files))


def validate_sample_id(sample_id: str) -> None:
    """Require an unchanged ASCII-safe machine identifier.

    This function validates but never strips, rewrites, or otherwise normalizes
    the supplied identifier.
    """
    if not isinstance(sample_id, str) or not sample_id:
        raise InputValidationError(
            "sample_id must be a non-empty ASCII-safe string."
        )
    if sample_id in {".", ".."} or _SAMPLE_ID_PATTERN.fullmatch(sample_id) is None:
        raise InputValidationError(
            f"Invalid sample_id {sample_id!r}. Use only ASCII letters, digits, "
            "underscore, hyphen, and dot; a leading dot is not allowed."
        )


def _coerce_files(files: Iterable[PathInput]) -> tuple[Path, ...]:
    """Copy a file iterable into ordered, expanded, unresolved Paths."""
    if isinstance(files, (str, Path)):
        raise InputValidationError(
            "Input files must be an iterable of string or Path values, not one "
            "bare path."
        )

    try:
        supplied_files = tuple(files)
    except TypeError as error:
        raise InputValidationError(
            "Input files must be an iterable of string or Path values."
        ) from error

    if not supplied_files:
        raise InputValidationError("Input dataset must contain at least one file.")

    normalized: list[Path] = []
    for index, path in enumerate(supplied_files):
        if not isinstance(path, (str, Path)):
            raise InputValidationError(
                f"Input file at index {index} must be a string or Path."
            )
        if isinstance(path, str) and not path.strip():
            raise InputValidationError(
                f"Input file at index {index} must not be an empty path."
            )
        normalized.append(Path(path).expanduser())
    return tuple(normalized)


def _validate_unique_files(files: tuple[Path, ...]) -> None:
    """Reject duplicate normalized path spellings without resolving links."""
    identities: set[str] = set()
    for path in files:
        identity = os.path.normcase(os.path.normpath(str(path)))
        if identity in identities:
            raise InputValidationError(
                f"Input dataset contains duplicate file path: '{path}'."
            )
        identities.add(identity)


def _infer_normalized_input_type(files: tuple[Path, ...]) -> InputType:
    """Infer exactly one supported suffix type from normalized paths."""
    inferred_types = {_input_type_for_path(path) for path in files}
    if len(inferred_types) != 1:
        types = ", ".join(sorted(input_type.value for input_type in inferred_types))
        raise InputValidationError(
            f"Mixed primary input types are not allowed: {types}."
        )
    return next(iter(inferred_types))


def _input_type_for_path(path: Path) -> InputType:
    """Return the supported type for one suffix without guessing."""
    lowered_name = path.name.lower()
    if any(lowered_name.endswith(suffix) for suffix in _FASTQ_SUFFIXES):
        return InputType.FASTQ
    if lowered_name.endswith(".bam"):
        return InputType.BAM
    if lowered_name.endswith(".cram"):
        return InputType.CRAM
    raise InputValidationError(
        f"Unable to infer input type from unsupported file suffix: '{path}'."
    )


def _validate_file_count(
    input_type: InputType,
    files: tuple[Path, ...],
) -> None:
    """Enforce HiFi FASTQ and single-alignment cardinality rules."""
    if input_type in {InputType.BAM, InputType.CRAM} and len(files) != 1:
        raise InputValidationError(
            f"{input_type.value.upper()} input requires exactly one file; "
            f"received {len(files)}."
        )


__all__ = [
    "InputDataset",
    "InputType",
    "Sample",
    "infer_input_type",
    "validate_sample_id",
]
