"""TSV sample-sheet, pedigree, and cohort metadata models."""

from __future__ import annotations

import csv
import os
from collections.abc import Iterator
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import TextIO

from hifivar import validation
from hifivar.exceptions import InputValidationError
from hifivar.logging_utils import get_logger
from hifivar.sample import InputDataset, InputType, Sample, validate_sample_id


PathInput = str | Path
Trio = tuple["SampleRecord", "SampleRecord", "SampleRecord"]

REQUIRED_COLUMNS = frozenset({"sample_id", "input"})
OPTIONAL_COLUMNS = frozenset(
    {"input_type", "sex", "father", "mother", "phenotype", "group"}
)
SUPPORTED_COLUMNS = REQUIRED_COLUMNS | OPTIONAL_COLUMNS

_SEX_ALIASES = {
    "m": "male",
    "male": "male",
    "f": "female",
    "female": "female",
    "unknown": "unknown",
    ".": "unknown",
}
_LOGGER = get_logger(__name__)


class Sex(str, Enum):
    """User-declared sample sex metadata; never inferred by HiFiVar."""

    MALE = "male"
    FEMALE = "female"
    UNKNOWN = "unknown"


@dataclass(frozen=True, slots=True)
class SampleRecord:
    """One validated sample plus optional pedigree/cohort metadata."""

    sample: Sample
    sex: Sex | None = None
    father: str | None = None
    mother: str | None = None
    phenotype: str | None = None
    group: str | None = None

    def __post_init__(self) -> None:
        """Normalize missing metadata and validate local relationships."""
        if not isinstance(self.sample, Sample):
            raise InputValidationError(
                "SampleRecord sample must be a Sample instance."
            )
        if self.sex is not None and not isinstance(self.sex, Sex):
            raise InputValidationError("SampleRecord sex must be a Sex value or None.")

        father = _normalize_optional_text(self.father)
        mother = _normalize_optional_text(self.mother)
        phenotype = _normalize_optional_text(self.phenotype)
        group = _normalize_optional_text(self.group)
        object.__setattr__(self, "father", father)
        object.__setattr__(self, "mother", mother)
        object.__setattr__(self, "phenotype", phenotype)
        object.__setattr__(self, "group", group)

        sample_id = self.sample.sample_id
        for role, parent_id in (("father", father), ("mother", mother)):
            if parent_id is None:
                continue
            try:
                validate_sample_id(parent_id)
            except InputValidationError as error:
                raise InputValidationError(
                    f"Invalid {role} sample ID {parent_id!r} for sample "
                    f"'{sample_id}': {error}"
                ) from error
            if parent_id == sample_id:
                raise InputValidationError(
                    f"Sample '{sample_id}' cannot be its own {role}."
                )
        if father is not None and father == mother:
            raise InputValidationError(
                f"Sample '{sample_id}' cannot have the same father and mother "
                f"'{father}'."
            )

    def to_dict(self) -> dict[str, object]:
        """Return flattened JSON/YAML-friendly sample metadata."""
        metadata = self.sample.to_dict()
        metadata.update(
            {
                "sex": self.sex.value if self.sex is not None else None,
                "father": self.father,
                "mother": self.mother,
                "phenotype": self.phenotype,
                "group": self.group,
            }
        )
        return metadata


@dataclass(frozen=True, slots=True)
class SampleSheet:
    """An ordered, validated collection of sample metadata records."""

    source_path: Path
    records: tuple[SampleRecord, ...]

    def __post_init__(self) -> None:
        """Normalize immutable storage and validate sheet-wide invariants."""
        if not isinstance(self.source_path, (str, Path)):
            raise InputValidationError(
                "SampleSheet source_path must be a string or Path."
            )
        source_path = Path(self.source_path).expanduser()
        records = tuple(self.records)
        if not records:
            raise InputValidationError("SampleSheet must contain at least one sample.")
        if any(not isinstance(record, SampleRecord) for record in records):
            raise InputValidationError(
                "SampleSheet records must contain only SampleRecord objects."
            )
        object.__setattr__(self, "source_path", source_path)
        object.__setattr__(self, "records", records)

        _validate_unique_samples(records)
        _validate_unique_inputs(records)
        records_by_id = {record.sample.sample_id: record for record in records}
        _validate_parent_relationships(records, records_by_id)
        _validate_pedigree_acyclic(records, records_by_id)

    @classmethod
    def from_tsv(cls, path: PathInput) -> SampleSheet:
        """Read and validate a UTF-8 or UTF-8-BOM TSV sample sheet."""
        sheet_path = validation.validate_file(path)
        if sheet_path.suffix.lower() != ".tsv":
            raise InputValidationError(
                f"Sample sheet must use the .tsv suffix: '{sheet_path}'."
            )

        try:
            with sheet_path.open(
                "rt",
                encoding="utf-8-sig",
                newline="",
            ) as handle:
                rows = list(
                    csv.reader(
                        _iter_content_lines(handle),
                        delimiter="\t",
                    )
                )
        except (OSError, UnicodeError, csv.Error) as error:
            raise InputValidationError(
                f"Unable to read sample sheet '{sheet_path}': {error}"
            ) from error

        if not rows:
            raise InputValidationError(
                f"Sample sheet '{sheet_path}' has no header or sample records."
            )
        header = tuple(column.strip() for column in rows[0])
        _validate_header(header, sheet_path)

        records: list[SampleRecord] = []
        for record_number, row in enumerate(rows[1:], start=1):
            if len(row) != len(header):
                raise InputValidationError(
                    f"Sample sheet '{sheet_path}' record {record_number} has "
                    f"{len(row)} columns; expected {len(header)}."
                )
            values = dict(zip(header, row, strict=True))
            try:
                records.append(_parse_record(values, sheet_path.parent))
            except InputValidationError as error:
                raise InputValidationError(
                    f"Invalid sample sheet '{sheet_path}' record "
                    f"{record_number}: {error}"
                ) from error

        sheet = cls(source_path=sheet_path, records=tuple(records))
        pedigree_count = sum(
            record.father is not None or record.mother is not None
            for record in sheet.records
        )
        _LOGGER.info(
            "Loaded %d samples from %s",
            len(sheet),
            sheet.source_path,
        )
        _LOGGER.debug(
            "Sample sheet metadata: sample_count=%d, pedigree_record_count=%d",
            len(sheet),
            pedigree_count,
        )
        return sheet

    def __len__(self) -> int:
        """Return the number of sample records."""
        return len(self.records)

    @property
    def sample_ids(self) -> tuple[str, ...]:
        """Return sample IDs in deterministic TSV order."""
        return tuple(record.sample.sample_id for record in self.records)

    def get_record(self, sample_id: str) -> SampleRecord:
        """Return one record by exact sample ID or raise ``KeyError``."""
        for record in self.records:
            if record.sample.sample_id == sample_id:
                return record
        raise KeyError(sample_id)

    def get_trios(self) -> tuple[Trio, ...]:
        """Return complete child/father/mother trios in child sheet order."""
        records_by_id = {
            record.sample.sample_id: record for record in self.records
        }
        trios: list[Trio] = []
        for child in self.records:
            if child.father is not None and child.mother is not None:
                trios.append(
                    (
                        child,
                        records_by_id[child.father],
                        records_by_id[child.mother],
                    )
                )
        return tuple(trios)

    def to_dict(self) -> dict[str, object]:
        """Return JSON/YAML-friendly source and ordered record metadata."""
        return {
            "source_path": str(self.source_path),
            "records": [record.to_dict() for record in self.records],
        }


def parse_sex(value: str | None) -> Sex | None:
    """Normalize declared sex aliases without inferring missing metadata."""
    if value is None or not value.strip():
        return None
    normalized = value.strip().casefold()
    try:
        return Sex(_SEX_ALIASES[normalized])
    except KeyError as error:
        allowed = "M, F, male, female, unknown, or ."
        raise InputValidationError(
            f"Invalid sex value {value!r}; expected {allowed}."
        ) from error


def _iter_content_lines(handle: TextIO) -> Iterator[str]:
    """Yield non-empty, non-comment lines to the TSV reader."""
    for line in handle:
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        yield line


def _validate_header(header: tuple[str, ...], sheet_path: Path) -> None:
    """Require unique canonical columns and reject unknown metadata."""
    if not header or any(not column for column in header):
        raise InputValidationError(
            f"Sample sheet '{sheet_path}' contains an empty column name."
        )
    if len(set(header)) != len(header):
        raise InputValidationError(
            f"Sample sheet '{sheet_path}' contains duplicate column names."
        )
    missing = sorted(REQUIRED_COLUMNS.difference(header))
    if missing:
        raise InputValidationError(
            f"Sample sheet '{sheet_path}' is missing required column(s): "
            f"{', '.join(missing)}."
        )
    unknown = sorted(set(header).difference(SUPPORTED_COLUMNS))
    if unknown:
        raise InputValidationError(
            f"Sample sheet '{sheet_path}' contains unknown column(s): "
            f"{', '.join(unknown)}."
        )


def _parse_record(values: dict[str, str], sheet_directory: Path) -> SampleRecord:
    """Convert one canonical TSV row through Phase 1.2 models."""
    sample_id = values["sample_id"].strip()
    input_paths = _parse_input_paths(values["input"], sheet_directory)
    explicit_type = _parse_input_type(values.get("input_type"))
    dataset = InputDataset.from_files(input_paths, input_type=explicit_type)
    sample = Sample(sample_id=sample_id, input=dataset)
    return SampleRecord(
        sample=sample,
        sex=parse_sex(values.get("sex")),
        father=values.get("father"),
        mother=values.get("mother"),
        phenotype=values.get("phenotype"),
        group=values.get("group"),
    )


def _parse_input_paths(value: str, sheet_directory: Path) -> tuple[Path, ...]:
    """Split semicolon paths and interpret relative paths beside the sheet."""
    if not value.strip():
        raise InputValidationError("input must contain at least one file path.")
    components = value.split(";")
    if any(not component.strip() for component in components):
        raise InputValidationError(
            "input contains an empty semicolon-separated file component."
        )

    parsed: list[Path] = []
    for component in components:
        path = Path(component.strip()).expanduser()
        if not path.is_absolute():
            path = sheet_directory / path
        parsed.append(path)
    return tuple(parsed)


def _parse_input_type(value: str | None) -> InputType | None:
    """Normalize an optional case-insensitive explicit input type."""
    if value is None or not value.strip():
        return None
    normalized = value.strip().casefold()
    try:
        return InputType(normalized)
    except ValueError as error:
        allowed = ", ".join(input_type.value for input_type in InputType)
        raise InputValidationError(
            f"Invalid input_type {value!r}; expected one of: {allowed}."
        ) from error


def _normalize_optional_text(value: str | None) -> str | None:
    """Convert blank/dot metadata to one internal missing representation."""
    if value is None:
        return None
    if not isinstance(value, str):
        raise InputValidationError("Optional sample metadata must be text or None.")
    normalized = value.strip()
    return None if not normalized or normalized == "." else normalized


def _validate_unique_samples(records: tuple[SampleRecord, ...]) -> None:
    """Reject duplicate sample IDs without overwriting earlier records."""
    seen: set[str] = set()
    for record in records:
        sample_id = record.sample.sample_id
        if sample_id in seen:
            raise InputValidationError(
                f"SampleSheet contains duplicate sample_id '{sample_id}'."
            )
        seen.add(sample_id)


def _validate_unique_inputs(records: tuple[SampleRecord, ...]) -> None:
    """Reject normalized input-path reuse across different samples."""
    owners: dict[str, str] = {}
    for record in records:
        sample_id = record.sample.sample_id
        for path in record.sample.input.files:
            identity = os.path.normcase(os.path.normpath(str(path)))
            previous_owner = owners.get(identity)
            if previous_owner is not None:
                raise InputValidationError(
                    f"Input file '{path}' is reused by samples "
                    f"'{previous_owner}' and '{sample_id}'."
                )
            owners[identity] = sample_id


def _validate_parent_relationships(
    records: tuple[SampleRecord, ...],
    records_by_id: dict[str, SampleRecord],
) -> None:
    """Require present parents and consistent declared parent sex metadata."""
    for child in records:
        child_id = child.sample.sample_id
        for role, parent_id, required_sex in (
            ("father", child.father, Sex.MALE),
            ("mother", child.mother, Sex.FEMALE),
        ):
            if parent_id is None:
                continue
            parent = records_by_id.get(parent_id)
            if parent is None:
                raise InputValidationError(
                    f"Sample '{child_id}' {role} '{parent_id}' is not present "
                    "in the SampleSheet."
                )
            if parent.sex not in {None, Sex.UNKNOWN, required_sex}:
                raise InputValidationError(
                    f"Sample '{child_id}' {role} '{parent_id}' has conflicting "
                    f"declared sex '{parent.sex.value}'."
                )


def _validate_pedigree_acyclic(
    records: tuple[SampleRecord, ...],
    records_by_id: dict[str, SampleRecord],
) -> None:
    """Detect cycles across father and mother links with a small DFS."""
    state: dict[str, int] = {}

    def visit(sample_id: str) -> None:
        current_state = state.get(sample_id, 0)
        if current_state == 1:
            raise InputValidationError(
                f"Pedigree cycle detected involving sample '{sample_id}'."
            )
        if current_state == 2:
            return

        state[sample_id] = 1
        record = records_by_id[sample_id]
        for parent_id in (record.father, record.mother):
            if parent_id is not None:
                visit(parent_id)
        state[sample_id] = 2

    for record in records:
        visit(record.sample.sample_id)


__all__ = [
    "OPTIONAL_COLUMNS",
    "REQUIRED_COLUMNS",
    "SUPPORTED_COLUMNS",
    "SampleRecord",
    "SampleSheet",
    "Sex",
    "Trio",
    "parse_sex",
]
