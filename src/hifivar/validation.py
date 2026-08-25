"""Lightweight, non-mutating input validation for HiFiVar.

This module deliberately checks only properties that can be established with
the Python standard library. It never creates indexes, normalizes records, or
renames contigs. Tool-aware structural validation belongs in later wrappers.
"""

from __future__ import annotations

import gzip
import hashlib
import os
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import TextIO

from hifivar.exceptions import (
    InputValidationError,
    OutputValidationError,
    ReferenceError,
)
from hifivar.logging_utils import get_logger


PathInput = str | Path
CHECKSUM_CHUNK_SIZE = 8 * 1024 * 1024

_FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")
_FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
_VCF_SUFFIXES = (".vcf", ".vcf.gz")
_BED_SUFFIXES = (".bed", ".bed.gz")
_ALIGNMENT_SUFFIXES = (".bam", ".cram")
_LOGGER = get_logger(__name__)


def validate_file(
    path: PathInput,
    *,
    require_nonempty: bool = True,
) -> Path:
    """Validate a readable input file and return its un-resolved path.

    Normal symbolic links are accepted. Broken symbolic links fail because
    their targets do not exist. The path is intentionally not resolved so that
    diagnostics retain the path supplied by the caller.
    """
    input_path = _coerce_path(path, "Input file")
    _LOGGER.debug("Validating input file: %s", input_path)
    _validate_existing_file(
        input_path,
        require_nonempty=require_nonempty,
        error_type=InputValidationError,
        label="Input file",
    )
    _LOGGER.debug("Input file validation succeeded: %s", input_path)
    return input_path


def validate_directory(
    path: PathInput,
    *,
    must_exist: bool = True,
    writable: bool = False,
) -> Path:
    """Validate a directory without creating it.

    ``must_exist=False`` permits a missing future output directory but still
    rejects an existing non-directory or a broken symbolic link.
    """
    directory = _coerce_path(path, "Directory")
    _LOGGER.debug("Validating directory: %s", directory)

    if directory.is_symlink() and not directory.exists():
        raise InputValidationError(
            f"Directory '{directory}' is a broken symbolic link."
        )
    if not directory.exists():
        if must_exist:
            raise InputValidationError(f"Directory '{directory}' is missing.")
        return directory
    if not directory.is_dir():
        raise InputValidationError(f"Path '{directory}' is not a directory.")
    if not os.access(directory, os.R_OK):
        raise InputValidationError(f"Directory '{directory}' is not readable.")
    if writable and not os.access(directory, os.W_OK):
        raise InputValidationError(f"Directory '{directory}' is not writable.")

    _LOGGER.debug("Directory validation succeeded: %s", directory)
    return directory


def validate_output_file(
    path: PathInput,
    *,
    require_nonempty: bool = True,
) -> Path:
    """Validate an expected output using ``OutputValidationError`` on failure."""
    output_path = _coerce_path(path, "Output file", OutputValidationError)
    _LOGGER.debug("Validating output file: %s", output_path)
    _validate_existing_file(
        output_path,
        require_nonempty=require_nonempty,
        error_type=OutputValidationError,
        label="Output file",
    )
    _LOGGER.debug("Output file validation succeeded: %s", output_path)
    return output_path


def validate_fasta(path: PathInput, *, require_fai: bool = False) -> Path:
    """Lightly validate the first FASTA record and optionally its FAI index."""
    fasta_path = validate_file(path)
    _require_suffix(fasta_path, _FASTA_SUFFIXES, "FASTA")
    _LOGGER.debug("Validating FASTA header and first record: %s", fasta_path)

    try:
        with _open_text_auto(fasta_path) as handle:
            first_content = _first_nonempty_line(handle)
            if first_content is None or not first_content.startswith(">"):
                raise InputValidationError(
                    f"FASTA file '{fasta_path}' is missing a valid first header."
                )

            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    raise InputValidationError(
                        f"FASTA file '{fasta_path}' has a header with no sequence."
                    )
                break
            else:
                raise InputValidationError(
                    f"FASTA file '{fasta_path}' has a header with no sequence."
                )
    except (OSError, UnicodeError) as error:
        raise InputValidationError(
            f"Unable to read FASTA file '{fasta_path}': {error}"
        ) from error

    if require_fai:
        validate_fasta_index(fasta_path)
    _LOGGER.debug("FASTA validation succeeded: %s", fasta_path)
    return fasta_path


def validate_fasta_index(fasta_path: PathInput) -> Path:
    """Validate the conventional ``<fasta>.fai`` path and basic FAI content."""
    reference_path = _coerce_path(fasta_path, "FASTA file", ReferenceError)
    index_path = Path(f"{reference_path}.fai")
    if not index_path.exists():
        raise ReferenceError(
            f"FASTA index missing for '{reference_path}': expected '{index_path}'."
        )
    read_fai_contigs(index_path)
    return index_path


def read_fai_contigs(path: PathInput) -> dict[str, int]:
    """Stream a FAI file and return unique contig names with their lengths."""
    index_path = _coerce_path(path, "FASTA index", ReferenceError)
    _validate_existing_file(
        index_path,
        require_nonempty=True,
        error_type=ReferenceError,
        label="FASTA index",
    )
    contigs: dict[str, int] = {}

    try:
        with index_path.open("rt", encoding="utf-8", newline="") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                fields = raw_line.rstrip("\r\n").split("\t")
                if len(fields) < 5:
                    raise ReferenceError(
                        f"FASTA index '{index_path}' line {line_number} has fewer "
                        "than 5 tab-separated columns."
                    )
                contig = fields[0]
                if not contig:
                    raise ReferenceError(
                        f"FASTA index '{index_path}' line {line_number} has an "
                        "empty contig name."
                    )
                if contig in contigs:
                    raise ReferenceError(
                        f"FASTA index '{index_path}' contains duplicate contig "
                        f"'{contig}'."
                    )
                try:
                    length, offset, line_bases, line_width = map(int, fields[1:5])
                except ValueError as error:
                    raise ReferenceError(
                        f"FASTA index '{index_path}' line {line_number} contains "
                        "non-integer coordinate fields."
                    ) from error
                if length <= 0 or offset < 0 or line_bases <= 0 or line_width <= 0:
                    raise ReferenceError(
                        f"FASTA index '{index_path}' line {line_number} contains "
                        "invalid coordinate values."
                    )
                contigs[contig] = length
    except ReferenceError:
        raise
    except (OSError, UnicodeError) as error:
        raise ReferenceError(
            f"Unable to read FASTA index '{index_path}': {error}"
        ) from error

    return contigs


def validate_fastq(path: PathInput) -> Path:
    """Lightly validate only the first four-line FASTQ record."""
    fastq_path = validate_file(path)
    _require_suffix(fastq_path, _FASTQ_SUFFIXES, "FASTQ")
    _LOGGER.debug("Validating first FASTQ record: %s", fastq_path)

    try:
        with _open_text_auto(fastq_path) as handle:
            record = [handle.readline() for _ in range(4)]
    except (OSError, UnicodeError) as error:
        raise InputValidationError(
            f"Unable to read FASTQ file '{fastq_path}': {error}"
        ) from error

    if any(line == "" for line in record):
        raise InputValidationError(
            f"FASTQ file '{fastq_path}' has an incomplete first record."
        )
    header, sequence, separator, quality = (
        line.rstrip("\r\n") for line in record
    )
    if not header.startswith("@"):
        raise InputValidationError(
            f"FASTQ file '{fastq_path}' first header must start with '@'."
        )
    if not separator.startswith("+"):
        raise InputValidationError(
            f"FASTQ file '{fastq_path}' first separator must start with '+'."
        )
    if not sequence:
        raise InputValidationError(
            f"FASTQ file '{fastq_path}' first sequence is empty."
        )
    if len(sequence) != len(quality):
        raise InputValidationError(
            f"FASTQ file '{fastq_path}' first sequence and quality lengths differ."
        )

    _LOGGER.debug("FASTQ validation succeeded: %s", fastq_path)
    return fastq_path


def validate_alignment_file(
    path: PathInput,
    *,
    require_index: bool = False,
) -> Path:
    """Validate a BAM/CRAM path and optionally check index presence.

    This function does not parse binary alignment structure or headers.
    """
    alignment_path = validate_file(path)
    _require_suffix(alignment_path, _ALIGNMENT_SUFFIXES, "alignment")

    if require_index:
        suffix = alignment_path.suffix.lower()
        index_suffix = ".bai" if suffix == ".bam" else ".crai"
        candidates = (
            Path(f"{alignment_path}{index_suffix}"),
            alignment_path.with_suffix(index_suffix),
        )
        _require_existing_index(alignment_path, candidates, "alignment")

    _LOGGER.debug("Alignment path validation succeeded: %s", alignment_path)
    return alignment_path


def validate_vcf(path: PathInput, *, require_index: bool = False) -> Path:
    """Stream and validate required VCF headers plus optional index presence."""
    vcf_path = validate_file(path)
    _require_suffix(vcf_path, _VCF_SUFFIXES, "VCF")
    _LOGGER.debug("Validating VCF headers: %s", vcf_path)
    has_fileformat = False
    has_columns = False

    try:
        with _open_text_auto(vcf_path) as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\r\n")
                if line.startswith("##fileformat=VCF"):
                    has_fileformat = True
                elif line == "#CHROM" or line.startswith("#CHROM\t"):
                    has_columns = True
                    break
                elif not line.startswith("#"):
                    break
    except (OSError, UnicodeError) as error:
        raise InputValidationError(
            f"Unable to read VCF file '{vcf_path}': {error}"
        ) from error

    if not has_fileformat:
        raise InputValidationError(
            f"VCF file '{vcf_path}' is missing the ##fileformat=VCF header."
        )
    if not has_columns:
        raise InputValidationError(
            f"VCF file '{vcf_path}' is missing the #CHROM column header."
        )

    if require_index:
        if not vcf_path.name.lower().endswith(".vcf.gz"):
            raise InputValidationError(
                f"VCF index checking requires a .vcf.gz file: '{vcf_path}'."
            )
        candidates = (Path(f"{vcf_path}.tbi"), Path(f"{vcf_path}.csi"))
        _require_existing_index(vcf_path, candidates, "VCF")

    _LOGGER.debug("VCF validation succeeded: %s", vcf_path)
    return vcf_path


def validate_bed(path: PathInput) -> Path:
    """Stream all BED records and validate the first three columns."""
    bed_path = validate_file(path)
    _require_suffix(bed_path, _BED_SUFFIXES, "BED")
    _LOGGER.debug("Validating BED records: %s", bed_path)
    data_lines = 0

    try:
        with _open_text_auto(bed_path) as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                lowered = line.lower()
                if lowered.startswith("track ") or lowered == "track":
                    continue
                if lowered.startswith("browser ") or lowered == "browser":
                    continue

                fields = line.split("\t")
                if len(fields) < 3:
                    raise InputValidationError(
                        f"BED file '{bed_path}' line {line_number} has fewer "
                        "than 3 tab-separated columns."
                    )
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                except ValueError as error:
                    raise InputValidationError(
                        f"BED file '{bed_path}' line {line_number} has invalid "
                        "integer coordinates."
                    ) from error
                if start < 0:
                    raise InputValidationError(
                        f"BED file '{bed_path}' line {line_number} has a negative "
                        "start coordinate."
                    )
                if end <= start:
                    raise InputValidationError(
                        f"BED file '{bed_path}' line {line_number} must have "
                        "end greater than start."
                    )
                data_lines += 1
    except InputValidationError:
        raise
    except (OSError, UnicodeError) as error:
        raise InputValidationError(
            f"Unable to read BED file '{bed_path}': {error}"
        ) from error

    if data_lines == 0:
        raise InputValidationError(
            f"BED file '{bed_path}' contains no data records."
        )
    _LOGGER.debug("BED validation succeeded: %s", bed_path)
    return bed_path


def validate_contig_compatibility(
    reference_contigs: Iterable[str],
    query_contigs: Iterable[str],
) -> None:
    """Require query contigs to be a subset of unique reference contigs."""
    reference_list = list(reference_contigs)
    query_list = list(query_contigs)
    _validate_contig_names(reference_list, "reference")
    _validate_contig_names(query_list, "query", allow_duplicates=True)

    if not reference_list:
        raise ReferenceError("Reference contig collection is empty.")
    if not query_list:
        raise ReferenceError("Query contig collection is empty.")
    if len(set(reference_list)) != len(reference_list):
        raise ReferenceError("Reference contig collection contains duplicates.")

    missing = sorted(set(query_list) - set(reference_list))
    if missing:
        missing_text = ", ".join(missing)
        raise ReferenceError(
            "REFERENCE_CONTIG_MISMATCH: query contigs are absent from the "
            f"reference: {missing_text}. Normalize inputs explicitly; HiFiVar "
            "does not rename contigs automatically."
        )


def compute_sha256(path: PathInput) -> str:
    """Compute a SHA256 digest using fixed-size streaming chunks."""
    input_path = validate_file(path, require_nonempty=False)
    digest = hashlib.sha256()
    _LOGGER.debug("Computing SHA256 checksum: %s", input_path)

    try:
        with input_path.open("rb") as handle:
            while chunk := handle.read(CHECKSUM_CHUNK_SIZE):
                digest.update(chunk)
    except OSError as error:
        raise InputValidationError(
            f"Unable to read input file '{input_path}' for checksum: {error}"
        ) from error

    checksum = digest.hexdigest()
    _LOGGER.debug("SHA256 checksum completed: %s", input_path)
    return checksum


def _coerce_path(
    path: PathInput,
    label: str,
    error_type: type[InputValidationError] | type[OutputValidationError]
    | type[ReferenceError] = InputValidationError,
) -> Path:
    """Convert supported path inputs without resolving symbolic links."""
    if not isinstance(path, (str, Path)):
        raise error_type(f"{label} path must be a string or Path.")
    if isinstance(path, str) and not path.strip():
        raise error_type(f"{label} path must not be empty.")
    return Path(path).expanduser()


def _validate_existing_file(
    path: Path,
    *,
    require_nonempty: bool,
    error_type: type[InputValidationError] | type[OutputValidationError]
    | type[ReferenceError],
    label: str,
) -> None:
    """Apply common existence, type, size, and readability checks."""
    if path.is_symlink() and not path.exists():
        raise error_type(f"{label} '{path}' is a broken symbolic link.")
    if not path.exists():
        raise error_type(f"{label} '{path}' is missing.")
    if not path.is_file():
        raise error_type(f"{label} '{path}' is not a file.")
    try:
        size = path.stat().st_size
    except OSError as error:
        raise error_type(f"Unable to inspect {label.lower()} '{path}': {error}") from error
    if require_nonempty and size == 0:
        raise error_type(f"{label} '{path}' is empty.")
    if not os.access(path, os.R_OK):
        raise error_type(f"{label} '{path}' is not readable.")


def _require_suffix(path: Path, supported: tuple[str, ...], label: str) -> None:
    """Require one of a fixed set of case-insensitive filename suffixes."""
    lowered_name = path.name.lower()
    if not any(lowered_name.endswith(suffix) for suffix in supported):
        suffix_text = ", ".join(supported)
        raise InputValidationError(
            f"Unsupported {label} suffix for '{path}'. Expected one of: "
            f"{suffix_text}."
        )


def _require_existing_index(
    data_path: Path,
    candidates: tuple[Path, ...],
    label: str,
) -> Path:
    """Return the first deterministic readable index candidate."""
    for candidate in candidates:
        if candidate.exists() and candidate.is_file() and os.access(candidate, os.R_OK):
            return candidate
    expected = " or ".join(f"'{candidate}'" for candidate in candidates)
    raise InputValidationError(
        f"{label} index missing for '{data_path}': expected {expected}."
    )


@contextmanager
def _open_text_auto(path: Path) -> Iterator[TextIO]:
    """Open plain or gzip-compressed UTF-8 text without reading it eagerly."""
    if path.name.lower().endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            yield handle
    else:
        with path.open("rt", encoding="utf-8", newline="") as handle:
            yield handle


def _first_nonempty_line(handle: TextIO) -> str | None:
    """Return the first non-empty stripped line from a text stream."""
    for raw_line in handle:
        line = raw_line.strip()
        if line:
            return line
    return None


def _validate_contig_names(
    contigs: list[str],
    label: str,
    *,
    allow_duplicates: bool = False,
) -> None:
    """Validate explicit contig collection element types and values."""
    for contig in contigs:
        if not isinstance(contig, str) or not contig:
            raise ReferenceError(
                f"The {label} contig collection contains an invalid name."
            )
    if not allow_duplicates and len(set(contigs)) != len(contigs):
        raise ReferenceError(
            f"The {label} contig collection contains duplicate contigs."
        )


__all__ = [
    "CHECKSUM_CHUNK_SIZE",
    "compute_sha256",
    "read_fai_contigs",
    "validate_alignment_file",
    "validate_bed",
    "validate_contig_compatibility",
    "validate_directory",
    "validate_fasta",
    "validate_fasta_index",
    "validate_fastq",
    "validate_file",
    "validate_output_file",
    "validate_vcf",
]
