"""Tests for the Phase 1.2 sample and primary-input data models."""

from __future__ import annotations

import gzip
import json
from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest
import yaml

from hifivar.exceptions import InputValidationError
from hifivar.sample import InputDataset, InputType, Sample, infer_input_type


VALID_FASTQ = "@read1\nACGT\n+\nIIII\n"


def write_text(path: Path, content: str = "placeholder") -> Path:
    """Write one UTF-8 fixture and return its path."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def write_fastq(path: Path, content: str = VALID_FASTQ) -> Path:
    """Write one plain or gzip-compressed tiny FASTQ fixture."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.name.lower().endswith(".gz"):
        with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
            handle.write(content)
    else:
        path.write_text(content, encoding="utf-8")
    return path


def test_input_type_contains_only_primary_sequence_inputs() -> None:
    """The enum should not grow derived or variant-file types prematurely."""
    assert [(item.name, item.value) for item in InputType] == [
        ("FASTQ", "fastq"),
        ("BAM", "bam"),
        ("CRAM", "cram"),
    ]


def test_single_fastq_dataset(tmp_path: Path) -> None:
    """One long-read FASTQ is a valid primary input dataset."""
    fastq = write_fastq(tmp_path / "movie1.fastq")

    dataset = InputDataset.from_files([fastq])

    assert dataset.input_type is InputType.FASTQ
    assert dataset.files == (fastq,)


@pytest.mark.parametrize("suffix", (".fastq.gz", ".fq.gz"))
def test_gzipped_fastq_dataset(tmp_path: Path, suffix: str) -> None:
    """Both supported gzip FASTQ suffixes reuse streaming validation."""
    fastq = write_fastq(tmp_path / f"movie1{suffix}")

    dataset = InputDataset.from_files([fastq])

    assert dataset.input_type is InputType.FASTQ


def test_multiple_fastq_files_preserve_input_order(tmp_path: Path) -> None:
    """Multiple HiFi movie files remain one ordered FASTQ dataset."""
    movie2 = write_fastq(tmp_path / "movie2.fastq.gz")
    movie1 = write_fastq(tmp_path / "movie1.fastq.gz")

    dataset = InputDataset.from_files([movie2, movie1])

    assert dataset.files == (movie2, movie1)


def test_empty_input_file_list_is_rejected() -> None:
    """A primary input dataset must contain at least one file."""
    with pytest.raises(InputValidationError, match=r"at least one"):
        InputDataset.from_files([])


def test_duplicate_fastq_path_is_rejected(tmp_path: Path) -> None:
    """Manifest mistakes must not be hidden by silent de-duplication."""
    fastq = write_fastq(tmp_path / "movie1.fastq")

    with pytest.raises(InputValidationError, match=r"duplicate.*movie1"):
        InputDataset.from_files([fastq, fastq])


def test_malformed_fastq_reuses_lightweight_validation(tmp_path: Path) -> None:
    """Invalid first FASTQ records fail through the Phase 0 API."""
    fastq = write_fastq(tmp_path / "bad.fastq", "not-a-header\nACGT\n+\nIIII\n")

    with pytest.raises(InputValidationError, match=r"header.*@"):
        InputDataset.from_files([fastq])


def test_bam_placeholder_uses_path_level_validation(tmp_path: Path) -> None:
    """A non-empty BAM placeholder proves only the documented path boundary."""
    bam = write_text(tmp_path / "sample.bam")

    dataset = InputDataset.from_files([bam])

    assert dataset.input_type is InputType.BAM
    assert dataset.files == (bam,)


def test_cram_placeholder_uses_path_level_validation(tmp_path: Path) -> None:
    """A non-empty CRAM placeholder is not a binary-integrity claim."""
    cram = write_text(tmp_path / "sample.cram")

    dataset = InputDataset.from_files([cram])

    assert dataset.input_type is InputType.CRAM
    assert dataset.files == (cram,)


@pytest.mark.parametrize("suffix", (".bam", ".cram"))
def test_multiple_alignment_inputs_are_rejected(
    tmp_path: Path,
    suffix: str,
) -> None:
    """Alignment merging must be an explicit future workflow operation."""
    first = write_text(tmp_path / f"first{suffix}")
    second = write_text(tmp_path / f"second{suffix}")

    with pytest.raises(InputValidationError, match=r"exactly one"):
        InputDataset.from_files([first, second])


def test_fastq_and_bam_mixed_primary_input_is_rejected(tmp_path: Path) -> None:
    """One dataset cannot silently combine raw and aligned inputs."""
    fastq = write_fastq(tmp_path / "sample.fastq")
    bam = write_text(tmp_path / "sample.bam")

    with pytest.raises(InputValidationError, match=r"Mixed primary"):
        InputDataset.from_files([fastq, bam])


def test_bam_and_cram_mixed_primary_input_is_rejected(tmp_path: Path) -> None:
    """BAM and CRAM are separate primary modes even though both are aligned."""
    bam = write_text(tmp_path / "sample.bam")
    cram = write_text(tmp_path / "sample.cram")

    with pytest.raises(InputValidationError, match=r"Mixed primary"):
        InputDataset.from_files([bam, cram])


@pytest.mark.parametrize("name", ("unknown.txt", "sample.bam.gz", "sample.cram.gz"))
def test_unknown_or_nonstandard_suffix_is_rejected(
    tmp_path: Path,
    name: str,
) -> None:
    """Inference should fail rather than guess unsupported formats."""
    path = write_text(tmp_path / name)

    with pytest.raises(InputValidationError, match=r"unsupported file suffix"):
        InputDataset.from_files([path])


def test_explicit_matching_input_type_is_accepted(tmp_path: Path) -> None:
    """An explicit enum is allowed when it agrees with every suffix."""
    fastq = write_fastq(tmp_path / "sample.fq")

    dataset = InputDataset.from_files(
        [fastq],
        input_type=InputType.FASTQ,
    )

    assert dataset.input_type is InputType.FASTQ


def test_explicit_input_type_mismatch_is_rejected(tmp_path: Path) -> None:
    """Explicit metadata cannot override the observed filename type."""
    bam = write_text(tmp_path / "sample.bam")

    with pytest.raises(InputValidationError, match=r"does not match"):
        InputDataset.from_files([bam], input_type=InputType.FASTQ)


def test_bare_string_input_type_is_rejected(tmp_path: Path) -> None:
    """Callers should use InputType rather than spreading raw type strings."""
    fastq = write_fastq(tmp_path / "sample.fastq")

    with pytest.raises(InputValidationError, match=r"InputType"):
        InputDataset.from_files(
            [fastq],
            input_type="fastq",  # type: ignore[arg-type]
        )


@pytest.mark.parametrize(
    "sample_id",
    ("HG002", "NA12878", "sample_001", "family1.child", "sample-1"),
)
def test_valid_ascii_safe_sample_ids(
    tmp_path: Path,
    sample_id: str,
) -> None:
    """Documented machine-safe identifiers should remain unchanged."""
    dataset = InputDataset.from_files([write_fastq(tmp_path / "sample.fastq")])

    sample = Sample(sample_id=sample_id, input=dataset)

    assert sample.sample_id == sample_id


@pytest.mark.parametrize(
    "sample_id",
    (
        "",
        "   ",
        "sample 1",
        "sample\t1",
        "sample\n1",
        "sample/1",
        "sample\\1",
        "..",
        "../HG002",
        ".sample",
        "sample:1",
        "样本001",
    ),
)
def test_invalid_sample_ids_are_rejected(
    tmp_path: Path,
    sample_id: str,
) -> None:
    """Unsafe IDs fail instead of being rewritten into a different identity."""
    dataset = InputDataset.from_files([write_fastq(tmp_path / "sample.fastq")])

    with pytest.raises(InputValidationError, match=r"sample_id"):
        Sample(sample_id=sample_id, input=dataset)


def test_sample_id_is_not_silently_normalized(tmp_path: Path) -> None:
    """Spaces are not automatically changed to underscores."""
    dataset = InputDataset.from_files([write_fastq(tmp_path / "sample.fastq")])

    with pytest.raises(InputValidationError):
        Sample(sample_id="Sample 1", input=dataset)


def test_unicode_input_path_with_ascii_sample_id(tmp_path: Path) -> None:
    """Unicode filesystem paths remain distinct from machine ID policy."""
    fastq = write_fastq(tmp_path / "测序数据" / "HG002.fastq.gz")

    sample = Sample(
        sample_id="HG002",
        input=InputDataset.from_files([fastq]),
    )

    assert sample.input.files == (fastq,)


def test_fastq_sample_serialization_is_json_and_yaml_friendly(
    tmp_path: Path,
) -> None:
    """Paths and enum values should become ordinary serialization types."""
    movie1 = write_fastq(tmp_path / "movie1.fastq.gz")
    movie2 = write_fastq(tmp_path / "movie2.fastq.gz")
    sample = Sample(
        sample_id="HG002",
        input=InputDataset.from_files([movie1, movie2]),
    )

    metadata = sample.to_dict()

    assert metadata == {
        "sample_id": "HG002",
        "input": {
            "type": "fastq",
            "files": [str(movie1), str(movie2)],
        },
    }
    assert json.loads(json.dumps(metadata)) == metadata
    assert yaml.safe_load(yaml.safe_dump(metadata)) == metadata


def test_bam_serialization_is_json_and_yaml_friendly(tmp_path: Path) -> None:
    """Single-alignment inputs share the same stable representation."""
    bam = write_text(tmp_path / "HG002.bam")
    dataset = InputDataset.from_files([bam])

    metadata = dataset.to_dict()

    assert metadata == {"type": "bam", "files": [str(bam)]}
    json.dumps(metadata)
    yaml.safe_dump(metadata)


def test_relative_input_path_is_not_made_absolute(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Serialization reflects stored paths rather than inventing provenance."""
    monkeypatch.chdir(tmp_path)
    relative = write_fastq(Path("data") / "sample.fastq")

    dataset = InputDataset.from_files([relative])

    assert dataset.files == (relative,)
    assert dataset.to_dict()["files"] == [str(relative)]


def test_caller_owned_file_list_cannot_mutate_dataset(tmp_path: Path) -> None:
    """Factory construction copies mutable input collections into a tuple."""
    first = write_fastq(tmp_path / "movie1.fastq")
    second = write_fastq(tmp_path / "movie2.fastq")
    supplied = [first]

    dataset = InputDataset.from_files(supplied)
    supplied.append(second)

    assert dataset.files == (first,)


def test_dataset_and_sample_are_frozen(tmp_path: Path) -> None:
    """Primary identity and input metadata cannot change during a run."""
    dataset = InputDataset.from_files([write_fastq(tmp_path / "sample.fastq")])
    sample = Sample(sample_id="HG002", input=dataset)

    with pytest.raises(FrozenInstanceError):
        dataset.input_type = InputType.BAM  # type: ignore[misc]
    with pytest.raises(FrozenInstanceError):
        sample.sample_id = "HG003"  # type: ignore[misc]


def test_bam_index_is_optional_during_dataset_creation(tmp_path: Path) -> None:
    """The data model can represent a BAM before an explicit indexing step."""
    bam = write_text(tmp_path / "sample.bam")

    dataset = InputDataset.from_files([bam])

    assert dataset.input_type is InputType.BAM
    assert not Path(f"{bam}.bai").exists()


@pytest.mark.parametrize(
    "data_name,index_name",
    (
        ("sample.bam", "sample.bam.bai"),
        ("sample.cram", "sample.cram.crai"),
    ),
)
def test_explicit_alignment_index_validation_succeeds(
    tmp_path: Path,
    data_name: str,
    index_name: str,
) -> None:
    """Future callers may explicitly require a present BAM/CRAM index."""
    data_path = write_text(tmp_path / data_name)
    write_text(tmp_path / index_name, "index-placeholder")
    dataset = InputDataset.from_files([data_path])

    dataset.validate_index()


def test_explicit_alignment_index_validation_reports_missing(
    tmp_path: Path,
) -> None:
    """Index validation reports absence and never creates an index."""
    bam = write_text(tmp_path / "sample.bam")
    dataset = InputDataset.from_files([bam])

    with pytest.raises(InputValidationError, match=r"index missing"):
        dataset.validate_index()

    assert not Path(f"{bam}.bai").exists()


def test_fastq_index_validation_is_not_applicable(tmp_path: Path) -> None:
    """FASTQ should not acquire an artificial index concept."""
    dataset = InputDataset.from_files([write_fastq(tmp_path / "sample.fastq")])

    with pytest.raises(InputValidationError, match=r"not applicable"):
        dataset.validate_index()


@pytest.mark.parametrize(
    "name,expected",
    (
        ("sample.fastq", InputType.FASTQ),
        ("sample.fq.gz", InputType.FASTQ),
        ("sample.bam", InputType.BAM),
        ("sample.cram", InputType.CRAM),
    ),
)
def test_input_type_inference(name: str, expected: InputType) -> None:
    """Inference uses only the documented unambiguous suffix table."""
    assert infer_input_type([Path(name)]) is expected


def test_bare_path_is_not_treated_as_iterable_files(tmp_path: Path) -> None:
    """One input still needs a collection to avoid string-character parsing."""
    fastq = write_fastq(tmp_path / "sample.fastq")

    with pytest.raises(InputValidationError, match=r"bare path"):
        InputDataset.from_files(fastq)
