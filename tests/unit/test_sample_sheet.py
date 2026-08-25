"""Tests for the Phase 1.3 TSV sample-sheet and pedigree models."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from hifivar.exceptions import InputValidationError
from hifivar.sample import InputType
from hifivar.sample_sheet import SampleSheet, Sex


VALID_FASTQ = "@read1\nACGT\n+\nIIII\n"
FULL_HEADER = (
    "sample_id\tinput\tinput_type\tsex\tfather\tmother\tphenotype\tgroup"
)


def write_fastq(root: Path, name: str) -> Path:
    """Write one tiny FASTQ relative to a fixture root."""
    path = root / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(VALID_FASTQ, encoding="utf-8")
    return path


def write_alignment(root: Path, name: str) -> Path:
    """Write one non-empty path-level BAM/CRAM placeholder."""
    path = root / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("placeholder", encoding="utf-8")
    return path


def write_sheet(
    root: Path,
    header: str,
    rows: list[str],
    *,
    name: str = "samples.tsv",
    bom: bool = False,
    prefix: str = "",
) -> Path:
    """Write a TSV fixture and return its path."""
    root.mkdir(parents=True, exist_ok=True)
    path = root / name
    text = prefix + header + "\n" + "\n".join(rows) + "\n"
    path.write_text(text, encoding="utf-8-sig" if bom else "utf-8")
    return path


def basic_row(
    sample_id: str,
    input_path: str,
    *,
    input_type: str = "",
    sex: str = "",
    father: str = "",
    mother: str = "",
    phenotype: str = "",
    group: str = "",
) -> str:
    """Return one row matching FULL_HEADER."""
    return "\t".join(
        (
            sample_id,
            input_path,
            input_type,
            sex,
            father,
            mother,
            phenotype,
            group,
        )
    )


def test_single_sample_tsv_loads(tmp_path: Path) -> None:
    """The minimal required schema should create one record."""
    write_fastq(tmp_path, "data/HG002.fastq")
    path = write_sheet(
        tmp_path,
        "sample_id\tinput",
        ["HG002\tdata/HG002.fastq"],
    )

    sheet = SampleSheet.from_tsv(path)

    assert len(sheet) == 1
    assert sheet.source_path == path
    assert sheet.sample_ids == ("HG002",)


def test_multiple_samples_load_in_tsv_order(tmp_path: Path) -> None:
    """Record and sample-ID order should remain deterministic."""
    write_fastq(tmp_path, "data/S2.fastq")
    write_fastq(tmp_path, "data/S1.fastq")
    path = write_sheet(
        tmp_path,
        "input\tsample_id",
        ["data/S2.fastq\tS2", "data/S1.fastq\tS1"],
    )

    sheet = SampleSheet.from_tsv(path)

    assert sheet.sample_ids == ("S2", "S1")
    assert [record.sample.sample_id for record in sheet.records] == ["S2", "S1"]


@pytest.mark.parametrize(
    "header,keyword",
    (
        ("input\tsex", "sample_id"),
        ("sample_id\tsex", "input"),
    ),
)
def test_missing_required_column_is_rejected(
    tmp_path: Path,
    header: str,
    keyword: str,
) -> None:
    """Both required columns must be explicitly named."""
    path = write_sheet(tmp_path, header, ["value\tvalue"])

    with pytest.raises(InputValidationError, match=keyword):
        SampleSheet.from_tsv(path)


def test_duplicate_header_is_rejected(tmp_path: Path) -> None:
    """Duplicate columns cannot silently overwrite metadata."""
    path = write_sheet(
        tmp_path,
        "sample_id\tinput\tinput",
        ["S1\ta.fastq\tb.fastq"],
    )

    with pytest.raises(InputValidationError, match=r"duplicate column"):
        SampleSheet.from_tsv(path)


def test_unknown_column_is_rejected(tmp_path: Path) -> None:
    """Typos and unmodeled metadata columns fail under the strict schema."""
    path = write_sheet(
        tmp_path,
        "sample_id\tinput\tphentype",
        ["S1\ta.fastq\taffected"],
    )

    with pytest.raises(InputValidationError, match=r"unknown.*phentype"):
        SampleSheet.from_tsv(path)


def test_duplicate_sample_id_is_rejected(tmp_path: Path) -> None:
    """Later rows must not overwrite an earlier sample."""
    write_fastq(tmp_path, "a.fastq")
    write_fastq(tmp_path, "b.fastq")
    path = write_sheet(
        tmp_path,
        "sample_id\tinput",
        ["S1\ta.fastq", "S1\tb.fastq"],
    )

    with pytest.raises(InputValidationError, match=r"duplicate sample_id.*S1"):
        SampleSheet.from_tsv(path)


def test_invalid_sample_id_reuses_sample_validation(tmp_path: Path) -> None:
    """The TSV layer must not create a second sample-ID policy."""
    write_fastq(tmp_path, "a.fastq")
    path = write_sheet(
        tmp_path,
        "sample_id\tinput",
        ["Sample 1\ta.fastq"],
    )

    with pytest.raises(InputValidationError, match=r"sample_id"):
        SampleSheet.from_tsv(path)


def test_fastq_input_type_is_inferred(tmp_path: Path) -> None:
    """A blank optional type delegates to InputDataset inference."""
    write_fastq(tmp_path, "sample.fastq")
    path = write_sheet(
        tmp_path,
        FULL_HEADER,
        [basic_row("S1", "sample.fastq")],
    )

    record = SampleSheet.from_tsv(path).get_record("S1")

    assert record.sample.input.input_type is InputType.FASTQ


def test_bam_input_type_is_inferred_at_path_level(tmp_path: Path) -> None:
    """BAM inference retains the Phase 1.2 path-only validation boundary."""
    bam = write_alignment(tmp_path, "sample.bam")
    path = write_sheet(
        tmp_path,
        "sample_id\tinput",
        ["S1\tsample.bam"],
    )

    record = SampleSheet.from_tsv(path).get_record("S1")

    assert record.sample.input.input_type is InputType.BAM
    assert record.sample.input.files == (bam,)


def test_explicit_input_type_is_case_insensitive(tmp_path: Path) -> None:
    """Canonical enum values may be written in uppercase for human TSVs."""
    write_alignment(tmp_path, "sample.cram")
    path = write_sheet(
        tmp_path,
        "sample_id\tinput\tinput_type",
        ["S1\tsample.cram\tCRAM"],
    )

    record = SampleSheet.from_tsv(path).get_record("S1")

    assert record.sample.input.input_type is InputType.CRAM


def test_explicit_input_type_mismatch_is_rejected(tmp_path: Path) -> None:
    """TSV declarations cannot override the observed suffix type."""
    write_alignment(tmp_path, "sample.bam")
    path = write_sheet(
        tmp_path,
        "sample_id\tinput\tinput_type",
        ["S1\tsample.bam\tfastq"],
    )

    with pytest.raises(InputValidationError, match=r"does not match"):
        SampleSheet.from_tsv(path)


def test_unknown_input_type_is_rejected(tmp_path: Path) -> None:
    """The parser should not guess unknown explicit primary-input modes."""
    write_fastq(tmp_path, "sample.fastq")
    path = write_sheet(
        tmp_path,
        "sample_id\tinput\tinput_type",
        ["S1\tsample.fastq\treads"],
    )

    with pytest.raises(InputValidationError, match=r"Invalid input_type"):
        SampleSheet.from_tsv(path)


def test_multiple_fastq_semicolon_paths_preserve_order(tmp_path: Path) -> None:
    """Semicolon-separated HiFi movie files are trimmed and ordered."""
    first = write_fastq(tmp_path, "movie1.fastq")
    second = write_fastq(tmp_path, "movie2.fastq")
    path = write_sheet(
        tmp_path,
        "sample_id\tinput",
        ["S1\t movie1.fastq ; movie2.fastq "],
    )

    dataset = SampleSheet.from_tsv(path).get_record("S1").sample.input

    assert dataset.files == (first, second)


def test_empty_semicolon_input_component_is_rejected(tmp_path: Path) -> None:
    """Empty file components must expose manifest mistakes."""
    write_fastq(tmp_path, "a.fastq")
    write_fastq(tmp_path, "b.fastq")
    path = write_sheet(
        tmp_path,
        "sample_id\tinput",
        ["S1\ta.fastq;;b.fastq"],
    )

    with pytest.raises(InputValidationError, match=r"empty semicolon"):
        SampleSheet.from_tsv(path)


@pytest.mark.parametrize(
    "value,keyword",
    (
        ("a.bam;b.bam", "exactly one"),
        ("a.fastq;a.bam", "Mixed primary"),
    ),
)
def test_invalid_multi_file_primary_inputs_reuse_dataset_rules(
    tmp_path: Path,
    value: str,
    keyword: str,
) -> None:
    """The TSV parser should delegate cardinality and mixing decisions."""
    write_alignment(tmp_path, "a.bam")
    write_alignment(tmp_path, "b.bam")
    write_fastq(tmp_path, "a.fastq")
    path = write_sheet(tmp_path, "sample_id\tinput", [f"S1\t{value}"])

    with pytest.raises(InputValidationError, match=keyword):
        SampleSheet.from_tsv(path)


def test_relative_input_path_is_based_on_sheet_directory(tmp_path: Path) -> None:
    """Sheet portability must not depend on the current process directory."""
    project = tmp_path / "project"
    fastq = write_fastq(project, "data/S1.fastq")
    path = write_sheet(
        project,
        "sample_id\tinput",
        ["S1\tdata/S1.fastq"],
    )

    record = SampleSheet.from_tsv(path).get_record("S1")

    assert record.sample.input.files == (fastq,)


def test_absolute_input_path_is_retained(tmp_path: Path) -> None:
    """Absolute TSV paths should not be prefixed with the sheet directory."""
    fastq = write_fastq(tmp_path / "inputs", "S1.fastq")
    path = write_sheet(
        tmp_path / "sheet",
        "sample_id\tinput",
        [f"S1\t{fastq}"],
    )

    record = SampleSheet.from_tsv(path).get_record("S1")

    assert record.sample.input.files == (fastq,)


def test_input_reuse_across_samples_is_rejected(tmp_path: Path) -> None:
    """One primary file cannot silently belong to two different samples."""
    write_fastq(tmp_path, "shared.fastq")
    path = write_sheet(
        tmp_path,
        "sample_id\tinput",
        ["S1\tshared.fastq", "S2\tshared.fastq"],
    )

    with pytest.raises(InputValidationError, match=r"reused.*S1.*S2"):
        SampleSheet.from_tsv(path)


def test_blank_and_comment_lines_are_ignored(tmp_path: Path) -> None:
    """Human-maintained sheets may contain comments between records."""
    write_fastq(tmp_path, "S1.fastq")
    write_fastq(tmp_path, "S2.fastq")
    path = tmp_path / "samples.tsv"
    path.write_text(
        "# leading comment\n\n"
        "sample_id\tinput\n"
        "S1\tS1.fastq\n"
        "  # middle comment\n\n"
        "S2\tS2.fastq\n",
        encoding="utf-8",
    )

    sheet = SampleSheet.from_tsv(path)

    assert sheet.sample_ids == ("S1", "S2")


def test_utf8_bom_sheet_and_unicode_path_load(tmp_path: Path) -> None:
    """Windows-authored BOM TSVs and Unicode data directories are supported."""
    fastq = write_fastq(tmp_path, "测序数据/S1.fastq")
    path = write_sheet(
        tmp_path,
        "sample_id\tinput",
        ["S1\t测序数据/S1.fastq"],
        bom=True,
    )

    record = SampleSheet.from_tsv(path).get_record("S1")

    assert record.sample.input.files == (fastq,)


@pytest.mark.parametrize(
    "value,expected",
    (
        ("M", Sex.MALE),
        ("male", Sex.MALE),
        ("F", Sex.FEMALE),
        ("female", Sex.FEMALE),
        ("unknown", Sex.UNKNOWN),
        (".", Sex.UNKNOWN),
    ),
)
def test_sex_aliases_are_normalized(
    tmp_path: Path,
    value: str,
    expected: Sex,
) -> None:
    """Declared aliases normalize without any biological inference."""
    write_fastq(tmp_path, "S1.fastq")
    path = write_sheet(
        tmp_path,
        "sample_id\tinput\tsex",
        [f"S1\tS1.fastq\t{value}"],
    )

    assert SampleSheet.from_tsv(path).get_record("S1").sex is expected


def test_invalid_sex_is_rejected(tmp_path: Path) -> None:
    """Sex metadata is a small explicit enum, not an arbitrary string."""
    write_fastq(tmp_path, "S1.fastq")
    path = write_sheet(
        tmp_path,
        "sample_id\tinput\tsex",
        ["S1\tS1.fastq\tpossibly"],
    )

    with pytest.raises(InputValidationError, match=r"Invalid sex"):
        SampleSheet.from_tsv(path)


def test_missing_parent_fields_become_none(tmp_path: Path) -> None:
    """Blank and dot parent markers share one internal representation."""
    write_fastq(tmp_path, "S1.fastq")
    path = write_sheet(
        tmp_path,
        FULL_HEADER,
        [basic_row("S1", "S1.fastq", father=".", mother="")],
    )

    record = SampleSheet.from_tsv(path).get_record("S1")

    assert record.father is None
    assert record.mother is None


@pytest.mark.parametrize("role", ("father", "mother"))
def test_self_parent_is_rejected(tmp_path: Path, role: str) -> None:
    """A sample cannot reference itself as either parent."""
    write_fastq(tmp_path, "S1.fastq")
    fields = {role: "S1"}
    path = write_sheet(
        tmp_path,
        FULL_HEADER,
        [basic_row("S1", "S1.fastq", **fields)],
    )

    with pytest.raises(InputValidationError, match=rf"own {role}"):
        SampleSheet.from_tsv(path)


def test_same_father_and_mother_is_rejected(tmp_path: Path) -> None:
    """The two declared parent roles must not point to one sample."""
    write_fastq(tmp_path, "child.fastq")
    path = write_sheet(
        tmp_path,
        FULL_HEADER,
        [basic_row("child", "child.fastq", father="P1", mother="P1")],
    )

    with pytest.raises(InputValidationError, match=r"same father and mother"):
        SampleSheet.from_tsv(path)


def test_missing_parent_sample_is_rejected(tmp_path: Path) -> None:
    """Non-empty parent IDs must resolve inside the current SampleSheet."""
    write_fastq(tmp_path, "child.fastq")
    path = write_sheet(
        tmp_path,
        FULL_HEADER,
        [basic_row("child", "child.fastq", father="absent")],
    )

    with pytest.raises(InputValidationError, match=r"absent.*not present"):
        SampleSheet.from_tsv(path)


@pytest.mark.parametrize(
    "role,parent_sex,keyword",
    (("father", "female", "father"), ("mother", "male", "mother")),
)
def test_parent_declared_sex_conflict_is_rejected(
    tmp_path: Path,
    role: str,
    parent_sex: str,
    keyword: str,
) -> None:
    """Declared parent roles must agree with non-unknown declared sex."""
    write_fastq(tmp_path, "parent.fastq")
    write_fastq(tmp_path, "child.fastq")
    parent_id = "P1"
    relationship = {role: parent_id}
    path = write_sheet(
        tmp_path,
        FULL_HEADER,
        [
            basic_row(parent_id, "parent.fastq", sex=parent_sex),
            basic_row("child", "child.fastq", **relationship),
        ],
    )

    with pytest.raises(InputValidationError, match=rf"{keyword}.*conflicting"):
        SampleSheet.from_tsv(path)


@pytest.mark.parametrize("parent_sex", ("", "unknown", "."))
def test_unknown_parent_sex_is_allowed(
    tmp_path: Path,
    parent_sex: str,
) -> None:
    """Missing/unknown metadata is not treated as an inferred conflict."""
    write_fastq(tmp_path, "parent.fastq")
    write_fastq(tmp_path, "child.fastq")
    path = write_sheet(
        tmp_path,
        FULL_HEADER,
        [
            basic_row("P1", "parent.fastq", sex=parent_sex),
            basic_row("child", "child.fastq", father="P1"),
        ],
    )

    assert len(SampleSheet.from_tsv(path)) == 2


def test_pedigree_cycle_is_rejected(tmp_path: Path) -> None:
    """A small DFS should detect cycles across parent relationships."""
    write_fastq(tmp_path, "A.fastq")
    write_fastq(tmp_path, "B.fastq")
    path = write_sheet(
        tmp_path,
        FULL_HEADER,
        [
            basic_row("A", "A.fastq", sex="male", father="B"),
            basic_row("B", "B.fastq", sex="male", father="A"),
        ],
    )

    with pytest.raises(InputValidationError, match=r"Pedigree cycle"):
        SampleSheet.from_tsv(path)


def test_valid_trio_is_extracted_in_child_father_mother_order(
    tmp_path: Path,
) -> None:
    """Only complete trios should be returned with explicit role ordering."""
    for sample_id in ("father", "mother", "child"):
        write_fastq(tmp_path, f"{sample_id}.fastq")
    path = write_sheet(
        tmp_path,
        FULL_HEADER,
        [
            basic_row("father", "father.fastq", sex="M"),
            basic_row("mother", "mother.fastq", sex="F"),
            basic_row(
                "child",
                "child.fastq",
                father="father",
                mother="mother",
            ),
        ],
    )

    sheet = SampleSheet.from_tsv(path)
    trio = sheet.get_trios()[0]

    assert tuple(record.sample.sample_id for record in trio) == (
        "child",
        "father",
        "mother",
    )


def test_partial_pedigree_is_valid_but_not_a_complete_trio(tmp_path: Path) -> None:
    """One known parent remains useful metadata without becoming a trio."""
    write_fastq(tmp_path, "father.fastq")
    write_fastq(tmp_path, "child.fastq")
    path = write_sheet(
        tmp_path,
        FULL_HEADER,
        [
            basic_row("father", "father.fastq", sex="male"),
            basic_row("child", "child.fastq", father="father"),
        ],
    )

    sheet = SampleSheet.from_tsv(path)

    assert sheet.get_record("child").father == "father"
    assert sheet.get_trios() == ()


def test_acyclic_multigeneration_pedigree_is_supported(tmp_path: Path) -> None:
    """Pedigrees are not restricted to one trio generation."""
    for sample_id in ("grandfather", "father", "child"):
        write_fastq(tmp_path, f"{sample_id}.fastq")
    path = write_sheet(
        tmp_path,
        FULL_HEADER,
        [
            basic_row("grandfather", "grandfather.fastq", sex="male"),
            basic_row(
                "father",
                "father.fastq",
                sex="male",
                father="grandfather",
            ),
            basic_row("child", "child.fastq", father="father"),
        ],
    )

    sheet = SampleSheet.from_tsv(path)

    assert sheet.get_record("father").father == "grandfather"
    assert sheet.get_record("child").father == "father"


def test_phenotype_and_group_are_open_text_metadata(tmp_path: Path) -> None:
    """Project labels are preserved without a clinical ontology."""
    write_fastq(tmp_path, "S1.fastq")
    path = write_sheet(
        tmp_path,
        FULL_HEADER,
        [
            basic_row(
                "S1",
                "S1.fastq",
                phenotype="rare-disease-candidate",
                group="batch_A",
            )
        ],
    )

    record = SampleSheet.from_tsv(path).get_record("S1")

    assert record.phenotype == "rare-disease-candidate"
    assert record.group == "batch_A"


def test_sheet_serialization_is_json_and_yaml_friendly(tmp_path: Path) -> None:
    """Source, input, pedigree, and cohort metadata use standard types."""
    write_fastq(tmp_path, "S1.fastq")
    path = write_sheet(
        tmp_path,
        FULL_HEADER,
        [basic_row("S1", "S1.fastq", sex="unknown", group="case")],
    )
    sheet = SampleSheet.from_tsv(path)

    metadata = sheet.to_dict()

    assert metadata["source_path"] == str(path)
    assert metadata["records"][0]["sample_id"] == "S1"  # type: ignore[index]
    assert metadata["records"][0]["sex"] == "unknown"  # type: ignore[index]
    assert json.loads(json.dumps(metadata)) == metadata
    assert yaml.safe_load(yaml.safe_dump(metadata)) == metadata


def test_get_record_missing_id_raises_key_error(tmp_path: Path) -> None:
    """Lookup failure has the conventional explicit mapping-style behavior."""
    write_fastq(tmp_path, "S1.fastq")
    path = write_sheet(
        tmp_path,
        "sample_id\tinput",
        ["S1\tS1.fastq"],
    )
    sheet = SampleSheet.from_tsv(path)

    with pytest.raises(KeyError, match="absent"):
        sheet.get_record("absent")


def test_unicode_machine_sample_id_remains_rejected(tmp_path: Path) -> None:
    """SampleSheet must retain the Phase 1.2 ASCII machine-ID policy."""
    write_fastq(tmp_path, "sample.fastq")
    path = write_sheet(
        tmp_path,
        "sample_id\tinput",
        ["样本001\tsample.fastq"],
    )

    with pytest.raises(InputValidationError, match=r"sample_id"):
        SampleSheet.from_tsv(path)


def test_header_only_sheet_is_rejected(tmp_path: Path) -> None:
    """A schema without records is not a usable SampleSheet."""
    path = write_sheet(tmp_path, "sample_id\tinput", [])

    with pytest.raises(InputValidationError, match=r"at least one sample"):
        SampleSheet.from_tsv(path)


def test_non_tsv_suffix_is_rejected(tmp_path: Path) -> None:
    """Phase 1.3 supports one manifest format only."""
    path = write_sheet(
        tmp_path,
        "sample_id\tinput",
        ["S1\tS1.fastq"],
        name="samples.csv",
    )

    with pytest.raises(InputValidationError, match=r"\.tsv suffix"):
        SampleSheet.from_tsv(path)
