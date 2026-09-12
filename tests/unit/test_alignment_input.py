"""Tests for alignment input classification, validation, and routing.

R2-MAJOR-001 fix: Tests all routing paths with parser-based classification.
R2-MAJOR-002 fix: Tests that non-coordinate-sorted BAM raises clear error.
R2-MAJOR-005 fix: Tests that CRAM is explicitly deferred.
R2-MAJOR-008 fix: Tests exercise all claimed execution paths.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from hifivar.alignment import AlignmentAction, AlignmentOutputFormat
from hifivar.alignment_input import (
    AlignmentHeaderMetadata,
    AlignedStatus,
    InputClassification,
    SortOrder,
    classify_input,
    classify_inputs_batch,
    validate_aligned_bam_for_reuse,
    validate_uBAM_for_alignment,
)
from hifivar.exceptions import InputValidationError, ReferenceError
from hifivar.sample import InputType


class FakeHeaderParser:
    """Deterministic BAM header parser test double."""

    def __init__(self) -> None:
        self.headers: dict[str, AlignmentHeaderMetadata] = {}
        self.default_header = AlignmentHeaderMetadata(
            sort_order=SortOrder.COORDINATE,
            has_sq_lines=True,
            sample_names=("SAMPLE1",),
            read_group_ids=("RG1",),
            has_hd_header=True,
            reference_contigs=("chr1", "chr2"),
        )

    def parse_header(self, path: Path) -> AlignmentHeaderMetadata:
        return self.headers.get(str(path), self.default_header)


def write_fastq(root: Path, name: str = "sample.fastq") -> Path:
    root.mkdir(parents=True, exist_ok=True)
    path = root / name
    path.write_text("@read1\nACGT\n+\nIIII\n", encoding="utf-8")
    return path


def write_bam(root: Path, name: str = "sample.bam") -> Path:
    root.mkdir(parents=True, exist_ok=True)
    path = root / name
    path.write_bytes(b"fake-bam-content")
    return path


def write_cram(root: Path, name: str = "sample.cram") -> Path:
    root.mkdir(parents=True, exist_ok=True)
    path = root / name
    path.write_bytes(b"fake-cram-content")
    return path


def make_header(
    *,
    sort_order: SortOrder = SortOrder.COORDINATE,
    has_sq_lines: bool = True,
    sample_names: tuple[str, ...] = ("SAMPLE1",),
    read_group_ids: tuple[str, ...] = ("RG1",),
    has_hd_header: bool = True,
    reference_contigs: tuple[str, ...] = ("chr1", "chr2"),
) -> AlignmentHeaderMetadata:
    return AlignmentHeaderMetadata(
        sort_order=sort_order,
        has_sq_lines=has_sq_lines,
        sample_names=sample_names,
        read_group_ids=read_group_ids,
        has_hd_header=has_hd_header,
        reference_contigs=reference_contigs,
    )


class TestAlignedStatus:
    def test_aligned_status_values_are_stable(self) -> None:
        assert [s.value for s in AlignedStatus] == [
            "unaligned",
            "aligned",
            "unknown",
        ]


class TestSortOrder:
    def test_sort_order_values_are_stable(self) -> None:
        assert [s.value for s in SortOrder] == [
            "queryname",
            "coordinate",
            "unordered",
            "unknown",
        ]


class TestAlignmentHeaderMetadata:
    def test_valid_metadata_creation(self) -> None:
        header = make_header()
        assert header.sort_order is SortOrder.COORDINATE
        assert header.has_sq_lines is True
        assert header.sample_names == ("SAMPLE1",)
        assert header.unique_sample_names == ("SAMPLE1",)

    def test_unique_sample_names_preserves_order(self) -> None:
        header = make_header(sample_names=("S2", "S1", "S2"))
        assert header.unique_sample_names == ("S2", "S1")

    def test_metadata_serialization(self) -> None:
        header = make_header()
        payload = header.to_dict()
        assert payload["sort_order"] == "coordinate"
        assert payload["has_sq_lines"] is True
        assert payload["sample_names"] == ["SAMPLE1"]
        assert payload["reference_contigs_count"] == 2

    def test_rejects_invalid_sort_order(self) -> None:
        with pytest.raises(InputValidationError, match="SortOrder"):
            AlignmentHeaderMetadata(
                sort_order="coordinate",
                has_sq_lines=True,
                sample_names=("S1",),
                read_group_ids=("RG1",),
                has_hd_header=True,
                reference_contigs=("chr1",),
            )

    def test_rejects_invalid_has_sq_lines(self) -> None:
        with pytest.raises(InputValidationError, match="has_sq_lines"):
            AlignmentHeaderMetadata(
                sort_order=SortOrder.COORDINATE,
                has_sq_lines="yes",
                sample_names=("S1",),
                read_group_ids=("RG1",),
                has_hd_header=True,
                reference_contigs=("chr1",),
            )


class TestClassifyInput:
    def test_fastq_routes_to_align(self, tmp_path: Path) -> None:
        fastq = write_fastq(tmp_path)
        classification = classify_input(fastq)

        assert classification.input_type is InputType.FASTQ
        assert classification.aligned_status is AlignedStatus.UNALIGNED
        assert classification.routing_action is AlignmentAction.ALIGN
        assert classification.header is None

    def test_coordinate_sorted_bam_routes_to_reuse(self, tmp_path: Path) -> None:
        bam = write_bam(tmp_path)
        parser = FakeHeaderParser()
        parser.default_header = make_header(sort_order=SortOrder.COORDINATE)

        classification = classify_input(bam, parser=parser)

        assert classification.input_type is InputType.BAM
        assert classification.aligned_status is AlignedStatus.ALIGNED
        assert classification.routing_action is AlignmentAction.REUSE
        assert classification.header is not None

    def test_ubam_routes_to_align(self, tmp_path: Path) -> None:
        bam = write_bam(tmp_path)
        parser = FakeHeaderParser()
        parser.default_header = make_header(
            sort_order=SortOrder.UNORDERED,
            has_sq_lines=False,
            read_group_ids=("RG1",),
            sample_names=("SAMPLE1",),
        )

        classification = classify_input(bam, parser=parser)

        assert classification.input_type is InputType.BAM
        assert classification.aligned_status is AlignedStatus.UNALIGNED
        assert classification.routing_action is AlignmentAction.ALIGN

    def test_queryname_sorted_bam_raises_error(self, tmp_path: Path) -> None:
        """R2-MAJOR-002: Queryname-sorted BAM must raise error, not route to ALIGN."""
        bam = write_bam(tmp_path)
        parser = FakeHeaderParser()
        parser.default_header = make_header(
            sort_order=SortOrder.QUERYNAME,
            has_sq_lines=True,
        )

        with pytest.raises(InputValidationError, match="sort order.*cannot.*directly reused"):
            classify_input(bam, parser=parser)

    def test_unordered_aligned_bam_raises_error(self, tmp_path: Path) -> None:
        """R2-MAJOR-002: Unordered aligned BAM must raise error."""
        bam = write_bam(tmp_path)
        parser = FakeHeaderParser()
        parser.default_header = make_header(
            sort_order=SortOrder.UNORDERED,
            has_sq_lines=True,
        )

        with pytest.raises(InputValidationError, match="sort order.*cannot.*directly reused"):
            classify_input(bam, parser=parser)

    def test_unknown_sort_with_sq_lines_raises_error(self, tmp_path: Path) -> None:
        """R2-MAJOR-002: Unknown sort with SQ lines must raise error."""
        bam = write_bam(tmp_path)
        parser = FakeHeaderParser()
        parser.default_header = make_header(
            sort_order=SortOrder.UNKNOWN,
            has_sq_lines=True,
        )

        with pytest.raises(InputValidationError, match="sort order.*cannot.*directly reused"):
            classify_input(bam, parser=parser)

    def test_cram_raises_deferred_error(self, tmp_path: Path) -> None:
        """R2-MAJOR-005: CRAM input must be explicitly deferred."""
        cram = write_cram(tmp_path)
        with pytest.raises(InputValidationError, match="CRAM.*not yet supported"):
            classify_input(cram)

    def test_cram_with_parser_still_deferred(self, tmp_path: Path) -> None:
        """R2-MAJOR-005: CRAM is deferred even with parser available."""
        cram = write_cram(tmp_path)
        parser = FakeHeaderParser()
        with pytest.raises(InputValidationError, match="CRAM.*not yet supported"):
            classify_input(cram, parser=parser, reference_contigs=("chr1", "chr2"))

    def test_missing_file_fails(self, tmp_path: Path) -> None:
        with pytest.raises(InputValidationError, match="does not exist"):
            classify_input(tmp_path / "missing.fastq")

    def test_non_file_fails(self, tmp_path: Path) -> None:
        tmp_path.mkdir(parents=True, exist_ok=True)
        with pytest.raises(InputValidationError, match="not a file"):
            classify_input(tmp_path)

    def test_unsupported_suffix_fails(self, tmp_path: Path) -> None:
        path = tmp_path / "sample.xyz"
        path.write_bytes(b"content")
        with pytest.raises(InputValidationError, match="unsupported file suffix"):
            classify_input(path)

    def test_bam_without_parser_fails(self, tmp_path: Path) -> None:
        bam = write_bam(tmp_path)
        with pytest.raises(InputValidationError, match="BamHeaderParser"):
            classify_input(bam)

    def test_sample_identity_validation(self, tmp_path: Path) -> None:
        bam = write_bam(tmp_path)
        parser = FakeHeaderParser()
        parser.default_header = make_header(sample_names=("WRONG_SAMPLE",))

        with pytest.raises(InputValidationError, match="does not match expected"):
            classify_input(bam, expected_sample_id="EXPECTED_SAMPLE", parser=parser)

    def test_matching_sample_identity(self, tmp_path: Path) -> None:
        bam = write_bam(tmp_path)
        parser = FakeHeaderParser()
        parser.default_header = make_header(sample_names=("SAMPLE1",))

        classification = classify_input(
            bam, expected_sample_id="SAMPLE1", parser=parser
        )
        assert classification is not None

    def test_reference_compatibility_check(self, tmp_path: Path) -> None:
        bam = write_bam(tmp_path)
        parser = FakeHeaderParser()
        parser.default_header = make_header(reference_contigs=("chr1", "chrX"))

        with pytest.raises(ReferenceError, match="REFERENCE_CONTIG_MISMATCH"):
            classify_input(
                bam,
                reference_contigs=("chr1", "chr2"),
                parser=parser,
            )

    def test_compatible_reference_contigs(self, tmp_path: Path) -> None:
        bam = write_bam(tmp_path)
        parser = FakeHeaderParser()
        parser.default_header = make_header(reference_contigs=("chr1", "chr2"))

        classification = classify_input(
            bam,
            reference_contigs=("chr1", "chr2", "chr3"),
            parser=parser,
        )
        assert classification is not None

    def test_unicode_path_in_fastq(self, tmp_path: Path) -> None:
        fastq = write_fastq(tmp_path, name="sample_unicode_测试.fastq")
        classification = classify_input(fastq)
        assert classification.input_type is InputType.FASTQ
        assert classification.routing_action is AlignmentAction.ALIGN

    def test_space_path_in_fastq(self, tmp_path: Path) -> None:
        fastq = write_fastq(tmp_path, name="sample with spaces.fastq")
        classification = classify_input(fastq)
        assert classification.input_type is InputType.FASTQ
        assert classification.routing_action is AlignmentAction.ALIGN


class TestClassifyInputsBatch:
    def test_batch_classification(self, tmp_path: Path) -> None:
        fastq1 = write_fastq(tmp_path, "sample1.fastq")
        fastq2 = write_fastq(tmp_path, "sample2.fastq")

        classifications = classify_inputs_batch([fastq1, fastq2])
        assert len(classifications) == 2
        assert all(c.routing_action is AlignmentAction.ALIGN for c in classifications)

    def test_duplicate_path_fails(self, tmp_path: Path) -> None:
        fastq = write_fastq(tmp_path)
        with pytest.raises(InputValidationError, match="Duplicate input path"):
            classify_inputs_batch([fastq, fastq])

    def test_cross_file_sample_identity_conflict(self, tmp_path: Path) -> None:
        bam1 = write_bam(tmp_path, "sample1.bam")
        bam2 = write_bam(tmp_path, "sample2.bam")
        parser = FakeHeaderParser()
        parser.default_header = make_header(sample_names=("SHARED_SAMPLE",))

        with pytest.raises(InputValidationError, match="Duplicate sample identity"):
            classify_inputs_batch(
                [bam1, bam2],
                parser=parser,
            )

    def test_expected_sample_ids_mapping(self, tmp_path: Path) -> None:
        bam = write_bam(tmp_path)
        parser = FakeHeaderParser()
        parser.default_header = make_header(sample_names=("SAMPLE1",))

        classifications = classify_inputs_batch(
            [bam],
            expected_sample_ids={str(bam): "SAMPLE1"},
            parser=parser,
        )
        assert len(classifications) == 1

    def test_batch_rejects_cram(self, tmp_path: Path) -> None:
        """R2-MAJOR-005: Batch classification rejects CRAM."""
        cram = write_cram(tmp_path)
        with pytest.raises(InputValidationError, match="CRAM.*not yet supported"):
            classify_inputs_batch([cram])


class TestValidateUBAMForAlignment:
    def test_valid_ubam(self) -> None:
        classification = InputClassification(
            input_path=Path("test.bam"),
            input_type=InputType.BAM,
            aligned_status=AlignedStatus.UNALIGNED,
            header=make_header(
                sort_order=SortOrder.UNORDERED,
                has_sq_lines=False,
                read_group_ids=("RG1",),
                sample_names=("SAMPLE1",),
            ),
            routing_action=AlignmentAction.ALIGN,
        )
        validate_uBAM_for_alignment(classification)

    def test_rejects_non_bam(self) -> None:
        classification = InputClassification(
            input_path=Path("test.fastq"),
            input_type=InputType.FASTQ,
            aligned_status=AlignedStatus.UNALIGNED,
            header=None,
            routing_action=AlignmentAction.ALIGN,
        )
        with pytest.raises(InputValidationError, match="requires BAM input"):
            validate_uBAM_for_alignment(classification)

    def test_rejects_aligned_status(self) -> None:
        classification = InputClassification(
            input_path=Path("test.bam"),
            input_type=InputType.BAM,
            aligned_status=AlignedStatus.ALIGNED,
            header=make_header(),
            routing_action=AlignmentAction.REUSE,
        )
        with pytest.raises(InputValidationError, match="requires unaligned status"):
            validate_uBAM_for_alignment(classification)

    def test_rejects_coordinate_sorted(self) -> None:
        classification = InputClassification(
            input_path=Path("test.bam"),
            input_type=InputType.BAM,
            aligned_status=AlignedStatus.UNALIGNED,
            header=make_header(sort_order=SortOrder.COORDINATE),
            routing_action=AlignmentAction.ALIGN,
        )
        with pytest.raises(InputValidationError, match="must not be coordinate-sorted"):
            validate_uBAM_for_alignment(classification)

    def test_rejects_no_read_groups(self) -> None:
        classification = InputClassification(
            input_path=Path("test.bam"),
            input_type=InputType.BAM,
            aligned_status=AlignedStatus.UNALIGNED,
            header=make_header(
                sort_order=SortOrder.UNORDERED,
                read_group_ids=(),
                sample_names=(),
            ),
            routing_action=AlignmentAction.ALIGN,
        )
        with pytest.raises(InputValidationError, match="at least one read group"):
            validate_uBAM_for_alignment(classification)

    def test_rejects_no_sample_names(self) -> None:
        classification = InputClassification(
            input_path=Path("test.bam"),
            input_type=InputType.BAM,
            aligned_status=AlignedStatus.UNALIGNED,
            header=make_header(
                sort_order=SortOrder.UNORDERED,
                read_group_ids=("RG1",),
                sample_names=(),
            ),
            routing_action=AlignmentAction.ALIGN,
        )
        with pytest.raises(InputValidationError, match="SM.*tags"):
            validate_uBAM_for_alignment(classification)


class TestValidateAlignedBamForReuse:
    def test_valid_aligned_bam(self) -> None:
        classification = InputClassification(
            input_path=Path("test.bam"),
            input_type=InputType.BAM,
            aligned_status=AlignedStatus.ALIGNED,
            header=make_header(sort_order=SortOrder.COORDINATE),
            routing_action=AlignmentAction.REUSE,
        )
        validate_aligned_bam_for_reuse(classification)

    def test_rejects_unaligned(self) -> None:
        classification = InputClassification(
            input_path=Path("test.bam"),
            input_type=InputType.BAM,
            aligned_status=AlignedStatus.UNALIGNED,
            header=make_header(sort_order=SortOrder.UNORDERED),
            routing_action=AlignmentAction.ALIGN,
        )
        with pytest.raises(InputValidationError, match="requires aligned status"):
            validate_aligned_bam_for_reuse(classification)

    def test_rejects_non_coordinate_sorted(self) -> None:
        classification = InputClassification(
            input_path=Path("test.bam"),
            input_type=InputType.BAM,
            aligned_status=AlignedStatus.ALIGNED,
            header=make_header(sort_order=SortOrder.QUERYNAME),
            routing_action=AlignmentAction.REUSE,
        )
        with pytest.raises(InputValidationError, match="must be coordinate-sorted"):
            validate_aligned_bam_for_reuse(classification)

    def test_rejects_no_sq_lines(self) -> None:
        classification = InputClassification(
            input_path=Path("test.bam"),
            input_type=InputType.BAM,
            aligned_status=AlignedStatus.ALIGNED,
            header=make_header(has_sq_lines=False),
            routing_action=AlignmentAction.REUSE,
        )
        with pytest.raises(InputValidationError, match="must contain SQ header"):
            validate_aligned_bam_for_reuse(classification)

    def test_rejects_no_sample_names(self) -> None:
        classification = InputClassification(
            input_path=Path("test.bam"),
            input_type=InputType.BAM,
            aligned_status=AlignedStatus.ALIGNED,
            header=make_header(sample_names=()),
            routing_action=AlignmentAction.REUSE,
        )
        with pytest.raises(InputValidationError, match="sample identity"):
            validate_aligned_bam_for_reuse(classification)


class TestInputClassificationSerialization:
    def test_classification_serialization(self) -> None:
        classification = InputClassification(
            input_path=Path("test.bam"),
            input_type=InputType.BAM,
            aligned_status=AlignedStatus.ALIGNED,
            header=make_header(),
            routing_action=AlignmentAction.REUSE,
        )
        payload = classification.to_dict()
        assert payload["input_path"] == "test.bam"
        assert payload["input_type"] == "bam"
        assert payload["aligned_status"] == "aligned"
        assert payload["routing_action"] == "reuse"
        assert payload["header"] is not None

    def test_fastq_classification_serialization(self) -> None:
        classification = InputClassification(
            input_path=Path("test.fastq"),
            input_type=InputType.FASTQ,
            aligned_status=AlignedStatus.UNALIGNED,
            header=None,
            routing_action=AlignmentAction.ALIGN,
        )
        payload = classification.to_dict()
        assert payload["header"] is None
