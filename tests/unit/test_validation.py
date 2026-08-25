"""Tests for the Phase 0.8 lightweight validation layer."""

from __future__ import annotations

import gzip
import hashlib
import logging
from pathlib import Path

import pytest

from hifivar.exceptions import (
    InputValidationError,
    OutputValidationError,
    ReferenceError,
)
from hifivar.validation import (
    CHECKSUM_CHUNK_SIZE,
    compute_sha256,
    read_fai_contigs,
    validate_alignment_file,
    validate_bed,
    validate_contig_compatibility,
    validate_directory,
    validate_fasta,
    validate_fasta_index,
    validate_fastq,
    validate_file,
    validate_output_file,
    validate_vcf,
)


VALID_FASTA = ">chr1\nACGTACGT\n>chr2\nAAAA\n"
VALID_FAI = "chr1\t8\t6\t8\t9\nchr2\t4\t21\t4\t5\n"
VALID_FASTQ = "@read1\nACGT\n+\nIIII\n"
VALID_VCF = (
    "##fileformat=VCFv4.3\n"
    "##contig=<ID=chr1,length=1000>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "chr1\t10\t.\tA\tG\t.\tPASS\t.\n"
)


def write_text(path: Path, content: str) -> Path:
    """Write a UTF-8 test fixture and return its path."""
    path.write_text(content, encoding="utf-8")
    return path


def write_gzip_text(path: Path, content: str) -> Path:
    """Write a gzip-compressed UTF-8 test fixture and return its path."""
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write(content)
    return path


@pytest.mark.parametrize("as_string", (False, True))
def test_validate_file_accepts_existing_file(
    tmp_path: Path,
    as_string: bool,
) -> None:
    """Both Path and string inputs should return a Path."""
    input_path = write_text(tmp_path / "input.txt", "content")
    value = str(input_path) if as_string else input_path

    assert validate_file(value) == input_path


def test_validate_file_rejects_missing_path(tmp_path: Path) -> None:
    """Missing inputs should identify the path and reason."""
    missing = tmp_path / "missing.txt"

    with pytest.raises(InputValidationError, match=r"missing\.txt.*missing"):
        validate_file(missing)


def test_validate_file_rejects_directory(tmp_path: Path) -> None:
    """A directory must not pass file validation."""
    with pytest.raises(InputValidationError, match=r"not a file"):
        validate_file(tmp_path)


def test_validate_file_rejects_empty_by_default(tmp_path: Path) -> None:
    """Content-bearing input validation should default to non-empty."""
    empty = write_text(tmp_path / "empty.txt", "")

    with pytest.raises(InputValidationError, match=r"empty\.txt.*empty"):
        validate_file(empty)


def test_validate_file_can_allow_empty_input(tmp_path: Path) -> None:
    """Callers such as checksum computation may explicitly allow empty files."""
    empty = write_text(tmp_path / "empty.txt", "")

    assert validate_file(empty, require_nonempty=False) == empty


def test_validate_file_supports_unicode_path(tmp_path: Path) -> None:
    """Unicode paths used on Windows and Linux should remain intact."""
    unicode_path = write_text(tmp_path / "测试文件.txt", "样本")

    assert validate_file(unicode_path) == unicode_path


def test_validate_file_accepts_working_symlink(tmp_path: Path) -> None:
    """A symbolic link should be validated through its existing target."""
    target = write_text(tmp_path / "target.txt", "data")
    link = tmp_path / "link.txt"
    try:
        link.symlink_to(target)
    except OSError as error:
        pytest.skip(f"Symlink creation is unavailable: {error}")

    assert validate_file(link) == link


def test_validate_file_rejects_broken_symlink(tmp_path: Path) -> None:
    """Broken links must fail without losing the user-visible link path."""
    link = tmp_path / "broken.txt"
    try:
        link.symlink_to(tmp_path / "absent.txt")
    except OSError as error:
        pytest.skip(f"Symlink creation is unavailable: {error}")

    with pytest.raises(InputValidationError, match=r"broken\.txt.*broken"):
        validate_file(link)


def test_validate_directory_accepts_existing_directory(tmp_path: Path) -> None:
    """Existing readable directories should pass."""
    assert validate_directory(tmp_path) == tmp_path


def test_validate_directory_rejects_missing_directory(tmp_path: Path) -> None:
    """Input directories are not created implicitly."""
    missing = tmp_path / "missing-directory"

    with pytest.raises(InputValidationError, match=r"missing-directory.*missing"):
        validate_directory(missing)
    assert not missing.exists()


def test_validate_directory_can_describe_future_output(tmp_path: Path) -> None:
    """must_exist=False returns a path but does not create it."""
    future = tmp_path / "future-output"

    assert validate_directory(future, must_exist=False) == future
    assert not future.exists()


def test_validate_directory_rejects_file(tmp_path: Path) -> None:
    """Files should not pass directory validation."""
    file_path = write_text(tmp_path / "input.txt", "data")

    with pytest.raises(InputValidationError, match=r"not a directory"):
        validate_directory(file_path)


def test_validate_output_file_accepts_nonempty_output(tmp_path: Path) -> None:
    """Expected non-empty outputs should pass."""
    output = write_text(tmp_path / "result.txt", "result")

    assert validate_output_file(output) == output


@pytest.mark.parametrize("content", (None, ""))
def test_validate_output_file_uses_output_error(
    tmp_path: Path,
    content: str | None,
) -> None:
    """Missing and empty outputs should use OutputValidationError."""
    output = tmp_path / "result.txt"
    if content is not None:
        write_text(output, content)

    with pytest.raises(OutputValidationError, match=r"result\.txt"):
        validate_output_file(output)


def test_validate_fasta_accepts_tiny_reference(tmp_path: Path) -> None:
    """A header followed by sequence is sufficient for lightweight checking."""
    fasta = write_text(tmp_path / "reference.fa", VALID_FASTA)

    assert validate_fasta(fasta) == fasta


def test_validate_fasta_accepts_gzip_reference(tmp_path: Path) -> None:
    """Supported gzip FASTA input should be streamed through gzip."""
    fasta = write_gzip_text(tmp_path / "reference.fasta.gz", VALID_FASTA)

    assert validate_fasta(fasta) == fasta


def test_validate_fasta_rejects_unsupported_suffix(tmp_path: Path) -> None:
    """A FASTA-like body should not bypass the explicit suffix contract."""
    fasta = write_text(tmp_path / "reference.txt", VALID_FASTA)

    with pytest.raises(InputValidationError, match=r"suffix.*reference\.txt"):
        validate_fasta(fasta)


def test_validate_fasta_rejects_invalid_first_record(tmp_path: Path) -> None:
    """The first non-empty line must be a FASTA header."""
    fasta = write_text(tmp_path / "reference.fa", "ACGT\n")

    with pytest.raises(InputValidationError, match=r"header"):
        validate_fasta(fasta)


def test_validate_fasta_rejects_header_without_sequence(tmp_path: Path) -> None:
    """A header-only first record is not useful reference content."""
    fasta = write_text(tmp_path / "reference.fa", ">chr1\n>chr2\nAAAA\n")

    with pytest.raises(InputValidationError, match=r"no sequence"):
        validate_fasta(fasta)


def test_validate_fasta_rejects_empty_file(tmp_path: Path) -> None:
    """Empty references should fail common file validation first."""
    fasta = write_text(tmp_path / "reference.fa", "")

    with pytest.raises(InputValidationError, match=r"empty"):
        validate_fasta(fasta)


def test_validate_fasta_with_required_fai(tmp_path: Path) -> None:
    """A conventional present and parseable FAI should pass."""
    fasta = write_text(tmp_path / "reference.fa", VALID_FASTA)
    index = write_text(Path(f"{fasta}.fai"), VALID_FAI)

    assert validate_fasta(fasta, require_fai=True) == fasta
    assert validate_fasta_index(fasta) == index


def test_validate_fasta_reports_missing_required_fai(tmp_path: Path) -> None:
    """FAI absence is a reference error and never triggers index creation."""
    fasta = write_text(tmp_path / "reference.fa", VALID_FASTA)

    with pytest.raises(ReferenceError, match=r"index missing.*reference\.fa"):
        validate_fasta(fasta, require_fai=True)
    assert not Path(f"{fasta}.fai").exists()


def test_read_fai_contigs_returns_names_and_lengths(tmp_path: Path) -> None:
    """FAI parsing should preserve contig names and integer lengths."""
    index = write_text(tmp_path / "reference.fa.fai", VALID_FAI)

    assert read_fai_contigs(index) == {"chr1": 8, "chr2": 4}


def test_read_fai_contigs_rejects_duplicate_contig(tmp_path: Path) -> None:
    """Reference contig names must be unique."""
    index = write_text(
        tmp_path / "reference.fa.fai",
        "chr1\t8\t6\t8\t9\nchr1\t4\t21\t4\t5\n",
    )

    with pytest.raises(ReferenceError, match=r"duplicate contig.*chr1"):
        read_fai_contigs(index)


@pytest.mark.parametrize(
    "content,keyword",
    (("chr1\t8\n", "5"), ("chr1\teight\t6\t8\t9\n", "integer")),
)
def test_read_fai_contigs_rejects_malformed_rows(
    tmp_path: Path,
    content: str,
    keyword: str,
) -> None:
    """Basic column count and numeric fields should be checked."""
    index = write_text(tmp_path / "reference.fa.fai", content)

    with pytest.raises(ReferenceError, match=keyword):
        read_fai_contigs(index)


def test_validate_fastq_accepts_first_record(tmp_path: Path) -> None:
    """Only the first complete FASTQ record is inspected."""
    fastq = write_text(tmp_path / "reads.fastq", VALID_FASTQ)

    assert validate_fastq(fastq) == fastq


def test_validate_fastq_accepts_gzip(tmp_path: Path) -> None:
    """Gzip FASTQ text should be decoded without full decompression."""
    fastq = write_gzip_text(tmp_path / "reads.fq.gz", VALID_FASTQ)

    assert validate_fastq(fastq) == fastq


@pytest.mark.parametrize(
    "content,keyword",
    (
        ("read1\nACGT\n+\nIIII\n", "@"),
        ("@read1\nACGT\nseparator\nIIII\n", r"\+"),
        ("@read1\nACGT\n+\n", "incomplete"),
        ("@read1\nACGT\n+\nIII\n", "lengths differ"),
    ),
)
def test_validate_fastq_rejects_invalid_first_record(
    tmp_path: Path,
    content: str,
    keyword: str,
) -> None:
    """Each required four-line FASTQ invariant should fail clearly."""
    fastq = write_text(tmp_path / "reads.fastq", content)

    with pytest.raises(InputValidationError, match=keyword):
        validate_fastq(fastq)


@pytest.mark.parametrize("suffix", (".bam", ".cram"))
def test_alignment_path_validation_without_binary_claims(
    tmp_path: Path,
    suffix: str,
) -> None:
    """A non-empty placeholder only proves path-level alignment validation."""
    alignment = write_text(tmp_path / f"sample{suffix}", "placeholder")

    assert validate_alignment_file(alignment) == alignment


@pytest.mark.parametrize(
    "alignment_name,index_name",
    (
        ("sample.bam", "sample.bam.bai"),
        ("sample.bam", "sample.bai"),
        ("sample.cram", "sample.cram.crai"),
        ("sample.cram", "sample.crai"),
    ),
)
def test_alignment_path_and_index_validation(
    tmp_path: Path,
    alignment_name: str,
    index_name: str,
) -> None:
    """Both conventional index naming forms should be recognized."""
    alignment = write_text(tmp_path / alignment_name, "placeholder")
    write_text(tmp_path / index_name, "index-placeholder")

    assert validate_alignment_file(alignment, require_index=True) == alignment


def test_alignment_required_index_missing(tmp_path: Path) -> None:
    """Index checking should never invoke an external indexer."""
    alignment = write_text(tmp_path / "sample.bam", "placeholder")

    with pytest.raises(InputValidationError, match=r"index missing.*sample\.bam"):
        validate_alignment_file(alignment, require_index=True)


def test_alignment_rejects_unsupported_suffix(tmp_path: Path) -> None:
    """Only BAM and CRAM paths are within this API."""
    alignment = write_text(tmp_path / "sample.sam", "@HD")

    with pytest.raises(InputValidationError, match=r"suffix.*sample\.sam"):
        validate_alignment_file(alignment)


def test_alignment_rejects_empty_file(tmp_path: Path) -> None:
    """An empty alignment path cannot represent usable input."""
    alignment = write_text(tmp_path / "sample.bam", "")

    with pytest.raises(InputValidationError, match=r"empty"):
        validate_alignment_file(alignment)


@pytest.mark.parametrize("compressed", (False, True))
def test_validate_vcf_streams_required_headers(
    tmp_path: Path,
    compressed: bool,
) -> None:
    """Plain VCF and gzip VCF should share lightweight header checks."""
    if compressed:
        vcf = write_gzip_text(tmp_path / "sample.vcf.gz", VALID_VCF)
    else:
        vcf = write_text(tmp_path / "sample.vcf", VALID_VCF)

    assert validate_vcf(vcf) == vcf


@pytest.mark.parametrize(
    "content,keyword",
    (
        (
            "##contig=<ID=chr1,length=1000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "fileformat",
        ),
        ("##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n", "CHROM"),
    ),
)
def test_validate_vcf_rejects_missing_required_header(
    tmp_path: Path,
    content: str,
    keyword: str,
) -> None:
    """Both required VCF header markers must be present."""
    vcf = write_text(tmp_path / "sample.vcf", content)

    with pytest.raises(InputValidationError, match=keyword):
        validate_vcf(vcf)


def test_validate_vcf_rejects_empty_file(tmp_path: Path) -> None:
    """Empty VCFs should fail before header inspection."""
    vcf = write_text(tmp_path / "sample.vcf", "")

    with pytest.raises(InputValidationError, match=r"empty"):
        validate_vcf(vcf)


@pytest.mark.parametrize("index_suffix", (".tbi", ".csi"))
def test_validate_vcf_accepts_index_presence(
    tmp_path: Path,
    index_suffix: str,
) -> None:
    """TBI or CSI path presence is accepted without claiming index integrity."""
    vcf = write_gzip_text(tmp_path / "sample.vcf.gz", VALID_VCF)
    write_text(Path(f"{vcf}{index_suffix}"), "index-placeholder")

    assert validate_vcf(vcf, require_index=True) == vcf


def test_validate_vcf_reports_missing_index(tmp_path: Path) -> None:
    """A requested compressed VCF index must already exist."""
    vcf = write_gzip_text(tmp_path / "sample.vcf.gz", VALID_VCF)

    with pytest.raises(InputValidationError, match=r"index missing.*sample\.vcf\.gz"):
        validate_vcf(vcf, require_index=True)


def test_validate_vcf_index_requires_compressed_vcf(tmp_path: Path) -> None:
    """This lightweight API does not pretend a plain VCF is tabix-indexed."""
    vcf = write_text(tmp_path / "sample.vcf", VALID_VCF)

    with pytest.raises(InputValidationError, match=r"requires a \.vcf\.gz"):
        validate_vcf(vcf, require_index=True)


def test_unicode_gzip_vcf_path_and_logging(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Unicode paths should survive gzip access and diagnostic logging."""
    vcf = write_gzip_text(tmp_path / "测试文件.vcf.gz", VALID_VCF)

    with caplog.at_level(logging.DEBUG, logger="hifivar"):
        assert validate_vcf(vcf) == vcf

    assert "测试文件.vcf.gz" in caplog.text


def test_validate_bed_accepts_comments_and_directives(tmp_path: Path) -> None:
    """Comments, track, and browser declarations are not data records."""
    bed = write_text(
        tmp_path / "regions.bed",
        "# comment\ntrack name=test\nbrowser position chr1\n"
        "chr1\t0\t100\nchr1\t200\t300\n",
    )

    assert validate_bed(bed) == bed


def test_validate_bed_accepts_gzip(tmp_path: Path) -> None:
    """BED gzip input should be validated as a stream."""
    bed = write_gzip_text(tmp_path / "regions.bed.gz", "chr1\t0\t100\n")

    assert validate_bed(bed) == bed


@pytest.mark.parametrize(
    "content,keyword",
    (
        ("chr1\t0\n", "fewer"),
        ("chr1\tstart\t100\n", "integer"),
        ("chr1\t-1\t100\n", "negative"),
        ("chr1\t100\t100\n", "greater than start"),
    ),
)
def test_validate_bed_rejects_invalid_record(
    tmp_path: Path,
    content: str,
    keyword: str,
) -> None:
    """BED coordinate failures should identify their semantic cause."""
    bed = write_text(tmp_path / "regions.bed", content)

    with pytest.raises(InputValidationError, match=keyword):
        validate_bed(bed)


def test_validate_bed_rejects_comment_only_file(tmp_path: Path) -> None:
    """A non-empty file with no region records is not a usable BED input."""
    bed = write_text(tmp_path / "regions.bed", "# no regions\n")

    with pytest.raises(InputValidationError, match=r"no data records"):
        validate_bed(bed)


def test_contig_compatibility_accepts_reference_subset() -> None:
    """Query files may legitimately contain only part of the reference."""
    validate_contig_compatibility(
        ["chr1", "chr2", "chr3"],
        ["chr1", "chr3"],
    )


def test_contig_compatibility_rejects_missing_query_contig() -> None:
    """Query-only names should produce the stable mismatch marker."""
    with pytest.raises(ReferenceError, match=r"REFERENCE_CONTIG_MISMATCH.*chr99"):
        validate_contig_compatibility(["chr1", "chr2"], ["chr99"])


def test_contig_compatibility_does_not_rename_chr_prefix() -> None:
    """chr1 and 1 remain distinct and require explicit upstream normalization."""
    with pytest.raises(ReferenceError, match=r"REFERENCE_CONTIG_MISMATCH.*\b1\b"):
        validate_contig_compatibility(["chr1"], ["1"])


def test_contig_compatibility_rejects_duplicate_reference() -> None:
    """Duplicate reference contigs make compatibility ambiguous."""
    with pytest.raises(ReferenceError, match=r"duplicate"):
        validate_contig_compatibility(["chr1", "chr1"], ["chr1"])


def test_contig_compatibility_rejects_empty_query() -> None:
    """The generic helper must not silently accept an unknown contig set."""
    with pytest.raises(ReferenceError, match=r"Query contig collection is empty"):
        validate_contig_compatibility(["chr1"], [])


def test_compute_sha256_matches_hashlib(tmp_path: Path) -> None:
    """The public digest should match hashlib for a tiny binary file."""
    input_path = tmp_path / "input.bin"
    input_path.write_bytes(b"abc")

    assert compute_sha256(input_path) == hashlib.sha256(b"abc").hexdigest()


def test_compute_sha256_accepts_empty_file(tmp_path: Path) -> None:
    """An empty file has a valid, reproducible SHA256 digest."""
    input_path = tmp_path / "empty.bin"
    input_path.touch()

    assert compute_sha256(input_path) == hashlib.sha256(b"").hexdigest()


def test_compute_sha256_handles_unicode_content(tmp_path: Path) -> None:
    """Checksums operate on bytes and preserve UTF-8 file content semantics."""
    content = "HiFiVar 样本".encode()
    input_path = tmp_path / "校验.txt"
    input_path.write_bytes(content)

    assert compute_sha256(input_path) == hashlib.sha256(content).hexdigest()


def test_compute_sha256_reads_more_than_one_chunk(tmp_path: Path) -> None:
    """A multi-chunk synthetic file should produce the complete digest."""
    content = b"A" * (CHECKSUM_CHUNK_SIZE + 257)
    input_path = tmp_path / "multi-chunk.bin"
    input_path.write_bytes(content)

    assert compute_sha256(input_path) == hashlib.sha256(content).hexdigest()
