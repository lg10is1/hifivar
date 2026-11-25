"""Fail-closed benchmark reference gates, without external tools or repairs."""

from __future__ import annotations

import gzip
import hashlib
import json
import struct
import zlib
from pathlib import Path

import pytest

from hifivar.benchmark_inputs import audit_benchmark_regions, audit_sv_benchmark_input
from hifivar.exceptions import InputValidationError, ReferenceError
from hifivar.reference import Contig, ReferenceGenome


def _bgzf(path: Path, payload: bytes) -> None:
    compressor = zlib.compressobj(level=6, wbits=-15)
    compressed = compressor.compress(payload) + compressor.flush()
    size = 18 + len(compressed) + 8
    header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
    path.write_bytes(header + struct.pack("<H", size - 1) + compressed + struct.pack("<II", zlib.crc32(payload), len(payload)))


@pytest.fixture
def reference(tmp_path: Path) -> ReferenceGenome:
    fasta = tmp_path / "参考.fa"
    fasta.write_text(">chr1\n" + "A" * 100 + "\n>chr2\n" + "C" * 80 + "\n", encoding="utf-8")
    fai = Path(f"{fasta}.fai")
    fai.write_text("chr1\t100\t6\t100\t101\nchr2\t80\t113\t80\t81\n", encoding="utf-8")
    return ReferenceGenome.from_fasta(fasta, build="GRCh38")


def _vcf(
    tmp_path: Path,
    *,
    sample: str = "HG002",
    contigs: str = "##contig=<ID=chr1,length=100>\n",
    metadata: str = "",
    svtype: bool = True,
    body: str = "chr1\t1\tid\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL\tGT\t0/1\n",
) -> Path:
    path = tmp_path / "真值.sv.vcf.gz"
    info = '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n' if svtype else ""
    text = "##fileformat=VCFv4.4\n" + contigs + metadata + info
    text += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n" + body
    _bgzf(path, text.encode("utf-8"))
    _bgzf(Path(f"{path}.tbi"), b"TBI\x01")
    return path


def _bed(tmp_path: Path, text: str, *, compressed: bool = False) -> Path:
    path = tmp_path / ("范围.bed.gz" if compressed else "范围.bed")
    if compressed:
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
    else:
        path.write_text(text, encoding="utf-8")
    return path


def _hashes(paths: list[Path]) -> dict[str, str]:
    return {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def test_mixed_stvar_header_audit_does_not_filter_small_or_star_records(tmp_path, reference):
    body = ("chr1\t1\tsnp\tA\tG\t30\tPASS\t.\tGT\t1|1\n"
            "chr1\t2\tshort\tA\tAT\t30\tPASS\t.\tGT\t0|1\n"
            "chr1\t10\tsv\tA\t" + "A" * 60 + "\t30\tPASS\tSVTYPE=INS\tGT\t0|1\n"
            "chr1\t80\tstar\tA\t*\t30\tPASS\t.\tGT\t1|0\n")
    path = _vcf(tmp_path, body=body)
    before = _hashes([path, Path(f"{path}.tbi")])
    result = audit_sv_benchmark_input(path, "HG002", reference)
    assert result["audit_status"] == "PASS"
    assert result["record_validation"] == "NOT_RUN_HEADER_ONLY"
    assert result["sequence_check"] == "NOT_RUN_REQUIRES_BCFTOOLS"
    assert result["reference_build_check"] == "UNVERIFIED_MISSING_OR_UNRECOGNIZED_REFERENCE_METADATA"
    assert _hashes([path, Path(f"{path}.tbi")]) == before
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        assert "".join(line for line in handle if not line.startswith("#")) == body


def test_vcf_exact_name_subset_and_header_scope(tmp_path, reference):
    path = _vcf(tmp_path)
    result = audit_sv_benchmark_input(path, "HG002", reference)
    assert result["audit_status"] == "PASS"
    assert result["declared_contig_count"] == 1
    assert result["reference_contig_count"] == 2
    assert result["contig_length_check"] == "PASS"
    assert result["reference_build_check"] == "UNVERIFIED_MISSING_OR_UNRECOGNIZED_REFERENCE_METADATA"
    assert result["sequence_check"] == "NOT_RUN_REQUIRES_BCFTOOLS"
    assert result["record_validation"] == "NOT_RUN_HEADER_ONLY"
    assert len(result["header_sha256"]) == 64
    json.dumps(result)


def test_truth_24_contigs_may_be_subset_of_reference_195(tmp_path):
    contigs = tuple(Contig(f"chr{i}", 100) for i in range(1, 23)) + (Contig("chrX", 100), Contig("chrY", 100))
    ref = ReferenceGenome(tmp_path / "ref.fa", tmp_path / "ref.fa.fai", "GRCh38", contigs + tuple(Contig(f"decoy_{i}", 50) for i in range(171)))
    path = _vcf(tmp_path, contigs="".join(f"##contig=<ID={item.name},length={item.length}>\n" for item in contigs))
    result = audit_sv_benchmark_input(path, "HG002", ref)
    assert (result["declared_contig_count"], result["reference_contig_count"]) == (24, 195)


@pytest.mark.parametrize("metadata", ["##reference=GRCh37\n", "##reference=hg19\n", "##reference=file:///data/hg19/ref.fa\n", "##reference=GRCh37.p13.fa\n"])
def test_explicit_grch37_vs_grch38_fails_even_if_names_lengths_match(tmp_path, reference, metadata):
    with pytest.raises(ReferenceError, match="REFERENCE_BUILD_MISMATCH"):
        audit_sv_benchmark_input(_vcf(tmp_path, metadata=metadata), "HG002", reference)


@pytest.mark.parametrize("metadata", ["##reference=GRCh38\n", "##reference=hg38\n", "##reference=GRCh38_full_analysis_set.fa\n"])
def test_matching_explicit_builds(tmp_path, reference, metadata):
    result = audit_sv_benchmark_input(_vcf(tmp_path, metadata=metadata), "HG002", reference)
    assert result["reference_build_check"] == "PASS"
    assert result["explicit_vcf_builds"] == ["GRCh38"]


def test_reverse_build_mismatch_and_conflicting_headers_fail(tmp_path, reference):
    ref37 = ReferenceGenome(reference.fasta, reference.fai, "GRCh37", reference.contigs)
    with pytest.raises(ReferenceError, match="REFERENCE_BUILD_MISMATCH"):
        audit_sv_benchmark_input(_vcf(tmp_path, metadata="##reference=GRCh38\n"), "HG002", ref37)
    with pytest.raises(ReferenceError, match="REFERENCE_BUILD_MISMATCH"):
        audit_sv_benchmark_input(_vcf(tmp_path, metadata="##reference=GRCh37\n##reference=GRCh38\n"), "HG002", reference)


def test_unknown_declared_build_remains_unverified(tmp_path, reference):
    ref = ReferenceGenome(reference.fasta, reference.fai, None, reference.contigs)
    result = audit_sv_benchmark_input(_vcf(tmp_path, metadata="##reference=GRCh38\n"), "HG002", ref)
    assert result["reference_build_check"] == "UNVERIFIED_DECLARED_REFERENCE_BUILD"


def test_contig_length_mismatch_blocks_renamed_grch37_truth(tmp_path, reference):
    path = _vcf(tmp_path, contigs="##contig=<ID=chr1,length=249250621>\n")
    with pytest.raises(ReferenceError, match="CONTIG_LENGTH_MISMATCH"):
        audit_sv_benchmark_input(path, "HG002", reference)


def test_missing_length_strict_failure_and_explicit_diagnostic_unverified(tmp_path, reference):
    path = _vcf(tmp_path, contigs="##contig=<ID=chr1>\n")
    with pytest.raises(InputValidationError, match="CONTIG_LENGTH_UNVERIFIED"):
        audit_sv_benchmark_input(path, "HG002", reference)
    result = audit_sv_benchmark_input(path, "HG002", reference, require_contig_lengths=False)
    assert result["audit_status"] == "UNVERIFIED"
    assert result["contig_length_check"] == "UNVERIFIED_MISSING_LENGTHS"
    assert result["missing_contig_lengths"] == ["chr1"]


@pytest.mark.parametrize("length", ["0", "-1", "1.0", "NaN", "True", "１２"])
def test_invalid_header_lengths_fail_even_in_diagnostic_mode(tmp_path, reference, length):
    with pytest.raises(InputValidationError, match="Invalid VCF contig length"):
        audit_sv_benchmark_input(_vcf(tmp_path, contigs=f"##contig=<ID=chr1,length={length}>\n"), "HG002", reference, require_contig_lengths=False)


@pytest.mark.parametrize("contigs", ["", "##contig=<ID=chr1,length=100>\n##contig=<ID=chr1,length=100>\n", "##contig=<ID=chr1,length=100,length=100>\n", "##contig=<ID=chr1,length=100\n", '##contig=<ID=chr1,length=100,Description="unfinished>\n'])
def test_empty_duplicate_or_malformed_contig_headers_fail(tmp_path, reference, contigs):
    with pytest.raises(InputValidationError):
        audit_sv_benchmark_input(_vcf(tmp_path, contigs=contigs), "HG002", reference)


def test_quoted_contig_description_can_contain_commas(tmp_path, reference):
    result = audit_sv_benchmark_input(_vcf(tmp_path, contigs='##contig=<ID=chr1,length=100,Description="a, b">\n'), "HG002", reference)
    assert result["contig_length_check"] == "PASS"


@pytest.mark.parametrize("sample_id", ["", "UnknownSampleName", "unknownsamplename", "unknown", "HG 002", " HG002", "HG002\t"])
def test_missing_or_placeholder_expected_sample_fails(tmp_path, reference, sample_id):
    with pytest.raises(InputValidationError, match="SAMPLE_IDENTITY_UNVERIFIED"):
        audit_sv_benchmark_input(_vcf(tmp_path, sample=sample_id), sample_id, reference)


@pytest.mark.parametrize("sample", ["UnknownSampleName", "wrong", "", "HG002\tEXTRA"])
def test_observed_sample_mismatch_or_extra_sample_fails(tmp_path, reference, sample):
    with pytest.raises(InputValidationError):
        audit_sv_benchmark_input(_vcf(tmp_path, sample=sample), "HG002", reference)


@pytest.mark.parametrize("failure", ["contig", "svtype", "plain_gzip", "bad_tbi", "missing_tbi"])
def test_existing_sv_binary_and_header_validation_reused(tmp_path, reference, failure):
    path = _vcf(tmp_path, contigs="##contig=<ID=1,length=100>\n" if failure == "contig" else "##contig=<ID=chr1,length=100>\n", svtype=failure != "svtype")
    index = Path(f"{path}.tbi")
    if failure == "plain_gzip":
        text = gzip.decompress(path.read_bytes())
        path.write_bytes(gzip.compress(text))
    elif failure == "bad_tbi":
        _bgzf(index, b"WRONG")
    elif failure == "missing_tbi":
        index.unlink()
    with pytest.raises(InputValidationError):
        audit_sv_benchmark_input(path, "HG002", reference)


def test_header_only_does_not_claim_record_or_sequence_validation(tmp_path, reference):
    path = _vcf(tmp_path, body="not a valid variant record\n")
    result = audit_sv_benchmark_input(path, "HG002", reference)
    assert result["record_validation"] == "NOT_RUN_HEADER_ONLY"
    assert result["sequence_check"] == "NOT_RUN_REQUIRES_BCFTOOLS"


def test_vcf_reference_index_unchanged_on_success_and_failure(tmp_path, reference):
    path = _vcf(tmp_path)
    paths = [path, Path(f"{path}.tbi"), reference.fasta, reference.fai]
    before = _hashes(paths)
    audit_sv_benchmark_input(path, "HG002", reference)
    with pytest.raises(InputValidationError):
        audit_sv_benchmark_input(path, "WRONG", reference)
    assert _hashes(paths) == before


def test_bed_union_counts_overlap_duplicates_and_adjacency(tmp_path, reference):
    path = _bed(tmp_path, '#comment\ntrack name="test"\nbrowser position chr1:1-50\nchr1\t0\t10\nchr1\t5\t15\nchr1\t5\t15\nchr1\t15\t20\nchr1\t30\t40\nchr2\t0\t80\n')
    before = path.read_bytes()
    result = audit_benchmark_regions(path, reference)
    assert result["interval_count"] == 6
    assert result["input_bases"] == 125
    assert result["union_bases"] == 110
    assert result["overlap_bases_counted_once"] == 15
    assert result["overlapping_interval_count"] == 2
    assert result["adjacent_interval_count"] == 1
    assert result["per_contig"]["chr1"]["union_bases"] == 30
    assert result["input_modified"] is False
    assert path.read_bytes() == before
    json.dumps(result)


def test_bed_gzip_unicode_path_and_any_group_order(tmp_path, reference):
    result = audit_benchmark_regions(_bed(tmp_path, "chr2\t0\t80\nchr1\t90\t100\textra\n", compressed=True), reference)
    assert result["union_bases"] == 90
    assert result["contig_count"] == 2


@pytest.mark.parametrize("text", ["chr1\t-1\t10\n", "chr1\t0\t101\n", "chr1\t10\t10\n", "chr1\t20\t10\n", "chr1\t1.5\t10\n", "chr1\t+1\t10\n", "chr1\t０\t10\n", "chr1\t0\n"])
def test_bed_bad_coordinates_fail(tmp_path, reference, text):
    with pytest.raises(InputValidationError):
        audit_benchmark_regions(_bed(tmp_path, text), reference)


@pytest.mark.parametrize("text", ["chr1\t10\t20\nchr1\t0\t5\n", "chr1\t0\t5\nchr2\t0\t5\nchr1\t10\t15\n"])
def test_bed_unsorted_fails_without_silent_sorting(tmp_path, reference, text):
    path = _bed(tmp_path, text)
    before = path.read_bytes()
    with pytest.raises(InputValidationError, match="BED_NOT_SORTED"):
        audit_benchmark_regions(path, reference)
    assert path.read_bytes() == before


def test_bed_contig_mismatch_never_renames(tmp_path, reference):
    path = _bed(tmp_path, "1\t0\t10\n")
    with pytest.raises(ReferenceError, match="REFERENCE_CONTIG_MISMATCH"):
        audit_benchmark_regions(path, reference)
    assert path.read_text(encoding="utf-8") == "1\t0\t10\n"


@pytest.mark.parametrize("text", ["", "#comment\n\n", "track name=test\nbrowser position chr1:1-10\n"])
def test_empty_bed_fails(tmp_path, reference, text):
    with pytest.raises(InputValidationError):
        audit_benchmark_regions(_bed(tmp_path, text), reference)


def test_bed_nested_intervals_union_is_not_sum_or_previous_end_only(tmp_path, reference):
    result = audit_benchmark_regions(_bed(tmp_path, "chr1\t0\t100\nchr1\t10\t20\nchr1\t30\t40\n"), reference)
    assert result["input_bases"] == 120
    assert result["union_bases"] == 100
    assert result["overlapping_interval_count"] == 2


def test_invalid_api_argument_types_fail_explicitly(tmp_path, reference):
    path = _vcf(tmp_path)
    with pytest.raises(InputValidationError, match="boolean"):
        audit_sv_benchmark_input(path, "HG002", reference, require_contig_lengths="yes")
    with pytest.raises(InputValidationError, match="ReferenceGenome"):
        audit_sv_benchmark_input(path, "HG002", None)
    with pytest.raises(InputValidationError, match="nonempty path"):
        audit_benchmark_regions("", reference)
