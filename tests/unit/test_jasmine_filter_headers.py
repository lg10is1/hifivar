"""Header-only Jasmine FILTER recovery; never repair native body records."""

from collections import Counter
import gzip
from pathlib import Path
import struct
import zlib

import pytest

from hifivar.assembly_sv import SVEvidenceSource
from hifivar.command import CommandResult
from hifivar.exceptions import OutputValidationError
from hifivar.harmonization import (
    EvidenceRunStatus,
    SVEvidenceSourceArtifact,
    SVHarmonizationRequest,
)
import hifivar.jasmine as jasmine_module
from hifivar.jasmine import JasmineWrapper
from hifivar.reference import ReferenceGenome


Q5 = '##FILTER=<ID=q5,Description="Quality below 5">'
NEAR_GAP = ('##FILTER=<ID=NearReferenceGap,Description="Variant is near '
            '(< 1000 bp) from a gap (run of >= 50 Ns) in the reference assembly">')
INFO_FORMAT = (
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">',
    '##INFO=<ID=SUPP,Number=1,Type=Integer,Description="Sources">',
    '##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Source vector">',
    '##INFO=<ID=IDLIST,Number=.,Type=String,Description="Source IDs">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
)


def payload(records=(), *, filters=(), ending=b"\n", declarations=INFO_FORMAT):
    headers = ("##fileformat=VCFv4.4", "##contig=<ID=chr1,length=1000>",
               "##contig=<ID=chr3,length=1000>", *filters, *declarations,
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002")
    return ending.join(line.encode() for line in headers) + ending + b"".join(records)


def row(identity="v1", position=10, filters="q5", *, ending=b"\n", kind="INS", alt="AT", end=None):
    finish = position if end is None else end
    text = (f"chr1\t{position}\t{identity}\tA\t{alt}\t.\t{filters}\t"
            f"SVTYPE={kind};END={finish};SVLEN=0;SUPP=1;SUPP_VEC=1;IDLIST={identity}\tGT\t0/0")
    return text.encode() + ending


def body(data):
    return [line for line in data.splitlines(keepends=True) if not line.startswith(b"#")]


def finalize(tmp_path, records, *, source_filters=(Q5, NEAR_GAP), raw_filters=(), ending=b"\n", sources=None):
    raw = tmp_path / "native.raw.vcf"
    raw.write_bytes(payload(records, filters=raw_filters, ending=ending))
    if sources is None:
        source = tmp_path / "source.vcf"
        source.write_bytes(payload(records, filters=source_filters, ending=ending))
        sources = (source,)
    destination = tmp_path / "final.sorted.vcf"
    jasmine_module._sort_jasmine_vcf(raw, destination, tuple(sources), overwrite=False)
    return raw, destination


@pytest.mark.parametrize("chunk_size", [1, 100_000])
@pytest.mark.parametrize("ending", [b"\n", b"\r\n"])
def test_recovers_only_used_filters_preserving_all_body_bytes_and_native_bnd_end(
    tmp_path, monkeypatch, chunk_size, ending,
):
    monkeypatch.setattr(jasmine_module, "_SORT_CHUNK_RECORDS", chunk_size)
    records = [row("bnd", 180, "PASS", kind="BND", alt="]chr3:10]A", end=10, ending=ending),
               row("boundary", 0, "q5", ending=ending),
               row("near", 12, "NearReferenceGap", ending=ending)]
    raw, destination = finalize(tmp_path, records, ending=ending)
    result = destination.read_bytes()
    assert raw.read_bytes() == payload(records, ending=ending)
    assert Counter(body(result)) == Counter(records)
    assert body(result) == [records[1], records[2], records[0]]
    assert result.count((Q5 + "\n").encode()) == 1
    assert result.count((NEAR_GAP + "\n").encode()) == 1
    assert b"\t180\tbnd\tA\t]chr3:10]A\t" in result
    assert b"SVTYPE=BND;END=10;SVLEN=0" in result
    assert gzip.decompress(gzip.compress(result)) == result
    assert not list(tmp_path.glob("*.sort.*.tmp"))


def test_multiple_filter_tokens_remain_verbatim(tmp_path):
    original = row(filters="NearReferenceGap;q5")
    _, destination = finalize(tmp_path, [original])
    assert body(destination.read_bytes()) == [original]
    assert Q5 in destination.read_text(encoding="utf-8")
    assert NEAR_GAP in destination.read_text(encoding="utf-8")


def test_unused_filter_definitions_are_not_added(tmp_path):
    _, destination = finalize(tmp_path, [row(filters="PASS")])
    text = destination.read_text(encoding="utf-8")
    assert Q5 not in text and NEAR_GAP not in text


@pytest.mark.parametrize("reserved", ["PASS", "."])
def test_reserved_pass_and_missing_do_not_require_or_restore_declarations(tmp_path, reserved):
    declarations = ('##FILTER=<ID=PASS,Description="Passed">',
                    '##FILTER=<ID=.,Description="Missing">')
    _, destination = finalize(tmp_path, [row(filters=reserved)], source_filters=declarations)
    assert not any(line.startswith("##FILTER=") for line in destination.read_text(encoding="utf-8").splitlines())


@pytest.mark.parametrize("filters", ["", "PASS;q5", ".;q5", "q5;", ";q5"])
def test_malformed_reserved_or_empty_filter_combinations_fail_without_rewriting(tmp_path, filters):
    raw = tmp_path / "raw.vcf"
    original = payload([row(filters=filters)])
    raw.write_bytes(original)
    source = tmp_path / "source.vcf"
    source.write_bytes(payload(filters=(Q5,)))
    with pytest.raises(OutputValidationError, match="Malformed Jasmine FILTER"):
        jasmine_module._sort_jasmine_vcf(raw, tmp_path / "sorted.vcf", (source,), overwrite=False)
    assert raw.read_bytes() == original
    assert not (tmp_path / "sorted.vcf").exists()


def test_missing_used_unknown_filter_fails_and_cleans_chunk_files(tmp_path, monkeypatch):
    monkeypatch.setattr(jasmine_module, "_SORT_CHUNK_RECORDS", 1)
    raw = tmp_path / "raw.vcf"
    original = payload([row(filters="unknown")])
    raw.write_bytes(original)
    source = tmp_path / "source.vcf"
    source.write_bytes(payload(filters=(Q5,)))
    with pytest.raises(OutputValidationError, match="FILTER identifiers without source definitions: unknown"):
        jasmine_module._sort_jasmine_vcf(raw, tmp_path / "sorted.vcf", (source,), overwrite=False)
    assert raw.read_bytes() == original
    assert not (tmp_path / "sorted.vcf").exists()
    assert not list(tmp_path.glob("*.tmp"))


def test_same_source_definitions_are_deduplicated_verbatim(tmp_path):
    sources = []
    for number in range(2):
        source = tmp_path / f"source{number}.vcf"
        source.write_bytes(payload(filters=(Q5, Q5)))
        sources.append(source)
    _, destination = finalize(tmp_path, [row()], sources=sources)
    assert destination.read_text(encoding="utf-8").count(Q5) == 1


@pytest.mark.parametrize("used", [True, False])
def test_conflicting_source_definition_only_blocks_missing_used_filter(tmp_path, used):
    sources = []
    for number, definition in enumerate((Q5, '##FILTER=<ID=q5,Description="Different meaning">')):
        source = tmp_path / f"source{number}.vcf"
        source.write_bytes(payload(filters=(definition,)))
        sources.append(source)
    if used:
        with pytest.raises(OutputValidationError, match="conflicting source definitions: q5"):
            finalize(tmp_path, [row()], sources=sources)
        assert not (tmp_path / "final.sorted.vcf").exists()
    else:
        _, destination = finalize(tmp_path, [row(filters="PASS")], sources=sources)
        assert b"##FILTER=" not in destination.read_bytes()


def test_native_existing_definition_is_not_overwritten_or_duplicated(tmp_path):
    native = '##FILTER=<ID=q5,Description="Native declaration already present">'
    source1, source2 = tmp_path / "source1.vcf", tmp_path / "source2.vcf"
    source1.write_bytes(payload(filters=(Q5,)))
    source2.write_bytes(payload(filters=('##FILTER=<ID=q5,Description="Different source">',)))
    _, destination = finalize(tmp_path, [row()], raw_filters=(native,), sources=(source1, source2))
    text = destination.read_text(encoding="utf-8")
    assert text.count(native) == 1
    assert Q5 not in text
    assert "Different source" not in text


def test_unterminated_native_body_is_rejected_instead_of_silently_changed(tmp_path):
    raw = tmp_path / "raw.vcf"
    original = payload([row(ending=b"")])
    raw.write_bytes(original)
    source = tmp_path / "source.vcf"
    source.write_bytes(payload(filters=(Q5,)))
    with pytest.raises(OutputValidationError, match="complete line ending"):
        jasmine_module._sort_jasmine_vcf(raw, tmp_path / "sorted.vcf", (source,), overwrite=False)
    assert raw.read_bytes() == original


def bgzf(path, data):
    compressor = zlib.compressobj(level=6, wbits=-15)
    compressed = compressor.compress(data) + compressor.flush()
    size = 18 + len(compressed) + 8
    header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", size - 1)
    path.write_bytes(header + compressed + struct.pack("<II", zlib.crc32(data), len(data)))


def test_wrapper_compressed_output_keeps_native_body_and_restored_header(tmp_path):
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\n" + "A" * 1000 + "\n>chr3\n" + "A" * 1000 + "\n", encoding="utf-8")
    Path(str(fasta) + ".fai").write_text("chr1\t1000\t6\t1000\t1001\nchr3\t1000\t1013\t1000\t1001\n", encoding="utf-8")
    reference = ReferenceGenome.from_fasta(fasta, build="GRCh38")
    original_records = [row("native", 0, "q5")]
    source_vcf = tmp_path / "source.vcf.gz"
    bgzf(source_vcf, payload(original_records, filters=(Q5,)))
    bgzf(Path(str(source_vcf) + ".tbi"), b"TBI\x01")
    source = SVEvidenceSourceArtifact("HG002", SVEvidenceSource.READ, "cutesv", source_vcf,
                                    Path(str(source_vcf) + ".tbi"), EvidenceRunStatus.COMPLETED)
    request = SVHarmonizationRequest("HG002", reference, (source,), tmp_path / "work",
        tmp_path / "out/HG002.harmonized.sv.vcf.gz", tmp_path / "out/HG002.sv.evidence.tsv")

    class NativeOutputFixtureRunner:
        def require_executable(self, executable):
            return Path(executable)

        def run(self, command, **kwargs):
            args = tuple(str(item) for item in command)
            if "--version" in args:
                version = "1.1.5" if args[0] == "jasmine" else "1.24"
                return CommandResult(args, 0, version, "", 0.0, None, True)
            if args[0] == "jasmine":
                (request.work_directory / "HG002.jasmine.raw.vcf").write_bytes(payload(original_records))
            elif args[0] == "bgzip":
                bgzf(Path(kwargs["stdout_path"]), Path(args[-1]).read_bytes())
            elif args[0] == "tabix":
                bgzf(request.output_index, b"TBI\x01")
            return CommandResult(args, 0, "", "", 0.0, None, True)

    result = JasmineWrapper(runner=NativeOutputFixtureRunner()).run(request)
    assert result.artifact is not None
    with gzip.open(request.output_vcf, "rb") as handle:
        final = handle.read()
    assert body(final) == original_records
    assert final.count((Q5 + "\n").encode()) == 1
    assert (request.work_directory / "HG002.jasmine.raw.vcf").read_bytes() == payload(original_records)
