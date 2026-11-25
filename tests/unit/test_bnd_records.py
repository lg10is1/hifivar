"""Byte preservation, mate QC, ledger tamper and header-budget regressions."""
import base64
import gzip
import hashlib
import json
from pathlib import Path

import pytest

from hifivar.bnd_records import derive_record_stream, verify_record_round_trip
from hifivar.breakends import BndPolicy, POLICY, POLICY_SHA256
from hifivar.exceptions import InputValidationError, OutputValidationError


HEADER = ("##fileformat=VCFv4.4\n"
          "##contig=<ID=chr1,length=1000>\n##contig=<ID=chr2,length=2000>\n"
          '##INFO=<ID=END,Number=1,Type=Integer,Description="Local end">\n'
          '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Length, native metadata">\n'
          '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence, not identity">\n'
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n")
POLICY_OBJECT = BndPolicy(POLICY["policy_id"], 1, POLICY_SHA256, POLICY["mode"])


def record(identifier="b1", *, pos=20, alt="A[chr2:10[", info="SVTYPE=BND;END=10;SVLEN=0", ref="A"):
    return f"chr1\t{pos}\t{identifier}\t{ref}\t{alt}\t7.1\tLowQual\t{info}\tGT:PS\t1|0:123\n"


def derive(tmp_path, text, *, compressed=False):
    native = tmp_path / "native.vcf.gz" if compressed else tmp_path / "native.vcf"
    data = text.encode("utf-8")
    native.write_bytes(gzip.compress(data, mtime=0) if compressed else data)
    paths = [tmp_path / name for name in ("derived.vcf", "records.jsonl", "junctions.tsv")]
    summary = derive_record_stream(native, *paths, reference_lengths={"chr1": 1000, "chr2": 2000},
        sample_id="HG002", native_sha256=hashlib.sha256(native.read_bytes()).hexdigest(),
        policy=POLICY_OBJECT, work_directory=tmp_path)
    return native, *paths, summary


@pytest.mark.parametrize("ending", ["\n", "\r\n"])
@pytest.mark.parametrize("compressed", [False, True])
def test_records_are_lossless_only_with_mandatory_ledger(tmp_path, ending, compressed):
    non_bnd = record("ins0", pos=0, ref="N", alt="<INS>", info="SVTYPE=INS;END=0;SVLEN=70;KEEP=x%3By")
    bnd = record(info="X=one%3Btwo;END=00010;SVTYPE=BND;SVLEN=+0")
    # '+0' is not accepted as an integer token; use a preserved valid lexical zero.
    bnd = bnd.replace("SVLEN=+0", "SVLEN=-0")
    original = (HEADER + non_bnd + bnd).replace("\n", ending)
    native, changed, ledger, junctions, summary = derive(tmp_path, original, compressed=compressed)
    rows = [json.loads(line) for line in ledger.read_text().splitlines()]
    records = [line for line in changed.read_bytes().splitlines(keepends=True) if not line.startswith(b"#")]
    assert records[0] == non_bnd.replace("\n", ending).encode()
    assert rows[0]["action"] == "UNCHANGED_NON_BND"
    assert rows[1]["removed_info_tokens"] == [{"token_index": 1, "raw_token": "END=00010"}, {"token_index": 3, "raw_token": "SVLEN=-0"}]
    assert b"END=" not in records[1] and b"SVLEN=" not in records[1]
    assert b"X=one%3Btwo;SVTYPE=BND;HIFIVAR_BND_REPR=" in records[1]
    assert records[1].split(b"\t")[:7] == bnd.replace("\n", ending).encode().split(b"\t")[:7]
    assert records[1].split(b"\t")[8:] == bnd.replace("\n", ending).encode().split(b"\t")[8:]
    assert summary["round_trip"]["status"] == "PASS"
    assert summary["records"] == 2 and summary["bnd_records"] == 1
    assert len(junctions.read_text().splitlines()) == 2
    assert base64.b64decode(rows[1]["native_info_b64"]) == b"X=one%3Btwo;END=00010;SVTYPE=BND;SVLEN=-0"
    assert verify_record_round_trip(native, changed, ledger)["records"] == 2


@pytest.mark.parametrize("info", ["SVTYPE=BND", "SVTYPE=BND;END=.;SVLEN=.", "SVTYPE=BND;END=10", "SVTYPE=BND;SVLEN=0"])
def test_absent_or_missing_removed_values_are_explicit(tmp_path, info):
    _, changed, ledger, _, _ = derive(tmp_path, HEADER + record(info=info))
    row = json.loads(ledger.read_text())
    assert len(row["removed_info_tokens"]) == info.count(";")
    assert b"END=" not in changed.read_bytes().splitlines()[-1]


@pytest.mark.parametrize("info", [
    "SVTYPE=BND;END=10;END=11", "SVTYPE=BND;;END=10", "SVTYPE=BND;END=", "SVTYPE=BND;END",
    "SVTYPE=BND;SVLEN=0,1", "SVTYPE=BND;CIPOS=one,two", "SVTYPE=BND;CIPOS=1,-1",
    "SVTYPE=BND;MATEID", "SVTYPE=BND;MATEID=a,", "SVTYPE;END=10",
    "SVTYPE=BND;HIFIVAR_BND_REPR=old", "SVTYPE=BND;HIFIVAR_BND_LEDGER=old",
    "SVTYPE=BND;X=unescaped space", "SVTYPE=BND;X=a=b",
])
def test_bad_or_already_derived_info_fails(tmp_path, info):
    with pytest.raises(InputValidationError):
        derive(tmp_path, HEADER + record(info=info))


@pytest.mark.parametrize("header", [
    HEADER.replace("VCFv4.4", "VCFv4.3"),
    HEADER.replace("ID=END,Number=1", "ID=END,Number=A"),
    HEADER.replace("ID=SVLEN,Number=A,Type=Integer", "ID=SVLEN,Number=A,Type=String"),
    HEADER.replace("ID=CIPOS,Number=2,Type=Integer", "ID=CIPOS,Number=2,Type=String"),
    HEADER.replace("chr1,length=1000", "chr1,length=1001"),
    HEADER.replace("##contig=<ID=chr1,length=1000>\n", ""),
    HEADER.replace("\tHG002\n", "\tOTHER\n"),
    HEADER.replace("##INFO=<ID=END", '##INFO=<ID=HIFIVAR_BND_REPR,Number=1,Type=String,Description="x">\n##INFO=<ID=END'),
    HEADER.replace("#CHROM", '##INFO=<ID=END,Number=1,Type=Integer,Description="duplicate">\n#CHROM'),
])
def test_header_contract_is_not_silently_repaired(tmp_path, header):
    with pytest.raises(InputValidationError):
        derive(tmp_path, header + record())


def test_duplicate_native_records_have_distinct_ordinal_keys_and_mate_ambiguity(tmp_path):
    text = HEADER + record("b1", info="SVTYPE=BND;MATEID=b1") * 2
    _, _, ledger, _, summary = derive(tmp_path, text)
    rows = [json.loads(line) for line in ledger.read_text().splitlines()]
    assert rows[0]["native_record_sha256"] == rows[1]["native_record_sha256"]
    assert rows[0]["ledger_key"] != rows[1]["ledger_key"]
    assert rows[0]["mate_qc"]["links"][0]["id_resolution"] == "AMBIGUOUS_DUPLICATE_ID"
    assert summary["records"] == 2


def test_microhomology_does_not_require_exact_mate_coordinate(tmp_path):
    one = record("b1", alt="A[chr2:12[", info="SVTYPE=BND;MATEID=b2;CIPOS=-1,2")
    two = "chr2\t10\tb2\tC\t]chr1:20]C\t.\tPASS\tSVTYPE=BND;MATEID=b1;CIPOS=-3,3\tGT\t0/1\n"
    _, _, ledger, _, _ = derive(tmp_path, HEADER + one + two)
    link = json.loads(ledger.read_text().splitlines()[0])["mate_qc"]["links"][0]
    assert link["id_resolution"] == "RESOLVED" and link["reciprocal_id_link"]
    assert link["orientation_consistent"] and not link["exact_mate_coordinate"]
    assert link["coordinate_consistency"] == "WITHIN_EXPLICIT_MATE_CIPOS"


def test_unresolved_mate_is_recorded_without_synthetic_link(tmp_path):
    _, _, ledger, _, summary = derive(tmp_path, HEADER + record(info="SVTYPE=BND;MATEID=original_unprefixed"))
    row = json.loads(ledger.read_text())
    assert row["mate_qc"]["links"][0]["id_resolution"] == "UNRESOLVED"
    assert row["mate_ids"] == ["original_unprefixed"]
    assert summary["unresolved_mate_links"] == 1


@pytest.mark.parametrize("change", ["genotype", "info", "ledger", "count"])
def test_round_trip_detects_tamper(tmp_path, change):
    native, changed, ledger, _, _ = derive(tmp_path, HEADER + record())
    if change == "genotype":
        changed.write_bytes(changed.read_bytes().replace(b"1|0:123", b"0|0:123"))
    elif change == "info":
        changed.write_bytes(changed.read_bytes().replace(b"SVTYPE=BND", b"SVTYPE=DEL"))
    elif change == "ledger":
        row = json.loads(ledger.read_text())
        row["native_info_b64"] = base64.b64encode(b"SVTYPE=BND;END=999").decode()
        ledger.write_text(json.dumps(row) + "\n")
    else:
        with changed.open("ab") as writer:
            writer.write(record().encode())
    with pytest.raises(OutputValidationError):
        verify_record_round_trip(native, changed, ledger)


def test_existing_output_and_missing_terminator_fail_without_overwrite(tmp_path):
    native, changed, ledger, junctions, _ = derive(tmp_path, HEADER + record())
    before = changed.read_bytes()
    with pytest.raises(OutputValidationError):
        derive_record_stream(native, changed, ledger, junctions, reference_lengths={"chr1": 1000, "chr2": 2000},
            sample_id="HG002", native_sha256="a" * 64, policy=POLICY_OBJECT, work_directory=tmp_path)
    assert changed.read_bytes() == before
    sub = tmp_path / "unterminated"
    sub.mkdir()
    with pytest.raises(InputValidationError, match="terminator"):
        derive(sub, HEADER + record().rstrip("\n"))
