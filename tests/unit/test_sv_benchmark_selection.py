"""Fake optional dependencies test contracts, not native matching accuracy."""

from __future__ import annotations

import gzip
import io
import json
from pathlib import Path
import re
import struct
from types import SimpleNamespace
import zlib

import pytest

from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar import sv_benchmark_selection as selection


def bgzf(path: Path, data: bytes) -> None:
    obj = zlib.compressobj(wbits=-15)
    body = obj.compress(data) + obj.flush()
    header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
    path.write_bytes(header + struct.pack("<H", 25 + len(body)) + body
                     + struct.pack("<II", zlib.crc32(data), len(data)))


def make_vcf(tmp_path: Path, name="input.vcf.gz", *, body="", sample="HG002", extra="") -> Path:
    path = tmp_path / name
    text = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n'
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '##FILTER=<ID=LowQual,Description="Low quality">\n' + extra
            + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + sample + '\n' + body)
    bgzf(path, text.encode())
    bgzf(Path(f"{path}.tbi"), b"TBI\x01")
    return path


def row(pos, vid, alt="<INS>", info="SVTYPE=INS;SVLEN=60", flt="PASS", fmt="GT", value="0|1"):
    return f"chr1\t{pos}\t{vid}\tA\t{alt}\t30.10\t{flt}\t{info}\t{fmt}\t{value}\n"


class FakeVariantFile:
    def __init__(self, path):
        self.path = Path(path)
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            self.lines = list(stream)
        self.header = SimpleNamespace(samples=[], contigs={}, info={}, formats={}, filters={})
        for line in self.lines:
            if line.startswith("##contig"):
                match = re.search(r"ID=([^,>]+),length=(\d+)", line)
                self.header.contigs[match[1]] = SimpleNamespace(length=int(match[2]))
            for prefix, target in (("INFO", self.header.info), ("FORMAT", self.header.formats)):
                if line.startswith("##" + prefix):
                    match = re.search(r"ID=([^,]+),Number=([^,]+),Type=([^,>]+)", line)
                    target[match[1]] = SimpleNamespace(number=int(match[2]) if match[2].isdigit() else match[2], type=match[3])
            if line.startswith("##FILTER"):
                self.header.filters[re.search(r"ID=([^,>]+)", line)[1]] = True
            if line.startswith("#CHROM"):
                self.header.samples = line.rstrip().split("\t")[9:]

    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False

    def __iter__(self):
        for line in self.lines:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info = {}
            for part in fields[7].split(";"):
                key, sep, value = part.partition("=")
                if key != ".":
                    info[key] = int(value) if key == "SVLEN" else (value if sep else True)
            formats, values = fields[8].split(":"), fields[9].split(":")
            values += ["."] * (len(formats) - len(values))
            yield SimpleNamespace(contig=fields[0], pos=int(fields[1]), id=fields[2] if fields[2] != "." else None,
                                  ref=fields[3], alts=tuple(fields[4].split(",")), info=info,
                                  samples={self.header.samples[0]: dict(zip(formats, values))})


class FakeBGZFile(io.BytesIO):
    def __init__(self, path, mode):
        super().__init__()
        self.path = Path(path)
        assert mode == "wb"

    def __exit__(self, *_):
        bgzf(self.path, self.getvalue())
        self.close()
        return False


class FakeTabix:
    def __init__(self, path, index):
        with gzip.open(index, "rb") as handle:
            assert handle.read(4) == b"TBI\x01"
        self.path = path
        self.contigs = sorted({record.contig for record in FakeVariantFile(path)})

    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False

    def fetch(self, *_):
        return iter(())


class FakeNativeRecord:
    def __init__(self, parsed, params):
        assert params.dup_to_ins is False
        self.parsed = parsed

    def var_type(self):
        parsed = self.parsed
        if any(mark in parsed.alts[0] for mark in "[]"):
            kind = "BND"
        elif "SVTYPE" in parsed.info:
            kind = parsed.info["SVTYPE"]
        elif parsed.alts[0] == "*":
            kind = "UNK"
        elif len(parsed.ref) < len(parsed.alts[0]):
            kind = "INS"
        elif len(parsed.ref) > len(parsed.alts[0]):
            kind = "DEL"
        else:
            kind = "SNP"
        return SimpleNamespace(name=kind)

    def var_size(self):
        return abs(self.parsed.info["SVLEN"]) if "SVLEN" in self.parsed.info else abs(len(self.parsed.ref) - len(self.parsed.alts[0]))


@pytest.fixture
def runtime(tmp_path, monkeypatch):
    source = tmp_path / "fake_native_source.py"
    source.write_text("# TEST DOUBLE: no claim of native biological matching\n", encoding="utf-8")
    def tabix_index(path, *, preset, force):
        assert preset == "vcf" and force is False
        bgzf(Path(f"{path}.tbi"), b"TBI\x01")
    result = selection._NativeRuntime(
        SimpleNamespace(VariantParams=lambda **kw: SimpleNamespace(**kw), VariantRecord=FakeNativeRecord),
        SimpleNamespace(VariantFile=FakeVariantFile, BGZFile=FakeBGZFile, TabixFile=FakeTabix, tabix_index=tabix_index),
        {"truvari_version": "5.4.0", "pysam_version": "TEST_DOUBLE", "sources": {"fake": {"path": str(source), "sha256": selection.compute_sha256(source)}},
         "test_double": True, "dup_to_ins": False},
    )
    monkeypatch.setattr(selection, "_load_native_runtime", lambda: result)
    return result


def test_selection_preserves_all_retained_fields_order_filter_and_short_records(tmp_path, runtime):
    kept = [row(10, "short", "A" * 31, "SVLEN=30", "LowQual"), row(20, "del", "<DEL>", "SVTYPE=DEL;SVLEN=-100", ".")]
    body = kept[0] + row(15, "dup", "A" * 100, "SVTYPE=DUP;SVLEN=99") + kept[1]
    body += row(30, "star", "*", ".") + row(40, "bnd", "A]chr1:50]", "SVTYPE=INS;SVLEN=60")
    source = make_vcf(tmp_path, body=body)
    old = source.read_bytes(), Path(f"{source}.tbi").read_bytes()
    result = selection.derive_sv_benchmark_query(source, tmp_path / "new", sample_id="HG002")
    assert result["counts"] == {"source_records": 5, "retained_INS": 1, "retained_DEL": 1}
    assert result["excluded_by_native_type"] == {"BND": 1, "DUP": 1, "UNK": 1}
    with gzip.open(result["output_paths"]["vcf"], "rt", newline="") as stream:
        text = stream.read()
    assert "".join(line for line in text.splitlines(keepends=True) if not line.startswith("#")) == "".join(kept)
    assert "Canonical=false" in text
    assert (source.read_bytes(), Path(f"{source}.tbi").read_bytes()) == old
    assert result["source_immutability"].startswith("PASS_FULL_SHA256")


def test_selection_star_with_native_del_annotation_is_not_custom_prefiltered(tmp_path, runtime):
    source = make_vcf(tmp_path, body=row(10, "annotated_star", "*", "SVTYPE=DEL;SVLEN=-80"))
    result = selection.derive_sv_benchmark_query(source, tmp_path / "new", sample_id="HG002")
    assert result["counts"]["retained_DEL"] == 1


def test_zero_eligible_query_is_valid_not_perfect_accuracy(tmp_path, runtime):
    source = make_vcf(tmp_path, body=row(10, "only_inv", "<INV>", "SVTYPE=INV;SVLEN=70"))
    result = selection.derive_sv_benchmark_query(source, tmp_path / "new", sample_id="HG002")
    assert result["counts"]["retained_INS"] == result["counts"]["retained_DEL"] == 0
    assert result["scientific_accuracy"] == "NOT_ASSESSED"


@pytest.mark.parametrize("body", [row(1, "multi", "A,AT"), row(2, "bad_gt", value="0/2"), row(2, "unknown_info", info="MADEUP=2"), "broken\n", "\n"])
def test_malformed_or_multiallelic_fails_without_published_pair(tmp_path, runtime, body):
    source = make_vcf(tmp_path, body=body)
    with pytest.raises(OutputValidationError, match="derivation failed"):
        selection.derive_sv_benchmark_query(source, tmp_path / "new", sample_id="HG002")
    manifest = json.loads((tmp_path / "new/selection_manifest.json").read_text())
    assert manifest["status"] == "FAILED"
    assert not (tmp_path / "new/query.ins_del.vcf.gz").exists()


@pytest.mark.parametrize("position", [0, 1, 100000, 100001])
def test_vcf_telomeric_and_ordinary_positions_preserved_verbatim(tmp_path, runtime, position):
    original_row = row(position, "position_boundary", flt="LowQual", value="0/0")
    source = make_vcf(tmp_path, body=original_row)
    source_before = source.read_bytes(), Path(f"{source}.tbi").read_bytes()
    result = selection.derive_sv_benchmark_query(source, tmp_path / "new", sample_id="HG002")
    with gzip.open(result["output_paths"]["vcf"], "rt", newline="") as stream:
        retained_rows = [line for line in stream if not line.startswith("#")]
    assert retained_rows == [original_row]
    assert result["counts"] == {"source_records": 1, "retained_INS": 1, "retained_DEL": 0}
    assert (source.read_bytes(), Path(f"{source}.tbi").read_bytes()) == source_before
    assert result["scientific_accuracy"] == "NOT_ASSESSED"


@pytest.mark.parametrize("position", [-1, 100002, "1.0", "1e1", "NaN", "", "٠", "+1"])
def test_illegal_or_malformed_positions_still_fail_closed(tmp_path, runtime, position):
    source = make_vcf(tmp_path, body=row(position, "invalid_position"))
    source_before = source.read_bytes(), Path(f"{source}.tbi").read_bytes()
    with pytest.raises(OutputValidationError, match="derivation failed"):
        selection.derive_sv_benchmark_query(source, tmp_path / "new", sample_id="HG002")
    manifest = json.loads((tmp_path / "new/selection_manifest.json").read_text())
    assert manifest["status"] == "FAILED"
    assert not (tmp_path / "new/query.ins_del.vcf.gz").exists()
    assert (source.read_bytes(), Path(f"{source}.tbi").read_bytes()) == source_before


def test_dry_run_has_no_dependency_import_or_writes(tmp_path, monkeypatch):
    source = make_vcf(tmp_path)
    monkeypatch.setattr(selection, "_load_native_runtime", lambda: pytest.fail("native imported during dry run"))
    result = selection.derive_sv_benchmark_query(source, tmp_path / "new", sample_id="HG002", dry_run=True)
    assert result["status"] == "PLANNED" and not (tmp_path / "new").exists()


@pytest.mark.parametrize("sample", ["UnknownSampleName", "", "different"])
def test_sample_identity_is_not_repaired(tmp_path, sample):
    source = make_vcf(tmp_path)
    with pytest.raises(InputValidationError):
        selection.derive_sv_benchmark_query(source, tmp_path / "new", sample_id=sample, dry_run=True)


def test_existing_output_is_never_overwritten(tmp_path, runtime):
    source = make_vcf(tmp_path)
    target = tmp_path / "old"
    target.mkdir()
    (target / "sentinel").write_text("keep")
    with pytest.raises(OutputValidationError, match="fresh"):
        selection.derive_sv_benchmark_query(source, target, sample_id="HG002")
    assert (target / "sentinel").read_text() == "keep"


def test_source_mutation_is_not_accepted(tmp_path, runtime):
    source = make_vcf(tmp_path, body=row(10, "a"))
    old = runtime.pysam.tabix_index
    def mutate(*args, **kwargs):
        old(*args, **kwargs)
        source.touch()
    runtime.pysam.tabix_index = mutate
    with pytest.raises(OutputValidationError, match="SOURCE_IMMUTABILITY_FAIL"):
        selection.derive_sv_benchmark_query(source, tmp_path / "new", sample_id="HG002")
    assert json.loads((tmp_path / "new/selection_manifest.json").read_text())["source_immutability"] == "FAIL"


def test_missing_optional_dependencies_explicit(monkeypatch):
    def missing(_):
        raise ImportError("not installed")
    monkeypatch.setattr(selection.importlib, "import_module", missing)
    with pytest.raises(ToolVersionError, match="optional Truvari"):
        selection._load_native_runtime()


def test_wrong_native_version_fails_before_source_lookup(monkeypatch):
    monkeypatch.setattr(selection.importlib, "import_module", lambda name: SimpleNamespace(__version__="5.5.0"))
    with pytest.raises(ToolVersionError, match="exactly Truvari 5.4.0"):
        selection._load_native_runtime()


def test_classifier_identity_must_match_actual_source_hash(tmp_path, runtime):
    source = make_vcf(tmp_path)
    runtime.identity["sources"]["fake"]["sha256"] = "0" * 64
    with pytest.raises(ToolVersionError, match="source changed"):
        selection.derive_sv_benchmark_query(source, tmp_path / "new", sample_id="HG002")
    assert not (tmp_path / "new").exists()


def test_library_index_failure_is_preserved_explicitly(tmp_path, runtime):
    source = make_vcf(tmp_path, body=row(10, "a"))
    def fail_index(*_args, **_kwargs):
        raise OSError("index failed")
    runtime.pysam.tabix_index = fail_index
    with pytest.raises(OutputValidationError, match="index failed"):
        selection.derive_sv_benchmark_query(source, tmp_path / "new", sample_id="HG002")
    assert json.loads((tmp_path / "new/selection_manifest.json").read_text())["status"] == "FAILED"


def test_vcf_legal_trailing_missing_format_field_is_preserved(tmp_path, runtime):
    source = make_vcf(tmp_path, body=row(10, "a", fmt="GT:DP", value="0/1"),
                      extra='##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
    result = selection.derive_sv_benchmark_query(source, tmp_path / "new", sample_id="HG002")
    with gzip.open(result["output_paths"]["vcf"], "rt") as handle:
        assert "GT:DP\t0/1\n" in handle.read()
