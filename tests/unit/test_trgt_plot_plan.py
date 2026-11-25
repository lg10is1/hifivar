"""Strict YAML tests with tiny synthetic files, not native TRGT evidence."""
from __future__ import annotations

import hashlib
from pathlib import Path
import struct
import zlib

import pytest
import yaml

from hifivar.exceptions import InputValidationError
from hifivar.trgt_plot_plan import read_trgt_plot_plan, run_trgt_plot_plan


def _bgzf(path: Path, raw: bytes) -> None:
    compressor = zlib.compressobj(wbits=-15)
    compressed = compressor.compress(raw) + compressor.flush()
    total = 18 + len(compressed) + 8
    header = b"\x1f\x8b\x08\x04\0\0\0\0\0\xff\x06\0BC\x02\0" + struct.pack("<H", total - 1)
    path.write_bytes(header + compressed + struct.pack("<II", zlib.crc32(raw), len(raw)))


@pytest.fixture
def plot_plan_case(tmp_path):
    parent = tmp_path / "计划 inputs"
    parent.mkdir()
    inputs = parent / "inputs"
    inputs.mkdir()
    (parent / "results").mkdir()
    fasta = inputs / "reference.fa"
    fasta.write_bytes(b">chr1\nAACGACGACGTT\n")
    fai = Path(str(fasta) + ".fai")
    fai.write_text("chr1\t12\t6\t12\t13\n", encoding="utf-8")
    catalog = inputs / "catalog.bed"
    catalog.write_text("chr1\t2\t8\tID=L1;MOTIFS=CGA;STRUC=(CGA)n\n", encoding="utf-8")
    vcf = inputs / "HG002.tr.vcf.gz"
    _bgzf(vcf, (
        "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=12>\n"
        '##INFO=<ID=TRID,Number=1,Type=String,Description="ID">\n'
        '##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motifs">\n'
        '##INFO=<ID=STRUC,Number=1,Type=String,Description="Structure">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t2\t.\tACGACGA\tACGA\t.\tPASS\tTRID=L1;END=8;MOTIFS=CGA;STRUC=(CGA)n\tGT\t1/1\n"
    ).encode())
    tbi = Path(str(vcf) + ".tbi")
    _bgzf(tbi, b"TBI\x01")
    bam = inputs / "HG002.tr.spanning.bam"
    _bgzf(bam, b"BAM\x01SYNTHETIC NON-NATIVE BAM FIXTURE")
    bai = Path(str(bam) + ".bai")
    bai.write_bytes(b"BAI\x01SYNTHETIC NON-NATIVE INDEX FIXTURE")
    resources = {"fasta": fasta, "fai": fai, "catalog": catalog, "vcf": vcf,
                 "tbi": tbi, "spanning_bam": bam, "bai": bai}
    document = {
        "schema_version": 1, "sample_id": "HG002", "selector": "TRID:L1",
        "reference": {"fasta": "inputs/reference.fa", "fai": "inputs/reference.fa.fai", "build": "GRCh38"},
        "catalog": {"path": "inputs/catalog.bed", "build": "GRCh38"},
        "source": {"mode": "EXTERNAL_FILES", **{role: str(resources[role].relative_to(parent))
                   for role in ("vcf", "tbi", "spanning_bam", "bai")}},
        "expected_sha256": {role: hashlib.sha256(path.read_bytes()).hexdigest() for role, path in resources.items()},
        "runtime": {"executable": "missing-trgt-plan-fixture", "supported_version": "5.1.0"},
        "plot": {"type": "allele", "show": "motifs", "format": "svg", "max_allele_reads": 20},
        "execution": {"timeout_seconds": 60}, "output": {"directory": "results/new_plot"},
    }
    plan = parent / "plot.yaml"
    plan.write_text(yaml.safe_dump(document, sort_keys=False, allow_unicode=True), encoding="utf-8")
    return plan, document, resources


def _write(plan, document):
    plan.write_text(yaml.safe_dump(document, sort_keys=False, allow_unicode=True), encoding="utf-8")


def test_plan_resolves_from_its_parent_and_binds_raw_bytes(plot_plan_case, monkeypatch, tmp_path):
    plan, _, resources = plot_plan_case
    monkeypatch.chdir(tmp_path)
    parsed = read_trgt_plot_plan(plan)
    assert parsed.path == plan.absolute()
    assert parsed.sha256 == hashlib.sha256(plan.read_bytes()).hexdigest()
    assert parsed.request.plan_path == parsed.path and parsed.request.plan_sha256 == parsed.sha256
    assert parsed.request.repeat_id == "L1" and parsed.request.sample_id == "HG002"
    assert parsed.source_mode == "EXTERNAL_FILES"
    for name in ("vcf", "tbi", "spanning_bam", "bai"):
        assert getattr(parsed.request, name) == resources[name].absolute()
    assert parsed.request.output_dir == plan.parent / "results/new_plot"
    assert not parsed.request.output_dir.exists()


def test_plan_accepts_bom_crlf_and_relative_executable_path(plot_plan_case):
    plan, document, _ = plot_plan_case
    document["runtime"]["executable"] = "tools/trgt"
    raw = yaml.safe_dump(document, sort_keys=False).replace("\n", "\r\n").encode("utf-8-sig")
    plan.write_bytes(raw)
    parsed = read_trgt_plot_plan(plan)
    assert parsed.sha256 == hashlib.sha256(raw).hexdigest()
    assert parsed.executable == str(plan.parent / "tools/trgt")


@pytest.mark.parametrize("section", [None, "reference", "catalog", "source", "expected_sha256", "runtime", "plot", "execution", "output"])
def test_unknown_fields_at_every_mapping_level_fail(plot_plan_case, section):
    plan, document, _ = plot_plan_case
    (document if section is None else document[section])["unexpected"] = "not accepted"
    _write(plan, document)
    with pytest.raises(InputValidationError, match="unknown or missing"):
        read_trgt_plot_plan(plan)


@pytest.mark.parametrize("name", ["schema_version", "sample_id", "selector", "reference", "catalog", "source", "expected_sha256", "runtime", "plot", "execution", "output"])
def test_missing_required_fields_fail(plot_plan_case, name):
    plan, document, _ = plot_plan_case
    del document[name]
    _write(plan, document)
    with pytest.raises(InputValidationError, match="unknown or missing"):
        read_trgt_plot_plan(plan)


@pytest.mark.parametrize("section,key,value", [
    (None, "schema_version", True), (None, "schema_version", 2),
    (None, "sample_id", ""), (None, "selector", "L1"),
    (None, "selector", "TRID:"), (None, "selector", "TRID:TRID:L1"),
    (None, "selector", "TRID:L1,L2"), (None, "selector", "TRID:."),
    ("source", "mode", "GENERATED"), ("source", "mode", "VERIFIED_TRGT_RECEIPT"),
    ("runtime", "executable", "trgt --threads 4"), ("runtime", "executable", ["trgt"]),
    ("runtime", "supported_version", "5.2.0"), ("plot", "type", "waterfall"),
    ("plot", "show", "meth"), ("plot", "format", "pdf"),
    ("plot", "max_allele_reads", 0), ("plot", "max_allele_reads", True),
    ("execution", "timeout_seconds", True), ("execution", "timeout_seconds", 0),
    ("execution", "timeout_seconds", 1201), ("execution", "timeout_seconds", float("nan")),
    ("execution", "timeout_seconds", float("inf")),
    ("execution", "timeout_seconds", 10 ** 1000),
    ("catalog", "build", "GRCh37"), ("reference", "fai", "inputs/other.fai"),
    ("source", "vcf", "inputs/invalid\npath"), ("output", "directory", ""),
])
def test_unsupported_values_fail_before_tools(plot_plan_case, section, key, value):
    plan, document, _ = plot_plan_case
    (document if section is None else document[section])[key] = value
    _write(plan, document)
    with pytest.raises(InputValidationError):
        read_trgt_plot_plan(plan)


@pytest.mark.parametrize("role", ["fasta", "fai", "catalog", "vcf", "tbi", "spanning_bam", "bai"])
def test_each_hash_requires_exact_lowercase_hex(plot_plan_case, role):
    plan, document, _ = plot_plan_case
    document["expected_sha256"][role] = "A" * 64
    _write(plan, document)
    with pytest.raises(InputValidationError, match="SHA256"):
        read_trgt_plot_plan(plan)


@pytest.mark.parametrize("tail,match", [
    ("\nsample_id: OTHER\n", "duplicate"),
    ("\nunexpected:\n  key: 1\n  key: 2\n", "duplicate"),
    ("\nunexpected:\n  7: value\n", "keys must"),
    ("\nunexpected: &x {a: 1}\nsecond: *x\n", "aliases"),
    ("\nunexpected:\n  <<: {a: 1}\n", "merge"),
])
def test_ambiguous_yaml_is_rejected_recursively(plot_plan_case, tail, match):
    plan, _, _ = plot_plan_case
    plan.write_text(plan.read_text(encoding="utf-8") + tail, encoding="utf-8")
    with pytest.raises(InputValidationError, match=match):
        read_trgt_plot_plan(plan)


@pytest.mark.parametrize("raw", [b"\xff", b"", b"#" * (256 * 1024 + 1),
                                  b"# line\n" * 4097, b"[]", b"---\n{}\n---\n{}\n",
                                  b"unexpected: " + b"[" * 25 + b"0" + b"]" * 25,
                                  b"unexpected: " + b"9" * 5000],
                         ids=["invalid_utf8", "empty", "too_many_bytes", "too_many_lines",
                              "not_mapping", "multiple_documents", "too_deep", "integer_parser_limit"])
def test_invalid_or_excessive_plan_is_bounded(plot_plan_case, raw):
    plan, _, _ = plot_plan_case
    plan.write_bytes(raw)
    with pytest.raises(InputValidationError):
        read_trgt_plot_plan(plan)


def test_plan_change_during_request_construction_is_detected(plot_plan_case, monkeypatch):
    from hifivar.reference import ReferenceGenome
    plan, _, _ = plot_plan_case
    original = ReferenceGenome.from_fasta
    def changed(*args, **kwargs):
        result = original(*args, **kwargs)
        plan.write_text(plan.read_text(encoding="utf-8") + "\n# changed while parsed\n", encoding="utf-8")
        return result
    monkeypatch.setattr(ReferenceGenome, "from_fasta", changed)
    with pytest.raises(InputValidationError, match="changed"):
        read_trgt_plot_plan(plan)


def test_run_passes_only_declared_executable_and_request(plot_plan_case, monkeypatch):
    plan, document, _ = plot_plan_case
    observed = []
    class Wrapper:
        def __init__(self, command_prefix):
            assert command_prefix == (document["runtime"]["executable"],)
        def run(self, request, *, dry_run):
            observed.append((request, dry_run))
            return {"status": "PLANNED", "source_origin": "EXTERNAL_FILES_NATIVE_PRODUCER_RECEIPT_NOT_VERIFIED"}
    monkeypatch.setattr("hifivar.trgt_plot.TrgtPlotWrapper", Wrapper)
    result = run_trgt_plot_plan(plan, dry_run=True)
    assert result["status"] == "PLANNED" and observed[0][1] is True
    assert observed[0][0].plan_sha256 == hashlib.sha256(plan.read_bytes()).hexdigest()
    assert not observed[0][0].output_dir.exists()


def test_dry_run_type_is_not_coerced(plot_plan_case):
    with pytest.raises(InputValidationError, match="boolean"):
        run_trgt_plot_plan(plot_plan_case[0], dry_run=1)


def test_actual_cli_dry_run_validates_inputs_without_commands_or_writes(plot_plan_case, monkeypatch, capsys):
    from hifivar.cli import main
    from hifivar.command import CommandRunner
    import json

    plan, _, resources = plot_plan_case
    before = {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
              for path in (*resources.values(), plan)}
    original = CommandRunner.run
    def no_native(self, args, **kwargs):
        assert kwargs.get("dry_run") is True, "No real command or version probe is permitted"
        return original(self, args, **kwargs)
    monkeypatch.setattr(CommandRunner, "run", no_native)
    monkeypatch.setattr(CommandRunner, "require_executable", lambda *a: pytest.fail("No executable/version probe"))
    exit_code = main(["--dry-run", "review", "tr-plot", "--plan", str(plan)])
    captured = capsys.readouterr()
    assert exit_code == 0, captured.err
    assert json.loads(captured.out)["status"] == "PLANNED"
    assert not (plan.parent / "results/new_plot").exists()
    assert {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (*resources.values(), plan)} == before


def test_actual_cli_fake_executable_completes_technical_svg(plot_plan_case, capsys):
    """Exercise real subprocess routing, not an invented wrapper return value."""
    import json
    import os
    import sys

    from hifivar.cli import main

    plan, document, resources = plot_plan_case
    script = plan.parent / "fake_trgt.py"
    script.write_text(
        "import pathlib,sys\n"
        "if sys.argv[1:] == ['--version']:\n"
        "    print('trgt 5.1.0')\n"
        "elif sys.argv[1] == 'plot':\n"
        "    output = pathlib.Path(sys.argv[sys.argv.index('--image') + 1])\n"
        "    assert not (output.parent / 'completion.json').exists()\n"
        "    output.write_text('<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"10\" height=\"10\"><rect x=\"1\" y=\"1\" width=\"5\" height=\"5\" fill=\"#0072B2\"/></svg>', encoding='utf-8')\n"
        "else:\n"
        "    raise SystemExit(7)\n",
        encoding="utf-8",
    )
    if os.name == "nt":
        executable = plan.parent / "fake_trgt.cmd"
        executable.write_bytes(f'@echo off\r\nchcp 65001 >nul\r\n"{sys.executable}" "{script}" %*\r\n'.encode("utf-8"))
    else:
        executable = plan.parent / "fake_trgt"
        executable.write_text(f'#!{sys.executable}\n' + script.read_text(encoding="utf-8"), encoding="utf-8")
        executable.chmod(0o700)
    document["runtime"]["executable"] = str(executable)
    _write(plan, document)
    before = {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
              for path in (*resources.values(), plan)}
    exit_code = main(["review", "tr-plot", "--plan", str(plan)])
    captured = capsys.readouterr()
    assert exit_code == 0, captured.err
    result = json.loads(captured.out)
    assert result["status"] == "COMPLETE_TECHNICAL_SVG"
    assert result["manual_review"] == "NOT_REVIEWED"
    assert result["accuracy"] == "NOT_ESTABLISHED"
    assert result["source_origin"] == "EXTERNAL_FILES_NATIVE_PRODUCER_RECEIPT_NOT_VERIFIED"
    output = plan.parent / "results/new_plot"
    receipt = json.loads((output / "completion.json").read_text(encoding="utf-8"))
    assert receipt["status"] == "COMPLETE_TECHNICAL_SVG"
    assert receipt["plan"]["sha256"] == before[str(plan)]
    assert (output / "plot.svg").is_file()
    assert {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (*resources.values(), plan)} == before
