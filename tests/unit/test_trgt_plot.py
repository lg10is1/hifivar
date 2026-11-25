"""Fresh-directory TR plotting contracts with real fake-executable processes."""
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import struct
import sys
import zlib

import pytest

from hifivar.command import CommandRunner
from hifivar.exceptions import CommandExecutionError, InputValidationError, OutputValidationError, ToolVersionError
from hifivar.reference import ReferenceGenome
from hifivar.tr import TandemRepeatCatalog
import hifivar.trgt_plot as plot


def bgzf(path, payload):
    compressor = zlib.compressobj(wbits=-15)
    compressed = compressor.compress(payload) + compressor.flush()
    total = 18 + len(compressed) + 8
    header = b"\x1f\x8b\x08\x04\0\0\0\0\0\xff\x06\0BC\x02\0" + struct.pack("<H", total - 1)
    path.write_bytes(header + compressed + struct.pack("<II", zlib.crc32(payload), len(payload)))


def vcf_bytes(*, sample="HG002", gt="1/1", duplicate=False, ref="ACGACGA", end=8, motifs="CGA", structure="(CGA)n", trid="L1"):
    row = f"chr1\t2\t.\t{ref}\tACGA\t.\t.\tTRID={trid};END={end};MOTIFS={motifs};STRUC={structure}\tGT\t{gt}\n"
    return ("##fileformat=VCFv4.3\n##contig=<ID=chr1,length=12>\n"
            '##INFO=<ID=TRID,Number=1,Type=String,Description="ID">\n'
            '##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motifs">\n'
            '##INFO=<ID=STRUC,Number=1,Type=String,Description="Structure">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n" + row * (2 if duplicate else 1)).encode()


@pytest.fixture
def case(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    fasta = inputs / "ref.fa"
    fasta.write_bytes(b">chr1\nAACGACGACGTT\n")
    Path(str(fasta) + ".fai").write_bytes(b"chr1\t12\t6\t12\t13\n")
    catalog = inputs / "repeats.bed"
    catalog.write_bytes(b"chr1\t2\t8\tID=L1;MOTIFS=CGA;STRUC=(CGA)n\n")
    vcf = inputs / "HG002.tr.vcf.gz"
    bgzf(vcf, vcf_bytes())
    bgzf(Path(str(vcf) + ".tbi"), b"TBI\x01")
    bam = inputs / "HG002.tr.spanning.bam"
    bgzf(bam, b"BAM\x01" + b"SYNTHETIC UNIT FIXTURE, NOT NATIVE TRGT EVIDENCE")
    Path(str(bam) + ".bai").write_bytes(b"BAI\x01\0\0\0\0")
    request = plot.TrgtPlotRequest("HG002", "L1", ReferenceGenome.from_fasta(fasta, build="GRCh38"),
        TandemRepeatCatalog(catalog, "GRCh38"), vcf, Path(str(vcf) + ".tbi"), bam, Path(str(bam) + ".bai"),
        {name: "0" * 64 for name in plot._ROLES}, tmp_path / "fresh-output")
    return rebind(request)


def hashes(request):
    return {role: hashlib.sha256(path.read_bytes()).hexdigest() for role, path in request.inputs.items()}


def rebind(request):
    return replace(request, expected_sha256=hashes(request))


def wrapper(tmp_path, *, mode="valid", version="5.1.0", version_stderr=None, runner=None):
    script = tmp_path / "fake native tool.py"
    script.write_text("""import json,sys,time
from pathlib import Path
args=sys.argv[1:]
mode=MODE
stderr_token=STDERR_TOKEN
if args==['--version']:
 print('trgt VERSION')
 if stderr_token is not None:print('trgt '+stderr_token,file=sys.stderr)
 if mode=='version-fail':
  print('fixture version failure',file=sys.stderr);sys.exit(7)
else:
 print('fixture plot stdout');print('fixture plot stderr',file=sys.stderr)
 image=Path(args[args.index('--image')+1])
 if mode=='fail':sys.exit(9)
 if mode=='timeout':time.sleep(3)
 if mode=='mutate':
  with Path(args[args.index('--vcf')+1]).open('ab') as out:out.write(b'MUTATION')
 if mode=='empty':image.write_bytes(b'')
 elif mode=='invalid':image.write_text('<not-svg/>')
 elif mode=='active':image.write_text('<svg xmlns="http://www.w3.org/2000/svg"><style>@import url(https://invalid.example/x)</style><rect width="1" height="1"/></svg>')
 else:image.write_text('<?xml version="1.0"?><svg xmlns="http://www.w3.org/2000/svg" width="10" height="10"><rect width="5" height="5" fill="#0072B2"/></svg>')
""".replace("MODE", repr(mode)).replace("STDERR_TOKEN", repr(version_stderr)).replace("VERSION", version), encoding="utf-8")
    return plot.TrgtPlotWrapper((sys.executable, str(script)), runner=runner)


def test_real_fake_process_success_keeps_source_and_seals_receipt_last(case, tmp_path):
    before = hashes(case)
    adapter = wrapper(tmp_path)
    result = adapter.run(case)
    assert result["status"] == "COMPLETE_TECHNICAL_SVG"
    assert result["tool_version"] == result["tool_version_reported"] == "5.1.0"
    assert result["source_origin"] == "EXTERNAL_FILES_NATIVE_PRODUCER_RECEIPT_NOT_VERIFIED"
    assert result["manual_review"] == "NOT_REVIEWED"
    assert result["accuracy"] == "NOT_ESTABLISHED"
    assert "NOT_PARENTAL" in result["allele_panel_order"]
    assert result["visual_validation"] == "STRUCTURAL_SVG_ONLY_NOT_PIXEL_OR_MANUAL_QA"
    assert result["inputs"]["vcf"]["sha256"] == before["vcf"]
    assert json.loads((case.output_dir / "completion.json").read_text()) == result
    assert (case.output_dir / "completion.json").read_bytes() == (case.output_dir / ".completion.pending.json").read_bytes()
    assert "fixture plot stdout" in (case.output_dir / "command.stdout.log").read_text()
    assert "fixture plot stderr" in (case.output_dir / "command.stderr.log").read_text()
    assert not (case.output_dir / "failure.json").exists()
    assert hashes(case) == before
    assert not any(flag in result["command"] for flag in ("--threads", "--sample-name", "genotype"))


def test_observed_trgt_build_is_normalized_without_losing_reported_identity(case, tmp_path):
    """A fake process reports the exact token observed in the retained native failure."""
    reported = "5.1.0-ec66463"
    before = hashes(case)
    result = wrapper(tmp_path, version=reported).run(case)
    receipt = json.loads((case.output_dir / "completion.json").read_text())
    assert result["tool_version"] == receipt["tool_version"] == "5.1.0"
    assert result["tool_version_reported"] == receipt["tool_version_reported"] == reported
    assert result["status"] == "COMPLETE_TECHNICAL_SVG"
    assert result["version_returncode"] == result["returncode"] == 0
    assert result["manual_review"] == "NOT_REVIEWED" and result["accuracy"] == "NOT_ESTABLISHED"
    raw_stdout = (case.output_dir / "version.stdout.log").read_bytes()
    raw_stderr = (case.output_dir / "version.stderr.log").read_bytes()
    assert raw_stdout == f"trgt {reported}{os.linesep}".encode()
    assert raw_stderr == b""
    assert result["logs"]["version.stdout.log"]["sha256"] == hashlib.sha256(raw_stdout).hexdigest()
    assert result["logs"]["version.stderr.log"]["sha256"] == hashlib.sha256(raw_stderr).hexdigest()
    assert (case.output_dir / "completion.json").read_bytes() == (case.output_dir / ".completion.pending.json").read_bytes()
    assert not (case.output_dir / "failure.json").exists()
    assert hashes(case) == before


@pytest.mark.parametrize("reported", ["5.1.0", "5.1.0-ec66463"])
def test_matching_stdout_stderr_version_reports_are_unambiguous(case, tmp_path, reported):
    result = wrapper(tmp_path, version=reported, version_stderr=reported).run(case)
    assert result["tool_version"] == "5.1.0"
    assert result["tool_version_reported"] == reported
    assert (case.output_dir / "version.stdout.log").read_bytes() == f"trgt {reported}{os.linesep}".encode()
    assert (case.output_dir / "version.stderr.log").read_bytes() == f"trgt {reported}{os.linesep}".encode()


@pytest.mark.parametrize("report,expected", [("version v5.1.0", "5.1.0"), ("VERSION V5.1.0-ec66463", "5.1.0-ec66463")])
def test_existing_version_prefix_conventions_are_preserved(case, tmp_path, report, expected):
    result = wrapper(tmp_path, version=report).run(case)
    assert result["tool_version"] == "5.1.0" and result["tool_version_reported"] == expected


def test_dry_run_needs_no_executable_or_directory(case):
    adapter = plot.TrgtPlotWrapper(("uninstalled-trgt-for-preview",))
    result = adapter.run(case, dry_run=True)
    assert result["status"] == "PLANNED" and result["executed"] is False
    assert not case.output_dir.exists()
    assert result["command"] == adapter.command(case)


@pytest.mark.parametrize("prefix", ["trgt", b"trgt", (), ("",), ("\x00",), 1, iter(["trgt"])])
def test_command_prefix_requires_an_actual_argument_sequence(prefix):
    with pytest.raises(InputValidationError):
        plot.TrgtPlotWrapper(prefix)


@pytest.mark.parametrize("field,value", [("sample_id", "bad sample"), ("repeat_id", "TRID:L1"), ("repeat_id", "L1;exit"),
    ("max_allele_reads", 0), ("max_allele_reads", -1), ("max_allele_reads", True), ("max_allele_reads", 1.5),
    ("timeout_seconds", 0), ("timeout_seconds", True), ("timeout_seconds", float("nan")),
    ("timeout_seconds", float("inf")), ("timeout_seconds", 1201), ("timeout_seconds", 10 ** 1000)])
def test_invalid_request_fields_fail_closed(case, field, value):
    with pytest.raises(InputValidationError):
        replace(case, **{field: value})


@pytest.mark.parametrize("change", ["missing", "extra", "upper", "not-mapping"])
def test_exact_seven_sha_roles(case, change):
    values = dict(case.expected_sha256)
    if change == "missing":
        values.pop("fai")
    elif change == "extra":
        values["attempt_id"] = "0" * 64
    elif change == "upper":
        values["vcf"] = "A" * 64
    else:
        values = None
    with pytest.raises(InputValidationError):
        replace(case, expected_sha256=values)


@pytest.mark.parametrize("sample,gt,duplicate,ref,end,motifs,structure,trid", [
    ("OTHER", "1/1", False, "ACGACGA", 8, "CGA", "(CGA)n", "L1"),
    ("HG002\tOTHER", "1/1\t1/1", False, "ACGACGA", 8, "CGA", "(CGA)n", "L1"),
    ("HG002", "./.", False, "ACGACGA", 8, "CGA", "(CGA)n", "L1"),
    ("HG002", "0/.", False, "ACGACGA", 8, "CGA", "(CGA)n", "L1"),
    ("HG002", "2/2", False, "ACGACGA", 8, "CGA", "(CGA)n", "L1"),
    ("HG002", "1/1/1", False, "ACGACGA", 8, "CGA", "(CGA)n", "L1"),
    ("HG002", "1/1", True, "ACGACGA", 8, "CGA", "(CGA)n", "L1"),
    ("HG002", "0/0", False, "TTTTTTT", 8, "CGA", "(CGA)n", "L1"),
    ("HG002", "1/1", False, "ACGACGA", 7, "CGA", "(CGA)n", "L1"),
    ("HG002", "1/1", False, "ACGACGA", 8, "TTT", "(CGA)n", "L1"),
    ("HG002", "1/1", False, "ACGACGA", 8, "CGA", "(TTT)n", "L1"),
    ("HG002", "1/1", False, "ACGACGA", 8, "CGA", "(CGA)n", "OTHER")])
def test_selected_vcf_contract_before_external_execution(case, sample, gt, duplicate, ref, end, motifs, structure, trid):
    bgzf(case.vcf, vcf_bytes(sample=sample, gt=gt, duplicate=duplicate, ref=ref, end=end, motifs=motifs, structure=structure, trid=trid))
    request = rebind(case)
    with pytest.raises((InputValidationError, OutputValidationError)):
        plot.TrgtPlotWrapper(("unavailable",)).run(request, dry_run=True)
    assert not request.output_dir.exists()


@pytest.mark.parametrize("gt", ["0", "1", "0/0", "1/1", "0|1", "1|0"])
def test_called_haploid_diploid_gt_context_is_retained(case, gt):
    bgzf(case.vcf, vcf_bytes(gt=gt))
    result = plot.TrgtPlotWrapper(("unavailable",)).run(rebind(case), dry_run=True)
    assert result["selection"]["GT"] == gt


@pytest.mark.parametrize("kind", ["bam", "bai", "tbi", "catalog-tags", "catalog-duplicate", "fai"])
def test_malformed_inputs_do_not_create_outputs(case, kind):
    if kind == "bam":
        bgzf(case.spanning_bam, b"NOTBAM")
    elif kind == "bai":
        case.bai.write_bytes(b"BAD INDEX")
    elif kind == "tbi":
        bgzf(case.tbi, b"NOTTBI")
    elif kind == "catalog-tags":
        case.catalog.path.write_bytes(b"chr1\t2\t8\tID=L1;MOTIFS=CGA;MOTIFS=TTT;STRUC=(CGA)n\n")
    elif kind == "catalog-duplicate":
        case.catalog.path.write_bytes(case.catalog.path.read_bytes() * 2)
    else:
        case.reference.fai.write_bytes(b"chr1\t12\t1000\t12\t13\n")
    with pytest.raises((InputValidationError, OutputValidationError)):
        plot.TrgtPlotWrapper(("unavailable",)).run(rebind(case), dry_run=True)
    assert not case.output_dir.exists()


def test_sha_mismatch_before_any_directory_or_probe(case):
    case.bai.write_bytes(case.bai.read_bytes() + b"drift")
    with pytest.raises(InputValidationError, match="SHA256"):
        plot.TrgtPlotWrapper(("unavailable",)).run(case)
    assert not case.output_dir.exists()


@pytest.mark.parametrize("kind", ["directory", "file", "input", "input-child", "tool"])
def test_existing_or_overlapping_output_is_never_modified(case, tmp_path, kind):
    adapter = wrapper(tmp_path)
    output = case.output_dir
    if kind == "directory":
        output.mkdir()
        (output / "OLD.txt").write_bytes(b"OLD")
    elif kind == "file":
        output.write_bytes(b"OLD")
    elif kind == "input":
        output = case.vcf
    elif kind == "input-child":
        output = case.vcf / "child"
    else:
        output = Path(adapter.prefix[1])
    before = hashes(case)
    with pytest.raises(OutputValidationError):
        adapter.run(replace(case, output_dir=output))
    assert hashes(case) == before
    if kind == "directory":
        assert (output / "OLD.txt").read_bytes() == b"OLD"
    elif kind == "file":
        assert output.read_bytes() == b"OLD"


def test_late_directory_creation_is_preserved(case, tmp_path, monkeypatch):
    adapter = wrapper(tmp_path)
    real = Path.mkdir
    def late(path, *args, **kwargs):
        if path == case.output_dir:
            real(path)
            (path / "foreign").write_bytes(b"FOREIGN")
        return real(path, *args, **kwargs)
    monkeypatch.setattr(Path, "mkdir", late)
    with pytest.raises(OutputValidationError):
        adapter.run(case)
    assert (case.output_dir / "foreign").read_bytes() == b"FOREIGN"
    assert not (case.output_dir / "completion.json").exists()


def test_output_junction_ancestry_rejected_before_any_tool_or_mkdir(case, monkeypatch):
    """Mock the portable path query, not native junction or TRGT execution."""
    before = hashes(case)
    monkeypatch.setattr(Path, "is_junction", lambda path: path == case.output_dir.parent, raising=False)
    monkeypatch.setattr(CommandRunner, "run", lambda *a, **kw: pytest.fail("No tool may run through a junction output"))
    with pytest.raises(OutputValidationError, match="junctions"):
        plot.TrgtPlotWrapper(("uninstalled-junction-fixture",)).run(case, dry_run=True)
    assert not case.output_dir.exists()
    assert hashes(case) == before


@pytest.mark.parametrize("name", ["plot.svg", "command.stdout.log", "command.stderr.log", "version.stdout.log", "version.stderr.log"])
def test_late_foreign_output_or_log_is_never_truncated(case, tmp_path, monkeypatch, name):
    adapter = wrapper(tmp_path)
    real = plot._reserve
    def late(path):
        if path.name == name:
            path.write_bytes(b"LATE FOREIGN BYTES")
        return real(path)
    monkeypatch.setattr(plot, "_reserve", late)
    with pytest.raises(OutputValidationError):
        adapter.run(case)
    assert (case.output_dir / name).read_bytes() == b"LATE FOREIGN BYTES"
    assert not (case.output_dir / "completion.json").exists()


def test_output_hardlink_alias_cannot_damage_an_input(case, tmp_path):
    os.link(case.spanning_bam, case.output_dir)
    before = hashes(case)
    with pytest.raises(OutputValidationError):
        wrapper(tmp_path).run(case)
    assert hashes(case) == before


@pytest.mark.parametrize("kind", ["direct", "ancestor", "dangling"])
def test_output_symlink_paths_are_rejected_before_execution(case, tmp_path, kind):
    link = tmp_path / "output-link"
    try:
        link.symlink_to(tmp_path / "missing" if kind == "dangling" else case.catalog.path.parent, target_is_directory=True)
    except OSError as error:
        pytest.skip(f"Local symlink privilege unavailable: {type(error).__name__}")
    output = link / "fresh" if kind == "ancestor" else link
    before = hashes(case)
    with pytest.raises(OutputValidationError):
        plot.TrgtPlotWrapper(("unavailable",)).run(replace(case, output_dir=output), dry_run=True)
    assert hashes(case) == before


@pytest.mark.parametrize("mode", ["fail", "empty", "invalid", "active", "version-fail", "timeout"])
def test_fake_tool_failure_keeps_evidence_without_completion(case, tmp_path, mode):
    before = hashes(case)
    request = replace(case, timeout_seconds=0.3) if mode == "timeout" else case
    with pytest.raises((CommandExecutionError, OutputValidationError)):
        wrapper(tmp_path, mode=mode).run(request)
    assert request.output_dir.exists()
    assert not (request.output_dir / "completion.json").exists()
    failure = json.loads((request.output_dir / "failure.json").read_text())
    assert failure["status"] == "FAILED_OR_INTERRUPTED_NO_NORMAL_RETURN"
    assert failure["protected_identity"] == "PASS"
    assert hashes(case) == before


@pytest.mark.parametrize("version", ["5.0.0", "5.1.1", "5.1.0-rc1", "5.1.0.dev0", "5.1.0+unreviewed", "unknown",
    "5.1.0-deadbee", "5.1.0-ec66463-rc1", "5.1.0-ec66463+unreviewed",
    "5.1.0(unreviewed)", "5.1.0-ec66463(unreviewed)", "5.1.0 -rc1",
    "5.1.0 (unreviewed)", "5.1.0-ec66463 -rc1", "5.1.0-ec66463 (unreviewed)"])
def test_exact_trgt_version_only(case, tmp_path, version):
    with pytest.raises(ToolVersionError):
        wrapper(tmp_path, version=version).run(case)
    assert not (case.output_dir / "completion.json").exists()
    assert not (case.output_dir / "command.stderr.log").exists()


@pytest.mark.parametrize("stdout_version,stderr_version", [
    ("5.1.0", "5.1.0-ec66463"),
    ("5.1.0-ec66463", "5.1.0"),
    ("5.1.0-ec66463", "5.1.0-deadbee"),
])
def test_conflicting_version_reports_fail_before_plot(case, tmp_path, stdout_version, stderr_version):
    before = hashes(case)
    with pytest.raises(ToolVersionError):
        wrapper(tmp_path, version=stdout_version, version_stderr=stderr_version).run(case)
    assert not (case.output_dir / "completion.json").exists()
    assert not (case.output_dir / "command.stdout.log").exists()
    assert not (case.output_dir / "command.stderr.log").exists()
    assert (case.output_dir / "version.stdout.log").read_bytes() == f"trgt {stdout_version}{os.linesep}".encode()
    assert (case.output_dir / "version.stderr.log").read_bytes() == f"trgt {stderr_version}{os.linesep}".encode()
    failure = json.loads((case.output_dir / "failure.json").read_text())
    assert failure["status"] == "FAILED_OR_INTERRUPTED_NO_NORMAL_RETURN"
    assert failure["protected_identity"] == "PASS"
    assert hashes(case) == before


def test_failed_command_records_actual_post_mutation(case, tmp_path):
    with pytest.raises(InputValidationError):
        wrapper(tmp_path, mode="mutate").run(case)
    failure = json.loads((case.output_dir / "failure.json").read_text())
    assert failure["protected_identity"] == "FAIL"
    assert failure["protected_after"]["inputs"]["vcf"]["sha256"] != case.expected_sha256["vcf"]
    assert not (case.output_dir / "completion.json").exists()


@pytest.mark.parametrize("stage", ["before-commit", "after-valid-commit"])
def test_interrupt_rethrows_same_object_and_retains_evidence(case, tmp_path, monkeypatch, stage):
    adapter = wrapper(tmp_path)
    original = KeyboardInterrupt("fixture interruption")
    real_link = os.link
    def interrupted(source, target, *args, **kwargs):
        if stage == "after-valid-commit":
            real_link(source, target, *args, **kwargs)
        raise original
    monkeypatch.setattr(plot.os, "link", interrupted)
    with pytest.raises(KeyboardInterrupt) as caught:
        adapter.run(case)
    assert caught.value is original
    assert (case.output_dir / "completion.json").exists() == (stage == "after-valid-commit")
    failure = json.loads((case.output_dir / "failure.json").read_text())
    assert failure["status"] == "FAILED_OR_INTERRUPTED_NO_NORMAL_RETURN" and failure["protected_identity"] == "PASS"


@pytest.mark.parametrize("kind", ["old-final", "pending-content", "pending-inode"])
def test_receipt_commit_faults_never_return_success_or_overwrite_foreign(case, tmp_path, monkeypatch, kind):
    adapter = wrapper(tmp_path)
    real_link = os.link
    def changed(source, target, *args, **kwargs):
        if kind == "old-final":
            Path(target).write_bytes(b"FOREIGN COMPLETION")
        elif kind == "pending-content":
            source = Path(source)
            source.write_bytes(source.read_bytes().replace(b"NOT_REVIEWED", b"YES_REVIEWED"))
        else:
            source = Path(source)
            other = source.with_name("foreign-receipt.json")
            other.write_bytes(source.read_bytes())
            source.unlink()
            real_link(other, source)
        return real_link(source, target, *args, **kwargs)
    monkeypatch.setattr(plot.os, "link", changed)
    with pytest.raises(OutputValidationError):
        adapter.run(case)
    assert (case.output_dir / "failure.json").exists()
    if kind == "old-final":
        assert (case.output_dir / "completion.json").read_bytes() == b"FOREIGN COMPLETION"


@pytest.mark.parametrize("field", ["plan_path", "plan_sha256"])
def test_plan_fields_must_be_paired(case, field):
    with pytest.raises(InputValidationError):
        replace(case, **{field: Path("plan.json") if field == "plan_path" else "0" * 64})


def test_plan_sha_bound_before_execution_and_in_receipt(case, tmp_path):
    plan = tmp_path / "plan.json"
    plan.write_bytes(b"EXPLICIT RAW PLAN BYTES")
    checksum = hashlib.sha256(plan.read_bytes()).hexdigest()
    request = replace(case, plan_path=plan, plan_sha256=checksum)
    result = wrapper(tmp_path).run(request)
    assert result["plan"]["sha256"] == checksum
    assert plan.read_bytes() == b"EXPLICIT RAW PLAN BYTES"


def test_plan_changed_after_loader_is_rejected_before_directory(case, tmp_path):
    plan = tmp_path / "plan.json"
    plan.write_bytes(b"ORIGINAL RAW PLAN")
    request = replace(case, plan_path=plan, plan_sha256=hashlib.sha256(plan.read_bytes()).hexdigest())
    plan.write_bytes(b"CHANGED RAW PLAN!")
    with pytest.raises(InputValidationError, match="plan bytes changed"):
        wrapper(tmp_path).run(request)
    assert not case.output_dir.exists()


@pytest.mark.parametrize("payload", [
    '<style>@import url(https://invalid.example/a)</style>',
    '<rect fill="url(https://invalid.example/a)"/>',
    '<g xml:base="https://invalid.example/"><use href="#x"/></g>',
    '<animate attributeName="href" to="https://invalid.example/a"/>',
    '<set attributeName="href" to="https://invalid.example/a"/>',
    '<image href="https://invalid.example/a"/>', '<script/>', '<foreignObject/>', '<rect onload="evil()"/>'])
def test_svg_rejects_css_remote_active_or_dynamic_content(tmp_path, payload):
    path = tmp_path / "unsafe.svg"
    path.write_text('<svg xmlns="http://www.w3.org/2000/svg">' + payload + '</svg>', encoding="utf-8")
    with pytest.raises(OutputValidationError):
        plot._validate_svg(path)


def test_svg_rejects_external_stylesheet_processing_instruction(tmp_path):
    path = tmp_path / "unsafe.svg"
    path.write_bytes(b'<?xml version="1.0"?><?xml-stylesheet href="https://invalid.example/x"?><svg xmlns="http://www.w3.org/2000/svg"><rect/></svg>')
    with pytest.raises(OutputValidationError):
        plot._validate_svg(path)


def test_same_request_cannot_be_reentered_or_reused(case, tmp_path):
    adapter = wrapper(tmp_path)
    adapter.run(case)
    before = {p.name: p.read_bytes() for p in case.output_dir.iterdir() if p.is_file()}
    with pytest.raises(OutputValidationError):
        adapter.run(case)
    assert {p.name: p.read_bytes() for p in case.output_dir.iterdir() if p.is_file()} == before
