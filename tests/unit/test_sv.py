from pathlib import Path
import struct
import zlib

import pytest

from hifivar.command import CommandResult
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolNotFoundError
from hifivar.reference import ReferenceGenome
from hifivar.sv import BgzipTabixWrapper, SvCaller, VcfFinalizeRequest, VcfFinalizeStatus, create_structural_variant_artifact


def write_bgzf(path: Path, payload: bytes) -> None:
    compressor = zlib.compressobj(level=6, wbits=-15)
    compressed = compressor.compress(payload) + compressor.flush()
    total_size = 18 + len(compressed) + 8
    header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", total_size - 1)
    footer = struct.pack("<II", zlib.crc32(payload), len(payload))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(header + compressed + footer)


def vcf_text(*, sample="S1", contig="chr1", svtype=True):
    info = '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n' if svtype else ""
    return f"##fileformat=VCFv4.3\n##contig=<ID={contig},length=4>\n{info}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"


@pytest.fixture
def reference(tmp_path):
    fasta = tmp_path / "参考.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    return ReferenceGenome.from_fasta(fasta, build="GRCh38")


def make_artifact(tmp_path, reference, *, caller=SvCaller.PBSV, sample="S1", contig="chr1", svtype=True):
    vcf = tmp_path / f"S1.{caller.value}.sv.vcf.gz"
    write_bgzf(vcf, vcf_text(sample=sample, contig=contig, svtype=svtype).encode())
    write_bgzf(Path(f"{vcf}.tbi"), b"TBI\x01")
    return create_structural_variant_artifact(caller=caller, sample_id="S1", reference=reference, vcf_path=vcf, caller_version="1.0", commands=((caller.value, "call"),))


def test_common_artifact_records_identity_and_is_not_harmonized(tmp_path, reference):
    artifact = make_artifact(tmp_path, reference)
    assert artifact.to_dict()["harmonized"] is False
    assert artifact.to_dict()["reference_build"] == "GRCh38"
    assert artifact.commands == (("pbsv", "call"),)


@pytest.mark.parametrize("failure", ["sample", "contig", "svtype", "bgzf", "index"])
def test_lightweight_validation_rejects_invalid_artifacts(tmp_path, reference, failure):
    if failure in {"sample", "contig", "svtype"}:
        with pytest.raises(OutputValidationError):
            make_artifact(tmp_path, reference, sample="wrong" if failure == "sample" else "S1", contig="1" if failure == "contig" else "chr1", svtype=failure != "svtype")
        return
    vcf = tmp_path / "S1.pbsv.sv.vcf.gz"
    if failure == "bgzf":
        vcf.write_text(vcf_text(), encoding="utf-8")
    else:
        write_bgzf(vcf, vcf_text().encode())
    index = Path(f"{vcf}.tbi")
    index.write_bytes(b"bad")
    with pytest.raises(OutputValidationError):
        create_structural_variant_artifact(caller=SvCaller.PBSV, sample_id="S1", reference=reference, vcf_path=vcf, caller_version=None, commands=())


class FakeRunner:
    def __init__(self, request=None, missing=False):
        self.request, self.missing, self.calls = request, missing, []

    def require_executable(self, executable):
        if self.missing:
            raise ToolNotFoundError("missing")
        return Path(executable)

    def run(self, command, **kwargs):
        args = tuple(str(item) for item in command)
        self.calls.append((args, kwargs))
        if "--version" in args:
            return CommandResult(args, 0, f"{args[0]} (htslib) 1.21", "", 0.1, None, True)
        if args[0] == "bgzip" and not kwargs.get("dry_run"):
            write_bgzf(Path(kwargs["stdout_path"]), self.request.source_vcf.read_bytes())
        if args[0] == "tabix" and not kwargs.get("dry_run"):
            write_bgzf(self.request.output_index, b"TBI\x01")
        return CommandResult(args, 0, "", "", 0.1, None, not kwargs.get("dry_run", False))


def test_plain_vcf_finalization_is_two_explicit_commands(tmp_path, reference):
    source = tmp_path / "S1.pbsv.raw.vcf"
    source.write_text(vcf_text(), encoding="utf-8")
    request = VcfFinalizeRequest(SvCaller.PBSV, "S1", reference, source, tmp_path / "S1.pbsv.sv.vcf.gz", "2.11", (("pbsv", "call"),))
    runner = FakeRunner(request)
    wrapper = BgzipTabixWrapper(runner=runner)
    result = wrapper.run(request)
    assert result.status is VcfFinalizeStatus.COMPLETED
    assert [command.args[0] for command in result.commands] == ["bgzip", "tabix"]
    assert result.artifact is not None and result.artifact.commands[-2][0] == "bgzip"


def test_finalizer_dry_run_needs_no_tools_and_writes_nothing(tmp_path, reference):
    source = tmp_path / "S1.cutesv.raw.vcf"
    source.write_text(vcf_text(), encoding="utf-8")
    request = VcfFinalizeRequest(SvCaller.CUTESV, "S1", reference, source, tmp_path / "S1.cutesv.sv.vcf.gz")
    result = BgzipTabixWrapper(runner=FakeRunner(request, missing=True)).run(request, dry_run=True)
    assert result.status is VcfFinalizeStatus.PLANNED
    assert not request.output_vcf.exists()


def test_finalizer_refuses_overwrite(tmp_path, reference):
    source = tmp_path / "S1.pbsv.raw.vcf"
    source.write_text(vcf_text(), encoding="utf-8")
    output = tmp_path / "S1.pbsv.sv.vcf.gz"
    output.write_bytes(b"old")
    request = VcfFinalizeRequest(SvCaller.PBSV, "S1", reference, source, output)
    with pytest.raises(OutputValidationError):
        BgzipTabixWrapper(runner=FakeRunner(request)).run(request)


def test_output_name_is_caller_specific(tmp_path, reference):
    source = tmp_path / "S1.pbsv.raw.vcf"
    source.write_text(vcf_text(), encoding="utf-8")
    with pytest.raises(InputValidationError):
        VcfFinalizeRequest(SvCaller.PBSV, "S1", reference, source, tmp_path / "S1.cutesv.sv.vcf.gz")
