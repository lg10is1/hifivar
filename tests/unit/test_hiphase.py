from __future__ import annotations

import struct
import zlib
from pathlib import Path

import pytest

from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import (
    AlignmentArtifact,
    AlignmentSortOrder,
    AlignmentSource,
)
from hifivar.command import CommandResult
from hifivar.exceptions import (
    CommandExecutionError,
    InputValidationError,
    OutputValidationError,
    ToolNotFoundError,
    ToolVersionError,
)
from hifivar.hiphase import HiPhaseWrapper, PhasingResultStatus
from hifivar.phasing import PhasingRequest, PhasingResources
from hifivar.reference import ReferenceGenome
from hifivar.small import SmallVariantArtifact


def write_bgzf(path: Path, payload: bytes) -> None:
    compressor = zlib.compressobj(level=6, wbits=-15)
    compressed = compressor.compress(payload) + compressor.flush()
    total_size = 18 + len(compressed) + 8
    header = (
        b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
        + struct.pack("<H", total_size - 1)
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(
        header + compressed + struct.pack("<II", zlib.crc32(payload), len(payload))
    )


def make_request(tmp_path: Path, *, output_sample: str = "S1") -> PhasingRequest:
    tmp_path.mkdir(parents=True, exist_ok=True)
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t8\t6\t8\t9\n", encoding="utf-8")
    bam = tmp_path / "S1.bam"
    bam.write_bytes(b"BAM")
    bai = Path(f"{bam}.bai")
    bai.write_bytes(b"BAI")
    reference = ReferenceGenome.from_fasta(fasta, build="GRCh38")
    alignment = AlignmentArtifact(
        "S1",
        bam,
        AlignmentOutputFormat.BAM,
        reference,
        AlignmentSource.EXISTING,
        AlignmentSortOrder.COORDINATE,
        bai,
    )
    small = tmp_path / "S1.small.vcf.gz"
    small.write_bytes(b"small")
    small_index = Path(f"{small}.tbi")
    small_index.write_bytes(b"TBI")
    gvcf = tmp_path / "S1.g.vcf.gz"
    gvcf.write_bytes(b"gvcf")
    gvcf_index = Path(f"{gvcf}.tbi")
    gvcf_index.write_bytes(b"TBI")
    artifact = SmallVariantArtifact(
        "S1",
        "GRCh38",
        small,
        gvcf,
        small_index,
        gvcf_index,
        tool_version="1.10.0",
    )
    return PhasingRequest(
        alignment,
        artifact,
        tmp_path / "phasing" / f"{output_sample}.phased.vcf.gz",
        PhasingResources(threads=6),
    )


class FakeRunner:
    def __init__(
        self,
        *,
        missing: bool = False,
        fail: bool = False,
        materialize: bool = False,
        bad_version: bool = False,
    ) -> None:
        self.missing = missing
        self.fail = fail
        self.materialize = materialize
        self.bad_version = bad_version
        self.calls: list[tuple[tuple[str, ...], dict[str, object]]] = []

    def require_executable(self, executable: str) -> Path:
        if self.missing:
            raise ToolNotFoundError(f"missing {executable}")
        return Path(executable)

    def run(self, command, **kwargs):
        args = tuple(str(value) for value in command)
        self.calls.append((args, kwargs))
        if "--version" in args:
            output = "unknown" if self.bad_version else "HiPhase v1.7.0"
            return CommandResult(args, 0, output, "", 0.01, None, True)
        if kwargs.get("dry_run"):
            return CommandResult(args, None, None, None, 0.0, None, False)
        if self.fail:
            raise CommandExecutionError("synthetic HiPhase failure")
        if self.materialize and args[0] == "hiphase":
            output = Path(args[args.index("--output-vcf") + 1])
            header = (
                "##fileformat=VCFv4.2\n"
                "##contig=<ID=chr1,length=8>\n"
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            ).encode()
            write_bgzf(output, header)
        if self.materialize and args[0] == "tabix":
            write_bgzf(Path(f"{args[-1]}.tbi"), b"TBI\x01")
        return CommandResult(args, 0, "", "", 0.2, None, True)


def test_phasing_contract_and_deterministic_official_command(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    commands = HiPhaseWrapper(executable="/opt/hiphase").plan_commands(request)
    args = commands[0].args
    assert args[0] == "/opt/hiphase"
    assert args[args.index("--bam") + 1] == str(request.alignment.path.absolute())
    assert args[args.index("--vcf") + 1] == str(request.small_variants.vcf_path.absolute())
    assert args[args.index("--reference") + 1] == str(request.alignment.reference.fasta.absolute())
    assert args[args.index("--sample-name") + 1] == "S1"
    assert args[args.index("--threads") + 1] == "6"
    assert "--disable-global-realignment" in args
    assert commands[1].args[:4] == ("tabix", "-f", "-p", "vcf")
    assert all(item.to_dict()["shell"] is False for item in commands)


def test_contract_rejects_wrong_output_name_and_missing_index(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="must follow"):
        make_request(tmp_path, output_sample="wrong")
    request = make_request(tmp_path / "second")
    object.__setattr__(request.alignment, "index_path", None)
    with pytest.raises(InputValidationError, match="indexed BAM"):
        PhasingRequest(
            request.alignment,
            request.small_variants,
            tmp_path / "second" / "S1.phased.vcf.gz",
        )


def test_version_missing_and_parse_failure(tmp_path: Path) -> None:
    with pytest.raises(ToolNotFoundError):
        HiPhaseWrapper(runner=FakeRunner(missing=True)).detect_version()  # type: ignore[arg-type]
    with pytest.raises(ToolVersionError):
        HiPhaseWrapper(runner=FakeRunner(bad_version=True)).detect_version()  # type: ignore[arg-type]


def test_dry_run_needs_no_installed_tool_and_writes_nothing(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    runner = FakeRunner(missing=True)
    result = HiPhaseWrapper(runner=runner).run(request, dry_run=True)  # type: ignore[arg-type]
    assert result.status is PhasingResultStatus.PLANNED
    assert len(runner.calls) == 2
    assert not request.output_vcf.parent.exists()


def test_fake_execution_validates_and_records_provenance(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    result = HiPhaseWrapper(runner=FakeRunner(materialize=True)).run(request)  # type: ignore[arg-type]
    assert result.status is PhasingResultStatus.COMPLETED
    assert result.hiphase_version == "1.7.0"
    assert result.artifact is not None
    assert result.artifact.source_bam == request.alignment.path
    assert result.artifact.source_small_vcf == request.small_variants.vcf_path
    assert result.artifact.execution_backend == "native"


def test_failure_and_no_silent_overwrite(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    with pytest.raises(CommandExecutionError):
        HiPhaseWrapper(runner=FakeRunner(fail=True)).run(request)  # type: ignore[arg-type]
    request.output_vcf.parent.mkdir(parents=True, exist_ok=True)
    request.output_vcf.write_bytes(b"old")
    with pytest.raises(OutputValidationError, match="already exists"):
        make_request(tmp_path)


def test_output_validation_requires_phase_set_header(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    runner = FakeRunner(materialize=True)
    original = runner.run

    def without_phase_set(command, **kwargs):
        result = original(command, **kwargs)
        args = tuple(str(value) for value in command)
        if args[0] == "hiphase" and "--version" not in args:
            output = Path(args[args.index("--output-vcf") + 1])
            write_bgzf(
                output,
                (
                    "##fileformat=VCFv4.2\n"
                    "##contig=<ID=chr1,length=8>\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                ).encode(),
            )
        return result

    runner.run = without_phase_set  # type: ignore[method-assign]
    with pytest.raises(OutputValidationError, match="FORMAT/PS"):
        HiPhaseWrapper(runner=runner).run(request)  # type: ignore[arg-type]
