from __future__ import annotations

import struct
import zlib
from pathlib import Path

import pytest

from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import AlignmentArtifact, AlignmentSortOrder, AlignmentSource
from hifivar.command import CommandResult
from hifivar.exceptions import CommandExecutionError, OutputValidationError, ToolNotFoundError, ToolVersionError
from hifivar.reference import ReferenceGenome
from hifivar.tr import TandemRepeatCatalog
from hifivar.trgt import TrgtPreset, TrgtRequest, TrgtResources, TrgtResultStatus, TrgtWrapper


def write_bgzf(path: Path, payload: bytes) -> None:
    compressor = zlib.compressobj(level=6, wbits=-15)
    compressed = compressor.compress(payload) + compressor.flush()
    total_size = 18 + len(compressed) + 8
    header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", total_size - 1)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(header + compressed + struct.pack("<II", zlib.crc32(payload), len(payload)))


def trgt_header(sample: str) -> bytes:
    return (
        "##fileformat=VCFv4.3\n"
        "##contig=<ID=chr1,length=8>\n"
        '##INFO=<ID=TRID,Number=1,Type=String,Description="ID">\n'
        '##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motifs">\n'
        '##INFO=<ID=STRUC,Number=1,Type=String,Description="Structure">\n'
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
    ).encode()


def make_request(tmp_path: Path) -> TrgtRequest:
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t8\t6\t8\t9\n", encoding="utf-8")
    bam = tmp_path / "S1.bam"
    bam.write_bytes(b"BAM")
    bai = Path(f"{bam}.bai")
    bai.write_bytes(b"BAI")
    catalog = tmp_path / "catalog.bed"
    catalog.write_text("chr1\t1\t7\tID=TR1;MOTIFS=CG;STRUC=<TR>\n", encoding="utf-8")
    artifact = AlignmentArtifact(
        "S1", bam, AlignmentOutputFormat.BAM, ReferenceGenome.from_fasta(fasta, build="GRCh38"),
        AlignmentSource.EXISTING, AlignmentSortOrder.UNKNOWN, bai,
    )
    return TrgtRequest(
        artifact,
        TandemRepeatCatalog(catalog, "GRCh38"),
        tmp_path / "work" / "S1.trgt",
        tmp_path / "results" / "S1.tr.vcf.gz",
        tmp_path / "results" / "S1.tr.spanning.bam",
        "XY",
        TrgtResources(threads=6),
        TrgtPreset.WGS,
    )


class FakeRunner:
    def __init__(self, *, missing: bool = False, fail: bool = False, materialize: bool = False) -> None:
        self.missing, self.fail, self.materialize = missing, fail, materialize
        self.calls: list[tuple[tuple[str, ...], dict[str, object]]] = []

    def require_executable(self, executable: str) -> Path:
        if self.missing:
            raise ToolNotFoundError(f"missing {executable}")
        return Path(executable)

    def run(self, command, **kwargs):
        args = tuple(str(value) for value in command)
        self.calls.append((args, kwargs))
        if "--version" in args:
            versions = {"trgt": "trgt 5.1.0", "bcftools": "bcftools 1.21", "samtools": "samtools 1.21"}
            return CommandResult(args, 0, versions[args[0]], "", 0.01, None, True)
        if kwargs.get("dry_run"):
            return CommandResult(args, None, None, None, 0.0, None, False)
        if self.fail:
            raise CommandExecutionError("synthetic TRGT failure")
        if self.materialize:
            if args[0] == "trgt":
                prefix = Path(args[args.index("--output-prefix") + 1])
                prefix.parent.mkdir(parents=True, exist_ok=True)
                Path(f"{prefix}.vcf.gz").write_bytes(b"raw-vcf")
                Path(f"{prefix}.spanning.bam").write_bytes(b"raw-bam")
            elif args[:2] == ("bcftools", "sort"):
                write_bgzf(Path(args[args.index("-o") + 1]), trgt_header("S1"))
            elif args[:2] == ("bcftools", "index"):
                write_bgzf(Path(f"{args[-1]}.tbi"), b"TBI\x01")
            elif args[:2] == ("samtools", "sort"):
                output = Path(args[args.index("-o") + 1])
                output.parent.mkdir(parents=True, exist_ok=True)
                output.write_bytes(b"sorted-bam")
            elif args[:2] == ("samtools", "index"):
                Path(args[args.index("-o") + 1]).write_bytes(b"BAI")
        return CommandResult(args, 0, "", "", 0.2, None, True)


def test_command_plan_matches_official_genotype_and_sort_contract(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    commands = TrgtWrapper(executable="/opt/trgt").plan_commands(request)
    genotype = commands[0].args
    assert genotype[:2] == ("/opt/trgt", "genotype")
    assert genotype[genotype.index("--sample-name") + 1] == "S1"
    assert genotype[genotype.index("--karyotype") + 1] == "XY"
    assert genotype[genotype.index("--preset") + 1] == "wgs"
    assert commands[1].args[:2] == ("bcftools", "sort")
    assert "--threads" not in commands[1].args
    assert commands[2].args[:4] == ("bcftools", "index", "--tbi", "--threads")
    assert commands[3].args[:2] == ("samtools", "sort")
    assert "-@" in commands[3].args and "-@" in commands[4].args
    assert all("|" not in command.display_command for command in commands)


def test_versions_dry_run_and_missing_tools(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    runner = FakeRunner(missing=True)
    result = TrgtWrapper(runner=runner).run(request, dry_run=True)  # type: ignore[arg-type]
    assert result.status is TrgtResultStatus.PLANNED
    assert len(runner.calls) == 5
    assert not request.final_vcf.parent.exists()
    with pytest.raises(ToolNotFoundError):
        TrgtWrapper(runner=runner).detect_versions()  # type: ignore[arg-type]


def test_version_parse_failure(tmp_path: Path) -> None:
    runner = FakeRunner()
    original = runner.run
    runner.run = lambda command, **kwargs: CommandResult(tuple(command), 0, "unknown", "", 0.0, None, True)  # type: ignore[method-assign]
    with pytest.raises(ToolVersionError):
        TrgtWrapper(runner=runner).detect_versions()  # type: ignore[arg-type]
    runner.run = original  # type: ignore[method-assign]


def test_real_fake_execution_validates_independent_tr_outputs(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    result = TrgtWrapper(runner=FakeRunner(materialize=True)).run(request)  # type: ignore[arg-type]
    assert result.status is TrgtResultStatus.COMPLETED
    assert result.tool_versions == {"trgt": "5.1.0", "bcftools": "1.21", "samtools": "1.21"}
    assert result.artifact is not None
    assert result.artifact.vcf_path.name == "S1.tr.vcf.gz"
    assert result.artifact.spanning_bam_path.name == "S1.tr.spanning.bam"


def test_failure_and_no_silent_overwrite(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    with pytest.raises(CommandExecutionError):
        TrgtWrapper(runner=FakeRunner(fail=True)).run(request)  # type: ignore[arg-type]
    request.final_vcf.parent.mkdir(parents=True, exist_ok=True)
    request.final_vcf.write_bytes(b"old")
    with pytest.raises(OutputValidationError, match="already exists"):
        TrgtWrapper(runner=FakeRunner(materialize=True)).run(request)  # type: ignore[arg-type]


class LocusErrorRunner(FakeRunner):
    def run(self, command, **kwargs):
        result = super().run(command, **kwargs)
        args = tuple(str(value) for value in command)
        if args[0] != "trgt" or "--version" in args or kwargs.get("dry_run"):
            return result
        message = (
            "2026-08-21 12:11:24 [ERROR] - Locus processing: "
            "Error at BED line 2\n"
        )
        stderr_path = kwargs.get("stderr_path")
        if stderr_path is not None:
            path = Path(stderr_path)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(message, encoding="utf-8")
            return CommandResult(args, 0, "", None, 0.2, None, True)
        return CommandResult(args, 0, "", message, 0.2, None, True)


@pytest.mark.parametrize("persist_log", (False, True))
def test_zero_exit_locus_processing_error_rejects_partial_outputs(
    tmp_path: Path,
    persist_log: bool,
) -> None:
    request = make_request(tmp_path)
    runner = LocusErrorRunner(materialize=True)
    stderr_path = tmp_path / "logs" / "S1.trgt.log" if persist_log else None

    with pytest.raises(OutputValidationError, match="refusing incomplete"):
        TrgtWrapper(runner=runner).run(request, stderr_path=stderr_path)  # type: ignore[arg-type]

    biological_commands = [args for args, _ in runner.calls if "--version" not in args]
    assert len(biological_commands) == 1
    assert biological_commands[0][:2] == ("trgt", "genotype")
