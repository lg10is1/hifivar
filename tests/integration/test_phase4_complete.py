from __future__ import annotations

import json
from pathlib import Path
import struct
import zlib

import yaml
import pytest

from hifivar.command import CommandResult
from hifivar.config import load_config
from hifivar.context import AnalysisContext
from hifivar.cutesv import CuteSvWrapper
from hifivar.exceptions import InputValidationError
from hifivar.pbsv import PbsvWrapper
from hifivar.phase4 import Phase4RunStatus, run_phase4
from hifivar.sawfish import SawfishWrapper
from hifivar.sniffles2 import Sniffles2Wrapper
from hifivar.sv import BgzipTabixWrapper


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "src" / "hifivar" / "resources" / "configs" / "default.yaml"


def write_bgzf(path: Path, payload: bytes) -> None:
    compressor = zlib.compressobj(level=6, wbits=-15)
    compressed = compressor.compress(payload) + compressor.flush()
    total_size = 18 + len(compressed) + 8
    header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", total_size - 1)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(header + compressed + struct.pack("<II", zlib.crc32(payload), len(payload)))


def vcf_payload(sample: str) -> bytes:
    return (
        "##fileformat=VCFv4.3\n"
        "##contig=<ID=chr1,length=4>\n"
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
    ).encode()


class Phase4Runner:
    def __init__(self) -> None:
        self.commands: list[tuple[str, ...]] = []

    def require_executable(self, executable: str) -> Path:
        return Path(f"/opt/{executable}")

    def run(self, command, **kwargs):
        args = tuple(str(item) for item in command)
        self.commands.append(args)
        if "--version" in args:
            versions = {
                "sawfish": "sawfish 2.2.1",
                "sniffles": "Sniffles2 Version 2.8.0",
                "pbsv": "pbsv 2.11.0",
                "cuteSV": "cuteSV 2.1.4",
                "bgzip": "bgzip (htslib) 1.21",
                "tabix": "tabix (htslib) 1.21",
            }
            return CommandResult(args, 0, versions[args[0]], "", 0.01, None, True)
        if kwargs.get("dry_run"):
            return CommandResult(args, None, None, None, 0.0, None, False)
        if args[0] == "sawfish" and args[1] == "joint-call":
            directory = Path(args[args.index("--output-dir") + 1])
            vcf = directory / "genotyped.sv.vcf.gz"
            write_bgzf(vcf, vcf_payload("S1"))
            write_bgzf(Path(f"{vcf}.tbi"), b"TBI\x01")
        elif args[0] == "sniffles":
            vcf = Path(args[args.index("--vcf") + 1])
            write_bgzf(vcf, vcf_payload("S1"))
            write_bgzf(Path(f"{vcf}.tbi"), b"TBI\x01")
        elif args[0] == "pbsv" and args[1] == "discover":
            Path(args[-1]).parent.mkdir(parents=True, exist_ok=True)
            Path(args[-1]).write_bytes(b"signatures")
        elif args[0] == "pbsv" and args[1] == "call":
            Path(args[-1]).parent.mkdir(parents=True, exist_ok=True)
            Path(args[-1]).write_bytes(vcf_payload("S1"))
        elif args[0] == "cuteSV":
            Path(args[3]).parent.mkdir(parents=True, exist_ok=True)
            Path(args[3]).write_bytes(vcf_payload("S1"))
        elif args[0] == "bgzip":
            write_bgzf(Path(kwargs["stdout_path"]), Path(args[-1]).read_bytes())
        elif args[0] == "tabix":
            write_bgzf(Path(f"{args[-1]}.tbi"), b"TBI\x01")
        elif args[0] == "sawfish" and args[1] == "discover":
            pass
        else:  # pragma: no cover
            raise AssertionError(args)
        return CommandResult(args, 0, "", "", 0.25, None, True)


def test_analysis_context_to_four_independent_sv_artifacts_and_manifest(tmp_path: Path) -> None:
    fasta = tmp_path / "参考.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    bam = tmp_path / "样本.bam"
    bam.write_bytes(b"BAM")
    Path(f"{bam}.bai").write_bytes(b"BAI")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(f"sample_id\tinput\tinput_type\nS1\t{bam}\tbam\n", encoding="utf-8")
    user = tmp_path / "phase4.yaml"
    user.write_text(yaml.safe_dump({"reference": {"fasta": str(fasta), "build": "GRCh38"}, "samples": {"sheet": str(sheet)}, "sv": {"enabled": True}}, sort_keys=False), encoding="utf-8")
    context = AnalysisContext.from_config(load_config(DEFAULT_CONFIG, user_config=user))
    runner = Phase4Runner()
    report = run_phase4(
        context,
        tmp_path / "results",
        sawfish_wrapper=SawfishWrapper(runner=runner),
        sniffles2_wrapper=Sniffles2Wrapper(runner=runner),
        pbsv_wrapper=PbsvWrapper(runner=runner),
        cutesv_wrapper=CuteSvWrapper(runner=runner),
        finalizer=BgzipTabixWrapper(runner=runner),
    )
    assert report.status is Phase4RunStatus.COMPLETED
    assert report.tool_versions == {"sawfish": "2.2.1", "sniffles2": "2.8.0", "pbsv": "2.11.0", "cutesv": "2.1.4", "bgzip": "1.21", "tabix": "1.21"}
    artifacts = report.sample_results[0].artifacts
    assert [artifact.vcf_path.name for artifact in artifacts] == [
        "S1.sawfish.sv.vcf.gz", "S1.sniffles2.sv.vcf.gz", "S1.pbsv.sv.vcf.gz", "S1.cutesv.sv.vcf.gz"
    ]
    assert all(not artifact.harmonized for artifact in artifacts)
    assert len(report.sample_results[0].caller_results) == 4
    assert len(report.sample_results[0].finalization_results) == 2
    command_text = "\n".join(" ".join(command) for command in runner.commands)
    assert "pbsv discover" in command_text and "pbsv call" in command_text
    assert "|" not in command_text
    json_path = report.write_json(tmp_path / "reports" / "phase4.json")
    yaml_path = report.write_yaml(tmp_path / "reports" / "phase4.yaml")
    assert json.loads(json_path.read_text(encoding="utf-8")) == yaml.safe_load(yaml_path.read_text(encoding="utf-8"))


def test_phase4_dry_run_plans_four_callers_and_finalization_without_tools_or_writes(tmp_path: Path) -> None:
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    bam = tmp_path / "S1.bam"
    bam.write_bytes(b"BAM")
    Path(f"{bam}.bai").write_bytes(b"BAI")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(f"sample_id\tinput\tinput_type\nS1\t{bam}\tbam\n", encoding="utf-8")
    user = tmp_path / "config.yaml"
    user.write_text(yaml.safe_dump({"reference": {"fasta": str(fasta), "build": "GRCh38"}, "samples": {"sheet": str(sheet)}, "sv": {"enabled": True}}), encoding="utf-8")
    context = AnalysisContext.from_config(load_config(DEFAULT_CONFIG, user_config=user))
    runner = Phase4Runner()
    output = tmp_path / "planned-results"
    report = run_phase4(
        context,
        output,
        sawfish_wrapper=SawfishWrapper(runner=runner),
        sniffles2_wrapper=Sniffles2Wrapper(runner=runner),
        pbsv_wrapper=PbsvWrapper(runner=runner),
        cutesv_wrapper=CuteSvWrapper(runner=runner),
        finalizer=BgzipTabixWrapper(runner=runner),
        dry_run=True,
    )
    assert report.status is Phase4RunStatus.PLANNED
    assert len(report.sample_results[0].caller_results) == 4
    assert len(report.sample_results[0].finalization_results) == 2
    assert report.sample_results[0].artifacts == ()
    assert not output.exists()


def test_cram_with_bam_only_callers_fails_before_any_external_command(tmp_path: Path) -> None:
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    cram = tmp_path / "S1.cram"
    cram.write_bytes(b"CRAM")
    Path(f"{cram}.crai").write_bytes(b"CRAI")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(f"sample_id\tinput\tinput_type\nS1\t{cram}\tcram\n", encoding="utf-8")
    user = tmp_path / "config.yaml"
    user.write_text(yaml.safe_dump({"reference": {"fasta": str(fasta), "build": "GRCh38"}, "samples": {"sheet": str(sheet)}, "sv": {"enabled": True}}), encoding="utf-8")
    context = AnalysisContext.from_config(load_config(DEFAULT_CONFIG, user_config=user))
    runner = Phase4Runner()
    with pytest.raises(InputValidationError, match="CRAM"):
        run_phase4(
            context,
            tmp_path / "results",
            sawfish_wrapper=SawfishWrapper(runner=runner),
            sniffles2_wrapper=Sniffles2Wrapper(runner=runner),
            pbsv_wrapper=PbsvWrapper(runner=runner),
            cutesv_wrapper=CuteSvWrapper(runner=runner),
            finalizer=BgzipTabixWrapper(runner=runner),
        )
    assert runner.commands == []
