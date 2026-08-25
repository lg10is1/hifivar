"""Tiny FASTQ-to-small-variant mock integration through Phases 2 and 3."""

from __future__ import annotations

import json
import struct
import zlib
from pathlib import Path

import yaml

from hifivar.command import CommandResult
from hifivar.config import load_config
from hifivar.context import AnalysisContext
from hifivar.deepvariant import DeepVariantWrapper
from hifivar.pbmm2 import Pbmm2Wrapper
from hifivar.phase2 import run_phase2
from hifivar.phase3 import (
    Phase3RunStatus,
    collect_phase2_alignment_artifacts,
    run_phase3,
)
from hifivar.samtools import SamtoolsWrapper


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "src" / "hifivar" / "resources" / "configs" / "default.yaml"


def _write_bgzf(path: Path, payload: bytes) -> None:
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


class Phase2And3Runner:
    def __init__(self) -> None:
        self.commands: list[tuple[str, ...]] = []

    def require_executable(self, executable: str) -> Path:
        return Path(f"/opt/{executable}")

    def run(self, command: object, **kwargs: object) -> CommandResult:
        args = tuple(str(value) for value in command)  # type: ignore[union-attr]
        self.commands.append(args)
        if args[-1] == "--version":
            versions = {
                "pbmm2": "pbmm2 1.17.0",
                "samtools": "samtools 1.22.1",
                "run_deepvariant": "DeepVariant version 1.10.0",
            }
            return CommandResult(args, 0, versions[args[0]], "", 0.01, None, True)
        if args[0] == "pbmm2":
            Path(args[4]).write_bytes(b"coordinate-sorted-bam")
        elif args[0] == "samtools":
            Path(args[-1]).write_bytes(b"bai")
        elif args[0] == "run_deepvariant":
            values = {
                item.split("=", 1)[0]: item.split("=", 1)[1]
                for item in args
                if item.startswith(("--output_vcf=", "--output_gvcf=", "--sample_name="))
            }
            sample = values["--sample_name"]
            base = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=4>\n"
            columns = f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
            vcf = Path(values["--output_vcf"])
            gvcf = Path(values["--output_gvcf"])
            _write_bgzf(vcf, (base + columns).encode())
            _write_bgzf(
                gvcf,
                (
                    base
                    + "##DeepVariant_version=1.10.0\n"
                    + "##FILTER=<ID=RefCall,Description=\"Reference call\">\n"
                    + "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP\">\n"
                    + columns
                ).encode(),
            )
            _write_bgzf(Path(f"{vcf}.tbi"), b"TBI\x01")
            _write_bgzf(Path(f"{gvcf}.tbi"), b"TBI\x01")
        else:  # pragma: no cover
            raise AssertionError(f"Unexpected tool: {args[0]}")
        return CommandResult(args, 0, "", "", 1.0, None, True)


def test_analysis_context_through_alignment_and_deepvariant_provenance(
    tmp_path: Path,
) -> None:
    data = tmp_path / "data"
    data.mkdir()
    fasta = data / "reference.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    fastq = data / "S1.fastq"
    fastq.write_text("@r\nACGT\n+\nIIII\n", encoding="utf-8")
    sheet = data / "samples.tsv"
    sheet.write_text("sample_id\tinput\nS1\tS1.fastq\n", encoding="utf-8")
    user = data / "analysis.yaml"
    user.write_text(
        yaml.safe_dump(
            {
                "reference": {"fasta": "reference.fa", "build": "GRCh38"},
                "samples": {"sheet": "samples.tsv"},
                "alignment": {"threads": 4, "index_threads": 2},
                "small": {"enabled": True, "threads": 6},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    context = AnalysisContext.from_config(load_config(DEFAULT_CONFIG, user_config=user))
    runner = Phase2And3Runner()
    phase2 = run_phase2(
        context,
        tmp_path / "results" / "alignment",
        pbmm2_wrapper=Pbmm2Wrapper(runner=runner),  # type: ignore[arg-type]
        samtools_wrapper=SamtoolsWrapper(runner=runner),  # type: ignore[arg-type]
    )
    phase3 = run_phase3(
        context,
        tmp_path / "results",
        alignment_artifacts=collect_phase2_alignment_artifacts(phase2),
        deepvariant_wrapper=DeepVariantWrapper(runner=runner),  # type: ignore[arg-type]
    )

    assert phase3.status is Phase3RunStatus.COMPLETED
    assert phase3.tool_versions == {"deepvariant": "1.10.0"}
    result = phase3.sample_results[0]
    assert result.alignment.path.name == "S1.aligned.bam"
    assert result.call.artifact is not None
    assert result.call.artifact.vcf_path.name == "S1.small.vcf.gz"
    assert result.call.artifact.gvcf_path.name == "S1.g.vcf.gz"
    executables = [command[0] for command in runner.commands]
    assert executables.count("pbmm2") == 2
    assert executables.count("samtools") == 2
    assert executables.count("run_deepvariant") == 2
    forbidden = ("sawfish", "sniffles", "pbsv", "cutesv", "trgt")
    command_text = "\n".join(" ".join(command) for command in runner.commands).lower()
    assert all(word not in command_text for word in forbidden)

    json_path = phase3.write_json(tmp_path / "reports" / "phase3.json")
    yaml_path = phase3.write_yaml(tmp_path / "reports" / "phase3.yaml")
    assert json.loads(json_path.read_text(encoding="utf-8")) == yaml.safe_load(
        yaml_path.read_text(encoding="utf-8")
    )
