from __future__ import annotations

import json
import struct
import zlib
from pathlib import Path

import yaml

from hifivar.command import CommandResult
from hifivar.config import load_config
from hifivar.context import AnalysisContext
from hifivar.phase5 import Phase5RunStatus, run_phase5
from hifivar.trgt import TrgtWrapper


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "src" / "hifivar" / "resources" / "configs" / "default.yaml"


def write_bgzf(path: Path, payload: bytes) -> None:
    compressor = zlib.compressobj(level=6, wbits=-15)
    compressed = compressor.compress(payload) + compressor.flush()
    total_size = 18 + len(compressed) + 8
    header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", total_size - 1)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(header + compressed + struct.pack("<II", zlib.crc32(payload), len(payload)))


class Phase5Runner:
    def require_executable(self, executable: str) -> Path:
        return Path(executable)

    def run(self, command, **kwargs):
        args = tuple(str(value) for value in command)
        if "--version" in args:
            versions = {"trgt": "trgt 5.1.0", "bcftools": "bcftools 1.21", "samtools": "samtools 1.21"}
            return CommandResult(args, 0, versions[args[0]], "", 0.01, None, True)
        if args[0] == "trgt":
            prefix = Path(args[args.index("--output-prefix") + 1])
            prefix.parent.mkdir(parents=True, exist_ok=True)
            Path(f"{prefix}.vcf.gz").write_bytes(b"raw")
            Path(f"{prefix}.spanning.bam").write_bytes(b"raw-bam")
        elif args[:2] == ("bcftools", "sort"):
            payload = (
                "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=8>\n"
                "##INFO=<ID=TRID,Number=1,Type=String,Description=\"ID\">\n"
                "##INFO=<ID=MOTIFS,Number=.,Type=String,Description=\"Motifs\">\n"
                "##INFO=<ID=STRUC,Number=1,Type=String,Description=\"Structure\">\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            ).encode()
            write_bgzf(Path(args[args.index("-o") + 1]), payload)
        elif args[:2] == ("bcftools", "index"):
            write_bgzf(Path(f"{args[-1]}.tbi"), b"TBI\x01")
        elif args[:2] == ("samtools", "sort"):
            output = Path(args[args.index("-o") + 1])
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_bytes(b"sorted-bam")
        elif args[:2] == ("samtools", "index"):
            Path(args[args.index("-o") + 1]).write_bytes(b"BAI")
        return CommandResult(args, 0, "", "", 0.2, None, True)


def test_context_to_trgt_artifact_and_provenance(tmp_path: Path) -> None:
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t8\t6\t8\t9\n", encoding="utf-8")
    bam = tmp_path / "S1.bam"
    bam.write_bytes(b"BAM")
    Path(f"{bam}.bai").write_bytes(b"BAI")
    catalog = tmp_path / "catalog.bed"
    catalog.write_text("chr1\t1\t7\tID=TR1;MOTIFS=CG;STRUC=<TR>\n", encoding="utf-8")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(f"sample_id\tinput\tinput_type\tsex\nS1\t{bam}\tbam\tmale\n", encoding="utf-8")
    user = tmp_path / "phase5.yaml"
    user.write_text(yaml.safe_dump({
        "reference": {"fasta": str(fasta), "build": "GRCh38"},
        "samples": {"sheet": str(sheet)},
        "tr": {"enabled": True, "catalog": str(catalog), "catalog_reference_build": "GRCh38"},
    }, sort_keys=False), encoding="utf-8")
    context = AnalysisContext.from_config(load_config(DEFAULT_CONFIG, user_config=user))
    report = run_phase5(
        context,
        tmp_path / "results",
        wrapper=TrgtWrapper(runner=Phase5Runner()),  # type: ignore[arg-type]
    )
    assert report.status is Phase5RunStatus.COMPLETED
    assert report.sample_results[0].karyotype == "XY"
    assert report.tool_versions == {"trgt": "5.1.0", "bcftools": "1.21", "samtools": "1.21"}
    artifact = report.sample_results[0].trgt_result.artifact
    assert artifact is not None and artifact.vcf_path.name == "S1.tr.vcf.gz"
    json_path = report.write_json(tmp_path / "report.json")
    yaml_path = report.write_yaml(tmp_path / "report.yaml")
    assert json.loads(json_path.read_text(encoding="utf-8")) == yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
