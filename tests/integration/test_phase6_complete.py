from __future__ import annotations

import json
import struct
import zlib
from pathlib import Path

import yaml

from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import (
    AlignmentArtifact,
    AlignmentSortOrder,
    AlignmentSource,
)
from hifivar.command import CommandResult
from hifivar.config import load_config
from hifivar.context import AnalysisContext
from hifivar.hiphase import HiPhaseWrapper
from hifivar.phase6 import Phase6RunStatus, run_phase6
from hifivar.small import SmallVariantArtifact


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "src" / "hifivar" / "resources" / "configs" / "default.yaml"


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


class Phase6Runner:
    def require_executable(self, executable: str) -> Path:
        return Path(executable)

    def run(self, command, **kwargs):
        args = tuple(str(value) for value in command)
        if "--version" in args:
            return CommandResult(args, 0, "HiPhase v1.7.0", "", 0.01, None, True)
        if args[0] == "hiphase":
            output = Path(args[args.index("--output-vcf") + 1])
            write_bgzf(
                output,
                (
                    "##fileformat=VCFv4.2\n"
                    "##contig=<ID=chr1,length=8>\n"
                    '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                ).encode(),
            )
        elif args[0] == "tabix":
            write_bgzf(Path(f"{args[-1]}.tbi"), b"TBI\x01")
        return CommandResult(args, 0, "", "", 0.2, None, True)


def test_context_to_hiphase_artifact_and_provenance(tmp_path: Path) -> None:
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t8\t6\t8\t9\n", encoding="utf-8")
    bam = tmp_path / "S1.bam"
    bam.write_bytes(b"BAM")
    bai = Path(f"{bam}.bai")
    bai.write_bytes(b"BAI")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        f"sample_id\tinput\tinput_type\nS1\t{bam}\tbam\n",
        encoding="utf-8",
    )
    user = tmp_path / "phase6.yaml"
    user.write_text(
        yaml.safe_dump(
            {
                "reference": {"fasta": str(fasta), "build": "GRCh38"},
                "samples": {"sheet": str(sheet)},
                "phasing": {"enabled": True},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    config = load_config(DEFAULT_CONFIG, user_config=user)
    context = AnalysisContext.from_config(config)
    reference = context.reference
    alignment = AlignmentArtifact(
        "S1",
        bam,
        AlignmentOutputFormat.BAM,
        reference,
        AlignmentSource.EXISTING,
        AlignmentSortOrder.COORDINATE,
        bai,
    )
    vcf = tmp_path / "S1.small.vcf.gz"
    gvcf = tmp_path / "S1.g.vcf.gz"
    for path in (vcf, gvcf):
        path.write_bytes(b"VCF")
        Path(f"{path}.tbi").write_bytes(b"TBI")
    small = SmallVariantArtifact(
        "S1",
        "GRCh38",
        vcf,
        gvcf,
        Path(f"{vcf}.tbi"),
        Path(f"{gvcf}.tbi"),
        tool_version="1.10.0",
    )
    report = run_phase6(
        context,
        alignment_artifacts={"S1": alignment},
        small_variant_artifacts={"S1": small},
        output_directory=tmp_path / "results" / "phasing",
        config=config,
        wrapper=HiPhaseWrapper(runner=Phase6Runner()),  # type: ignore[arg-type]
    )
    assert report.status is Phase6RunStatus.COMPLETED
    artifact = report.sample_results[0].phasing.artifact
    assert artifact is not None
    assert artifact.vcf_path.name == "S1.phased.vcf.gz"
    assert artifact.hiphase_version == "1.7.0"
    json_path = report.write_json(tmp_path / "phase6.json")
    yaml_path = report.write_yaml(tmp_path / "phase6-report.yaml")
    assert json.loads(json_path.read_text(encoding="utf-8")) == yaml.safe_load(
        yaml_path.read_text(encoding="utf-8")
    )
