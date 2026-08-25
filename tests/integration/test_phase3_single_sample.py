"""Tiny Phase 3.2 single-sample execution integration without DeepVariant."""

from __future__ import annotations

import struct
import zlib
from pathlib import Path

from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import (
    AlignmentArtifact,
    AlignmentSortOrder,
    AlignmentSource,
)
from hifivar.command import CommandResult
from hifivar.deepvariant import DeepVariantWrapper
from hifivar.reference import ReferenceGenome
from hifivar.small import DeepVariantRequest, SmallVariantResultStatus


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
        header
        + compressed
        + struct.pack("<II", zlib.crc32(payload), len(payload))
    )


class TinyDeepVariantRunner:
    def require_executable(self, executable: str) -> Path:
        return Path("/opt/deepvariant/bin/run_deepvariant")

    def run(self, command: object, **kwargs: object) -> CommandResult:
        args = tuple(str(value) for value in command)  # type: ignore[union-attr]
        if args[-1] == "--version":
            return CommandResult(args, 0, "DeepVariant version 1.10.0", "", 0.1, None, True)
        if kwargs.get("dry_run") is True:
            return CommandResult(args, None, None, None, 0.0, None, False)
        values = {
            argument.split("=", 1)[0]: argument.split("=", 1)[1]
            for argument in args
            if argument.startswith(("--output_vcf=", "--output_gvcf=", "--sample_name="))
        }
        sample = values["--sample_name"]
        header = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=4>\n"
        columns = f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
        vcf = Path(values["--output_vcf"])
        gvcf = Path(values["--output_gvcf"])
        _write_bgzf(vcf, (header + columns).encode())
        _write_bgzf(
            gvcf,
            (
                header
                + "##DeepVariant_version=1.10.0\n"
                + "##FILTER=<ID=RefCall,Description=\"Reference call\">\n"
                + "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP\">\n"
                + columns
            ).encode(),
        )
        _write_bgzf(Path(f"{vcf}.tbi"), b"TBI\x01")
        _write_bgzf(Path(f"{gvcf}.tbi"), b"TBI\x01")
        return CommandResult(args, 0, "", "", 4.0, None, True)


def test_indexed_alignment_to_validated_small_variant_artifact(tmp_path: Path) -> None:
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    bam = tmp_path / "S1.bam"
    bam.write_bytes(b"coordinate-sorted-bam")
    bai = Path(f"{bam}.bai")
    bai.write_bytes(b"bai")
    artifact = AlignmentArtifact(
        sample_id="S1",
        path=bam,
        output_format=AlignmentOutputFormat.BAM,
        reference=ReferenceGenome.from_fasta(fasta, build="GRCh38"),
        source=AlignmentSource.EXISTING,
        sort_order=AlignmentSortOrder.UNKNOWN,
        index_path=bai,
    )
    request = DeepVariantRequest.create(artifact, tmp_path / "results" / "small")

    result = DeepVariantWrapper(runner=TinyDeepVariantRunner()).run(request)  # type: ignore[arg-type]

    assert result.status is SmallVariantResultStatus.COMPLETED
    assert result.artifact is not None
    assert result.artifact.sample_id == "S1"
    assert result.artifact.vcf_path.name == "S1.small.vcf.gz"
    assert result.artifact.gvcf_path.name == "S1.g.vcf.gz"
    assert result.artifact.vcf_index_path.is_file()
    assert result.artifact.gvcf_index_path.is_file()
