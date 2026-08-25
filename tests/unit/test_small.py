"""Tests for Phase 3 small-variant request and result models."""

from __future__ import annotations

from pathlib import Path
import struct
import zlib

import pytest

import hifivar.small as small_module

from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import (
    AlignmentArtifact,
    AlignmentSortOrder,
    AlignmentSource,
)
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.reference import ReferenceGenome
from hifivar.small import (
    DeepVariantRequest,
    SmallVariantCommandPlan,
    SmallVariantResources,
    SmallVariantResult,
    SmallVariantResultStatus,
    validate_small_variant_outputs,
)


def write_bgzf(path: Path, payload: bytes) -> None:
    compressor = zlib.compressobj(level=6, wbits=-15)
    compressed = compressor.compress(payload) + compressor.flush()
    total_size = 18 + len(compressed) + 8
    header = (
        b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
        + struct.pack("<H", total_size - 1)
    )
    footer = struct.pack("<II", zlib.crc32(payload), len(payload))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(header + compressed + footer)


def write_outputs(request: DeepVariantRequest, *, sample: str = "S1", contig: str = "chr1") -> None:
    base = f"##fileformat=VCFv4.2\n##contig=<ID={contig},length=4>\n"
    columns = f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
    gvcf_markers = (
        "##DeepVariant_version=1.10.0\n"
        "##FILTER=<ID=RefCall,Description=\"Reference call\">\n"
        "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP\">\n"
    )
    write_bgzf(request.output_vcf, (base + columns).encode())
    write_bgzf(
        request.output_gvcf,
        (base + gvcf_markers + columns).encode(),
    )
    write_bgzf(request.output_vcf_index, b"TBI\x01")
    write_bgzf(request.output_gvcf_index, b"TBI\x01")


def make_artifact(tmp_path: Path, *, cram: bool = False) -> AlignmentArtifact:
    tmp_path.mkdir(parents=True, exist_ok=True)
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    alignment = tmp_path / ("S1.cram" if cram else "S1.bam")
    alignment.write_bytes(b"alignment")
    index = Path(f"{alignment}.{'crai' if cram else 'bai'}")
    index.write_bytes(b"index")
    return AlignmentArtifact(
        sample_id="S1",
        path=alignment,
        output_format=(AlignmentOutputFormat.CRAM if cram else AlignmentOutputFormat.BAM),
        reference=ReferenceGenome.from_fasta(fasta, build="GRCh38"),
        source=AlignmentSource.EXISTING,
        sort_order=AlignmentSortOrder.UNKNOWN,
        index_path=index,
    )


def test_request_create_has_deterministic_separate_outputs(tmp_path: Path) -> None:
    request = DeepVariantRequest.create(
        make_artifact(tmp_path),
        tmp_path / "results" / "small",
        resources=SmallVariantResources(threads=12, memory_mb=48_000, runtime_minutes=600),
    )
    assert request.output_vcf.name == "S1.small.vcf.gz"
    assert request.output_gvcf.name == "S1.g.vcf.gz"
    assert request.output_vcf_index.name == "S1.small.vcf.gz.tbi"
    assert request.output_gvcf_index.name == "S1.g.vcf.gz.tbi"
    assert request.resources.threads == 12
    assert request.to_dict()["model_type"] == "PACBIO"


def test_request_supports_indexed_bam_and_cram(tmp_path: Path) -> None:
    assert DeepVariantRequest.create(make_artifact(tmp_path / "bam"), tmp_path / "out1")
    assert DeepVariantRequest.create(
        make_artifact(tmp_path / "cram", cram=True), tmp_path / "out2"
    )


@pytest.mark.parametrize("field", ("threads", "memory_mb", "runtime_minutes"))
def test_resources_require_positive_integers(field: str) -> None:
    values = {"threads": 1, "memory_mb": 1, "runtime_minutes": 1}
    values[field] = 0
    with pytest.raises(InputValidationError, match=field):
        SmallVariantResources(**values)


def test_request_rejects_wrong_output_policy_and_collisions(tmp_path: Path) -> None:
    artifact = make_artifact(tmp_path)
    with pytest.raises(InputValidationError, match="small.vcf.gz"):
        DeepVariantRequest(
            artifact=artifact,
            output_vcf=tmp_path / "S1.vcf.gz",
            output_gvcf=tmp_path / "S1.g.vcf.gz",
        )


def test_path_identity_respects_case_sensitive_platform_policy(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(small_module.os.path, "normcase", lambda value: value)
    upper = small_module._path_identity(tmp_path / "Sample.small.vcf.gz")
    lower = small_module._path_identity(tmp_path / "sample.small.vcf.gz")
    assert upper != lower


def test_request_default_refuses_existing_output(tmp_path: Path) -> None:
    artifact = make_artifact(tmp_path)
    output = tmp_path / "small" / "S1.small.vcf.gz"
    output.parent.mkdir()
    output.write_bytes(b"existing")
    with pytest.raises(OutputValidationError, match="already exists"):
        DeepVariantRequest(
            artifact=artifact,
            output_vcf=output,
            output_gvcf=output.parent / "S1.g.vcf.gz",
        )


def test_result_serialization_distinguishes_planned_and_completed(tmp_path: Path) -> None:
    request = DeepVariantRequest.create(make_artifact(tmp_path), tmp_path / "small")
    command = SmallVariantCommandPlan(("run_deepvariant",), "run_deepvariant")
    planned = SmallVariantResult(request, SmallVariantResultStatus.PLANNED, command)
    completed = SmallVariantResult(
        request,
        SmallVariantResultStatus.COMPLETED,
        command,
        tool_version="1.10.0",
        duration_seconds=3.5,
    )
    assert planned.executed is False
    assert completed.executed is True
    assert completed.to_dict()["tool_version"] == "1.10.0"


def test_validate_small_variant_outputs_checks_pair_and_indexes(tmp_path: Path) -> None:
    request = DeepVariantRequest.create(make_artifact(tmp_path), tmp_path / "small")
    write_outputs(request)
    artifact = validate_small_variant_outputs(request, tool_version="1.10.0")
    assert artifact.sample_id == "S1"
    assert artifact.tool_version == "1.10.0"
    assert artifact.reference_compatibility == "declared_not_header_verified"


def test_deepvariant_1_10_gvcf_header_without_non_ref_is_accepted(
    tmp_path: Path,
) -> None:
    request = DeepVariantRequest.create(make_artifact(tmp_path), tmp_path / "small")
    fixture = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "deepvariant_1_10_gvcf_header.txt"
    ).read_text(encoding="utf-8")
    assert "##ALT=<ID=NON_REF" not in fixture
    write_bgzf(request.output_gvcf, fixture.encode())
    base = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=4>\n"
    columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    write_bgzf(request.output_vcf, (base + columns).encode())
    write_bgzf(request.output_vcf_index, b"TBI\x01")
    write_bgzf(request.output_gvcf_index, b"TBI\x01")

    artifact = validate_small_variant_outputs(request, tool_version="1.10.0")

    assert artifact.gvcf_path == request.output_gvcf


@pytest.mark.parametrize("failure", ("sample", "contig", "index", "gvcf"))
def test_validate_small_variant_outputs_rejects_invalid_products(
    tmp_path: Path,
    failure: str,
) -> None:
    request = DeepVariantRequest.create(make_artifact(tmp_path), tmp_path / "small")
    write_outputs(
        request,
        sample="wrong" if failure == "sample" else "S1",
        contig="1" if failure == "contig" else "chr1",
    )
    if failure == "index":
        request.output_vcf_index.write_bytes(b"not-tabix")
    elif failure == "gvcf":
        base = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=4>\n"
        columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        write_bgzf(request.output_gvcf, (base + columns).encode())
    with pytest.raises(OutputValidationError):
        validate_small_variant_outputs(request)
