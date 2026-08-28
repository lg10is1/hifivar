from pathlib import Path

import pytest

from hifivar.assembly_sv import SVEvidenceSource
from hifivar.exceptions import OutputValidationError
from hifivar.harmonization import (
    EvidenceRunStatus,
    SVEvidenceSourceArtifact,
    SVHarmonizationRequest,
    write_evidence_table,
)
from hifivar.reference import ReferenceGenome


@pytest.fixture
def reference(tmp_path: Path) -> ReferenceGenome:
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\n" + "A" * 1000 + "\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text(
        "chr1\t1000\t6\t1000\t1001\n", encoding="utf-8"
    )
    return ReferenceGenome.from_fasta(fasta, build="GRCh38")


def _sources(tmp_path: Path) -> tuple[SVEvidenceSourceArtifact, ...]:
    result = []
    for caller in ("sawfish", "sniffles2", "pbsv", "cutesv"):
        vcf_path = tmp_path / f"{caller}.vcf.gz"
        index_path = Path(f"{vcf_path}.tbi")
        vcf_path.write_bytes(b"source-vcf")
        index_path.write_bytes(b"source-index")
        result.append(
            SVEvidenceSourceArtifact(
                sample_id="S1",
                source=SVEvidenceSource.READ,
                caller=caller,
                vcf_path=vcf_path,
                index_path=index_path,
                status=EvidenceRunStatus.COMPLETED,
            )
        )
    return tuple(result)


def _request(
    tmp_path: Path,
    reference: ReferenceGenome,
    sources: tuple[SVEvidenceSourceArtifact, ...],
) -> SVHarmonizationRequest:
    return SVHarmonizationRequest(
        sample_id="S1",
        reference=reference,
        sources=sources,
        work_directory=tmp_path / "work",
        output_vcf=tmp_path / "out" / "S1.harmonized.sv.vcf.gz",
        evidence_table=tmp_path / "out" / "S1.sv.evidence.tsv",
    )


def _merged_vcf(
    tmp_path: Path,
    identifiers: str,
    support_vector: str,
) -> Path:
    path = tmp_path / "merged.vcf"
    path.write_text(
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t10\tJ1\tA\t<DEL>\t.\tPASS\t"
        f"SVTYPE=DEL;END=20;IDLIST={identifiers};SUPP_VEC={support_vector}"
        "\tGT\t0/1\n",
        encoding="utf-8",
        newline="\n",
    )
    return path


@pytest.mark.parametrize(
    ("identifiers", "support_vector", "expected_callers"),
    [
        ("sawfish:0:0:0:0,pbsv.DEL.0", "1010", ("sawfish", "pbsv")),
        (
            "sawfish:0:1:0:0,Sniffles2.INS.1S0",
            "1100",
            ("sawfish", "sniffles2"),
        ),
        ("pbsv.INS.1,cuteSV.INS.0", "0011", ("pbsv", "cutesv")),
        (
            "sawfish:r1,Sniffles2.DEL.1,pbsv.DEL.1,cuteSV.DEL.1",
            "1111",
            ("sawfish", "sniffles2", "pbsv", "cutesv"),
        ),
        ("cuteSV.INS.0", "0001", ("cutesv",)),
    ],
)
def test_support_vector_is_authoritative_for_native_caller_ids(
    tmp_path: Path,
    reference: ReferenceGenome,
    identifiers: str,
    support_vector: str,
    expected_callers: tuple[str, ...],
) -> None:
    sources = _sources(tmp_path)
    request = _request(tmp_path, reference, sources)
    merged = _merged_vcf(tmp_path, identifiers, support_vector)
    original = merged.read_bytes()
    source_originals = {
        item.vcf_path: item.vcf_path.read_bytes() for item in sources if item.vcf_path
    }

    write_evidence_table(request, merged)

    fields = request.evidence_table.read_text(encoding="utf-8").splitlines()[1].split("\t")
    assert fields[7] == str(len(expected_callers))
    assert fields[8] == str(len(expected_callers))
    assert fields[10] == ",".join(expected_callers)
    assert fields[11] == identifiers
    assert fields[12].split(",") == [
        str(next(item.vcf_path for item in sources if item.caller == caller))
        for caller in expected_callers
    ]
    assert merged.read_bytes() == original
    assert {
        path: path.read_bytes() for path in source_originals
    } == source_originals


@pytest.mark.parametrize("support_vector", ["101", "10x0", "0000"])
def test_invalid_support_vector_is_rejected(
    tmp_path: Path,
    reference: ReferenceGenome,
    support_vector: str,
) -> None:
    sources = _sources(tmp_path)
    request = _request(tmp_path, reference, sources)
    merged = _merged_vcf(tmp_path, "native1", support_vector)

    with pytest.raises(OutputValidationError, match="SUPP_VEC"):
        write_evidence_table(request, merged)


def test_legacy_partial_id_mapping_is_unresolved(
    tmp_path: Path,
    reference: ReferenceGenome,
) -> None:
    sources = _sources(tmp_path)
    request = _request(tmp_path, reference, sources)
    merged = tmp_path / "legacy.vcf"
    merged.write_text(
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t10\tJ1\tA\t<DEL>\t.\tPASS\t"
        "SVTYPE=DEL;END=20;IDLIST=sawfish:r1,pbsv.DEL.0\tGT\t0/1\n",
        encoding="utf-8",
        newline="\n",
    )

    write_evidence_table(request, merged)

    fields = request.evidence_table.read_text(encoding="utf-8").splitlines()[1].split("\t")
    assert fields[6] == "UNRESOLVED"
    assert fields[7] == "0"
    assert fields[10] == ""
