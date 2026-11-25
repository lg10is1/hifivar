from pathlib import Path
import hashlib

import pytest

from hifivar.exceptions import InputValidationError, ReferenceError
from hifivar.reference import ReferenceGenome
from hifivar.tr import TandemRepeatCatalog


def make_reference(tmp_path: Path) -> ReferenceGenome:
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t8\t6\t8\t9\n", encoding="utf-8")
    return ReferenceGenome.from_fasta(fasta, build="GRCh38")


def test_catalog_stream_validation_accepts_required_trgt_tags(tmp_path: Path) -> None:
    reference = make_reference(tmp_path)
    catalog_path = tmp_path / "catalog.bed"
    catalog_path.write_text("chr1\t1\t7\tID=TR1;MOTIFS=CG;STRUC=<TR>\n", encoding="utf-8")
    catalog = TandemRepeatCatalog(catalog_path, "GRCh38")
    assert catalog.validate(reference) is catalog


@pytest.mark.parametrize(
    "line,pattern",
    [
        ("chr1\t1\t7\tID=TR1;MOTIFS=CG\n", "missing annotation"),
        ("chr1\t7\t1\tID=TR1;MOTIFS=CG;STRUC=<TR>\n", "invalid BED"),
        ("chr1\t1\t7\tID=TR1;MOTIFS=CG;STRUC=<TR>\nchr1\t1\t7\tID=TR1;MOTIFS=CG;STRUC=<TR>\n", "duplicate ID"),
    ],
)
def test_catalog_rejects_invalid_rows(tmp_path: Path, line: str, pattern: str) -> None:
    reference = make_reference(tmp_path)
    path = tmp_path / "catalog.bed"
    path.write_text(line, encoding="utf-8")
    with pytest.raises(InputValidationError, match=pattern):
        TandemRepeatCatalog(path).validate(reference)


def test_catalog_rejects_contig_and_build_mismatch(tmp_path: Path) -> None:
    reference = make_reference(tmp_path)
    path = tmp_path / "catalog.bed"
    path.write_text("1\t1\t7\tID=TR1;MOTIFS=CG;STRUC=<TR>\n", encoding="utf-8")
    with pytest.raises(ReferenceError, match="REFERENCE_CONTIG_MISMATCH"):
        TandemRepeatCatalog(path).validate(reference)
    path.write_text("chr1\t1\t7\tID=TR1;MOTIFS=CG;STRUC=<TR>\n", encoding="utf-8")
    with pytest.raises(ReferenceError, match="catalog build"):
        TandemRepeatCatalog(path, "GRCh37").validate(reference)


@pytest.mark.parametrize("end", [7, 8])
def test_catalog_accepts_bed_end_up_to_exact_reference_length(tmp_path: Path, end: int) -> None:
    reference = make_reference(tmp_path)
    path = tmp_path / "catalog.bed"
    path.write_text(f"chr1\t1\t{end}\tID=L1;MOTIFS=CG;STRUC=(CG)n\n", encoding="utf-8")
    before = path.read_bytes()
    catalog = TandemRepeatCatalog(path, "GRCh38")
    assert catalog.validate(reference) is catalog
    assert path.read_bytes() == before


@pytest.mark.parametrize("start,end", [(1, 9), (8, 9), (9, 10)])
def test_catalog_rejects_reference_out_of_bounds_without_rewriting(tmp_path: Path, start: int, end: int) -> None:
    reference = make_reference(tmp_path)
    path = tmp_path / "catalog.bed"
    path.write_text(f"chr1\t{start}\t{end}\tID=L1;MOTIFS=CG;STRUC=(CG)n\n", encoding="utf-8")
    protected = (path, reference.fasta, reference.fai)
    before = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in protected}
    with pytest.raises(InputValidationError, match="exceeds reference contig 'chr1' length 8"):
        TandemRepeatCatalog(path, "GRCh38").validate(reference)
    assert {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in protected} == before


def test_catalog_uses_each_contigs_own_length(tmp_path: Path) -> None:
    from hifivar.reference import Contig

    reference = ReferenceGenome(tmp_path / "ref.fa", tmp_path / "ref.fa.fai", "GRCh38",
                                (Contig("chr1", 8), Contig("chr2", 100)))
    path = tmp_path / "catalog.bed"
    path.write_text("chr2\t1\t100\tID=L1;MOTIFS=CG;STRUC=(CG)n\n", encoding="utf-8")
    assert TandemRepeatCatalog(path, "GRCh38").validate(reference)
    path.write_text("chr1\t1\t100\tID=L1;MOTIFS=CG;STRUC=(CG)n\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="exceeds reference contig 'chr1' length 8"):
        TandemRepeatCatalog(path, "GRCh38").validate(reference)


def test_out_of_bounds_catalog_fails_wrapper_before_version_probe(tmp_path: Path, monkeypatch) -> None:
    from hifivar.alignment import AlignmentOutputFormat
    from hifivar.alignment_postprocess import AlignmentArtifact, AlignmentSortOrder, AlignmentSource
    from hifivar.trgt import TrgtRequest, TrgtWrapper

    reference = make_reference(tmp_path)
    bam = tmp_path / "input.bam"
    bam.write_bytes(b"immutable BAM")
    bai = Path(str(bam) + ".bai")
    bai.write_bytes(b"immutable BAI")
    catalog = tmp_path / "catalog.bed"
    catalog.write_text("chr1\t1\t9\tID=L1;MOTIFS=CG;STRUC=(CG)n\n", encoding="utf-8")
    protected = (bam, bai, catalog, reference.fasta, reference.fai)
    before = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in protected}
    artifact = AlignmentArtifact("S1", bam, AlignmentOutputFormat.BAM, reference,
                                 AlignmentSource.EXISTING, AlignmentSortOrder.UNKNOWN, bai)
    request = TrgtRequest(artifact, TandemRepeatCatalog(catalog, "GRCh38"),
                          tmp_path / "work" / "S1.trgt", tmp_path / "out" / "S1.tr.vcf.gz",
                          tmp_path / "out" / "S1.tr.spanning.bam", "XY")
    monkeypatch.setattr(TrgtWrapper, "detect_versions", lambda self: pytest.fail("Must reject catalog before external execution"))
    with pytest.raises(InputValidationError, match="exceeds reference contig"):
        TrgtWrapper().run(request)
    assert {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in protected} == before
    assert not request.raw_output_prefix.parent.exists()
    assert not request.final_vcf.parent.exists()
