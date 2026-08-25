from pathlib import Path

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
