from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from hifivar.cohort import CohortDefinition, CohortSampleInput, SampleCallState
from hifivar.cohort_tracks import build_sv_cohort_tables, build_tr_cohort_tables
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.reference import Contig, ReferenceGenome


def cohort(tmp_path: Path) -> CohortDefinition:
    return CohortDefinition("C1", ("S1", "S2", "S3"), ReferenceGenome(tmp_path / "ref.fa", tmp_path / "ref.fa.fai", "GRCh38", (Contig("chr1", 1000),)))


def vcf(path: Path, sample: str, record: str, *, tr: bool = False) -> None:
    info_headers = "##INFO=<ID=TRID,Number=1,Type=String,Description=x>\n##INFO=<ID=MOTIFS,Number=1,Type=String,Description=x>\n" if tr else "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=x>\n"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(f"##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n{info_headers}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n{record}")
    Path(f"{path}.tbi").write_bytes(b"index")


def test_sv_native_representation_and_missing_not_homref(tmp_path: Path) -> None:
    source = tmp_path / "S1.sv.vcf.gz"
    vcf(source, "S1", "chr1\t10\td1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=20\tGT\t0/1\n")
    empty = tmp_path / "S2.empty.vcf.gz"
    vcf(empty, "S2", "")
    inputs = (
        CohortSampleInput("S1", SampleCallState.CALLED, source, Path(f"{source}.tbi"), "jasmine", "1.1.5", "GRCh38"),
        CohortSampleInput("S2", SampleCallState.NO_CALLS, empty, Path(f"{empty}.tbi"), "jasmine", "1.1.5", "GRCh38"),
        CohortSampleInput("S3", SampleCallState.FAILED),
    )
    result = build_sv_cohort_tables(cohort(tmp_path), inputs, site_table=tmp_path / "sites.tsv", sample_matrix=tmp_path / "matrix.tsv")
    matrix = (tmp_path / "matrix.tsv").read_text(encoding="utf-8")
    assert "d1\tS2\tNOT_OBSERVED\t" in matrix
    assert "d1\tS3\tFAILED\t" in matrix
    assert "0/0" not in matrix
    assert result.metrics["cross_sample_clustering"] is False
    assert "allele_frequency" not in result.metrics
    assert result.metrics["native_site_count"] == 1


def test_sv_overwrite_and_contig_safety(tmp_path: Path) -> None:
    source = tmp_path / "S1.vcf.gz"
    vcf(source, "S1", "chr2\t10\tx\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=20\tGT\t0/1\n")
    inputs = (CohortSampleInput("S1", SampleCallState.CALLED, source, Path(f"{source}.tbi"), reference_build="GRCh38"), CohortSampleInput("S2", SampleCallState.FAILED), CohortSampleInput("S3", SampleCallState.FAILED))
    with pytest.raises(Exception, match="CONTIG_MISMATCH"):
        build_sv_cohort_tables(cohort(tmp_path), inputs, site_table=tmp_path / "sites.tsv", sample_matrix=tmp_path / "matrix.tsv")


def test_tr_locus_matrix_catalog_and_partial_states(tmp_path: Path) -> None:
    one, two = tmp_path / "S1.tr.vcf.gz", tmp_path / "S2.tr.vcf.gz"
    vcf(one, "S1", "chr1\t10\t.\tA\t<STR10>\t.\tPASS\tTRID=L1;MOTIFS=CAG;END=20\tGT:AL\t0/1:3,10\n", tr=True)
    vcf(two, "S2", "chr1\t10\t.\tA\t<STR11>\t.\tPASS\tTRID=L1;MOTIFS=CAG;END=20\tGT:AL\t1/1:11,11\n", tr=True)
    inputs = (
        CohortSampleInput("S1", SampleCallState.CALLED, one, Path(f"{one}.tbi"), "trgt", "5.1.0", "GRCh38", "catalog-sha256"),
        CohortSampleInput("S2", SampleCallState.CALLED, two, Path(f"{two}.tbi"), "trgt", "5.1.0", "GRCh38", "catalog-sha256"),
        CohortSampleInput("S3", SampleCallState.NOT_RUN),
    )
    result = build_tr_cohort_tables(cohort(tmp_path), inputs, locus_table=tmp_path / "loci.tsv", sample_matrix=tmp_path / "matrix.tsv", scratch_database=tmp_path / "scratch.sqlite")
    assert result.metrics["locus_count"] == 1
    text = (tmp_path / "matrix.tsv").read_text(encoding="utf-8")
    assert "L1\tS3\tNOT_RUN\t" in text
    assert "L1\tS1\tCALLED\t0/1" in text


def test_tr_requires_identical_catalog_and_representation(tmp_path: Path) -> None:
    one, two = tmp_path / "S1.tr.vcf.gz", tmp_path / "S2.tr.vcf.gz"
    vcf(one, "S1", "chr1\t10\t.\tA\t<T>\t.\tPASS\tTRID=L1;MOTIFS=CAG;END=20\tGT\t0/1\n", tr=True)
    vcf(two, "S2", "chr1\t11\t.\tA\t<T>\t.\tPASS\tTRID=L1;MOTIFS=CAG;END=21\tGT\t0/1\n", tr=True)
    base = lambda sample, path, catalog: CohortSampleInput(sample, SampleCallState.CALLED, path, Path(f"{path}.tbi"), "trgt", "5.1.0", "GRCh38", catalog)
    with pytest.raises(InputValidationError, match="catalog_id"):
        build_tr_cohort_tables(cohort(tmp_path), (base("S1", one, "A"), base("S2", two, "B"), CohortSampleInput("S3", SampleCallState.FAILED)), locus_table=tmp_path / "l1.tsv", sample_matrix=tmp_path / "m1.tsv", scratch_database=tmp_path / "s1.sqlite")
    with pytest.raises(OutputValidationError, match="inconsistent"):
        build_tr_cohort_tables(cohort(tmp_path), (base("S1", one, "A"), base("S2", two, "A"), CohortSampleInput("S3", SampleCallState.FAILED)), locus_table=tmp_path / "l2.tsv", sample_matrix=tmp_path / "m2.tsv", scratch_database=tmp_path / "s2.sqlite")
