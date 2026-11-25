"""Source-record identity is distinct from caller IDs and cross-sample loci."""

from __future__ import annotations

import csv
import gzip
import hashlib
import sqlite3
from pathlib import Path

import pytest

from hifivar.cohort import CohortDefinition, CohortSampleInput, SampleCallState
from hifivar.cohort import CohortTrack, CohortTrackResult
import hifivar.cohort_tracks as tracks
from hifivar.cohort_tracks import build_sv_cohort_tables, build_tr_cohort_tables
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.reference import Contig, ReferenceGenome
from hifivar.phase12 import run_phase12


def _row(record_id: str, *, pos: int = 10, end: int = 20, gt: str = "0/1",
         alt: str = "<DEL>") -> str:
    return f"chr1\t{pos}\t{record_id}\tN\t{alt}\t.\tPASS\tSVTYPE=DEL;END={end}\tGT\t{gt}\n"


def _source(root: Path, name: str, body: str) -> CohortSampleInput:
    path = root / f"{name}.sv.vcf.gz"
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
                     "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=x>\n"
                     f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{name}\n{body}")
    index = Path(f"{path}.tbi")
    index.write_bytes(b"fixture-index")
    return CohortSampleInput(name, SampleCallState.CALLED, path, index,
                             reference_build="GRCh38")


def _cohort(root: Path, names: tuple[str, ...]) -> CohortDefinition:
    return CohortDefinition("C1", names, ReferenceGenome(
        root / "ref.fa", root / "ref.fa.fai", "GRCh38", (Contig("chr1", 1000),)))


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _build(root: Path, inputs: tuple[CohortSampleInput, ...], name: str = "out"):
    sites, matrix = root / name / "sites.tsv", root / name / "matrix.tsv"
    result = build_sv_cohort_tables(_cohort(root, tuple(x.sample_id for x in inputs)),
                                   inputs, site_table=sites, sample_matrix=matrix)
    return result, sites, matrix


@pytest.mark.parametrize("body,expected_ids", [
    (_row("DUP") + _row("DUP", pos=30, end=40, gt="1/1"), ["DUP", "DUP"]),
    (_row("DUP") + _row("DUP", gt="1/1"), ["DUP", "DUP"]),
    (_row(".") + _row(".", end=40, gt="1/1"), ["chr1:10:N:<DEL>"] * 2),
    (_row("chr1:10:N:<DEL>", pos=30) + _row("."), ["chr1:10:N:<DEL>"] * 2),
    (_row("SAME") + _row("SAME"), ["SAME", "SAME"]),
])
def test_colliding_legacy_ids_keep_all_records_and_unique_cells(tmp_path, body, expected_ids):
    item = _source(tmp_path, "S1", body)
    absent = CohortSampleInput("S2", SampleCallState.FAILED)
    before = {p: hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (item.source_path, item.index_path)}
    result, sites, matrix = _build(tmp_path, (item, absent))
    rows, cells = _read(sites), _read(matrix)
    assert len(rows) == 2 and len(cells) == 4
    assert len({row["cohort_site_id"] for row in rows}) == 2
    assert len({(row["cohort_site_id"], row["sample"]) for row in cells}) == 4
    assert [row["source_variant_id"] for row in rows] == expected_ids
    assert [row["start"] for row in rows] == [x.split("\t")[1] for x in body.splitlines()]
    assert [row["genotype"] for row in cells if row["sample"] == "S1"] == [
        x.split("\t")[9] for x in body.splitlines()]
    assert all(row["state"] == "FAILED" and row["genotype"] == ""
               for row in cells if row["sample"] == "S2")
    assert result.metrics["disambiguated_record_count"] == 2
    assert result.metrics["site_key_uniqueness_verified"] is True
    assert result.metrics["cross_sample_clustering"] is False
    assert {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in before} == before


def test_reserved_keys_cannot_collide_with_a_native_id_or_legacy_suffix(tmp_path):
    body = (_row("DUP") + _row("DUP", pos=30) + _row("DUP:record=1", pos=50)
            + _row("@hifivar-sv-record-v2:S1:1", pos=60))
    result, sites, matrix = _build(tmp_path, (_source(tmp_path, "S1", body),))
    rows = _read(sites)
    assert len({row["cohort_site_id"] for row in rows}) == 4
    assert [row["cohort_site_id"] for row in rows[:2]] == [
        "@hifivar-sv-record-v2:S1:1", "@hifivar-sv-record-v2:S1:2"]
    assert rows[2]["cohort_site_id"] == "S1:DUP:record=1"
    assert rows[3]["cohort_site_id"] == "S1:@hifivar-sv-record-v2:S1:1"
    assert len(_read(matrix)) == result.metrics["native_site_count"] == 4


def test_ordinal_counts_records_not_headers_comments_or_blank_lines(tmp_path):
    body = _row("DUP") + "\n# retained comment\n" + _row("DUP", pos=30)
    _, sites, _ = _build(tmp_path, (_source(tmp_path, "S1", body),))
    assert [row["cohort_site_id"] for row in _read(sites)] == [
        "@hifivar-sv-record-v2:S1:1", "@hifivar-sv-record-v2:S1:2"]


def test_unique_native_and_missing_ids_keep_legacy_keys_and_columns(tmp_path):
    item = _source(tmp_path, "S1", _row("d1") + _row(".", pos=30))
    result, sites, matrix = _build(tmp_path, (item,))
    assert [row["cohort_site_id"] for row in _read(sites)] == ["S1:d1", "S1:chr1:30:N:<DEL>"]
    assert [row["source_variant_id"] for row in _read(sites)] == ["d1", "chr1:30:N:<DEL>"]
    assert result.metrics["disambiguated_record_count"] == 0
    assert result.metrics["site_key_policy"] == "legacy_unique_else_reserved_record_ordinal_v2"
    assert matrix.read_text(encoding="utf-8").splitlines()[0] == (
        "cohort_site_id\tsample\tstate\tgenotype\tsource_variant_id\tsource_vcf")


def test_same_native_id_across_samples_is_not_clustered_or_homref(tmp_path):
    inputs = (_source(tmp_path, "S1", _row("d1")),
              _source(tmp_path, "S2", _row("d1", gt="1/1")))
    result, sites, matrix = _build(tmp_path, inputs)
    assert [row["cohort_site_id"] for row in _read(sites)] == ["S1:d1", "S2:d1"]
    cells = _read(matrix)
    assert [(x["sample"], x["state"], x["genotype"]) for x in cells] == [
        ("S1", "CALLED", "0/1"), ("S2", "NOT_OBSERVED", ""),
        ("S1", "NOT_OBSERVED", ""), ("S2", "CALLED", "1/1")]
    assert result.metrics["disambiguated_record_count"] == 0


def test_multiallelic_record_and_genotype_are_not_split(tmp_path):
    body = _row("a;b", alt="<DEL>,<DUP>", gt="1/2")
    result, sites, matrix = _build(tmp_path, (_source(tmp_path, "S1", body),))
    assert len(_read(sites)) == result.metrics["native_site_count"] == 1
    assert _read(sites)[0]["alt"] == "<DEL>,<DUP>"
    assert _read(sites)[0]["source_variant_id"] == "a;b"
    assert _read(matrix)[0]["genotype"] == "1/2"


def test_collision_keys_are_deterministic_across_fresh_output_directories(tmp_path):
    inputs = (_source(tmp_path, "S1", _row("DUP") + _row("DUP", pos=30)),)
    one, s1, m1 = _build(tmp_path, inputs, "one")
    two, s2, m2 = _build(tmp_path, inputs, "two")
    assert s1.read_bytes() == s2.read_bytes() and m1.read_bytes() == m2.read_bytes()
    assert one.metrics == two.metrics
    assert not list(tmp_path.rglob(".hifivar-cohort-*"))


@pytest.mark.parametrize("tr", [False, True])
@pytest.mark.parametrize("tail", ["\textra-sample", "\textra-sample\tunexpected"])
def test_single_sample_records_reject_unadvertised_sample_columns(tmp_path, tr, tail):
    body = _row("d1").rstrip("\n") + tail + "\n"
    if tr:
        body = body.replace("SVTYPE=DEL;END=20", "TRID=L1;MOTIFS=CAG;END=20")
    item = _source(tmp_path, "S1", body)
    if tr:
        item = CohortSampleInput("S1", SampleCallState.CALLED, item.source_path,
                                 item.index_path, reference_build="GRCh38", catalog_id="catalog")
    before = item.source_path.read_bytes()
    sites, matrix = tmp_path / "sites.tsv", tmp_path / "matrix.tsv"
    with pytest.raises(OutputValidationError, match="Malformed single-sample"):
        if tr:
            build_tr_cohort_tables(_cohort(tmp_path, ("S1",)), (item,),
                locus_table=sites, sample_matrix=matrix, scratch_database=tmp_path / "scratch.sqlite")
        else:
            build_sv_cohort_tables(_cohort(tmp_path, ("S1",)), (item,),
                site_table=sites, sample_matrix=matrix)
    assert not sites.exists() and not matrix.exists()
    assert item.source_path.read_bytes() == before
    assert list(tmp_path.rglob(".hifivar-cohort-*"))


def test_same_file_cannot_be_relabelled_as_another_sample(tmp_path):
    one = _source(tmp_path, "S1", _row("d1"))
    two = CohortSampleInput("S2", SampleCallState.CALLED, one.source_path, one.index_path,
                            reference_build="GRCh38")
    with pytest.raises(OutputValidationError, match="sample mismatch"):
        _build(tmp_path, (one, two))
    assert not (tmp_path / "out/sites.tsv").exists()
    assert not (tmp_path / "out/matrix.tsv").exists()


def test_duplicate_input_sample_is_rejected_before_staging(tmp_path):
    one = _source(tmp_path, "S1", _row("d1"))
    with pytest.raises(InputValidationError):
        build_sv_cohort_tables(_cohort(tmp_path, ("S1",)), (one, one),
            site_table=tmp_path / "sites.tsv", sample_matrix=tmp_path / "matrix.tsv")
    assert not list(tmp_path.rglob(".hifivar-cohort-*"))


def test_interleaved_collision_groups_keep_complete_source_order(tmp_path):
    body = "".join(_row(key, pos=10 + i, gt=gt) for i, (key, gt) in enumerate([
        ("A", "0/1"), ("B", "1/1"), ("A", "1|0"), ("unique", "./."), ("B", "0|1")]))
    result, sites, matrix = _build(tmp_path, (_source(tmp_path, "S1", body),))
    assert [x["source_variant_id"] for x in _read(sites)] == ["A", "B", "A", "unique", "B"]
    assert [x["cohort_site_id"] for x in _read(sites)] == [
        "@hifivar-sv-record-v2:S1:1", "@hifivar-sv-record-v2:S1:2",
        "@hifivar-sv-record-v2:S1:3", "S1:unique", "@hifivar-sv-record-v2:S1:5"]
    assert [x["genotype"] for x in _read(matrix)] == ["0/1", "1/1", "1|0", "./.", "0|1"]
    assert result.metrics["disambiguated_record_count"] == 4


def test_two_samples_with_internal_duplicates_have_distinct_key_namespaces(tmp_path):
    inputs = tuple(_source(tmp_path, sample, _row("DUP") + _row("DUP", pos=30))
                   for sample in ("S1", "S2"))
    result, sites, matrix = _build(tmp_path, inputs)
    assert len({x["cohort_site_id"] for x in _read(sites)}) == 4
    assert len({(x["cohort_site_id"], x["sample"]) for x in _read(matrix)}) == 8
    assert result.metrics["disambiguated_record_count"] == 4
    assert all(x["genotype"] == "" for x in _read(matrix) if x["state"] == "NOT_OBSERVED")


def test_empty_source_produces_only_headers_and_cleans_staging(tmp_path):
    item = _source(tmp_path, "S1", "")
    result, sites, matrix = _build(tmp_path, (item,))
    assert _read(sites) == [] and _read(matrix) == []
    assert len(sites.read_text().splitlines()) == len(matrix.read_text().splitlines()) == 1
    assert result.metrics["native_site_count"] == result.metrics["disambiguated_record_count"] == 0
    assert result.metrics["site_key_uniqueness_verified"] is True
    assert not list(tmp_path.rglob(".hifivar-cohort-*"))


@pytest.mark.parametrize("stage", ["connect", "create", "insert", "count", "commit", "assigned"])
def test_sqlite_failures_preserve_sources_and_publish_neither_table(tmp_path, monkeypatch, stage):
    item = _source(tmp_path, "S1", _row("DUP") + _row("DUP", pos=30))
    before = item.source_path.read_bytes()
    connect = tracks.sqlite3.connect

    class FailingConnection(sqlite3.Connection):
        def execute(self, sql, parameters=()):
            if ((stage == "create" and sql.startswith("CREATE TABLE"))
                or (stage == "insert" and sql.startswith("INSERT INTO sv_records"))
                or (stage == "count" and sql.startswith("INSERT INTO sv_key_counts"))
                or (stage == "assigned" and sql.startswith("INSERT INTO sv_assigned_keys"))):
                raise sqlite3.OperationalError("injected disk/index failure")
            return super().execute(sql, parameters)

        def commit(self):
            if stage == "commit":
                raise sqlite3.OperationalError("injected disk commit failure")
            return super().commit()

    def fail_connect(*args, **kwargs):
        if stage == "connect":
            raise sqlite3.OperationalError("injected database open failure")
        return connect(*args, factory=FailingConnection, **kwargs)

    monkeypatch.setattr(tracks.sqlite3, "connect", fail_connect)
    with pytest.raises(OutputValidationError, match="SV cohort source keys"):
        _build(tmp_path, (item,))
    assert item.source_path.read_bytes() == before
    assert not (tmp_path / "out/sites.tsv").exists()
    assert not (tmp_path / "out/matrix.tsv").exists()
    assert list(tmp_path.rglob(".hifivar-cohort-*"))


def test_database_open_failure_is_isolated_to_the_sv_track(tmp_path, monkeypatch):
    item = _source(tmp_path, "S1", _row("d1"))
    definition = _cohort(tmp_path, ("S1",))
    completed = []

    def fail_connect(*args, **kwargs):
        raise sqlite3.OperationalError("injected database open failure")

    def independent_tr():
        completed.append("tr")
        return CohortTrackResult(CohortTrack.TR, True, SampleCallState.CALLED, "synthetic", None)

    monkeypatch.setattr(tracks.sqlite3, "connect", fail_connect)
    report = run_phase12(definition,
        enabled_tracks={CohortTrack.SV: True, CohortTrack.TR: True},
        operations={CohortTrack.SV: lambda: _build(tmp_path, (item,))[0], CohortTrack.TR: independent_tr})
    assert report.failed_tracks == (CohortTrack.SV,)
    assert report.successful_tracks == (CohortTrack.TR,) and completed == ["tr"]
    assert not (tmp_path / "out/sites.tsv").exists()
    assert not (tmp_path / "out/matrix.tsv").exists()
