"""Real-file failure isolation without changing cohort representation semantics."""

from __future__ import annotations

import gzip
import os
from pathlib import Path

import pytest

import hifivar.cohort_tracks as tracks
from hifivar.cohort import CohortDefinition, CohortSampleInput, SampleCallState
from hifivar.exceptions import InputValidationError, OutputValidationError, ReferenceError
from hifivar.reference import Contig, ReferenceGenome


def context(tmp_path: Path, names=("S1",)) -> CohortDefinition:
    return CohortDefinition("C1", names, ReferenceGenome(
        tmp_path / "ref.fa", tmp_path / "ref.fa.fai", "GRCh38", (Contig("chr1", 1000),)))


def source(tmp_path: Path, name="S1", body=None, *, tr=False) -> CohortSampleInput:
    path = tmp_path / f"{name}.{'tr' if tr else 'sv'}.vcf.gz"
    if body is None:
        body = "chr1\t10\td1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=20\tGT\t0/1\n"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
                     f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{name}\n" + body)
    return CohortSampleInput(name, SampleCallState.CALLED, path,
                            reference_build="GRCh38", catalog_id="catalog" if tr else None)


def destinations(tmp_path: Path):
    return tmp_path / "out/sites.tsv", tmp_path / "out/matrix.tsv"


def stages(tmp_path: Path):
    return tuple(tmp_path.rglob(".hifivar-cohort-*"))


@pytest.mark.parametrize("bad_tail,exception", [
    ("chr2\t30\td2\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=40\tGT\t0/1\n", ReferenceError),
    ("chr1\t30\tmalformed\n", OutputValidationError),
])
def test_sv_late_validation_retains_only_diagnostics_and_can_retry(tmp_path, caplog, bad_tail, exception):
    good = "chr1\t10\td1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=20\tGT\t0/1\n"
    item = source(tmp_path, body=good + bad_tail)
    before = item.source_path.read_bytes()
    sites, matrix = destinations(tmp_path)
    with pytest.raises(exception):
        tracks.build_sv_cohort_tables(context(tmp_path), (item,), site_table=sites, sample_matrix=matrix)
    assert item.source_path.read_bytes() == before
    assert not sites.exists() and not matrix.exists()
    assert stages(tmp_path) and "COHORT_STAGING_RETAINED" in caplog.text
    old_diagnostics = {p: p.read_bytes() for d in stages(tmp_path) for p in d.iterdir() if p.is_file()}
    source(tmp_path, body=good)  # Explicit fixture correction, not production repair.
    result = tracks.build_sv_cohort_tables(context(tmp_path), (item,), site_table=sites, sample_matrix=matrix)
    assert result.outputs == (sites, matrix)
    assert all(p.read_bytes() == data for p, data in old_diagnostics.items())


def test_tr_late_conflict_keeps_unique_sql_and_formal_names_free(tmp_path, caplog):
    def row(pos, trid, gt):
        return f"chr1\t{pos}\t.\tA\t<T>\t.\tPASS\tTRID={trid};END={pos+10};MOTIFS=CAG\tGT\t{gt}\n"
    one = source(tmp_path, "S1", row(10, "L1", "0/1") + row(30, "L2", "1/1"), tr=True)
    two = source(tmp_path, "S2", row(10, "L1", "1/1") + row(31, "L2", "0/1"), tr=True)
    sites, matrix = destinations(tmp_path)
    scratch = tmp_path / "scratch/requested.sqlite"
    scratch.parent.mkdir()
    scratch.write_bytes(b"pre-existing diagnostic, not owned by this call")
    before = scratch.read_bytes()
    with pytest.raises(OutputValidationError, match="inconsistent"):
        tracks.build_tr_cohort_tables(context(tmp_path, ("S1", "S2")), (one, two),
            locus_table=sites, sample_matrix=matrix, scratch_database=scratch)
    assert not sites.exists() and not matrix.exists() and scratch.read_bytes() == before
    assert "COHORT_STAGING_RETAINED" in caplog.text
    retained = {p: p.read_bytes() for d in stages(tmp_path) for p in d.iterdir() if p.is_file()}
    assert any(p.suffix == ".sqlite" for p in retained)
    source(tmp_path, "S2", row(10, "L1", "1/1") + row(30, "L2", "0/1"), tr=True)
    result = tracks.build_tr_cohort_tables(context(tmp_path, ("S1", "S2")), (one, two),
        locus_table=sites, sample_matrix=matrix, scratch_database=scratch)
    assert result.outputs == (sites, matrix) and result.metrics["locus_count"] == 2
    assert scratch.read_bytes() == before
    assert all(p.read_bytes() == data for p, data in retained.items())


@pytest.mark.parametrize("arrival", ["io_error", "late_second", "replace_first_and_late_second", "edit_first_and_late_second", "edit_same_size_and_mtime"])
def test_second_publication_failure_rolls_back_only_owned_first(tmp_path, monkeypatch, arrival):
    item = source(tmp_path)
    before = item.source_path.read_bytes()
    sites, matrix = destinations(tmp_path)
    link = os.link
    calls = []
    edited_first = []
    def fail_second(staged, target, *args, **kwargs):
        calls.append(Path(target))
        if len(calls) == 2:
            if arrival == "replace_first_and_late_second":
                sites.unlink()  # This test explicitly simulates a late external replacement.
                sites.write_bytes(b"external first")
            if arrival == "edit_first_and_late_second":
                sites.write_bytes(b"external edit on the existing first inode")
            if arrival == "edit_same_size_and_mtime":
                metadata = sites.stat()
                edited_first.append(b"x" * len(sites.read_bytes()))
                sites.write_bytes(edited_first[0])
                os.utime(sites, ns=(metadata.st_atime_ns, metadata.st_mtime_ns))
                assert sites.stat().st_mtime_ns == metadata.st_mtime_ns
            if arrival != "io_error":
                matrix.write_bytes(b"external second")
                return link(staged, target, *args, **kwargs)
            raise OSError("injected second-publication I/O failure")
        return link(staged, target, *args, **kwargs)
    monkeypatch.setattr(os, "link", fail_second)
    with pytest.raises(OutputValidationError, match="publication"):
        tracks.build_sv_cohort_tables(context(tmp_path), (item,), site_table=sites, sample_matrix=matrix)
    assert calls == [sites, matrix]
    assert item.source_path.read_bytes() == before
    if arrival == "replace_first_and_late_second":
        assert sites.read_bytes() == b"external first"
    elif arrival == "edit_first_and_late_second":
        assert sites.read_bytes() == b"external edit on the existing first inode"
    elif arrival == "edit_same_size_and_mtime":
        assert sites.read_bytes() == edited_first[0]
    else:
        assert not sites.exists()
    if arrival == "io_error":
        assert not matrix.exists()
    else:
        assert matrix.read_bytes() == b"external second"
    assert stages(tmp_path)
    if arrival == "io_error":
        retained = {p: p.read_bytes() for d in stages(tmp_path) for p in d.iterdir() if p.is_file()}
        monkeypatch.undo()
        tracks.build_sv_cohort_tables(context(tmp_path), (item,), site_table=sites, sample_matrix=matrix)
        assert all(p.read_bytes() == data for p, data in retained.items())


@pytest.mark.parametrize("relation", ["same", "ancestor"])
def test_output_path_conflicts_fail_before_staging(tmp_path, relation):
    item = source(tmp_path)
    sites = tmp_path / "out/site.tsv"
    matrix = sites if relation == "same" else sites / "nested.tsv"
    with pytest.raises((InputValidationError, OutputValidationError)):
        tracks.build_sv_cohort_tables(context(tmp_path), (item,), site_table=sites, sample_matrix=matrix)
    assert not sites.exists() and not matrix.exists() and not stages(tmp_path)


def test_existing_output_and_input_are_not_overwritten(tmp_path):
    item = source(tmp_path)
    before = item.source_path.read_bytes()
    sites, matrix = destinations(tmp_path)
    sites.parent.mkdir()
    sites.write_bytes(b"existing result")
    with pytest.raises(OutputValidationError):
        tracks.build_sv_cohort_tables(context(tmp_path), (item,), site_table=sites, sample_matrix=matrix)
    assert sites.read_bytes() == b"existing result" and not matrix.exists()
    assert item.source_path.read_bytes() == before and not stages(tmp_path)
    with pytest.raises((InputValidationError, OutputValidationError)):
        tracks.build_sv_cohort_tables(context(tmp_path), (item,), site_table=item.source_path, sample_matrix=matrix)
    assert item.source_path.read_bytes() == before


def test_all_rows_validate_before_either_formal_path_is_published(tmp_path, monkeypatch):
    item = source(tmp_path)
    sites, matrix = destinations(tmp_path)
    records = tracks._vcf_records
    observed = []
    def inspect(path):
        for row in records(path):
            observed.append((sites.exists(), matrix.exists()))
            yield row
        observed.append((sites.exists(), matrix.exists()))
    monkeypatch.setattr(tracks, "_vcf_records", inspect)
    tracks.build_sv_cohort_tables(context(tmp_path), (item,), site_table=sites, sample_matrix=matrix)
    assert observed == [(False, False), (False, False)]


def test_success_preserves_exact_bytes_paths_and_distinct_parent_directories(tmp_path):
    item = source(tmp_path)
    sites, matrix = tmp_path / "one/sites.tsv", tmp_path / "two/matrix.tsv"
    result = tracks.build_sv_cohort_tables(context(tmp_path), (item,), site_table=sites, sample_matrix=matrix)
    assert result.outputs == (sites, matrix) and result.metrics["native_site_count"] == 1
    assert sites.read_bytes() == (
        "cohort_site_id\tsource_sample\tsource_variant_id\tcontig\tstart\tend\tsvtype\tref\talt\tsample_support_count\tcallable_sample_count\tsample_support_fraction\tsource_vcf\n"
        f"S1:d1\tS1\td1\tchr1\t10\t20\tDEL\tN\t<DEL>\t1\t1\t1\t{item.source_path}\n"
    ).encode()
    assert matrix.read_bytes() == (
        "cohort_site_id\tsample\tstate\tgenotype\tsource_variant_id\tsource_vcf\n"
        f"S1:d1\tS1\tCALLED\t0/1\td1\t{item.source_path}\n"
    ).encode()
    assert not stages(tmp_path)


def test_success_does_not_clean_foreign_file_inserted_into_staging(tmp_path, monkeypatch):
    item = source(tmp_path)
    sites, matrix = destinations(tmp_path)
    link = os.link
    foreign = []
    def publish(staged, target, *args, **kwargs):
        result = link(staged, target, *args, **kwargs)
        if Path(target) == matrix:
            path = Path(staged).parent / "foreign.keep"
            path.write_bytes(b"external diagnostic")
            foreign.append(path)
        return result
    monkeypatch.setattr(os, "link", publish)
    tracks.build_sv_cohort_tables(context(tmp_path), (item,), site_table=sites, sample_matrix=matrix)
    assert len(foreign) == 1 and foreign[0].read_bytes() == b"external diagnostic"


def test_success_does_not_remove_replaced_staging_file(tmp_path, monkeypatch):
    item = source(tmp_path)
    sites, matrix = destinations(tmp_path)
    link = os.link
    first_stage = []
    def publish(staged, target, *args, **kwargs):
        result = link(staged, target, *args, **kwargs)
        if Path(target) == sites:
            first_stage.append(Path(staged))
        else:
            first_stage[0].unlink()
            first_stage[0].write_bytes(b"external staging replacement")
        return result
    monkeypatch.setattr(os, "link", publish)
    tracks.build_sv_cohort_tables(context(tmp_path), (item,), site_table=sites, sample_matrix=matrix)
    assert first_stage[0].read_bytes() == b"external staging replacement"
    assert b"S1:d1" in sites.read_bytes() and b"external" not in sites.read_bytes()


def test_stage_open_failure_never_creates_formal_tables(tmp_path, monkeypatch, caplog):
    item = source(tmp_path)
    sites, matrix = destinations(tmp_path)
    real_open = tracks._new_text
    calls = []
    def fail_second(path):
        calls.append(path)
        if len(calls) == 2:
            raise OSError("injected second staged-file open failure")
        return real_open(path)
    monkeypatch.setattr(tracks, "_new_text", fail_second)
    with pytest.raises(OSError, match="staged-file open failure"):
        tracks.build_sv_cohort_tables(context(tmp_path), (item,), site_table=sites, sample_matrix=matrix)
    assert not sites.exists() and not matrix.exists()
    assert "COHORT_STAGING_RETAINED" in caplog.text


def test_tr_success_bytes_and_missing_states_unchanged(tmp_path, caplog):
    caplog.set_level("INFO", logger="hifivar.cohort_tracks")
    one = source(tmp_path, "S1", "chr1\t10\t.\tA\t<T>\t.\tPASS\tTRID=L1;END=20;MOTIFS=CAG\tGT\t0/1\n", tr=True)
    two = CohortSampleInput("S2", SampleCallState.NOT_RUN)
    sites, matrix = tmp_path / "a/loci.tsv", tmp_path / "b/matrix.tsv"
    scratch = tmp_path / "c/database.sqlite"
    result = tracks.build_tr_cohort_tables(context(tmp_path, ("S1", "S2")), (one, two),
        locus_table=sites, sample_matrix=matrix, scratch_database=scratch)
    assert result.outputs == (sites, matrix) and result.metrics["locus_count"] == 1
    assert sites.read_bytes() == (
        "trid\tcontig\tstart\tend\tmotifs\tsample_support_count\tcallable_sample_count\tsample_support_fraction\tcatalog_id\n"
        "L1\tchr1\t10\t20\tCAG\t1\t1\t1\tcatalog\n").encode()
    assert matrix.read_bytes() == (
        "trid\tsample\tstate\tgenotype\tsource_vcf\n"
        f"L1\tS1\tCALLED\t0/1\t{one.source_path}\nL1\tS2\tNOT_RUN\t\t\n").encode()
    assert not scratch.exists() and not stages(tmp_path)
    assert "Cohort scratch database:" in caplog.text and ".hifivar-cohort-" in caplog.text


def test_source_key_disambiguation_preserves_both_native_records(tmp_path):
    body = ("chr1\t10\tDUP\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=20\tGT\t0/1\n"
            "chr1\t30\tDUP\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=40\tGT\t1/1\n")
    item = source(tmp_path, body=body)
    sites, matrix = destinations(tmp_path)
    result = tracks.build_sv_cohort_tables(context(tmp_path), (item,), site_table=sites, sample_matrix=matrix)
    rows = [row.split("\t") for row in matrix.read_text().splitlines()[1:]]
    assert result.metrics["native_site_count"] == 2
    assert [row[0] for row in rows] == ["@hifivar-sv-record-v2:S1:1", "@hifivar-sv-record-v2:S1:2"]
    assert [row[4] for row in rows] == ["DUP", "DUP"]
    assert result.metrics["disambiguated_record_count"] == 2
    assert [row[3] for row in rows] == ["0/1", "1/1"]


@pytest.mark.parametrize("protected", ["final", "missing_index"])
def test_scratch_parent_cannot_create_a_directory_at_a_protected_file_path(tmp_path, protected):
    item = source(tmp_path, body="chr1\t10\t.\tA\t<T>\t.\tPASS\tTRID=L1;END=20;MOTIFS=CAG\tGT\t0/1\n", tr=True)
    sites, matrix = destinations(tmp_path)
    absent_index = tmp_path / "absent.tbi"
    item = CohortSampleInput("S1", SampleCallState.CALLED, item.source_path, absent_index,
                            reference_build="GRCh38", catalog_id="catalog")
    parent = sites if protected == "final" else absent_index
    with pytest.raises(InputValidationError, match="scratch parent"):
        tracks.build_tr_cohort_tables(context(tmp_path), (item,), locus_table=sites,
            sample_matrix=matrix, scratch_database=parent / "scratch.sqlite")
    assert not sites.exists() and not matrix.exists() and not absent_index.exists()


@pytest.mark.parametrize("changed_index", [0, 1])
@pytest.mark.parametrize("change", ["edit", "replace"])
def test_seal_to_link_staging_changes_never_return_complete(tmp_path, monkeypatch, changed_index, change):
    item = source(tmp_path)
    before = item.source_path.read_bytes()
    targets = destinations(tmp_path)
    link = os.link
    observed = []
    def alter_before_link(staged, target, *args, **kwargs):
        index = len(observed)
        observed.append(Path(staged))
        if index == changed_index:
            if change == "replace":
                Path(staged).unlink()
            Path(staged).write_bytes(b"external staging bytes after seal")
        return link(staged, target, *args, **kwargs)
    monkeypatch.setattr(os, "link", alter_before_link)
    with pytest.raises(OutputValidationError, match="publication"):
        tracks.build_sv_cohort_tables(context(tmp_path), (item,), site_table=targets[0], sample_matrix=targets[1])
    assert item.source_path.read_bytes() == before
    assert targets[changed_index].read_bytes() == b"external staging bytes after seal"
    assert not targets[1 - changed_index].exists()
    assert stages(tmp_path)


def test_published_first_change_before_second_finishes_is_an_explicit_failure(tmp_path, monkeypatch):
    item = source(tmp_path)
    sites, matrix = destinations(tmp_path)
    link = os.link
    def edit_published(staged, target, *args, **kwargs):
        result = link(staged, target, *args, **kwargs)
        if Path(target) == matrix:
            sites.write_bytes(b"external final bytes")
        return result
    monkeypatch.setattr(os, "link", edit_published)
    with pytest.raises(OutputValidationError, match="publication"):
        tracks.build_sv_cohort_tables(context(tmp_path), (item,), site_table=sites, sample_matrix=matrix)
    assert sites.read_bytes() == b"external final bytes" and not matrix.exists()


@pytest.mark.parametrize("error", [OutputValidationError("injected validation failure"), KeyboardInterrupt("cancelled")])
def test_recoverable_exception_rolls_back_owned_first_without_retyping_interrupt(tmp_path, monkeypatch, error):
    item = source(tmp_path)
    sites, matrix = destinations(tmp_path)
    link = os.link
    def fail_second(staged, target, *args, **kwargs):
        if Path(target) == matrix:
            raise error
        return link(staged, target, *args, **kwargs)
    monkeypatch.setattr(os, "link", fail_second)
    with pytest.raises(type(error)) as caught:
        tracks.build_sv_cohort_tables(context(tmp_path), (item,), site_table=sites, sample_matrix=matrix)
    if isinstance(error, KeyboardInterrupt):
        assert caught.value is error
    assert not sites.exists() and not matrix.exists() and stages(tmp_path)


@pytest.mark.parametrize("interrupted_index", [0, 1])
def test_interrupt_after_real_link_rolls_back_unregistered_owned_output(tmp_path, monkeypatch, interrupted_index):
    item = source(tmp_path)
    before = item.source_path.read_bytes()
    targets = destinations(tmp_path)
    link = os.link
    calls = []
    interrupt = KeyboardInterrupt("real link completed, before Python call returned")
    def link_then_interrupt(staged, target, *args, **kwargs):
        index = len(calls)
        calls.append(Path(target))
        result = link(staged, target, *args, **kwargs)
        if index == interrupted_index:
            assert Path(target).is_file()  # Real filesystem link, not a mocked successful write.
            raise interrupt
        return result
    monkeypatch.setattr(os, "link", link_then_interrupt)
    with pytest.raises(KeyboardInterrupt) as caught:
        tracks.build_sv_cohort_tables(context(tmp_path), (item,), site_table=targets[0], sample_matrix=targets[1])
    assert caught.value is interrupt
    assert calls == list(targets[:interrupted_index + 1])
    assert not targets[0].exists() and not targets[1].exists()
    assert item.source_path.read_bytes() == before and stages(tmp_path)


@pytest.mark.parametrize("arrival_index", [0, 1])
def test_failed_link_candidate_does_not_own_late_external_file(tmp_path, monkeypatch, arrival_index):
    item = source(tmp_path)
    targets = destinations(tmp_path)
    link = os.link
    calls = []
    def late_file_before_link(staged, target, *args, **kwargs):
        index = len(calls)
        calls.append(Path(target))
        if index == arrival_index:
            Path(target).write_bytes(b"late external output, not linked from staging")
        return link(staged, target, *args, **kwargs)
    monkeypatch.setattr(os, "link", late_file_before_link)
    with pytest.raises(OutputValidationError, match="publication"):
        tracks.build_sv_cohort_tables(context(tmp_path), (item,), site_table=targets[0], sample_matrix=targets[1])
    assert targets[arrival_index].read_bytes() == b"late external output, not linked from staging"
    assert not targets[1 - arrival_index].exists() and stages(tmp_path)
