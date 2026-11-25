"""TRGT locus identity must be explicit; storage failures stay in their track."""
from __future__ import annotations

import csv
import gzip
import hashlib
from pathlib import Path
import sqlite3
from unittest import mock

import pytest

from hifivar.cohort import CohortDefinition, CohortSampleInput, CohortTrack, SampleCallState
import hifivar.cohort_tracks as tracks
from hifivar.exceptions import OutputValidationError
from hifivar.phase12 import run_phase12
from hifivar.reference import Contig, ReferenceGenome


def _fixture(root: Path, *, pos='10', end='20', trid='L1', info=None, gt='0/1'):
    fasta = root / 'ref.fa'
    fasta.write_text('>chr1\n' + 'A' * 1000 + '\n', encoding='ascii')
    fai = root / 'ref.fa.fai'
    fai.write_text('chr1\t1000\t6\t1000\t1001\n', encoding='ascii')
    reference = ReferenceGenome(fasta, fai, 'GRCh38', (Contig('chr1', 1000),))
    cohort = CohortDefinition('C1', ('S1',), reference)
    fields = ['MOTIFS=CAG']
    if trid is not None:
        fields.append('TRID=' + trid)
    if end is not None:
        fields.append('END=' + end)
    body_info = info if info is not None else ';'.join(fields)
    source = root / 'S1.tr.vcf.gz'
    with gzip.open(source, 'wt', encoding='utf-8', newline='') as out:
        out.write('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
                  '##INFO=<ID=TRID,Number=1,Type=String,Description="Locus">\n'
                  '##INFO=<ID=END,Number=1,Type=Integer,Description="Endpoint">\n'
                  '##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motifs">\n'
                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
                  f'chr1\t{pos}\t.\tA\tC,CC\t.\tPASS\t{body_info}\tGT\t{gt}\n')
    index = Path(str(source) + '.tbi')
    index.write_bytes(b'SYNTHETIC_NONFUNCTIONAL_INDEX_NOT_A_TABIX_TEST')
    item = CohortSampleInput('S1', SampleCallState.CALLED, source, index,
                            reference_build='GRCh38', catalog_id='catalog')
    return cohort, item


def _hashes(cohort, item):
    return {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in
            (cohort.reference.fasta, cohort.reference.fai, item.source_path, item.index_path)}


def _build(root, cohort, item):
    return tracks.build_tr_cohort_tables(cohort, (item,), locus_table=root / 'out/loci.tsv',
        sample_matrix=root / 'out/matrix.tsv', scratch_database=root / 'out/scratch.sqlite')


def _rows(path):
    with path.open(encoding='utf-8', newline='') as source:
        return list(csv.DictReader(source, delimiter='\t'))


def _sv_input(root):
    source = root / 'S1.sv.vcf.gz'
    with gzip.open(source, 'wt', encoding='utf-8', newline='') as out:
        out.write('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
                  '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV class">\n'
                  '##INFO=<ID=END,Number=1,Type=Integer,Description="Endpoint">\n'
                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
                  'chr1\t10\td1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=20\tGT\t0/1\n')
    index = Path(str(source) + '.tbi')
    index.write_bytes(b'SYNTHETIC_NONFUNCTIONAL_INDEX_NOT_A_TABIX_TEST')
    return CohortSampleInput('S1', SampleCallState.CALLED, source, index, reference_build='GRCh38')


@pytest.mark.parametrize('end', [None, '', '.', 'abc', '20.0', '10,20', ' 20', '20 ',
                                      '-1', '+20', '0', '9', '9223372036854775808'])
def test_invalid_or_missing_end_fails_without_coordinate_imputation(tmp_path, end):
    cohort, item = _fixture(tmp_path, end=end)
    before = _hashes(cohort, item)
    with pytest.raises(OutputValidationError, match='END'):
        _build(tmp_path, cohort, item)
    assert not (tmp_path / 'out/loci.tsv').exists()
    assert not (tmp_path / 'out/matrix.tsv').exists()
    assert before == _hashes(cohort, item)
    assert list(tmp_path.rglob('.hifivar-cohort-*'))


@pytest.mark.parametrize('pos', ['', '.', 'bad', '10.0', '-1', '0', '+10', ' 10',
                                     '\u0661\u0660', '9223372036854775808'])
def test_invalid_pos_uses_typed_output_error(tmp_path, pos):
    cohort, item = _fixture(tmp_path, pos=pos)
    before = _hashes(cohort, item)
    with pytest.raises(OutputValidationError, match='POS'):
        _build(tmp_path, cohort, item)
    assert not (tmp_path / 'out/loci.tsv').exists()
    assert not (tmp_path / 'out/matrix.tsv').exists()
    assert before == _hashes(cohort, item)


@pytest.mark.parametrize('trid', [None, '', '.', ' ', ' . ', 'L1,L2', '.,.'])
def test_missing_trid_cannot_be_a_real_locus_key(tmp_path, trid):
    cohort, item = _fixture(tmp_path, trid=trid)
    with pytest.raises(OutputValidationError, match='TRID'):
        _build(tmp_path, cohort, item)
    assert not (tmp_path / 'out/loci.tsv').exists()
    assert not (tmp_path / 'out/matrix.tsv').exists()


@pytest.mark.parametrize('info', [
    'TRID=L1;MOTIFS=CAG;END=20;END=30',
    'TRID=L1;TRID=L2;MOTIFS=CAG;END=20',
    'TRID=L1;MOTIFS=CAG;END',
    'TRID;MOTIFS=CAG;END=20',
])
def test_ambiguous_or_flag_locus_identity_is_rejected(tmp_path, info):
    cohort, item = _fixture(tmp_path, info=info)
    with pytest.raises(OutputValidationError, match='TRID|END'):
        _build(tmp_path, cohort, item)
    assert not (tmp_path / 'out/loci.tsv').exists()
    assert not (tmp_path / 'out/matrix.tsv').exists()


@pytest.mark.parametrize('gt', ['0/1', '1|0', '1/2', './.'])
def test_native_missing_id_column_and_gt_text_remain_valid(tmp_path, gt):
    cohort, item = _fixture(tmp_path, gt=gt)
    before = _hashes(cohort, item)
    result = _build(tmp_path, cohort, item)
    assert result.state is SampleCallState.CALLED
    locus = _rows(tmp_path / 'out/loci.tsv')[0]
    assert (locus['trid'], locus['start'], locus['end'], locus['motifs']) == ('L1', '10', '20', 'CAG')
    assert _rows(tmp_path / 'out/matrix.tsv')[0]['genotype'] == gt
    assert before == _hashes(cohort, item)
    assert not list(tmp_path.rglob('.hifivar-cohort-*'))


def _failing_connect(stage):
    original = tracks.sqlite3.connect
    class Failure(sqlite3.Connection):
        def execute(self, sql, parameters=()):
            if ((stage == 'create' and sql.startswith('CREATE TABLE calls'))
                or (stage == 'insert' and sql.startswith('INSERT INTO calls'))
                or (stage == 'query' and sql.startswith('SELECT trid,'))
                or (stage == 'sample_query' and sql.startswith('SELECT sample,'))):
                raise sqlite3.OperationalError('synthetic TR storage failure')
            return super().execute(sql, parameters)
        def commit(self):
            if stage == 'commit':
                raise sqlite3.OperationalError('synthetic TR commit failure')
            return super().commit()
        def close(self):
            super().close()
            if stage == 'close':
                raise sqlite3.OperationalError('synthetic TR close failure')
    def connect(*args, **kwargs):
        if stage == 'connect':
            raise sqlite3.OperationalError('synthetic TR open failure')
        return original(*args, factory=Failure, **kwargs)
    return connect


@pytest.mark.parametrize('stage', ['connect', 'create', 'insert', 'commit', 'query', 'sample_query', 'close'])
def test_sqlite_failures_are_typed_and_never_publish_half_tables(tmp_path, stage):
    cohort, item = _fixture(tmp_path)
    before = _hashes(cohort, item)
    with mock.patch.object(tracks.sqlite3, 'connect', _failing_connect(stage)):
        with pytest.raises(OutputValidationError, match='TR cohort'):
            _build(tmp_path, cohort, item)
    assert not (tmp_path / 'out/loci.tsv').exists()
    assert not (tmp_path / 'out/matrix.tsv').exists()
    assert before == _hashes(cohort, item)
    assert list(tmp_path.rglob('.hifivar-cohort-*'))


@pytest.mark.parametrize('stage', ['connect', 'query', 'close'])
def test_phase12_returns_prior_sv_success_and_failed_tr(tmp_path, stage):
    cohort, item = _fixture(tmp_path)
    before = _hashes(cohort, item)
    sv = _sv_input(tmp_path)
    sv_before = {path: hashlib.sha256(path.read_bytes()).hexdigest()
                 for path in (sv.source_path, sv.index_path)}
    def healthy_sv():
        return tracks.build_sv_cohort_tables(cohort, (sv,),
            site_table=tmp_path / 'sv/sites.tsv', sample_matrix=tmp_path / 'sv/matrix.tsv')
    def failed_tr():
        with mock.patch.object(tracks.sqlite3, 'connect', _failing_connect(stage)):
            return _build(tmp_path, cohort, item)
    report = run_phase12(cohort,
        enabled_tracks={CohortTrack.SV: True, CohortTrack.TR: True},
        operations={CohortTrack.SV: healthy_sv, CohortTrack.TR: failed_tr})
    assert report.successful_tracks == (CohortTrack.SV,)
    assert report.failed_tracks == (CohortTrack.TR,)
    assert (tmp_path / 'sv/sites.tsv').is_file() and (tmp_path / 'sv/matrix.tsv').is_file()
    assert not (tmp_path / 'out/loci.tsv').exists() and not (tmp_path / 'out/matrix.tsv').exists()
    assert before == _hashes(cohort, item)
    assert sv_before == {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in sv_before}


def test_percent_encoded_identifier_is_preserved_without_decoding(tmp_path):
    cohort, item = _fixture(tmp_path, trid='L1%2CL2')
    before = _hashes(cohort, item)
    _build(tmp_path, cohort, item)
    assert _rows(tmp_path / 'out/loci.tsv')[0]['trid'] == 'L1%2CL2'
    assert _rows(tmp_path / 'out/matrix.tsv')[0]['trid'] == 'L1%2CL2'
    assert before == _hashes(cohort, item)


def test_duplicate_real_trid_retains_sqlite_cause_and_no_formal_outputs(tmp_path):
    cohort, item = _fixture(tmp_path)
    with gzip.open(item.source_path, 'at', encoding='utf-8') as out:
        out.write('chr1\t10\t.\tA\tC,CC\t.\tPASS\tMOTIFS=CAG;TRID=L1;END=20\tGT\t0/1\n')
    before = _hashes(cohort, item)
    with pytest.raises(OutputValidationError, match='Duplicate TR locus') as raised:
        _build(tmp_path, cohort, item)
    assert isinstance(raised.value.__cause__, sqlite3.IntegrityError)
    assert not (tmp_path / 'out/loci.tsv').exists()
    assert not (tmp_path / 'out/matrix.tsv').exists()
    assert before == _hashes(cohort, item)


def test_close_failure_preserves_primary_validation_error(tmp_path, caplog):
    cohort, item = _fixture(tmp_path, end='.')
    before = _hashes(cohort, item)
    with mock.patch.object(tracks.sqlite3, 'connect', _failing_connect('close')):
        with pytest.raises(OutputValidationError, match='END') as raised:
            _build(tmp_path, cohort, item)
    assert raised.value.__cause__ is None
    assert 'also failed to close' in caplog.text
    assert not (tmp_path / 'out/loci.tsv').exists()
    assert not (tmp_path / 'out/matrix.tsv').exists()
    assert before == _hashes(cohort, item)


def test_close_only_failure_is_not_hidden_by_callers_active_exception(tmp_path):
    cohort, item = _fixture(tmp_path)
    before = _hashes(cohort, item)
    try:
        raise ValueError('unrelated caller exception')
    except ValueError:
        with mock.patch.object(tracks.sqlite3, 'connect', _failing_connect('close')):
            with pytest.raises(OutputValidationError, match='Unable to close TR cohort') as raised:
                _build(tmp_path, cohort, item)
    assert isinstance(raised.value.__cause__, sqlite3.OperationalError)
    assert not (tmp_path / 'out/loci.tsv').exists()
    assert not (tmp_path / 'out/matrix.tsv').exists()
    assert before == _hashes(cohort, item)


@pytest.mark.parametrize('error_type', [sqlite3.OperationalError, KeyboardInterrupt, SystemExit])
def test_double_failure_preserves_primary_cause_or_process_control(tmp_path, error_type, caplog):
    cohort, item = _fixture(tmp_path)
    before = _hashes(cohort, item)
    original_connect = tracks.sqlite3.connect
    primary = error_type('primary interrupted operation')
    class Broken(sqlite3.Connection):
        def execute(self, sql, parameters=()):
            if sql.startswith('INSERT INTO calls'):
                raise primary
            return super().execute(sql, parameters)
        def close(self):
            super().close()
            raise sqlite3.OperationalError('secondary close failure')
    def connect(*args, **kwargs):
        return original_connect(*args, factory=Broken, **kwargs)
    def operation():
        return _build(tmp_path, cohort, item)
    with mock.patch.object(tracks.sqlite3, 'connect', connect):
        if error_type is sqlite3.OperationalError:
            with pytest.raises(OutputValidationError, match='Unable to build TR cohort') as raised:
                operation()
            assert raised.value.__cause__ is primary
        else:
            with pytest.raises(error_type) as raised:
                run_phase12(cohort, enabled_tracks={CohortTrack.TR: True},
                            operations={CohortTrack.TR: operation})
            assert raised.value is primary
    assert 'secondary close failure' in caplog.text
    assert not (tmp_path / 'out/loci.tsv').exists()
    assert not (tmp_path / 'out/matrix.tsv').exists()
    assert before == _hashes(cohort, item)
