"""Test archival figure boundaries without requiring plotting dependencies."""

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'project_visuals/scripts/benchmark_qualification.py'
SPEC = importlib.util.spec_from_file_location('figure_qualification', SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_archival_sawfish_excluded_and_data_immutable():
    source = ROOT / 'project_visuals/data/hg002_cmrg_sv_benchmark.tsv'
    before = source.read_bytes()
    callers, manifest = MODULE.select_cmrg_callers(
        source, source.with_name('cmrg_figure_qualifications.json'),
    )
    assert callers == ['Sniffles2', 'pbsv', 'cuteSV']
    assert 'UnknownSampleName' in manifest['excluded_callers']['Sawfish']
    assert source.read_bytes() == before
    assert 'not whole-genome' in manifest['scope']


@pytest.fixture
def review_files(tmp_path):
    source = tmp_path / 'source.tsv'
    source.write_text('caller\tversion\nS1\t1\nS2\t1\n', encoding='utf-8')
    policy = tmp_path / 'policy.json'
    payload = {
        'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'scope': 'bounded', 'interpretation': 'qualified historical',
        'callers': {c: {'figure_disposition': 'INCLUDE_QUALIFIED_HISTORICAL'}
                    for c in ['S1', 'S2']},
    }
    policy.write_text(json.dumps(payload), encoding='utf-8')
    return source, policy, payload


def test_changed_source_fails(review_files):
    source, policy, _ = review_files
    source.write_bytes(source.read_bytes() + b'\n')
    with pytest.raises(ValueError, match='source changed'):
        MODULE.select_cmrg_callers(source, policy)


def test_unknown_caller_fails(review_files):
    source, policy, payload = review_files
    payload['callers'].pop('S2')
    policy.write_text(json.dumps(payload), encoding='utf-8')
    with pytest.raises(ValueError, match='membership'):
        MODULE.select_cmrg_callers(source, policy)


def test_unreviewed_disposition_fails(review_files):
    source, policy, payload = review_files
    payload['callers']['S1']['figure_disposition'] = 'ASSUME_VALID'
    policy.write_text(json.dumps(payload), encoding='utf-8')
    with pytest.raises(ValueError, match='Unreviewed'):
        MODULE.select_cmrg_callers(source, policy)
