"""Public batch template contracts; parsing is not execution or accuracy evidence."""
from copy import deepcopy
import csv
import hashlib
from pathlib import Path
import re

import pytest
import yaml

from hifivar.command import CommandRunner
from hifivar.config import load_config
from hifivar.downstream_plan import read_downstream_plan
from hifivar.exceptions import InputValidationError

ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = ROOT / 'examples/batch_multivariant'
DEFAULT = ROOT / 'src/hifivar/resources/configs/default.yaml'
SAMPLES = ('SAMPLE_A', 'SAMPLE_B')


def test_example_merges_with_current_schema_and_preserves_safe_defaults():
    c = load_config(DEFAULT, user_config=EXAMPLE / 'config.yaml').to_dict()
    assert c['runtime']['threads'] == 8 and c['workflow']['memory_mb'] == 262144
    assert c['small']['region'] is None and c['small']['max_concurrent_samples'] == 1
    assert c['small']['execution_mode'] == 'apptainer' and c['small']['model_type'] == 'PACBIO'
    assert c['sv']['sawfish']['disable_cnv'] is False
    assert all(c['sv'][n]['enabled'] for n in ['sawfish', 'sniffles2', 'pbsv', 'cutesv'])
    assert c['sv']['harmonization']['enabled']
    assert c['phasing']['sv_caller'] == 'pbsv' and c['phasing']['include_tr'] is True
    assert c['tr']['karyotype'] == 'auto'
    assert c['workflow']['downstream_source_plan'] == 'source_plan.yaml'
    assert c['cohort']['input_manifest'] is None
    assert c['runtime']['tmpdir'] is c['logging']['file'] is None
    assert c['paths']['workdir'] is c['paths']['outdir'] is None
    assert all(c[n]['enabled'] is False for n in ['alignment', 'assembly', 'assembly_sv', 'annotation', 'benchmark', 'review', 'low_depth'])
    tracks = [c['small'], c['tr'], c['phasing'], c['sv']['harmonization']]
    tracks += [c['sv'][n] for n in ['sawfish', 'sniffles2', 'pbsv', 'cutesv']]
    assert all(t['threads'] <= c['runtime']['threads'] and t['memory_mb'] <= c['workflow']['memory_mb'] for t in tracks)
    assert c['cohort']['small_variants']['memory_gb'] * 1024 <= c['workflow']['memory_mb']
    assert all(not t['overwrite'] for t in [c['small'], c['sv'], c['sv']['harmonization'], c['tr'], c['phasing'], c['cohort']])


def test_sample_sheet_has_real_tabs_and_declared_sex_not_inferred():
    text = (EXAMPLE / 'samples.tsv').read_text(encoding='utf-8')
    assert '\\t' not in text
    rows = list(csv.DictReader(text.splitlines(), delimiter='\t'))
    assert list(rows[0]) == ['sample_id', 'input', 'input_type', 'sex']
    assert tuple(r['sample_id'] for r in rows) == SAMPLES
    assert [r['sex'] for r in rows] == ['female', 'male']
    assert all(r['input_type'] == 'bam' and r['input'].startswith('/work/inputs/') for r in rows)


@pytest.fixture
def configured_plan(tmp_path):
    c = load_config(DEFAULT, user_config=EXAMPLE / 'config.yaml').to_dict()
    catalog = tmp_path / 'catalog.bed'
    catalog.write_text('chr1\t0\t6\tID=L1;MOTIFS=AC;STRUC=(AC)n\n', encoding='utf-8')
    c['tr']['catalog'] = str(catalog)
    document = yaml.safe_load((EXAMPLE / 'source_plan.yaml').read_text(encoding='utf-8'))
    path = tmp_path / 'source_plan.yaml'
    path.write_text(yaml.safe_dump(document, sort_keys=False), encoding='utf-8')
    return c, document, path, catalog


def test_generated_template_uses_native_roles_and_actual_catalog_without_commands(configured_plan, monkeypatch):
    c, _, path, catalog = configured_plan
    original = deepcopy(c)
    before = path.read_bytes()
    monkeypatch.setattr(CommandRunner, 'run', lambda *a, **k: pytest.fail('Parsing cannot run tools'))
    p = read_downstream_plan(path, sample_ids=SAMPLES, reference_build='GRCh38', config=c)
    assert len(p.cohort) == 6 and p.annotation == ()
    for sample in SAMPLES:
        rows = {r.track: r for r in p.cohort if r.sample_id == sample}
        assert set(rows) == {'small_variants', 'sv', 'tr'}
        assert (rows['small_variants'].producer, rows['small_variants'].output_role) == ('small.deepvariant', 'gvcf')
        assert rows['sv'].producer == 'sv.harmonization'
        assert rows['tr'].catalog_id == 'sha256:' + hashlib.sha256(catalog.read_bytes()).hexdigest()
        assert all(r.origin == 'GENERATED' and r.path is r.index_path is r.source_version is None for r in rows.values())
    assert path.read_bytes() == before and c == original
    assert p.sha256 == hashlib.sha256(before).hexdigest()


@pytest.mark.parametrize('mutation', ['small_vcf', 'future_path', 'unknown_sample', 'missing_row', 'duplicate_row', 'catalog_mismatch', 'disabled_tr', 'wrong_build'])
def test_bad_template_substitutions_fail_closed(configured_plan, mutation):
    c, d, path, _ = configured_plan
    if mutation == 'small_vcf': d['cohort'][0]['output_role'] = 'vcf'
    if mutation == 'future_path': d['cohort'][0]['path'] = 'future.g.vcf.gz'
    if mutation == 'unknown_sample': d['cohort'][0]['sample_id'] = 'UNDECLARED'
    if mutation == 'missing_row': d['cohort'].pop()
    if mutation == 'duplicate_row': d['cohort'].append(deepcopy(d['cohort'][0]))
    if mutation == 'catalog_mismatch': d['cohort'][2]['catalog_id'] = 'sha256:' + '0' * 64
    if mutation == 'disabled_tr': c['tr']['enabled'] = False
    if mutation == 'wrong_build': d['cohort'][0]['reference_build'] = 'GRCh37'
    path.write_text(yaml.safe_dump(d, sort_keys=False), encoding='utf-8')
    with pytest.raises(InputValidationError):
        read_downstream_plan(path, sample_ids=SAMPLES, reference_build='GRCh38', config=c)


def test_guide_distinguishes_template_execution_accuracy_and_publication():
    text = (ROOT / 'docs/batch_multivariant_run.md').read_text(encoding='utf-8')
    for required in ['RC5 release source', 'biological N=1', 'not two', 'not one all-module DAG',
                     'separate explicit command', 'not new cross-sample calling', 'GENERATED_SOURCE_PLAN_ONLY',
                     'not permission', 'tabix >=1.23.1', '192 GiB', 'not measured requirements', 'NOT_REVIEWED']:
        assert required in text, required
    assert 'hifivar --config config.yaml --dry-run run' in text
    assert 'hifivar --config config.yaml run' in text
    assert 'hifivar report --run-dir' in text
    assert '--rerun-incomplete' not in text
    assert '中文' in text


@pytest.mark.parametrize('relative', ['README.md', 'README.zh-CN.md', 'ARCHITECTURE.md', 'docs/capabilities.md', 'docs/joint_phasing.md', 'docs/indexed_bam_run.md'])
def test_current_entry_docs_link_complete_batch_guide(relative):
    assert 'batch_multivariant_run.md' in (ROOT / relative).read_text(encoding='utf-8')


def test_new_guide_and_example_local_links_resolve_without_private_paths():
    for path in [ROOT/'docs/batch_multivariant_run.md', EXAMPLE/'README.md']:
        text = path.read_text(encoding='utf-8')
        for target in re.findall(r'\]\(([^)]+)\)', text):
            if '://' not in target and not target.startswith('#'):
                assert (path.parent / target.split('#', 1)[0]).resolve().is_file(), target
        assert not re.search(r'(?i)(?:[a-z]:[\\/]|/home/|/lustre/|/data/[^/ \n]+/1_projects/)', text)
