"""Launcher routing contracts; synthetic inputs are not accuracy evidence."""
from pathlib import Path
import runpy
from types import SimpleNamespace

import pytest
import yaml

from hifivar import launcher, preflight
from hifivar.config import HiFiVarConfig, load_config
from hifivar.exceptions import ConfigurationError, InputValidationError
from hifivar.harmonization import EvidenceRunStatus

ROOT = Path(__file__).parents[2]
DEFAULT = ROOT / 'src/hifivar/resources/configs/default.yaml'


@pytest.fixture
def case(tmp_path):
    (tmp_path / 'ref.fa').write_text('>chr1\n' + 'A' * 1000 + '\n', encoding='utf8')
    (tmp_path / 'ref.fa.fai').write_text('chr1\t1000\t6\t1000\t1001\n', encoding='utf8')
    for sample in ('S2', 'S1'):
        (tmp_path / f'{sample}.bam').write_bytes(b'SYNTHETIC_NOT_BAM')
        (tmp_path / f'{sample}.bam.bai').write_bytes(b'SYNTHETIC_NOT_INDEX')
    (tmp_path / 'samples.tsv').write_text('sample_id\tinput\tinput_type\nS2\tS2.bam\tbam\nS1\tS1.bam\tbam\n', encoding='utf8')
    path = tmp_path / 'config.yaml'
    path.write_text(yaml.safe_dump({
        'reference': {'fasta': 'ref.fa', 'build': 'CUSTOM_TEST'},
        'samples': {'sheet': 'samples.tsv'}, 'runtime': {'threads': 4},
        'workflow': {'memory_mb': 2000}, 'sv': {'enabled': True,
            'sawfish': {'enabled': True, 'threads': 2, 'memory_mb': 1000},
            'sniffles2': {'enabled': True, 'threads': 2, 'memory_mb': 1000},
            'pbsv': {'enabled': False}, 'cutesv': {'enabled': False},
            'harmonization': {'enabled': True, 'threads': 2, 'memory_mb': 1000,
                             'jasmine_executable': 'tools/jasmine'}}}), encoding='utf8')
    return load_config(DEFAULT, user_config=path), tmp_path


def changed(cfg, h=None, sv=None):
    data = cfg.to_dict()
    data['sv'].update(sv or {})
    data['sv']['harmonization'].update(h or {})
    return HiFiVarConfig(data, cfg.sources)


def test_launch_generated_per_sample_harmonization_no_side_effects(case):
    cfg, root = case
    before = cfg.to_dict()
    effective, context, sheet, _ = launcher.prepare_launch(cfg, root / 'planned')
    assert context.sample_ids == ('S2', 'S1')
    assert sheet.index('S2\t') < sheet.index('S1\t')
    assert effective['sv']['harmonization']['jasmine_executable'] == str(root / 'tools/jasmine')
    assert [name for name, _ in launcher._enabled_tracks(effective)] == ['sv.sawfish', 'sv.sniffles2', 'sv.harmonization']
    assert cfg.to_dict() == before and not (root / 'planned').exists()


@pytest.mark.parametrize('h,sv,match', [
    ({}, {'enabled': False}, 'requires generated or external'),
    ({'input_vcfs': {'S2': []}}, {}, 'must be a mapping'),
    ({'input_vcfs': {'UNKNOWN': {'sawfish': 'elsewhere.vcf.gz'}}}, {}, 'undeclared sample'),
    ({'overwrite': True}, {}, 'overwrite=false'),
    ({'threads': 5}, {}, 'threads exceeds'),
    ({'memory_mb': 2001}, {}, 'memory_mb exceeds'),
    ({'distance_type': 'unknown'}, {}, 'linear or nonlinear'),
    ({}, {name: {'enabled': False} for name in ('sawfish', 'sniffles2', 'pbsv', 'cutesv')}, 'at least one'),
])
def test_invalid_harmonization_fails_before_run_creation(case, h, sv, match):
    cfg, root = case
    with pytest.raises(ConfigurationError, match=match):
        launcher.prepare_launch(changed(cfg, h, sv), root / 'planned')
    assert not (root / 'planned').exists()


def test_disabled_harmonization_has_no_dependency_group(case):
    cfg, _ = case
    groups = preflight._tool_groups(changed(cfg, {'enabled': False}).to_dict(), object())
    assert not {'jasmine', 'truvari'} & {name for name, _, _ in groups}


def test_enabled_harmonization_dependency_groups(case):
    cfg, _ = case
    groups = preflight._tool_groups(cfg.to_dict(), object())
    assert {'jasmine', 'truvari'} <= {name for name, _, _ in groups}
    jasmine = next(values for name, values, _ in groups if name == 'jasmine')
    assert set(jasmine) == {'jasmine', 'bgzip', 'tabix'}


@pytest.mark.parametrize('sample', ['S2', 'S1'])
def test_phase9_bridge_owns_exactly_one_sample(case, monkeypatch, sample):
    from hifivar import phase9
    cfg, root = case
    data, _, _, _ = launcher.prepare_launch(cfg, root / 'planned')
    data['samples']['sheet'] = str(root / 'samples.tsv')
    received = []

    def fake_phase9(context, evidence, **kwargs):
        assert context.sample_ids == (sample,)
        assert set(evidence) == {sample}
        assert context.source_sample_sheet == root / 'samples.tsv'
        sources = evidence[sample]
        assert {s.caller for s in sources if s.status is EvidenceRunStatus.COMPLETED} == {'sawfish', 'sniffles2'}
        assert all(s.sample_id == sample for s in sources)
        received.append(context)
        # This test checks sample ownership only, not native execution receipts.
        return SimpleNamespace(
            write_json=lambda *args, **kwargs: None,
            sample_results=(SimpleNamespace(
                jasmine=SimpleNamespace(tool_versions={'jasmine': 'test-only'}),
                truvari=tuple(SimpleNamespace(version='test-only') for _ in range(2)),
            ),),
        )

    monkeypatch.setattr(phase9, 'run_phase9', fake_phase9)
    paths = []
    for caller in ('sawfish', 'sniffles2'):
        path = root / 'results/sv' / f'{sample}.{caller}.sv.vcf.gz'
        paths.extend((str(path), str(path) + '.tbi'))
    sm = SimpleNamespace(wildcards=SimpleNamespace(sample=sample), config=data,
        input=SimpleNamespace(sources=paths),
        output=SimpleNamespace(vcf=str(root / 'results/sv_harmonized' / sample / f'{sample}.harmonized.sv.vcf.gz'),
                               index=str(root / 'results/sv_harmonized' / sample / f'{sample}.harmonized.sv.vcf.gz.tbi'),
                               evidence=str(root / 'results/sv_harmonized' / sample / f'{sample}.sv.evidence.tsv'),
                               provenance=str(root / 'never_written.json')),
        params=SimpleNamespace(work_root=str(root / 'work')))
    runpy.run_path(str(ROOT / 'workflow/scripts/run_harmonization.py'), init_globals={'snakemake': sm})
    assert len(received) == 1 and not (root / 'never_written.json').exists()


def test_bridge_unknown_sample_rejected_before_tools(case):
    cfg, root = case
    data, _, _, _ = launcher.prepare_launch(cfg, root / 'planned')
    data['samples']['sheet'] = str(root / 'samples.tsv')
    with pytest.raises(InputValidationError, match='not unique'):
        runpy.run_path(str(ROOT / 'workflow/scripts/run_harmonization.py'), init_globals={
            'snakemake': SimpleNamespace(wildcards=SimpleNamespace(sample='UNKNOWN'), config=data)})


def test_external_assembly_never_infers_haplotype_roles_from_filename(case, monkeypatch):
    from hifivar import phase9
    cfg, root=case
    data=cfg.to_dict()
    data['reference']['fasta']=str(root/'ref.fa')
    data['samples']['sheet']=str(root/'samples.tsv')
    data['paths']['outdir']=str(root/'results')
    data['sv']['enabled']=False
    source=root/'S1.pav.assembly.sv.vcf.gz'
    data['sv']['harmonization']['input_vcfs']={'S1':{'pav':str(source)}}
    received=[]
    def inspect(context,evidence,**kwargs):
        actual=next(x for x in evidence['S1'] if x.caller=='pav')
        assert actual.haplotypes==(), 'A filename is not evidence for two haplotype inputs'
        received.append(actual)
        return SimpleNamespace(write_json=lambda *a,**k:None,
            sample_results=(SimpleNamespace(jasmine=SimpleNamespace(tool_versions={'jasmine':'fixture-only'}),
                truvari=(SimpleNamespace(version='fixture-only'),)),))
    monkeypatch.setattr(phase9,'run_phase9',inspect)
    monkeypatch.delenv('HIFIVAR_ATTEMPT_CONTEXT',raising=False)
    output=root/'results/sv_harmonized/S1'
    sm=SimpleNamespace(wildcards=SimpleNamespace(sample='S1'),config=data,
        input=SimpleNamespace(sources=[str(source),str(source)+'.tbi']),
        output=SimpleNamespace(vcf=str(output/'S1.harmonized.sv.vcf.gz'),index=str(output/'S1.harmonized.sv.vcf.gz.tbi'),
            evidence=str(output/'S1.sv.evidence.tsv'),provenance=str(output/'S1.phase9.provenance.json')),
        params=SimpleNamespace(work_root=str(root/'work')))
    runpy.run_path(str(ROOT/'workflow/scripts/run_harmonization.py'),init_globals={'snakemake':sm})
    assert len(received)==1 and not output.exists()
