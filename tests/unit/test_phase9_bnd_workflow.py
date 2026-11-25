"""Opt-in orchestration tests, not native indexing or truth benchmarks."""
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock
import gzip

import pytest

from hifivar.assembly_sv import SVEvidenceSource
from hifivar.bnd_representation import BndRepresentationResult
from hifivar.breakends import BndPolicy, POLICY, POLICY_SHA256
from hifivar.config import load_config
from hifivar.context import AnalysisContext
from hifivar.exceptions import ConfigurationError, OutputValidationError, InputValidationError
from hifivar.harmonization import EvidenceRunStatus, HarmonizedSvArtifact, SVEvidenceSourceArtifact
from hifivar.jasmine import JasmineCommandPlan, JasmineResult, JasmineResultStatus
from hifivar.phase9 import Phase9Settings, phase9_bnd_paths, run_phase9
from hifivar.reference import ReferenceGenome
from hifivar.sample import InputDataset, Sample
from hifivar.sample_sheet import SampleRecord

DEFAULT = Path(__file__).resolve().parents[2] / 'src/hifivar/resources/configs/default.yaml'


def policy_config():
    return BndPolicy(POLICY['policy_id'], 1, POLICY_SHA256, POLICY['mode']).to_dict()


@pytest.fixture
def fixture(tmp_path):
    fasta = tmp_path / 'ref.fa'
    fasta.write_text('>chr1\n' + 'A' * 1000 + '\n')
    Path(str(fasta) + '.fai').write_text('chr1\t1000\t6\t1000\t1001\n')
    bam = tmp_path / 'S1.bam'
    bam.write_bytes(b'SYNTHETIC')
    source = tmp_path / 'S1.sawfish.sv.vcf.gz'
    source.write_bytes(gzip.compress(b'# synthetic source\n'))
    Path(str(source) + '.tbi').write_bytes(b'SYNTHETIC INDEX')
    context = AnalysisContext(ReferenceGenome.from_fasta(fasta, build='GRCh38'),
        (SampleRecord(Sample('S1', InputDataset.from_files([bam]))),), {})
    evidence = SVEvidenceSourceArtifact('S1', SVEvidenceSource.READ, 'sawfish', source, Path(str(source) + '.tbi'), EvidenceRunStatus.COMPLETED)
    cfg = deepcopy(load_config(DEFAULT).to_dict())
    cfg['sv']['harmonization']['enabled'] = True
    out = tmp_path / 'results'
    calls = []

    def native(request, *, dry_run=False, **kwargs):
        calls.append(('jasmine', request))
        commands = (JasmineCommandPlan('jasmine', ('jasmine', 'file_list=declared', 'out_file=declared')),)
        if dry_run:
            return JasmineResult(request, JasmineResultStatus.PLANNED, commands, {}, 0, None)
        request.output_vcf.parent.mkdir(parents=True, exist_ok=True)
        request.work_directory.mkdir(parents=True, exist_ok=True)
        for p in (request.output_vcf, request.output_index, request.evidence_table):
            p.write_bytes(b'SYNTHETIC NATIVE OUTPUT')
        raw = request.work_directory / 'S1.jasmine.raw.vcf'
        raw.write_bytes(b'SYNTHETIC RAW')
        artifact = HarmonizedSvArtifact('S1', request.reference.fasta, 'GRCh38', request.output_vcf,
            request.output_index, request.evidence_table, request.sources, '1.1.5', tuple(x.args for x in commands), intermediate_files=(raw,))
        return JasmineResult(request, JasmineResultStatus.COMPLETED, commands, {'jasmine': '1.1.5'}, .1, artifact)

    def compare(request, **kwargs):
        calls.append(('truvari', request))
        return SimpleNamespace(to_dict=lambda: {'scope': 'SYNTHETIC_CONCORDANCE_NOT_TRUTH'})

    def derive(request):
        calls.append(('bnd', request))
        paths = phase9_bnd_paths(out, 'S1')
        paths['directory'].mkdir(parents=True)
        return BndRepresentationResult(paths['directory'], paths['manifest'], paths['vcf'], paths['index'], paths['ledger'], paths['junctions'], {})

    return SimpleNamespace(context=context, evidence={'S1': (evidence,)}, config=cfg, output=out,
        work=tmp_path / 'work', calls=calls, jasmine=SimpleNamespace(run=Mock(side_effect=native)),
        truvari=SimpleNamespace(run=Mock(side_effect=compare)), bnd=SimpleNamespace(run=Mock(side_effect=derive)))


def execute(f, **kwargs):
    return run_phase9(f.context, f.evidence, output_directory=f.output, work_directory=f.work,
        config=f.config, jasmine_wrapper=f.jasmine, truvari_wrapper=f.truvari, bnd_wrapper=f.bnd, **kwargs)


def test_default_disabled_keeps_native_comparisons_and_creates_no_bundle(fixture):
    f = fixture
    report = execute(f)
    assert [name for name, _ in f.calls] == ['jasmine', 'truvari']
    assert f.calls[-1][1].base_vcf == f.calls[0][1].output_vcf
    assert report.sample_results[0].bnd_representation is None
    assert report.sample_results[0].to_dict()['bnd_representation'] is None
    assert not phase9_bnd_paths(f.output, 'S1')['directory'].exists()


def test_enabled_binds_completed_native_result_and_preserves_native_role(fixture):
    f = fixture
    f.config['sv']['harmonization']['bnd_representation'] = policy_config()
    report = execute(f)
    assert [name for name, _ in f.calls] == ['jasmine', 'bnd', 'truvari']
    request = f.bnd.run.call_args.args[0]
    result = report.sample_results[0]
    assert request.native_result is result.jasmine
    assert request.output_directory == phase9_bnd_paths(f.output, 'S1')['directory']
    assert request.policy.policy_sha256 == POLICY_SHA256
    assert f.calls[-1][1].base_vcf == result.jasmine.request.output_vcf
    assert f.calls[-1][1].base_vcf != result.bnd_representation.derived_vcf
    assert result.to_dict()['bnd_representation']['native_output_role_unchanged'] is True
    assert report.settings.to_dict()['bnd_representation']['policy_sha256'] == POLICY_SHA256


def test_enabled_dry_run_is_only_an_explicit_plan(fixture):
    f = fixture
    f.config['sv']['harmonization']['bnd_representation'] = policy_config()
    report = execute(f, dry_run=True)
    f.bnd.run.assert_not_called()
    item = report.sample_results[0].to_dict()['bnd_representation']
    assert item['status'] == 'planned' and item['executed'] is False
    assert item['consumer_opt_in_required'] is True
    assert not f.output.exists()


@pytest.mark.parametrize('overwrite', [True, False])
def test_existing_bundle_fails_before_any_native_tool(fixture, overwrite):
    f = fixture
    f.config['sv']['harmonization'].update(bnd_representation=policy_config(), overwrite=overwrite)
    directory = phase9_bnd_paths(f.output, 'S1')['directory']
    directory.mkdir(parents=True)
    marker = directory / 'protected'
    marker.write_bytes(b'PRESERVE')
    with pytest.raises(OutputValidationError, match='fresh per-sample'):
        execute(f)
    assert f.calls == [] and marker.read_bytes() == b'PRESERVE'


def test_bnd_failure_does_not_schedule_comparisons_or_remove_native_outputs(fixture):
    f = fixture
    f.config['sv']['harmonization']['bnd_representation'] = policy_config()
    f.bnd.run.side_effect = OutputValidationError('synthetic derivation failure')
    with pytest.raises(OutputValidationError, match='synthetic'):
        execute(f)
    f.truvari.run.assert_not_called()
    assert (f.output / 'S1/S1.harmonized.sv.vcf.gz').read_bytes() == b'SYNTHETIC NATIVE OUTPUT'


def test_later_sample_existing_bundle_blocks_all_samples_before_tools(fixture):
    f = fixture
    bam2 = f.work.parent / 'S2.bam'
    bam2.write_bytes(b'SYNTHETIC SECOND INPUT')
    f.context = AnalysisContext(f.context.reference,
        (*f.context.samples, SampleRecord(Sample('S2', InputDataset.from_files([bam2])))), {})
    f.config['sv']['harmonization']['bnd_representation'] = policy_config()
    phase9_bnd_paths(f.output, 'S2')['directory'].mkdir(parents=True)
    with pytest.raises(OutputValidationError, match='fresh per-sample'):
        execute(f)
    assert f.calls == [] and not (f.output / 'S1').exists()


@pytest.mark.parametrize('value', [None, True, 'on', {'enabled': 'true'}, {'enabled': False, 'extra': 1},
                                  {'enabled': True}, {**policy_config(), 'policy_sha256': 'wrong'}])
def test_strict_nested_policy_validation(fixture, value):
    f = fixture
    f.config['sv']['harmonization']['bnd_representation'] = value
    with pytest.raises(ConfigurationError):
        Phase9Settings.from_config(f.config)
    assert f.calls == []


def test_config_loader_requires_parent_and_exact_policy(tmp_path):
    import yaml
    path = tmp_path / 'config.yaml'
    path.write_text(yaml.safe_dump({'sv': {'harmonization': {'enabled': True, 'bnd_representation': policy_config()}}}))
    assert Phase9Settings.from_config(load_config(DEFAULT, user_config=path)).bnd_policy is not None
    path.write_text(yaml.safe_dump({'sv': {'harmonization': {'enabled': False, 'bnd_representation': policy_config()}}}))
    with pytest.raises(ConfigurationError, match='requires'):
        load_config(DEFAULT, user_config=path)


@pytest.mark.parametrize('sample', ['../escape', '..', '', 'sample/name'])
def test_bundle_path_rejects_unsafe_sample_identifiers(tmp_path, sample):
    with pytest.raises(InputValidationError):
        phase9_bnd_paths(tmp_path, sample)
