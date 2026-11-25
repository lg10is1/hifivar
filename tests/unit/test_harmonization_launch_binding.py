"""Policy fingerprint must include executable scientific policy resources."""
from hifivar.launcher import runtime_identity

def test_workflow_policy_bytes_change_runtime_identity(tmp_path):
    path = tmp_path / 'representation_policy.json'
    path.write_text('{"policy": 1}', encoding='utf-8')
    before = runtime_identity(tmp_path)
    path.write_text('{"policy": 2}', encoding='utf-8')
    assert runtime_identity(tmp_path)['payload_sha256'] != before['payload_sha256']

import json
import os
import sys
from pathlib import Path
import pytest
import yaml
from hifivar import launcher
from hifivar.attempt_receipts import (ATTEMPT_CONTEXT_ENV, artifact_identity,
    attempt_step, create_attempt_context, read_step_receipt)
from hifivar.command import CommandRunner
from hifivar.exceptions import ConfigurationError, InputValidationError
from hifivar.harmonization_inputs import harmonization_sources


def test_external_only_is_a_real_enabled_track():
    config = {'sv': {'enabled': False, 'harmonization': {'enabled': True,
        'input_vcfs': {'S1': {'pav': '/ready/S1.pav.assembly.sv.vcf.gz'}},
        'threads': 1, 'memory_mb': 32, 'distance_type': 'linear'}},
        'runtime': {'threads': 1}, 'workflow': {'memory_mb': 32}}
    assert [name for name, _ in launcher._enabled_tracks(config)] == ['sv.harmonization']
    assert not launcher._calling_enabled(config)


def test_selection_order_and_generated_paths_are_explicit(tmp_path):
    cfg = {'sv': {'enabled': True, 'sawfish': {'enabled': True},
        'harmonization': {'input_vcfs': {'S1': {
            'svim_asm': str(tmp_path / 'S1.svim_asm.assembly.sv.vcf.gz'),
            'pbsv': str(tmp_path / 'S1.pbsv.sv.vcf.gz')}}}}}
    sources = harmonization_sources(cfg, 'S1', tmp_path / 'run')
    assert [x['caller'] for x in sources] == ['sawfish', 'pbsv', 'svim_asm']
    assert sources[0]['vcf'] == tmp_path / 'run/results/sv/S1.sawfish.sv.vcf.gz'
    assert [x['origin'] for x in sources] == ['GENERATED', 'EXTERNAL', 'EXTERNAL']


@pytest.mark.parametrize('external,match', [
    ({'other': '/ready/S1.other.sv.vcf.gz'}, 'unsupported'),
    ({'pav': '/ready/wrong.vcf.gz'}, 'must be named'),
    ({'pav': 'S1.pav.assembly.sv.vcf.gz'}, 'absolute'),
    ({'pav': ''}, 'nonempty'),
])
def test_source_selection_rejects_ambiguous_sources(tmp_path, external, match):
    with pytest.raises(ConfigurationError, match=match):
        harmonization_sources({'sv': {'harmonization': {'input_vcfs': {'S1': external}}}}, 'S1', tmp_path)


def test_caller_cannot_be_generated_and_external(tmp_path):
    with pytest.raises(ConfigurationError, match='both generated and external'):
        harmonization_sources({'sv': {'enabled': True, 'sawfish': {'enabled': True},
            'harmonization': {'input_vcfs': {'S1': {'sawfish': str(tmp_path / 'S1.sawfish.sv.vcf.gz')}}}}}, 'S1', tmp_path)


@pytest.fixture
def owned_harmonization(tmp_path, monkeypatch, request):
    # Synthetic control fixture for binding tests, not a native workflow claim.
    root = tmp_path / 'run'
    root.mkdir()
    generated = getattr(request, 'param', None) == 'generated'
    vcf = root / 'results/sv/S1.sawfish.sv.vcf.gz' if generated else tmp_path / 'S1.pav.assembly.sv.vcf.gz'
    index = Path(str(vcf) + '.tbi')
    if not generated:
        vcf.write_bytes(b'synthetic canonical PAV bytes')
        index.write_bytes(b'synthetic index')
    config = {'sv': {'enabled': False, 'harmonization': {'enabled': True,
        'input_vcfs': {'S1': {'pav': str(vcf)}}}}}
    if generated:
        config['sv'].update(enabled=True, sawfish={'enabled': True})
        config['sv']['harmonization']['input_vcfs'] = {}
    sheet = 'sample_id\tinput\nS1\t/context.fastq\n'
    (root / 'samples.tsv').write_text(sheet, encoding='utf-8')
    (root / 'effective.yaml').write_text(yaml.safe_dump(config), encoding='utf-8')
    runtime = {'payload_sha256': 'a' * 64, 'hifivar_version': 'fixture-only'}
    monkeypatch.setattr(launcher, 'runtime_identity', lambda _: runtime)
    import hifivar.package_resources
    monkeypatch.setattr(hifivar.package_resources, 'installed_workflow_root', lambda: root)
    identity = {'schema': 2, 'run_dir': str(root), 'config': config,
        'sample_sheet': sheet, 'runtime': runtime,
        'inputs': [] if generated else [artifact_identity(vcf), artifact_identity(index)]}
    manifest = {'identity': identity, 'control_sha256': {
        name: artifact_identity(root / name)['sha256'] for name in ('effective.yaml', 'samples.tsv')}}
    (root / 'launch_manifest.json').write_text(json.dumps(manifest), encoding='utf-8')
    attempt = root / 'attempts/attempt_test'
    attempt.mkdir(parents=True)
    (attempt / 'status.json').write_text(json.dumps({'status': 'DAG_RUNNING', 'dry_run': False}), encoding='utf-8')
    context = create_attempt_context(root, attempt, dry_run=False)
    monkeypatch.setenv(ATTEMPT_CONTEXT_ENV, str(context))
    return root, attempt, vcf, index


def complete_synthetic(receipt, outputs):
    # A real bounded Python process writes fixture bytes; no fake caller exit.
    environment = {k: v for k, v in os.environ.items() if k.upper() == 'SYSTEMROOT'}
    CommandRunner().run([sys.executable, '-I', '-B', '-c',
        'import pathlib,sys;[(pathlib.Path(p).parent.mkdir(parents=True,exist_ok=True), pathlib.Path(p).write_bytes(b"synthetic-output")) for p in sys.argv[1:]]',
        *map(str, outputs.values())], env=environment, inherit_env=False)
    receipt.complete(outputs=outputs, tool_versions={'synthetic-process': '1'}, validation={'validated': True})


def invocation(data):
    root, _, vcf, index = data
    return {'track': 'sv.harmonization', 'samples': ('S1',),
        'outputs': {role: root / 'results' / role for role in ('vcf', 'index', 'evidence', 'provenance')},
        'inputs': {'pav': {'vcf': vcf, 'index': index}}, 'input_manifest': root / 'effective.yaml'}


def test_external_attempt_receipt_binds_actual_sources_and_readback(owned_harmonization):
    root, attempt, vcf, index = owned_harmonization
    args = invocation(owned_harmonization)
    with attempt_step(**args) as receipt:
        complete_synthetic(receipt, args['outputs'])
    record = read_step_receipt(root, attempt.name, 'sv.harmonization', ('S1',))
    binding = record['consumer_binding']
    assert binding['source_order'] == ['pav']
    assert binding['source_origins'] == {'pav': 'EXTERNAL'}
    assert binding['source_inputs']['pav'] == {'vcf': artifact_identity(vcf), 'index': artifact_identity(index)}
    assert record['scientific_accuracy'] == 'NOT_ASSESSED'
    vcf.write_bytes(b'modified after completion')
    with pytest.raises(InputValidationError, match='source bytes differ'):
        read_step_receipt(root, attempt.name, 'sv.harmonization', ('S1',))


@pytest.mark.parametrize('change', ['no_inputs', 'wrong_caller', 'wrong_path', 'wrong_index', 'wrong_control', 'changed_bytes'])
def test_actual_source_substitution_rejected_before_execution(owned_harmonization, change):
    root, _, vcf, index = owned_harmonization
    args = invocation(owned_harmonization)
    if change == 'no_inputs': args.pop('inputs')
    elif change == 'wrong_caller': args['inputs'] = {'sawfish': args['inputs']['pav']}
    elif change == 'wrong_path': args['inputs']['pav']['vcf'] = index
    elif change == 'wrong_index': args['inputs']['pav']['index'] = vcf
    elif change == 'wrong_control': args['input_manifest'] = root / 'samples.tsv'
    else: vcf.write_bytes(b'modified before execution')
    with pytest.raises(InputValidationError):
        with attempt_step(**args): pytest.fail('must reject before commands')


def test_source_changes_during_execution_cannot_seal(owned_harmonization):
    _, _, vcf, _ = owned_harmonization
    args = invocation(owned_harmonization)
    with pytest.raises(InputValidationError, match='source bytes differ'):
        with attempt_step(**args) as receipt:
            vcf.write_bytes(b'modified during execution')
            complete_synthetic(receipt, args['outputs'])


@pytest.mark.parametrize('owned_harmonization', ['generated'], indirect=True)
def test_generated_source_requires_current_attempt_producer(owned_harmonization):
    root, attempt, vcf, index = owned_harmonization
    args = invocation(owned_harmonization)
    args['inputs'] = {'sawfish': {'vcf': vcf, 'index': index}}
    pair = {'vcf': vcf, 'index': index}
    with attempt_step('sv.sawfish', ('S1',), outputs=pair) as producer:
        complete_synthetic(producer, pair)
    with attempt_step(**args) as consumer:
        complete_synthetic(consumer, args['outputs'])
    record = read_step_receipt(root, attempt.name, 'sv.harmonization', ('S1',))
    assert record['consumer_binding']['source_origins'] == {'sawfish': 'GENERATED'}
    assert record['consumer_binding']['producer_receipts'][0]['track'] == 'sv.sawfish'
    # Loss of the precise producer receipt invalidates the consumer as well.
    producer_path = Path(record['consumer_binding']['producer_receipts'][0]['receipt']['path'])
    producer_path.write_text('{}', encoding='utf-8')
    with pytest.raises(InputValidationError):
        read_step_receipt(root, attempt.name, 'sv.harmonization', ('S1',))
