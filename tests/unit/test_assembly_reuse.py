"""Ready input validation and native producer identity are separate facts."""
import hashlib
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest
import yaml

from hifivar import launcher
from hifivar.assembly import HaplotypeAssemblyArtifact
from hifivar.assembly_reuse import HANDOFF_NAME, read_ready_handoff, workflow_assemblies, ready_work_directory
from hifivar.assembly_source_roles import assembly_role_binding, source_haplotype_binding
from hifivar.attempt_receipts import ATTEMPT_CONTEXT_ENV, attempt_step, read_step_receipt
from hifivar.command import CommandRunner
from hifivar.config import HiFiVarConfig, load_config, validate_config
from hifivar.exceptions import ConfigurationError, InputValidationError
from hifivar.workflow_report import build_launcher_report

ROOT = Path(__file__).parents[2]
DEFAULT = ROOT / 'src/hifivar/resources/configs/default.yaml'


@pytest.fixture
def ready_case(tmp_path, monkeypatch):
    (tmp_path/'ref.fa').write_text('>chr1\nACGT\n', encoding='utf-8')
    (tmp_path/'ref.fa.fai').write_text('chr1\t4\t6\t4\t5\n', encoding='utf-8')
    (tmp_path/'reads.fastq').write_text('@m/1/ccs\nACGT\n+\nIIII\n', encoding='utf-8')
    (tmp_path/'samples.tsv').write_text('sample_id\tinput\nS1\treads.fastq\n', encoding='utf-8')
    paths = []
    for role in (1, 2):
        path = tmp_path/f'original_h{role}.fa';path.write_text(f'>h{role}\nACGT\n', encoding='utf-8');paths.append(path)
    rows = ['sample_id\trole\tfasta\tsha256']
    for role in (2, 1):
        path = paths[role-1];rows.append(f'S1\thaplotype{role}\t{path.name}\t{hashlib.sha256(path.read_bytes()).hexdigest()}')
    sheet = tmp_path/'ready.tsv';sheet.write_text('\n'.join(rows)+'\n', encoding='utf-8')
    settings = {'reference':{'fasta':'ref.fa','build':'custom'}, 'samples':{'sheet':'samples.tsv'},
                'runtime':{'threads':2},'workflow':{'memory_mb':4000},
                'assembly_sv':{'enabled':True,'input_manifest':'ready.tsv',
                               'pav':{'enabled':False},'svim_asm':{'enabled':True,'threads':2,'memory_mb':1000}}}
    user = tmp_path/'user.yaml';user.write_text(yaml.safe_dump(settings), encoding='utf-8')
    config = load_config(DEFAULT, user_config=user)
    monkeypatch.setattr(launcher.platform,'system',lambda:'Linux')
    # A scheduler simulation creates no scientific caller outputs/receipts.
    class Scheduler:
        def __init__(self):self.commands=[]
        def run(self, command, **kwargs):
            self.commands.append(command)
            for key in ('stdout_path','stderr_path'):
                if kwargs.get(key):Path(kwargs[key]).write_text('', encoding='utf-8')
            return SimpleNamespace(stdout='9.25.2\n' if '--version' in command else '',returncode=0,duration_seconds=.001)
    return SimpleNamespace(config=config, paths=paths, manifest=sheet, directory=tmp_path, scheduler=Scheduler())


def changed(config, section, values):
    value=config.to_dict();value[section].update(values);return HiFiVarConfig(value,config.sources)


def launched(case):
    root=case.directory/'run'
    assert launcher.WorkflowLauncher(case.scheduler).run(case.config,run_dir=root)=='WORKFLOW_SUCCEEDED'
    attempts=list((root/'attempts').iterdir());assert len(attempts)==1
    attempt=attempts[0];context=attempt/'execution_context.json'
    binding=json.loads(context.read_text(encoding='utf-8'))['ready_assembly_handoff']
    return SimpleNamespace(root=root,attempt=attempt,context=context,binding=binding)


def test_ready_manifest_is_frozen_without_enabling_assembly_or_creating_output(ready_case):
    case=ready_case;root=case.directory/'run'
    effective,context,_,inputs=launcher.prepare_launch(case.config,root)
    assert not effective['assembly']['enabled'] and not effective['alignment']['enabled']
    assert effective['assembly_sv']['input_manifest']==str(case.manifest.resolve())
    assert all(path.resolve() in inputs for path in (*case.paths,case.manifest))
    assert context.sample_ids==('S1',) and not root.exists()


@pytest.mark.parametrize('caller',['pav','svim_asm'])
def test_ready_generated_annotation_and_cohort_intent_keeps_frozen_fasta_inputs(ready_case,caller):
    case=ready_case;root=case.directory/'run'
    values=case.config.to_dict()
    values['reference']['build']='GRCh38'  # Synthetic four-base fixture, not a real human benchmark.
    values['assembly_sv']['pav'].update(enabled=caller=='pav',version='2.4.6',threads=1,memory_mb=1000)
    pav_snakefile=case.directory/'PAV.Snakefile'
    pav_snakefile.write_text('# Synthetic deployment dependency; never executed here\n',encoding='utf-8')
    values['assembly_sv']['pav']['snakefile']=str(pav_snakefile)
    values['assembly_sv']['svim_asm']['enabled']=caller=='svim_asm'
    values['annotation'].update(enabled=True,vep_enabled=True,vep_cache_directory='cache',vep_cache_version='116',threads=1,memory_mb=1000)
    values['cohort'].update(enabled=True,cohort_id='C1')
    values['cohort']['sv'].update(enabled=True,threads=1,memory_mb=1000)
    selector={'sample_id':'S1','track':'sv','origin':'GENERATED','producer':'assembly_sv.'+caller,'output_role':'vcf','reference_build':'GRCh38'}
    plan=case.directory/'sources.yaml'
    plan.write_text(yaml.safe_dump({'schema_version':1,'annotation':[dict(selector)],'cohort':[dict(selector)]}),encoding='utf-8')
    values['workflow']['downstream_source_plan']=str(plan)
    config=HiFiVarConfig(values,case.config.sources)
    effective,context,_,inputs=launcher.prepare_launch(config,root)
    assert context.sample_ids==('S1',) and not root.exists()
    assert not effective['assembly']['enabled'] and effective['assembly_sv'][caller]['enabled']
    assert effective['workflow']['downstream_source_plan']==str(root/'downstream_plan.yaml')
    assert all(p.resolve() in inputs for p in (*case.paths,case.manifest,plan))
    assert not (root/'results').exists() and not case.scheduler.commands
    assert config.to_dict()==values


@pytest.mark.parametrize('section,values',[('assembly',{'enabled':True}),
    ('assembly_sv',{'enabled':False}),('assembly_sv',{'input_manifest':''}),
    ('assembly_sv',{'input_manifest':42})])
def test_ready_mode_rejects_ambiguous_or_invalid_config(ready_case,section,values):
    with pytest.raises(ConfigurationError):validate_config(changed(ready_case.config,section,values))


def test_ready_source_digest_change_fails_before_output_or_probe(ready_case):
    ready_case.paths[0].write_text('>h1\nTGCA\n',encoding='utf-8')
    with pytest.raises(InputValidationError,match='SHA256'):
        launcher.WorkflowLauncher(ready_case.scheduler).run(ready_case.config,run_dir=ready_case.directory/'run')
    assert not ready_case.scheduler.commands and not (ready_case.directory/'run').exists()


def test_input_preparation_has_no_native_hifiasm_receipt(ready_case):
    case=ready_case;o=launched(case)
    record,items=read_ready_handoff(o.root,o.attempt,o.binding,verify_inputs=True)
    assert record['execution_kind']=='VALIDATION_ONLY' and record['native_commands']==[]
    assert record['new_hifiasm_execution'] is False
    assert [item.role.value for item in items]==['haplotype1','haplotype2']
    assert not (o.attempt/'receipts').exists()
    assert read_step_receipt(o.root,o.attempt.name,'assembly',('S1',)) is None
    assert not (o.root/'work/hifiasm').exists()


def test_same_byte_changed_validation_record_cannot_borrow_context(ready_case):
    o=launched(ready_case);path=o.attempt/HANDOFF_NAME
    record=json.loads(path.read_text());record['new_hifiasm_execution']=True
    path.write_text(json.dumps(record),encoding='utf-8')
    with pytest.raises(InputValidationError,match='sealed context'):
        read_ready_handoff(o.root,o.attempt,o.binding)


def test_changed_original_manifest_is_not_accepted_as_ready(ready_case):
    o=launched(ready_case);ready_case.manifest.write_text('changed',encoding='utf-8')
    with pytest.raises(InputValidationError,match='manifest differs'):
        read_ready_handoff(o.root,o.attempt,o.binding)


def test_selected_fasta_content_is_revalidated(ready_case):
    o=launched(ready_case);ready_case.paths[0].write_text('>h1\nTGCA\n',encoding='utf-8')
    with pytest.raises(InputValidationError,match='differs'):
        read_ready_handoff(o.root,o.attempt,o.binding,verify_inputs=True,sample_id='S1')


def test_ready_workflow_input_paths_cannot_be_swapped(ready_case,monkeypatch):
    o=launched(ready_case);monkeypatch.setenv(ATTEMPT_CONTEXT_ENV,str(o.context))
    status=json.loads((o.attempt/'status.json').read_text());status['status']='DAG_RUNNING'
    (o.attempt/'status.json').write_text(json.dumps(status),encoding='utf-8')
    with pytest.raises(InputValidationError,match='paths/roles'):
        workflow_assemblies('S1',tuple(reversed(ready_case.paths)),ready_case.config)


def test_attempt_scoped_workdir_and_explicit_ready_roles(ready_case,monkeypatch):
    o=launched(ready_case);monkeypatch.setenv(ATTEMPT_CONTEXT_ENV,str(o.context))
    # Read-only binding does not certify the scheduler simulation as a caller.
    _,items=read_ready_handoff(o.root,o.attempt,o.binding,verify_inputs=True)
    binding=assembly_role_binding('S1',items,o.context)
    assert binding['new_hifiasm_execution'] is False and 'assembly_receipt' not in binding
    assert binding['input_validation']==o.binding
    assert ready_work_directory(o.root/'work','svim_asm','S1').endswith(o.attempt.name)
    with pytest.raises(InputValidationError,match='work root'):
        ready_work_directory(ready_case.directory/'other','svim_asm','S1')
    wrong=tuple(HaplotypeAssemblyArtifact(x.sample_id,x.role,x.path,None,None,x.file_size) for x in items)
    with pytest.raises(InputValidationError,match='origin'):
        assembly_role_binding('S1',wrong,o.context)


def synthetic_source(ready_case,monkeypatch):
    """Actual tiny process witnesses; deliberately not biological caller validation."""
    o=launched(ready_case);monkeypatch.setenv(ATTEMPT_CONTEXT_ENV,str(o.context))
    import hifivar.package_resources
    runtime=json.loads((o.root/'launch_manifest.json').read_text())['identity']['runtime']
    monkeypatch.setattr(launcher,'runtime_identity',lambda _:runtime)
    status=json.loads((o.attempt/'status.json').read_text());status['status']='DAG_RUNNING'
    (o.attempt/'status.json').write_text(json.dumps(status),encoding='utf-8')
    _,items=read_ready_handoff(o.root,o.attempt,o.binding,verify_inputs=True)
    outputs={role:o.root/'results'/name for role,name in [('vcf','synthetic.vcf.gz'),('index','synthetic.vcf.gz.tbi')]}
    with attempt_step('assembly_sv.svim_asm',('S1',),outputs=outputs) as receipt:
        env={k:v for k,v in os.environ.items() if k.upper()=='SYSTEMROOT'}
        CommandRunner().run([sys.executable,'-I','-B','-c',
            'import pathlib,sys;[(pathlib.Path(p).parent.mkdir(parents=True,exist_ok=True),pathlib.Path(p).write_bytes(b"synthetic")) for p in sys.argv[1:]]',
            *map(str,outputs.values())],env=env,inherit_env=False)
        receipt.complete(outputs=outputs,tool_versions={'synthetic-process':'1'},
            validation={'validated':True,'scope':'SYNTHETIC_PROCESS_AND_BYTES_NOT_CALLER_VALIDATION'},
            metadata={'assembly_role_binding':assembly_role_binding('S1',items,o.context)})
    status['status']='WORKFLOW_SUCCEEDED';(o.attempt/'status.json').write_text(json.dumps(status),encoding='utf-8')
    return o


def test_ready_native_receipt_and_report_keep_assembly_disabled(ready_case,monkeypatch):
    o=synthetic_source(ready_case,monkeypatch)
    source=read_step_receipt(o.root,o.attempt.name,'assembly_sv.svim_asm',('S1',))
    assert source['execution_kind']=='NATIVE' and source['commands']
    value=source_haplotype_binding('S1','svim_asm',external=False,context_path=o.context)
    assert value['roles']==['haplotype1','haplotype2'] and 'READY_SOURCE' in value['status']
    result=build_launcher_report(o.root,ready_case.directory/'report',o.attempt.name)
    report=json.loads(Path(result['output_paths']['json']).read_text(encoding='utf-8'))
    tracks={x['name']:x for x in report['tracks']}
    assert tracks['assembly']['status']=='DISABLED'
    assert tracks['assembly_sv.svim_asm']['status']=='COMPLETE'
    assert tracks['assembly_input_reuse']['status']=='COMPLETE'
    assert tracks['assembly_input_reuse']['qc']['new_hifiasm_execution'] is False


def test_changed_ready_input_invalidates_selected_native_receipt(ready_case,monkeypatch):
    o=synthetic_source(ready_case,monkeypatch)
    ready_case.paths[0].write_text('>changed\nACGT\n',encoding='utf-8')
    with pytest.raises(InputValidationError):
        read_step_receipt(o.root,o.attempt.name,'assembly_sv.svim_asm',('S1',))


def test_scheduler_success_without_caller_receipt_remains_partial(ready_case):
    o=launched(ready_case)
    result=build_launcher_report(o.root,ready_case.directory/'report',o.attempt.name)
    report=json.loads(Path(result['output_paths']['json']).read_text(encoding='utf-8'))
    tracks={x['name']:x for x in report['tracks']}
    assert tracks['assembly_sv.svim_asm']['status']=='PARTIAL'
    assert tracks['assembly_input_reuse']['status']=='COMPLETE'
    assert tracks['assembly']['status']=='DISABLED'


def test_consumer_config_does_not_reactivate_ready_producer(ready_case):
    from hifivar.downstream_handoff import consumer_config, producer_config
    original=ready_case.config.to_dict()
    consumer=consumer_config(original)
    assert consumer['assembly_sv']['input_manifest'] is None
    assert not consumer['assembly_sv']['enabled'] and not consumer['assembly']['enabled']
    validate_config(consumer)
    assert producer_config(original)['assembly_sv']['input_manifest']=='ready.tsv'
    assert original['assembly_sv']['input_manifest']=='ready.tsv'


def test_unlinked_generated_bridge_preserves_unknown_gfa_and_version(ready_case,monkeypatch):
    monkeypatch.delenv(ATTEMPT_CONTEXT_ENV,raising=False)
    config=ready_case.config.to_dict();config['assembly_sv']['input_manifest']=None
    values=workflow_assemblies('S1',ready_case.paths,config)
    for item in values:
        data=item.to_dict()
        assert data['source_gfa'] is data['hifiasm_version'] is data['reference_independent'] is None
