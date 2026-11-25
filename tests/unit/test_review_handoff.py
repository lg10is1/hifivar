"""Owned review scheduling/receipt contracts with explicit synthetic tools.

BAM, callers, Snakemake and IGV are synthetic in these tests. Isolated Python
children provide command observations; this is not native caller/image accuracy.
"""
import base64
import copy
import json
import os
from pathlib import Path
from types import SimpleNamespace
import sys

import pytest
import yaml

import hifivar.attempt_receipts as receipts
from hifivar.command import CommandRunner
from hifivar.config import load_config, validate_config
from hifivar.exceptions import CommandExecutionError, ConfigurationError, InputValidationError, OutputValidationError, WorkflowError
from hifivar.igv import IgvWrapper
from hifivar.launcher import WorkflowLauncher
from hifivar.review import read_review_selection
from hifivar.review_handoff import run_review_step, validate_review_materialization
from hifivar.workflow_report import build_launcher_report
from test_launcher_generated_sources import DEFAULT, PortableScheduler, _attempt, _bgzf, _status, generated_config

from test_review import PNG


@pytest.fixture
def review_case(generated_config):
    case=generated_config
    user=yaml.safe_load((case.base/'user.yaml').read_text())
    # Calling plus explicit review, without unrelated annotation/cohort.
    user['annotation']['enabled']=False
    user['cohort']['enabled']=False
    user['workflow']['downstream_source_plan']=None
    user['review']={'enabled':True,'source_plan':'review_sources.yaml','threads':1,'memory_mb':1000}
    rows=[{'sample_id':sample,'track':'sv','origin':'GENERATED','producer':'sv.sawfish','output_role':'vcf',
        'candidates':[{'review_id':sample+'_DEL','variant_id':'v1','variant_type':'DEL'}]} for sample in ('S1','S2')]
    (case.base/'review_sources.yaml').write_text(yaml.safe_dump({'schema_version':1,'review':rows}),encoding='utf-8')
    (case.base/'user.yaml').write_text(yaml.safe_dump(user),encoding='utf-8')
    case.config=load_config(DEFAULT,user_config=case.base/'user.yaml')
    return case


class FixtureIgvRunner:
    """Translate simulated tool calls into witnessed isolated Python children."""
    def __init__(self,targets,fail=False,bad_png=False):
        self.targets=targets;self.fail=fail;self.bad_png=bad_png
    def require_executable(self,*args,**kwargs):
        return Path(sys.executable)
    def run(self,argv,**kwargs):
        if '--version' in argv:
            code="print('IGV 9.9.0 synthetic fixture')"
        else:
            payloads={str(p):base64.b64encode(b'NOT_A_PNG' if self.bad_png else PNG).decode() for t in self.targets for p in t.screenshot_paths}
            code="import base64,json,sys;from pathlib import Path;d=json.loads(sys.argv[1]);[(Path(p).parent.mkdir(parents=True,exist_ok=True),Path(p).write_bytes(base64.b64decode(b))) for p,b in d.items()];sys.exit("+('7' if self.fail else '0')+")"
        args=[sys.executable,'-I','-B','-c',code]
        if '--version' not in argv:args.append(json.dumps(payloads))
        return CommandRunner().run(args,inherit_env=False,env={'SystemRoot':r'C:\Windows'} if os.name=='nt' else {},**kwargs)


class ReviewScheduler(PortableScheduler):
    def __init__(self,case,*,failure=None):
        super().__init__(case);self.failure=failure
    def _produce(self,root,environment):
        if not self.fixture.config["sv"]["enabled"]:
            return
        self.fixture.receipt_proxy.environ={receipts.ATTEMPT_CONTEXT_ENV:environment[receipts.ATTEMPT_CONTEXT_ENV]}
        if self.failure=='missing_producer':return
        for sample in ('S1','S2'):
            out=root/'results'/sample/'sawfish.vcf.gz';index=Path(str(out)+'.tbi')
            text=('##fileformat=VCFv4.3\n##reference=GRCh38\n##contig=<ID=chr1,length=1000>\n'
                '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
                '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\n'
                'chr1\t100\tv1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=150\tGT\t0/1\n')
            outputs={'vcf':out,'index':index}
            payloads={str(out):base64.b64encode(_bgzf(text.encode())).decode(),str(index):base64.b64encode(_bgzf(b'TBI\x01synthetic')).decode()}
            with receipts.attempt_step('sv.sawfish',(sample,),outputs=outputs) as step:
                CommandRunner().run([sys.executable,'-I','-B','-c',
                    "import json,base64,sys;from pathlib import Path;[(Path(p).parent.mkdir(parents=True,exist_ok=True),Path(p).write_bytes(base64.b64decode(b))) for p,b in json.loads(sys.argv[1]).items()]",json.dumps(payloads)],
                    inherit_env=False,env={'SystemRoot':r'C:\Windows'} if os.name=='nt' else {})
                step.complete(outputs=outputs,tool_versions={'sawfish':'2.2.1-synthetic'},validation={'validated':True,'scope':'SYNTHETIC'},metadata={'reference_build':'GRCh38'})
            self.native_outputs.extend(outputs.values())
    def run(self,command,**kwargs):
        if '-m' not in command or 'snakemake' not in command or '--version' in command:
            return super().run(command,**kwargs)
        cfg_path=Path(command[command.index('--configfile')+1]);cfg=yaml.safe_load(cfg_path.read_text())
        if not cfg['review'].get('enabled'):
            return super().run(command,**kwargs)
        self.calls.append((list(command),dict(kwargs)))
        self.dags.append({'stage':'review','config':cfg,'path':cfg_path,'command':command})
        assert not any(cfg[name]['enabled'] for name in ('small','sv','tr','phasing','assembly','assembly_sv','annotation','cohort','alignment'))
        assert not cfg['sv']['harmonization']['enabled']
        self.fixture.receipt_proxy.environ={receipts.ATTEMPT_CONTEXT_ENV:kwargs['env'][receipts.ATTEMPT_CONTEXT_ENV]}
        from hifivar.context import AnalysisContext
        context=AnalysisContext.from_config(cfg);selection=Path(cfg['review']['selection_file'])
        output=Path(cfg['paths']['outdir'])/'review'
        targets=read_review_selection(selection,alignments={r.sample.sample_id:r.sample.input.files[0] for r in context.samples},
            reference_fasta=context.reference.fasta,output_directory=output,default_flank_bp=cfg['review']['flank_bp'],reference_build='GRCh38')
        attempt=cfg_path.parent.parent
        if self.failure=='selection_tamper':selection.write_text(selection.read_text()+'\n')
        if self.failure=='alignment_log_tamper':(attempt/'S1.header.stdout.log').write_text('@HD\tSO:unsorted\n')
        if self.failure!='missing_review_receipt':
            run_review_step(cfg,context,targets,output=output,selection=selection,
                igv_wrapper=IgvWrapper(executable='synthetic-igv',runner=FixtureIgvRunner(targets,fail=self.failure=='igv_exit',bad_png=self.failure=='invalid_png')))
        if self.failure=='output_tamper':targets[0].screenshot_paths[0].write_bytes(b'changed')
        for key in ('stdout_path','stderr_path'):
            if kwargs.get(key):Path(kwargs[key]).write_text('')
        return SimpleNamespace(args=tuple(command),stdout='',returncode=0,duration_seconds=0.01)


def test_review_source_config_modes_are_exclusive(review_case):
    cfg=review_case.config.to_dict()
    cfg['review']['selection_file']='also.tsv'
    with pytest.raises(ConfigurationError,match='mutually exclusive'):validate_config(cfg)


def test_generated_calling_then_review_receipt_and_report(review_case):
    case=review_case;original=copy.deepcopy(case.config.to_dict());runner=ReviewScheduler(case)
    assert WorkflowLauncher(runner).run(case.config,run_dir=case.root)=='WORKFLOW_SUCCEEDED'
    assert [d['stage'] for d in runner.dags]==['producer','review']
    assert not runner.dags[0]['config']['review']['enabled']
    attempt=_attempt(case)
    record=receipts.read_step_receipt(case.root,attempt.name,'review',('S2','S1'))
    assert record['metadata']['human_review']=='NOT_REVIEWED'
    assert record['consumer_binding']['selection']['review_ids']==['S1_DEL','S2_DEL']
    assert len(record['commands'])==2 and all(x['returncode']==0 for x in record['commands'])
    assert len(record['outputs'])==7 and record['tool_versions']=={'igv':'9.9.0'}
    assert record['consumer_binding']['alignment']['mode'].startswith('THIS_ATTEMPT_NATIVE_HEADER')
    assert len(record['consumer_binding']['producer_receipts'])==2
    assert case.config.to_dict()==original
    summary=build_launcher_report(case.root,case.base/'report',attempt.name)
    report=json.loads(Path(summary['output_paths']['json']).read_text())
    review=next(t for t in report['tracks'] if t['name']=='review')
    assert review['status']=='COMPLETE'
    assert review['qc']['accuracy_status']=='NOT_ASSESSED'
    assert review['qc']['source_rows'][0]['human_review']=='NOT_REVIEWED'


@pytest.mark.parametrize('failure',['missing_producer','missing_review_receipt','selection_tamper','alignment_log_tamper','igv_exit','invalid_png','output_tamper'])
def test_review_failure_is_not_promoted_by_scheduler_success(review_case,failure):
    case=review_case;runner=ReviewScheduler(case,failure=failure)
    with pytest.raises((CommandExecutionError,InputValidationError,OutputValidationError,WorkflowError)):
        WorkflowLauncher(runner).run(case.config,run_dir=case.root)
    assert _status(case)['status']=='FAILED'
    assert all(p.exists() for p in runner.native_outputs)
    if failure=='missing_producer':assert len(runner.dags)==1
    if failure in {'missing_review_receipt','selection_tamper','alignment_log_tamper','igv_exit','invalid_png'}:
        assert not (_attempt(case)/'receipts'/receipts.receipt_name('review',('S2','S1'))).exists()


def test_review_dry_run_does_not_run_igv_or_materialize_called_sources(review_case):
    case=review_case;runner=ReviewScheduler(case)
    assert WorkflowLauncher(runner).run(case.config,run_dir=case.root,dry_run=True)=='REVIEW_SOURCE_PLAN_ONLY'
    assert len(runner.dags)==1 and '--dry-run' in runner.dags[0]['command']
    assert not (_attempt(case)/'review').exists() and not runner.native_outputs


def test_review_preserves_existing_annotation_cohort_stage(review_case):
    case=review_case;user=yaml.safe_load((case.base/'user.yaml').read_text())
    user['annotation']['enabled']=user['cohort']['enabled']=True
    user['workflow']['downstream_source_plan']='sources.yaml'
    (case.base/'user.yaml').write_text(yaml.safe_dump(user))
    case.config=load_config(DEFAULT,user_config=case.base/'user.yaml')
    runner=ReviewScheduler(case)
    assert WorkflowLauncher(runner).run(case.config,run_dir=case.root)=='WORKFLOW_SUCCEEDED'
    assert [d['stage'] for d in runner.dags]==['producer','downstream','review']
    assert not runner.dags[1]['config']['review']['enabled']


def test_missing_review_plan_is_rejected_before_output_directory(review_case):
    case=review_case;user=yaml.safe_load((case.base/'user.yaml').read_text())
    user['review']['source_plan']='absent.yaml'
    (case.base/'user.yaml').write_text(yaml.safe_dump(user));case.config=load_config(DEFAULT,user_config=case.base/'user.yaml')
    with pytest.raises(InputValidationError):WorkflowLauncher(ReviewScheduler(case)).run(case.config,run_dir=case.root)
    assert not case.root.exists()



def test_external_only_review_does_not_invent_new_caller_execution(review_case):
    case=review_case
    user=yaml.safe_load((case.base/'user.yaml').read_text());user['sv']['enabled']=False
    rows=[];before={}
    for sample in ('S1','S2'):
        source=case.base/(sample+'.external.vcf.gz')
        text=('##fileformat=VCFv4.3\n##reference=GRCh38\n##contig=<ID=chr1,length=1000>\n'
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\n'
            'chr1\t100\tv1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=150\tGT\t0/1\n')
        source.write_bytes(_bgzf(text.encode()));index=Path(str(source)+'.tbi');index.write_bytes(_bgzf(b'TBI\x01synthetic'))
        before[str(source)]=source.read_bytes();before[str(index)]=index.read_bytes()
        rows.append({'sample_id':sample,'track':'sv','origin':'EXTERNAL','path':str(source),'index_path':str(index),
            'source_tool':'sawfish','source_version':'DECLARED_UNPROBED','reference_build':'GRCh38',
            'candidates':[{'review_id':sample+'_DEL','variant_id':'v1','variant_type':'DEL'}]})
    (case.base/'user.yaml').write_text(yaml.safe_dump(user))
    (case.base/'review_sources.yaml').write_text(yaml.safe_dump({'schema_version':1,'review':rows}))
    case.config=load_config(DEFAULT,user_config=case.base/'user.yaml');runner=ReviewScheduler(case)
    assert WorkflowLauncher(runner).run(case.config,run_dir=case.root)=='WORKFLOW_SUCCEEDED'
    record=receipts.read_step_receipt(case.root,_attempt(case).name,'review',('S2','S1'))
    assert record['consumer_binding']['producer_receipts']==[] and runner.native_outputs==[]
    bindings=json.loads((_attempt(case)/'review/source_bindings.json').read_text())
    assert all(b['version_source']=='USER_DECLARED_NOT_REPROBED' for b in bindings['source_bindings'])
    assert all(Path(p).read_bytes()==v for p,v in before.items())


def test_review_receipt_cannot_be_borrowed_by_a_new_attempt(review_case):
    case=review_case;runner=ReviewScheduler(case)
    WorkflowLauncher(runner).run(case.config,run_dir=case.root)
    attempt=_attempt(case)
    other=case.root/'attempts'/'other_attempt';other.mkdir()
    (other/'execution_context.json').write_bytes((attempt/'execution_context.json').read_bytes())
    with pytest.raises(InputValidationError,match='path/identity'):
        receipts.read_step_receipt(case.root,other.name,'review',('S2','S1'))


@pytest.mark.parametrize('field',['source_bindings','effective','materialization','plan','preflight'])
def test_postrun_source_control_changes_cannot_be_reported_as_complete(review_case,field):
    case=review_case;runner=ReviewScheduler(case)
    WorkflowLauncher(runner).run(case.config,run_dir=case.root)
    attempt=_attempt(case)
    path={'source_bindings':attempt/'review/source_bindings.json','effective':attempt/'review/effective.yaml',
          'materialization':attempt/'review/materialization.json','plan':case.root/'review_source_plan.yaml',
          'preflight':attempt/'input_preflight.json'}[field]
    path.write_bytes(path.read_bytes()+b'\n ')
    with pytest.raises(InputValidationError):
        receipts.read_step_receipt(case.root,attempt.name,'review',('S2','S1'))
    with pytest.raises(InputValidationError):build_launcher_report(case.root,case.base/'changed_report',attempt.name)


def test_invalid_stage_cannot_issue_review_receipt(review_case):
    case=review_case;runner=ReviewScheduler(case)
    WorkflowLauncher(runner).run(case.config,run_dir=case.root)
    attempt=_attempt(case);status=json.loads((attempt/'status.json').read_text());status['status']='DAG_RUNNING'
    (attempt/'status.json').write_text(json.dumps(status))
    with pytest.raises(InputValidationError,match='stage'):
        with receipts.attempt_step('review',('S2','S1'),outputs={'batch':case.root/'results/other.batch'}):
            pytest.fail('Wrong stage must be rejected before step execution')
