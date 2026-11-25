"""Owned benchmark contracts with synthetic comparator children, not native accuracy."""
from pathlib import Path
from types import SimpleNamespace
import base64,copy,gzip,hashlib,json,os,sys
import pytest,yaml
import hifivar.attempt_receipts as receipts
import hifivar.benchmark_handoff as bridge
from hifivar.command import CommandRunner
from hifivar.config import load_config,validate_config
from hifivar.exceptions import ConfigurationError,InputValidationError,OutputValidationError,WorkflowError,ReferenceError
from hifivar.launcher import WorkflowLauncher
from hifivar.workflow_report import build_launcher_report
from hifivar.benchmark import BenchmarkMetric,BenchmarkVariantClass
from hifivar.happy import HappyResult,HappyResultStatus
from hifivar.truvari import TruvariResult,TruvariResultStatus
from test_launcher_generated_sources import DEFAULT,PortableScheduler,_bgzf,_attempt,_status,generated_config


def vcf(sample='S1',kind='sv'):
 header=('##fileformat=VCFv4.3\n##reference=GRCh38\n##contig=<ID=chr1,length=1000>\n'
  '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
  '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
  '##INFO=<ID=TRID,Number=1,Type=String,Description="Locus">\n'
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
  '##FORMAT=<ID=AL,Number=.,Type=String,Description="Allele length">\n'
  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\n')
 row=('chr1\t100\tv1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=150;SVLEN=-50\tGT\t0/1\n' if kind=='sv' else
      'chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n' if kind=='small' else
      'chr1\t99\tL\tAA\tA\t.\tPASS\tTRID=L;END=100\tGT:AL\t0/1:1,2\n')
 return (header+row).encode()


def write_vcf(path,kind='sv'):
 path.write_bytes(_bgzf(vcf(kind=kind)));Path(str(path)+'.tbi').write_bytes(_bgzf(b'TBI\x01synthetic'))


@pytest.fixture
def benchmark_case(generated_config,monkeypatch):
 case=generated_config
 user=yaml.safe_load((case.base/'user.yaml').read_text());user['annotation']['enabled']=user['cohort']['enabled']=False
 user['workflow']['downstream_source_plan']=None
 (case.base/'catalog.bed').write_text('chr1\t99\t100\tID=L;MOTIFS=A;STRUC=(A)n\n')
 catalog='sha256:'+hashlib.sha256((case.base/'catalog.bed').read_bytes()).hexdigest()
 user['tr']={'enabled':False,'catalog':'catalog.bed','catalog_reference_build':'GRCh38'}
 user['benchmark']={'enabled':True,'source_plan':'benchmark_sources.yaml','benchmark_id':'B1','sample_id':'S1',
  'sv':{'enabled':True,'truth_vcf':'truth.sv.vcf.gz','truth_version':'fixture1','truth_source':'SYNTHETIC_NOT_REAL_TRUTH',
        'confident_bed':'regions.bed','confident_bed_version':'fixture1','threads':1,'memory_mb':1000,'size_bins':[50,100]}}
 write_vcf(case.base/'truth.sv.vcf.gz')
 row={'sample_id':'S1','track':'sv','origin':'GENERATED','producer':'sv.sawfish','output_role':'vcf'}
 (case.base/'benchmark_sources.yaml').write_text(yaml.safe_dump({'schema_version':1,'benchmark':{'sv':row}}))
 (case.base/'user.yaml').write_text(yaml.safe_dump(user));case.config=load_config(DEFAULT,user_config=case.base/'user.yaml')
 monkeypatch.setattr(bridge,'os',SimpleNamespace(environ=case.receipt_proxy.environ,path=os.path))
 case.catalog=catalog;case.monkeypatch=monkeypatch
 return case


class FakeTruvari:
 def run(self,request,**kwargs):
  header=vcf();empty=header[:header.index(b'chr1\t100')]
  payload={str(request.output_directory/n):base64.b64encode(gzip.compress(header if n.startswith('tp') else empty)).decode()
           for n in ('tp-base.vcf.gz','tp-comp.vcf.gz','fp.vcf.gz','fn.vcf.gz')}
  payload[str(request.output_directory/'summary.json')]=base64.b64encode(json.dumps({'TP-base':1,'TP-comp':1,'FP':0,'FN':0,'base cnt':1,'comp cnt':1,'precision':1.,'recall':1.,'f1':1.}).encode()).decode()
  cmd=[sys.executable,'-I','-B','-c','import base64,json,sys;from pathlib import Path;[(Path(p).parent.mkdir(parents=True,exist_ok=True),Path(p).write_bytes(base64.b64decode(v))) for p,v in json.loads(sys.argv[1]).items()]',json.dumps(payload)]
  CommandRunner().run(cmd,inherit_env=False,env={'SystemRoot':r'C:\Windows'} if os.name=='nt' else {},**kwargs)
  return TruvariResult(request,TruvariResultStatus.COMPLETED,tuple(cmd),'5.4.0-synthetic',.01,request.output_directory/'summary.json')


class FakeHappy:
 def run(self,request,**kwargs):
  payload={str(request.summary_csv):base64.b64encode(b'Type,Filter\nSNP,PASS\n').decode(),
           str(request.metrics_json_gz):base64.b64encode(gzip.compress(b'{"synthetic":true}')).decode()}
  cmd=[sys.executable,'-I','-B','-c','import base64,json,sys;from pathlib import Path;[(Path(p).parent.mkdir(parents=True,exist_ok=True),Path(p).write_bytes(base64.b64decode(v))) for p,v in json.loads(sys.argv[1]).items()]',json.dumps(payload)]
  CommandRunner().run(cmd,inherit_env=False,env={'SystemRoot':r'C:\Windows'} if os.name=='nt' else {},**kwargs)
  return HappyResult(request,HappyResultStatus.COMPLETED,tuple(cmd),'0.3.15','synthetic',.01,
   (BenchmarkMetric('f1',1.,BenchmarkVariantClass.SMALL_VARIANT,'SNP','SYNTHETIC'),),request.metrics_json_gz)


class BenchmarkScheduler(PortableScheduler):
 def __init__(self,case,failure=None):super().__init__(case);self.failure=failure
 def _produce(self,root,environment):
  if not self.fixture.config['sv']['enabled']:return
  if self.failure=='missing_producer':return
  super()._produce(root,environment)
 def run(self,command,**kwargs):
  if '-m' not in command or 'snakemake' not in command or '--version' in command:return super().run(command,**kwargs)
  path=Path(command[command.index('--configfile')+1]);cfg=yaml.safe_load(path.read_text())
  if not cfg['benchmark']['enabled']:return super().run(command,**kwargs)
  self.dags.append({'stage':'benchmark','config':cfg,'path':path,'command':command})
  self.fixture.receipt_proxy.environ={receipts.ATTEMPT_CONTEXT_ENV:kwargs['env'][receipts.ATTEMPT_CONTEXT_ENV]}
  bridge.os.environ=self.fixture.receipt_proxy.environ
  from hifivar.context import AnalysisContext
  context=AnalysisContext.from_config(cfg)
  assert not any(cfg[n]['enabled'] for n in ('alignment','small','sv','tr','phasing','assembly','assembly_sv','annotation','cohort','review'))
  if self.failure=='materialization':path.write_text(path.read_text()+'\n#changed')
  if self.failure!='missing_receipt':bridge.execute(cfg,context,happy_wrapper=FakeHappy(),truvari_wrapper=FakeTruvari())
  if self.failure=='output':(Path(cfg['paths']['outdir'])/'benchmark/B1/benchmark_metrics.tsv').write_text('changed')
  for key in ('stdout_path','stderr_path'):Path(kwargs[key]).write_text('')
  return SimpleNamespace(returncode=0,duration_seconds=.01)


def test_generated_sv_benchmark_and_report(benchmark_case):
 case=benchmark_case;runner=BenchmarkScheduler(case)
 assert WorkflowLauncher(runner).run(case.config,run_dir=case.root)=='WORKFLOW_SUCCEEDED'
 assert [d['stage'] for d in runner.dags]==['producer','benchmark']
 record=receipts.read_step_receipt(case.root,_attempt(case).name,'benchmark',('S1',))
 assert record['execution_kind']=='NATIVE' and len(record['commands'])==1
 assert record['consumer_binding']['source_receipts'][0]['producer_attempt_id']==_attempt(case).name
 assert record['metadata']['universal_accuracy'] is False
 report=build_launcher_report(case.root,case.base/'report',_attempt(case).name)
 b=json.loads(Path(report['output_paths']['json']).read_text());track=next(t for t in b['tracks'] if t['name']=='benchmark')
 assert track['status']=='COMPLETE' and track['samples']==['S1']
 assert track['qc']['truth_scope'][0]['origin']=='GENERATED'


@pytest.mark.parametrize('failure',['missing_producer','missing_receipt','materialization','output'])
def test_benchmark_fail_closed(benchmark_case,failure):
 case=benchmark_case
 with pytest.raises((InputValidationError,OutputValidationError,WorkflowError)):
  WorkflowLauncher(BenchmarkScheduler(case,failure),).run(case.config,run_dir=case.root)
 assert _status(case)['status']=='FAILED'


def test_benchmark_dry_run_is_plan_only(benchmark_case):
 case=benchmark_case;runner=BenchmarkScheduler(case)
 assert WorkflowLauncher(runner).run(case.config,run_dir=case.root,dry_run=True)=='BENCHMARK_SOURCE_PLAN_ONLY'
 assert len(runner.dags)==1 and not (_attempt(case)/'benchmark').exists()


def external_all(case,names=('small_variants','sv','assembly_sv','tr')):
 user=yaml.safe_load((case.base/'user.yaml').read_text());user['sv']['enabled']=False;rows={}
 for name in names:
  kind='small' if name=='small_variants' else 'tr' if name=='tr' else 'sv'
  q=case.base/(name+'.query.vcf.gz');t=case.base/(name+'.truth.vcf.gz');write_vcf(q,kind);write_vcf(t,kind)
  settings={'enabled':True,'truth_vcf':str(t),'truth_version':'fixture1','truth_source':'SYNTHETIC_NOT_REAL_TRUTH','threads':1,'memory_mb':1000}
  if name!='tr':settings.update(confident_bed='regions.bed',confident_bed_version='fixture1')
  else:settings['catalog_id']=case.catalog
  user['benchmark'][name]=settings
  rows[name]={'sample_id':'S1','track':kind,'origin':'EXTERNAL','path':str(q),'index_path':str(q)+'.tbi','source_tool':'fixture','source_version':'SYNTHETIC','reference_build':'GRCh38'}
  if name=='tr':rows[name]['catalog_id']=case.catalog
 for name in bridge.TRACKS:
  if name not in names:user['benchmark'][name]={'enabled':False}
 (case.base/'user.yaml').write_text(yaml.safe_dump(user));(case.base/'benchmark_sources.yaml').write_text(yaml.safe_dump({'schema_version':1,'benchmark':rows}))
 case.config=load_config(DEFAULT,user_config=case.base/'user.yaml')
 return rows


def test_all_four_tracks_external_chain(benchmark_case):
 case=benchmark_case;external_all(case)
 assert WorkflowLauncher(BenchmarkScheduler(case)).run(case.config,run_dir=case.root)=='WORKFLOW_SUCCEEDED'
 rec=receipts.read_step_receipt(case.root,_attempt(case).name,'benchmark',('S1',))
 assert len(rec['commands'])==3 and rec['consumer_binding']['source_receipts']==[]
 result=json.loads(Path(rec['outputs']['manifest']['path']).read_text())
 assert len(result['results'])==4 and {r['variant_class'] for r in result['results']}=={'small_variant','sv','assembly_sv','tr'}
 assert rec['execution_kind']=='NATIVE'


def test_tr_only_is_in_process_not_native(benchmark_case):
 case=benchmark_case;external_all(case,('tr',))
 WorkflowLauncher(BenchmarkScheduler(case)).run(case.config,run_dir=case.root)
 rec=receipts.read_step_receipt(case.root,_attempt(case).name,'benchmark',('S1',))
 assert rec['execution_kind']=='IN_PROCESS' and rec['commands']==[]


@pytest.mark.parametrize('change',['unknown_track','missing_track','wrong_sample','gvcf','truth_build','truth_sample','truth_length','truth_alias','bad_bed','unsafe_id'])
def test_plan_and_truth_rejections_precede_mkdir(benchmark_case,change):
 case=benchmark_case;user=yaml.safe_load((case.base/'user.yaml').read_text());doc=yaml.safe_load((case.base/'benchmark_sources.yaml').read_text())
 if change=='unknown_track':doc['benchmark']['bogus']=doc['benchmark']['sv']
 if change=='missing_track':doc['benchmark'].clear()
 if change=='wrong_sample':doc['benchmark']['sv']['sample_id']='UNKNOWN'
 if change=='gvcf':doc['benchmark']['sv'].update(producer='small.deepvariant',output_role='gvcf',track='small')
 if change=='unsafe_id':user['benchmark']['benchmark_id']='../../outside'
 if change in ('truth_build','truth_sample','truth_length'):
  text=vcf().replace(b'GRCh38',b'GRCh37') if change=='truth_build' else vcf().replace(b'S1',b'S9') if change=='truth_sample' else vcf().replace(b'length=1000',b'length=999')
  (case.base/'truth.sv.vcf.gz').write_bytes(_bgzf(text))
 if change=='bad_bed':(case.base/'regions.bed').write_text('chr1\t10\t1001\n')
 if change=='truth_alias':
  doc['benchmark']['sv']={'sample_id':'S1','track':'sv','origin':'EXTERNAL','path':str(case.base/'truth.sv.vcf.gz'),'index_path':str(case.base/'truth.sv.vcf.gz')+'.tbi','source_tool':'fixture','source_version':'synthetic','reference_build':'GRCh38'}
 (case.base/'user.yaml').write_text(yaml.safe_dump(user));(case.base/'benchmark_sources.yaml').write_text(yaml.safe_dump(doc));case.config=load_config(DEFAULT,user_config=case.base/'user.yaml')
 with pytest.raises((InputValidationError,ConfigurationError,OutputValidationError,ReferenceError)):
  WorkflowLauncher(BenchmarkScheduler(case)).run(case.config,run_dir=case.root)
 assert not case.root.exists()


@pytest.mark.parametrize('name',['effective.yaml','source_bindings.json','materialization.json'])
def test_report_rejects_changed_benchmark_controls(benchmark_case,name):
 case=benchmark_case;WorkflowLauncher(BenchmarkScheduler(case)).run(case.config,run_dir=case.root)
 path=_attempt(case)/'benchmark'/name;path.write_text(path.read_text()+'\n ')
 with pytest.raises((InputValidationError,OutputValidationError)):
  build_launcher_report(case.root,case.base/'report',_attempt(case).name)


def test_benchmark_config_conflicting_modes_rejected(benchmark_case):
 cfg=benchmark_case.config.to_dict();cfg['benchmark']['sv']['query_vcf']='old.vcf.gz'
 with pytest.raises(ConfigurationError,match='mutually exclusive'):validate_config(cfg)


@pytest.mark.parametrize('role',['query','truth','regions','reference'])
def test_report_rejects_changed_source_bytes(benchmark_case,role):
 case=benchmark_case;external_all(case,('sv',))
 WorkflowLauncher(BenchmarkScheduler(case)).run(case.config,run_dir=case.root)
 path={'query':case.base/'sv.query.vcf.gz','truth':case.base/'sv.truth.vcf.gz',
       'regions':case.base/'regions.bed','reference':case.base/'ref.fa'}[role]
 with path.open('ab') as f:f.write(b'changed')
 with pytest.raises((InputValidationError,OutputValidationError,ReferenceError)):
  build_launcher_report(case.root,case.base/'report',_attempt(case).name)


def test_benchmark_receipt_cannot_move_to_another_attempt(benchmark_case):
 case=benchmark_case;WorkflowLauncher(BenchmarkScheduler(case)).run(case.config,run_dir=case.root)
 attempt=_attempt(case);other=case.root/'attempts/another';other.mkdir()
 (other/'execution_context.json').write_bytes((attempt/'execution_context.json').read_bytes())
 with pytest.raises(InputValidationError,match='path/identity'):
  receipts.read_step_receipt(case.root,other.name,'benchmark',('S1',))


def test_benchmark_receipt_native_kind_cannot_be_relabelled(benchmark_case):
 case=benchmark_case;WorkflowLauncher(BenchmarkScheduler(case)).run(case.config,run_dir=case.root)
 attempt=_attempt(case);path=attempt/'receipts'/receipts.receipt_name('benchmark',('S1',))
 data=json.loads(path.read_text());data['execution_kind']='IN_PROCESS';path.write_text(json.dumps(data))
 with pytest.raises(InputValidationError,match='native/in-process'):
  receipts.read_step_receipt(case.root,attempt.name,'benchmark',('S1',))


@pytest.mark.parametrize('role',['query','truth'])
@pytest.mark.parametrize('allele',[b'<NON_REF>',b'<*>'])
def test_external_gvcf_blocks_rejected_even_under_vcf_name(benchmark_case,role,allele):
 case=benchmark_case;external_all(case,('sv',))
 path=case.base/('sv.'+role+'.vcf.gz')
 path.write_bytes(_bgzf(vcf().replace(b'<DEL>',allele)))
 with pytest.raises(InputValidationError,match='gVCF'):
  WorkflowLauncher(BenchmarkScheduler(case)).run(case.config,run_dir=case.root)
 assert not case.root.exists()


def test_spanning_deletion_star_is_not_a_gvcf_block(benchmark_case):
 case=benchmark_case;external_all(case,('sv',))
 (case.base/'sv.truth.vcf.gz').write_bytes(_bgzf(vcf().replace(b'<DEL>',b'*')))
 assert WorkflowLauncher(BenchmarkScheduler(case)).run(case.config,run_dir=case.root)=='WORKFLOW_SUCCEEDED'
