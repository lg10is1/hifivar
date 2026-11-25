"""Real owned benchmark DAG planning; structural fixtures are not accuracy data."""
from __future__ import annotations

import hashlib
import os
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

from hifivar import launcher
from hifivar.benchmark_handoff import materialize
from hifivar.config import load_config
from hifivar.context import AnalysisContext
from test_snakemake_launcher_downstream import _bgzf, _SchedulerBoundary

ROOT = Path(__file__).parents[2]
DEFAULT = ROOT / 'src/hifivar/resources/configs/default.yaml'


def _exercise_owned_benchmark_plan(tmp_path, monkeypatch, names):
    pytest.importorskip('snakemake')
    (tmp_path/'ref.fa').write_text('>chr1\n'+'A'*1000+'\n')
    (tmp_path/'ref.fa.fai').write_text('chr1\t1000\t6\t1000\t1001\n')
    (tmp_path/'HG002.bam').write_bytes(b'UNUSED_CONTEXT_NOT_A_REAL_BAM')
    (tmp_path/'samples.tsv').write_text('sample_id\tinput\nHG002\tHG002.bam\n')
    (tmp_path/'regions.bed').write_text('chr1\t90\t160\n')
    catalog='sha256:'+hashlib.sha256(b'SYNTHETIC_CATALOG_IDENTITY').hexdigest()
    rows={}
    benchmark={'enabled': True, 'source_plan': 'sources.yaml', 'benchmark_id': 'B1', 'sample_id': 'HG002'}
    for name in names:
        for role in ('query','truth'):
            path=tmp_path/(name+'.'+role+'.vcf.gz')
            text=('##fileformat=VCFv4.3\n##reference=GRCh38\n##contig=<ID=chr1,length=1000>\n'
                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n'
                  'chr1\t100\tv1\tA\tC\t.\tPASS\t.\tGT\t0/1\n')
            path.write_bytes(_bgzf(text.encode()))
            Path(str(path)+'.tbi').write_bytes(_bgzf(b'TBI\x01STRUCTURAL_FIXTURE_ONLY'))
        benchmark[name]={'enabled': True, 'truth_vcf': name+'.truth.vcf.gz',
            'truth_version': 'synthetic1', 'truth_source': 'SYNTHETIC_NOT_REAL_TRUTH',
            'threads': 1, 'memory_mb': 1000}
        if name!='tr':
            benchmark[name].update(confident_bed='regions.bed', confident_bed_version='synthetic1')
        else:
            benchmark[name]['catalog_id']=catalog
        rows[name]={'sample_id': 'HG002', 'track': 'small' if name=='small_variants' else 'tr' if name=='tr' else 'sv',
            'origin': 'EXTERNAL', 'path': name+'.query.vcf.gz', 'index_path': name+'.query.vcf.gz.tbi',
            'source_tool': 'synthetic_fixture', 'source_version': 'synthetic1', 'reference_build': 'GRCh38'}
        if name=='tr': rows[name]['catalog_id']=catalog
    (tmp_path/'sources.yaml').write_text(yaml.safe_dump({'schema_version':1,'benchmark':rows}))
    (tmp_path/'user.yaml').write_text(yaml.safe_dump({
        'reference': {'fasta':'ref.fa','build':'GRCh38'}, 'samples':{'sheet':'samples.tsv'},
        'runtime':{'threads':2},'workflow':{'memory_mb':4096}, 'benchmark':benchmark}))
    monkeypatch.setattr(launcher.platform,'system',lambda:'Linux')
    monkeypatch.setattr(launcher,'installed_workflow_root',lambda:ROOT/'workflow')
    run=tmp_path/'run'
    assert launcher.WorkflowLauncher(_SchedulerBoundary()).run(
        load_config(DEFAULT,user_config=tmp_path/'user.yaml'),run_dir=run,dry_run=True)=='BENCHMARK_SOURCE_PLAN_ONLY'
    # Only construct consumer controls for real Snakemake planning. No receipt or
    # native execution is claimed by this test, and no structural VCF is queried.
    attempt=next((run/'attempts').iterdir())
    config=yaml.safe_load((run/'effective.yaml').read_text())
    context=AnalysisContext.from_config(config)
    control,_=materialize(config,context,run,attempt,calling_config=config)
    env={**os.environ,'PYTHONPATH':str(ROOT/'src'),'HIFIVAR_ATTEMPT_CONTEXT':str(attempt/'execution_context.json'),
        'TMPDIR':str(tmp_path),'TMP':str(tmp_path),'TEMP':str(tmp_path),'XDG_CACHE_HOME':str(tmp_path/'cache'),
        'XDG_CONFIG_HOME':str(tmp_path/'config'),'PYTHONDONTWRITEBYTECODE':'1'}
    result=subprocess.run([sys.executable,'-B','-m','snakemake','--snakefile',str(ROOT/'workflow/Snakefile'),
        '--configfile',str(control),'--cores','2','--dry-run','--printshellcmds','--shared-fs-usage',
        'input-output','persistence','software-deployment','software-deployment-cache','sources','storage-local-copies','--','all'],
        cwd=run,env=env,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=60)
    output=result.stdout+result.stderr
    assert result.returncode==0,output
    assert 'benchmark_owned' in output
    for absent in ('deepvariant_small','sv_sawfish','benchmark_small_variants','benchmark_assembly_sv','benchmark_manifest'):
        assert 'rule '+absent+':' not in output,output
    assert not (run/'results/benchmark').exists()
    assert not (attempt/'receipts').exists() or not list((attempt/'receipts').glob('*benchmark*'))


@pytest.mark.parametrize('names', [
    ('small_variants',), ('sv',), ('assembly_sv',), ('tr',),
    ('small_variants', 'sv', 'assembly_sv', 'tr'),
])
def test_owned_benchmark_snakemake_dag(tmp_path, monkeypatch, names):
    _exercise_owned_benchmark_plan(tmp_path, monkeypatch, names)


@pytest.mark.skipif(os.name!='posix',reason='Real Snakemake subprocess execution requires Linux')
def test_owned_truvari_directory_is_not_precreated(tmp_path, monkeypatch):
    """Real scheduler + fake strict CLI, not a native biological benchmark."""
    import json
    import shlex
    import shutil
    from hifivar import __version__

    _exercise_owned_benchmark_plan(tmp_path,monkeypatch,('sv','assembly_sv'))
    program=tmp_path/'strict_truvari.py'
    program.write_text('''import gzip,json,sys
from pathlib import Path
if sys.argv[1]=='version':
 print('Truvari v5.4.0');raise SystemExit(0)
assert sys.argv[1]=='bench'
target=Path(sys.argv[sys.argv.index('-o')+1])
assert not target.exists(), 'Scheduler must not pre-create native output directory'
target.mkdir()
header=('##fileformat=VCFv4.3\\n##reference=GRCh38\\n##contig=<ID=chr1,length=1000>\\n'
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\\n'
        '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tHG002\\n')
record='chr1\\t100\\tv1\\tA\\t<DEL>\\t.\\tPASS\\tSVTYPE=DEL;END=150;SVLEN=-50\\tGT\\t0/1\\n'
for name in ('tp-base.vcf.gz','tp-comp.vcf.gz','fp.vcf.gz','fn.vcf.gz'):
 with gzip.open(target/name,'wt') as f:f.write(header+(record if name.startswith('tp') else ''))
(target/'summary.json').write_text(json.dumps({'TP-base':1,'TP-comp':1,'FP':0,'FN':0,
 'base cnt':1,'comp cnt':1,'precision':1.,'recall':1.,'f1':1.}))
''',encoding='utf-8')
    executable=tmp_path/'strict-truvari'
    executable.write_text('#!/bin/sh\nexec '+shlex.quote(sys.executable)+' -I -B '+shlex.quote(str(program))+' "$@"\n')
    executable.chmod(0o700)
    config_path=tmp_path/'user.yaml';cfg=yaml.safe_load(config_path.read_text())
    for name in ('sv','assembly_sv'):cfg['benchmark'][name]['truvari_executable']=str(executable)
    config_path.write_text(yaml.safe_dump(cfg))
    # The launcher intentionally uses -I. A PYTHONPATH-only source test cannot
    # supply that isolated subprocess. Materialize a test-only source fixture
    # in its own venv; actual wheel/sdist installation is a separate gate.
    runtime=tmp_path/'fixture_runtime'
    env={key:os.environ[key] for key in ('PATH','HOME','LANG','LC_ALL') if key in os.environ}
    env.update(TMPDIR=str(tmp_path),TMP=str(tmp_path),TEMP=str(tmp_path),PYTHONDONTWRITEBYTECODE='1')
    created=subprocess.run([sys.executable,'-B','-m','venv','--system-site-packages',str(runtime)],
        env=env,capture_output=True,text=True,timeout=60)
    assert created.returncode==0,created.stdout+created.stderr
    py=runtime/'bin/python'
    location=subprocess.run([str(py),'-I','-B','-c','import sysconfig;print(sysconfig.get_path("purelib"))'],
        env=env,capture_output=True,text=True,timeout=15)
    assert location.returncode==0,location.stderr
    site=Path(location.stdout.strip())
    assert site.resolve().is_relative_to(runtime.resolve())
    shutil.copytree(ROOT/'src/hifivar',site/'hifivar',ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copytree(ROOT/'workflow',runtime/'hifivar/workflow',ignore=shutil.ignore_patterns('__pycache__'))
    metadata=site/f'hifivar-{__version__}.dist-info';metadata.mkdir()
    (metadata/'METADATA').write_text(f'Metadata-Version: 2.1\nName: hifivar\nVersion: {__version__}\n')
    relative=os.path.relpath(runtime/'hifivar/workflow/Snakefile',site)
    (metadata/'RECORD').write_text(f'{relative},,\nhifivar-{__version__}.dist-info/METADATA,,\n')
    root=tmp_path/'executed'
    result=subprocess.run([str(py),'-I','-B','-m','hifivar','--config',str(config_path),'run','--run-dir',str(root)],
        env=env,cwd=tmp_path,capture_output=True,text=True,timeout=120)
    assert result.returncode==0,result.stdout+result.stderr
    attempt=next((root/'attempts').iterdir())
    for name in ('sv','assembly_sv'):
        output=root/'results/benchmark/B1'/name/'truvari'
        assert (output/'summary.json').is_file() and (output/'.snakemake_timestamp').is_file()
    verify='''import json,sys
from pathlib import Path
from hifivar.attempt_receipts import read_step_receipt
from hifivar.workflow_report import build_launcher_report
root=Path(sys.argv[1]);attempt=sys.argv[2]
receipt=read_step_receipt(root,attempt,'benchmark',('HG002',))
assert receipt['execution_kind']=='NATIVE' and len(receipt['commands'])==4
assert set(receipt['metadata']['benchmark_statuses'])=={'sv','assembly_sv'}
report=build_launcher_report(root,Path(sys.argv[3]),attempt)
data=json.loads(Path(report['output_paths']['json']).read_text())
assert next(t for t in data['tracks'] if t['name']=='benchmark')['status']=='COMPLETE'
'''
    checked=subprocess.run([str(py),'-I','-B','-c',verify,str(root),attempt.name,str(tmp_path/'report')],
        env=env,cwd=tmp_path,capture_output=True,text=True,timeout=60)
    assert checked.returncode==0,checked.stdout+checked.stderr
