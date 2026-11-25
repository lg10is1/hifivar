"""Actual Snakemake, synthetic tool executables: routing, not biological proof."""
import hashlib
import os
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

from hifivar.config import load_config
from hifivar.launcher import prepare_launch

ROOT = Path(__file__).parents[2]
DEFAULT = ROOT/'src/hifivar/resources/configs/default.yaml'

FAKE_TOOL = r'''import gzip,json,pathlib,struct,sys,zlib
name=pathlib.Path(sys.argv[0]).name;args=sys.argv[1:]
if '--version' in args:
 print('synthetic-contract-fixture 1.24.0');sys.exit(0)
def bgzf(data):
 c=zlib.compressobj(level=6,wbits=-15);d=c.compress(data)+c.flush();size=18+len(d)+8
 return b'\x1f\x8b\x08\x04\0\0\0\0\0\xff\x06\0BC\x02\0'+struct.pack('<H',size-1)+d+struct.pack('<II',zlib.crc32(data),len(data))
def vcf(sample):
 return ('##fileformat=VCFv4.3\n##source=PAV 2.4.6.0\n##contig=<ID=chr1,length=4>\n'
  '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">\n'
  '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="length">\n'
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="genotype">\n'
  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\n'
  'chr1\t2\tsynthetic1\tC\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=50\tGT\t0/1\n').encode()
if name=='minimap2':sys.stdout.write('@HD\tVN:1.6\n')
elif name=='samtools':
 if args[0]=='sort':pathlib.Path(args[args.index('-o')+1]).write_bytes(b'synthetic-bam-not-real')
 elif args[0]=='index':pathlib.Path(args[-1]+'.bai').write_bytes(b'synthetic-index-not-real')
 else:raise RuntimeError('Unexpected synthetic samtools command')
elif name=='svim-asm':
 sample=args[args.index('--sample')+1];native=pathlib.Path(args[args.index('--sample')+2]);native.mkdir(parents=True,exist_ok=True);(native/'variants.vcf').write_bytes(vcf(sample))
elif name=='pav-fixture':
 sample=pathlib.Path('assemblies.tsv').read_text().splitlines()[1].split('\t')[0]
 pathlib.Path(sample+'.vcf.gz').write_bytes(bgzf(vcf(sample)));pathlib.Path(sample+'.vcf.gz.tbi').write_bytes(bgzf(b'TBI\x01'))
elif name=='bgzip':sys.stdout.buffer.write(bgzf(pathlib.Path(args[-1]).read_bytes()))
elif name=='tabix':pathlib.Path(args[-1]+'.tbi').write_bytes(bgzf(b'TBI\x01'))
elif name=='hifiasm':raise RuntimeError('Disabled hifiasm MUST NOT execute in ready mode')
else:raise RuntimeError('Unexpected synthetic tool')
'''


def case(tmp_path, haplotypes):
    (tmp_path/'ref.fa').write_text('>chr1\nACGT\n',encoding='utf-8')
    (tmp_path/'ref.fa.fai').write_text('chr1\t4\t6\t4\t5\n',encoding='utf-8')
    (tmp_path/'reads.fastq').write_text('@m/1/ccs\nACGT\n+\nIIII\n',encoding='utf-8')
    (tmp_path/'samples.tsv').write_text('sample_id\tinput\nS1\treads.fastq\n',encoding='utf-8')
    (tmp_path/'PAV.Snakefile').write_text('# Synthetic adapter dependency, not real PAV\n',encoding='utf-8')
    tools=tmp_path/'tools';tools.mkdir()
    for name in ('minimap2','samtools','svim-asm','pav-fixture','bgzip','tabix','hifiasm'):
        p=tools/name;p.write_text('#!'+sys.executable+'\n'+FAKE_TOOL,encoding='utf-8');p.chmod(0o755)
    rows=['sample_id\trole\tfasta\tsha256'];identities={}
    for index in reversed(range(1,haplotypes+1)):
        p=tmp_path/f'original_h{index}.fa';p.write_text(f'>h{index}\nACGT\n',encoding='utf-8')
        digest=hashlib.sha256(p.read_bytes()).hexdigest();identities[p]=digest
        rows.append(f'S1\thaplotype{index}\t{p.name}\t{digest}')
    (tmp_path/'ready.tsv').write_text('\n'.join(rows)+'\n',encoding='utf-8')
    settings={'reference':{'fasta':'ref.fa','build':'custom'},'samples':{'sheet':'samples.tsv'},
        'runtime':{'threads':2},'workflow':{'memory_mb':4000},
        'assembly':{'enabled':False,'executable':str(tools/'hifiasm')},
        'assembly_sv':{'enabled':True,'input_manifest':'ready.tsv',
            'pav':{'enabled':True,'executable':str(tools/'pav-fixture'),'snakefile':'PAV.Snakefile','version':'2.4.6',
                   'bgzip_executable':str(tools/'bgzip'),'tabix_executable':str(tools/'tabix'),'threads':1,'memory_mb':1000},
            'svim_asm':{'enabled':True,'executable':str(tools/'svim-asm'),
                 'minimap2_executable':str(tools/'minimap2'),'samtools_executable':str(tools/'samtools'),
                 'bgzip_executable':str(tools/'bgzip'),'tabix_executable':str(tools/'tabix'),'threads':1,'memory_mb':1000}}}
    user=tmp_path/'user.yaml';user.write_text(yaml.safe_dump(settings),encoding='utf-8')
    run=tmp_path/'run';effective,_,sheet,_=prepare_launch(load_config(DEFAULT,user_config=user),run)
    run.mkdir();(run/'effective.yaml').write_text(yaml.safe_dump(effective),encoding='utf-8');(run/'samples.tsv').write_text(sheet,encoding='utf-8')
    return run,identities


@pytest.mark.parametrize('haplotypes',[1,2])
@pytest.mark.parametrize('dry_run',[True,False])
def test_snakemake_ready_inputs_do_not_register_or_execute_hifiasm(tmp_path,haplotypes,dry_run):
    pytest.importorskip('snakemake')
    if os.name=='nt':pytest.skip('Native Snakemake script execution is a Linux gate')
    run,identities=case(tmp_path,haplotypes)
    env={'PATH':str(Path(sys.executable).parent)+':/usr/bin:/bin','HOME':str(tmp_path),
         'TMPDIR':str(tmp_path),'TMP':str(tmp_path),'TEMP':str(tmp_path),'LC_ALL':'C','LANG':'C',
         'XDG_CACHE_HOME':str(tmp_path/'cache'),'XDG_CONFIG_HOME':str(tmp_path/'config'),
         'PYTHONDONTWRITEBYTECODE':'1','PYTHONPATH':str(ROOT/'src'),
         'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
    argv=[sys.executable,'-B','-m','snakemake','--snakefile',str(ROOT/'workflow/Snakefile'),
          '--configfile',str(run/'effective.yaml'),'--cores','2','--printshellcmds',
          '--shared-fs-usage','input-output','persistence','software-deployment','software-deployment-cache','sources','storage-local-copies']
    if dry_run:argv.append('--dry-run')
    result=subprocess.run([*argv,'--','all'],cwd=run,env=env,capture_output=True,text=True,
                          encoding='utf-8',errors='replace',timeout=90)
    output=result.stdout+result.stderr
    assert result.returncode==0,output
    assert 'hifiasm_assemble' not in output and not (run/'work/hifiasm').exists()
    for caller,rule in [('pav','pav_assembly_sv'),('svim_asm','svim_asm_call')]:
        assert rule in output
        assert (run/'results/assembly_sv/S1'/f'S1.{caller}.assembly.sv.vcf.gz').exists() is not dry_run
    assert not (run/'results/assembly').exists()
    assert not (run/'attempts').exists(), 'Standalone DAG must not invent launcher receipts'
    assert all(hashlib.sha256(path.read_bytes()).hexdigest()==digest for path,digest in identities.items())
