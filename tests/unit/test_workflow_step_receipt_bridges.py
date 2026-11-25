"""Source bridges run real tiny Python stand-ins via the real CommandRunner.

These are engineering fixtures, not native biological-tool or accuracy claims.
No CommandResult, receipt completion, wrapper validator, or runtime identity is
mocked. Only executable dispatch is adapted to a portable Python stand-in.
"""
from __future__ import annotations

import ast
import gzip
import hashlib
import importlib
import json
import os
from pathlib import Path
import runpy
import struct
import sys
from types import SimpleNamespace
import zlib

import pytest
import yaml

from hifivar import launcher
from hifivar.attempt_receipts import (
    ATTEMPT_CONTEXT_ENV, artifact_identity, create_attempt_context, read_step_receipt,
)
from hifivar.command import CommandRunner
from hifivar.exceptions import CommandExecutionError, InputValidationError, OutputValidationError
from hifivar.package_resources import installed_workflow_root


ROOT = Path(__file__).parents[2]
SCRIPTS = ROOT / "workflow/scripts"
BRIDGES = (
    "run_deepvariant.py", "run_sv_caller.py", "run_trgt.py", "run_hiphase.py",
    "run_hifiasm.py", "run_pav.py", "run_svim_asm.py", "run_harmonization.py",
    "run_annotation.py", "run_cohort_small.py", "run_cohort_sv.py", "run_cohort_tr.py",
)
VCF_HEADER = (
    '##fileformat=VCFv4.3\n##reference=GRCh38\n##contig=<ID=chr1,length=100>\n'
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n'
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
)


FAKE_TOOL = r'''
import gzip,json,pathlib,shutil,struct,sys,zlib
tool=pathlib.Path(sys.argv[1]).stem
args=sys.argv[2:]
base='##fileformat=VCFv4.3\n##reference=GRCh38\n##contig=<ID=chr1,length=100>\n'
base+='##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
base+='##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">\n'
base+='##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
base+='##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
def column(samples=('S1',)):
 return '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(samples)+'\n'
def value(name):return args[args.index(name)+1]
def eq(name):return next(a.partition('=')[2] for a in args if a.startswith(name+'='))
def write(path,data):
 p=pathlib.Path(path);p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data if isinstance(data,bytes) else data.encode())
def bgzf(path,data):
 c=zlib.compressobj(wbits=-15);x=c.compress(data)+c.flush();size=18+len(x)+8
 write(path,b'\x1f\x8b\x08\x04\0\0\0\0\0\xff\x06\0BC\x02\0'+struct.pack('<H',size-1)+x+struct.pack('<II',zlib.crc32(data),len(data)))
def sv(sample='S1',extra=''):
 return base+extra+column((sample,))+'chr1\t10\tv1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=70;SVLEN=-60\tGT\t0/1\n'
versions={'deepvariant':'DeepVariant version 1.10.0','sawfish':'sawfish 2.2.1','sniffles2':'Sniffles2 2.8.0','pbsv':'pbsv 2.11.0','cutesv':'cuteSV 2.1.4','trgt':'trgt 5.1.0','bcftools':'bcftools 1.21','samtools':'samtools 1.23.1','hiphase':'HiPhase 1.5.0','hifiasm':'hifiasm 0.25.0-r726','pav':'9.26.1','svim_asm':'svim-asm 1.0.3','minimap2':'2.31-r1302','bgzip':'bgzip (htslib) 1.23.1','tabix':'tabix (htslib) 1.23.1','jasmine':'Jasmine version 1.1.5','truvari':'Truvari v5.4.0','glnexus':'glnexus_cli release v1.4.1','vep':'ensembl-vep 116.0'}
if '--version' in args or '--help' in args or (tool=='truvari' and args==['version']):
 print(versions[tool]);sys.exit(0)
if pathlib.Path(sys.argv[1]).with_suffix('.fail').exists():
 print('intentional synthetic executable failure',file=sys.stderr);sys.exit(7)
if tool=='deepvariant':
 sample=eq('--sample_name');v=eq('--output_vcf');g=eq('--output_gvcf')
 bgzf(v,(base+column((sample,))+'chr1\t10\tv1\tA\tT\t.\tPASS\t.\tGT\t0/1\n').encode())
 markers='##DeepVariant_version=1.10.0\n##FILTER=<ID=RefCall,Description="Reference">\n##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="DP">\n'
 bgzf(g,(base+markers+column((sample,))).encode());bgzf(v+'.tbi',b'TBI\x01');bgzf(g+'.tbi',b'TBI\x01')
elif tool=='sawfish':
 d=pathlib.Path(value('--output-dir'));bgzf(d/'genotyped.sv.vcf.gz',sv().encode());bgzf(d/'genotyped.sv.vcf.gz.tbi',b'TBI\x01')
elif tool=='sniffles2':
 p=value('--vcf');bgzf(p,sv(value('--sample-id')).encode());bgzf(p+'.tbi',b'TBI\x01')
elif tool=='pbsv':
 if args[0]=='discover':write(args[-1],b'signatures')
 else:write(args[-1],sv())
elif tool=='cutesv':write(args[2],sv())
elif tool=='bgzip':
 data=pathlib.Path(args[-1]).read_bytes();c=zlib.compressobj(wbits=-15);x=c.compress(data)+c.flush();size=18+len(x)+8
 sys.stdout.buffer.write(b'\x1f\x8b\x08\x04\0\0\0\0\0\xff\x06\0BC\x02\0'+struct.pack('<H',size-1)+x+struct.pack('<II',zlib.crc32(data),len(data)))
elif tool=='tabix':
 if '-p' in args:bgzf(args[-1]+'.tbi',b'TBI\x01')
 elif '-l' in args:
  with gzip.open(args[-1],'rt') as f:print('\n'.join(dict.fromkeys(line.split('\t')[0] for line in f if not line.startswith('#'))))
 else:
  with gzip.open(args[0],'rt') as f:sys.stdout.write(''.join(line for line in f if not line.startswith('#')))
elif tool=='trgt':
 prefix=value('--output-prefix');extra='##INFO=<ID=TRID,Number=1,Type=String,Description="ID">\n##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motifs">\n##INFO=<ID=STRUC,Number=1,Type=String,Description="Structure">\n'
 bgzf(prefix+'.vcf.gz',(base+extra+column()+ 'chr1\t10\tTR1\tA\tAC\t.\tPASS\tTRID=TR1;MOTIFS=C;STRUC=<TR>;END=20\tGT\t0/1\n').encode());write(prefix+'.spanning.bam',b'fixture-bam')
elif tool=='bcftools':
 if args[0]=='index':
  suffix='.csi' if '--csi' in args else '.tbi';p=args[-1]+suffix
  if suffix=='.tbi':bgzf(p,b'TBI\x01')
  else:write(p,b'CSI\x01')
 elif args[0]=='sort':shutil.copyfile(args[-1],value('-o'))
 elif args[0]=='view':
  sample_names=tuple(pathlib.Path(args[-1]).read_text(encoding='utf-8').splitlines()[0].split('\t'))
  bgzf(value('-o'),(base+column(sample_names)+'chr1\t10\tv1\tA\tT\t.\tPASS\t.\tGT\t'+'\t'.join(['0/1']+['0/0']*(len(sample_names)-1))+'\n').encode())
elif tool=='samtools':
 if args[0]=='sort':write(value('-o'),b'fixture-sorted-bam')
 elif args[0]=='index':write(value('-o') if '-o' in args else args[-1]+'.bai',b'BAI\x01')
elif tool=='hiphase':bgzf(value('--output-vcf'),sv(extra='##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n').encode())
elif tool=='hifiasm':
 prefix=value('-o')
 for ending in ('.bp.p_ctg.gfa','.bp.hap1.p_ctg.gfa','.bp.hap2.p_ctg.gfa'):write(prefix+ending,'H\tVN:Z:1.0\nS\tctg1\tACGT\n')
elif tool=='pav':
 d=pathlib.Path(value('--directory'));sample=(d/'assemblies.tsv').read_text(encoding='utf-8').splitlines()[1].split('\t')[0]
 extra='##source=PAV 2.4.6.0\n##INFO=<ID=VARTYPE,Number=1,Type=String,Description="Class">\n'
 text=sv(sample,extra).replace('SVTYPE=DEL;','VARTYPE=SV;SVTYPE=DEL;');bgzf(d/(sample+'.vcf.gz'),text.encode());bgzf(d/(sample+'.vcf.gz.tbi'),b'TBI\x01')
elif tool=='minimap2':print('@HD\tSO:unsorted\n@SQ\tSN:chr1\tLN:100')
elif tool=='svim_asm':
 d=pathlib.Path(args[args.index('--sample')+2]);write(d/'variants.vcf',sv(value('--sample')))
elif tool=='jasmine':
 inputs=pathlib.Path(eq('file_list')).read_text(encoding='utf-8').splitlines();ids=[]
 for p in inputs:
  with gzip.open(p,'rt',encoding='utf-8') if p.endswith('.gz') else open(p,encoding='utf-8') as f:
   ids.append(next(line.split('\t')[2] for line in f if not line.startswith('#')))
 extra='##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Membership">\n##INFO=<ID=IDLIST,Number=.,Type=String,Description="IDs">\n'
 text=sv(extra=extra).replace('##fileformat=VCFv4.3','##fileformat=VCFv4.4').replace('ID=SVLEN,Number=1','ID=SVLEN,Number=A').replace('SVTYPE=DEL;','SUPP_VEC='+'1'*len(ids)+';IDLIST='+','.join(ids)+';SVTYPE=DEL;');write(eq('out_file'),text)
elif tool=='truvari':
 d=pathlib.Path(value('-o'));d.mkdir(parents=True)
 write(d/'summary.json',json.dumps({'TP-base':1,'TP-comp':1,'FP':0,'FN':0,'precision':1.0,'recall':1.0,'f1':1.0,'base cnt':1,'comp cnt':1}))
 for name in ('tp-base','tp-comp','fp','fn'):bgzf(d/(name+'.vcf.gz'),sv().encode())
elif tool=='glnexus':
 names=[]
 for p in args[args.index('--mem-gbytes')+2:]:
  with gzip.open(p,'rt',encoding='utf-8') as f:names.extend(next(l.rstrip().split('\t')[9:] for l in f if l.startswith('#CHROM')))
 print('\t'.join(sorted(names)))
elif tool=='annovar':
 prefix=value('-out');build=value('-buildver');source=args[0]
 write(prefix+'.'+build+'_multianno.txt','Chr\tStart\tEnd\tRef\tAlt\nchr1\t10\t10\tA\tT\n')
 shutil.copyfile(source,prefix+'.'+build+'_multianno.vcf')
elif tool=='vep':write(value('--output_file'),'#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\nv1\tchr1:10\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\n')
else:raise RuntimeError('Unknown synthetic tool '+tool)
'''


@pytest.fixture
def bridge(tmp_path, monkeypatch):
    root, external = tmp_path / "run", tmp_path / "input"
    root.mkdir()
    external.mkdir()
    reference = external / "reference.fa"
    reference.write_text(">chr1\n" + "A" * 100 + "\n", encoding="utf-8")
    Path(str(reference) + ".fai").write_text("chr1\t100\t6\t100\t101\n", encoding="utf-8")
    sheet = root / "samples.tsv"
    for sample in ("S2", "S1"):
        (external / (sample + ".bam")).write_bytes(b"synthetic-input-" + sample.encode())
        (external / (sample + ".bam.bai")).write_bytes(b"synthetic-index")
    sheet.write_text("sample_id\tinput\tinput_type\tsex\n" + "".join(
        f"{sample}\t{external / (sample + '.bam')}\tbam\tfemale\n" for sample in ("S2", "S1")), encoding="utf-8")
    fastq = external / "reads.fastq"
    fastq.write_text("@m/1/ccs\nACGT\n+\nIIII\n", encoding="utf-8")
    for name in ("hap1.fa", "hap2.fa"):
        (external / name).write_text(">ctg1\nACGT\n", encoding="utf-8")
    catalog = external / "catalog.bed"
    catalog.write_text("chr1\t9\t20\tID=TR1;MOTIFS=C;STRUC=<TR>\n", encoding="utf-8")
    database, cache = external / "db", external / "cache"
    database.mkdir()
    cache.mkdir()
    pav_snakefile = external / "PAV.Snakefile"
    pav_snakefile.write_text("# synthetic PAV deployment\n", encoding="utf-8")
    tools = root / "tools"
    tools.mkdir()
    names = ("deepvariant", "sawfish", "sniffles2", "pbsv", "cutesv", "trgt", "bcftools", "samtools", "hiphase", "hifiasm", "pav", "svim_asm", "minimap2", "bgzip", "tabix", "jasmine", "truvari", "glnexus", "annovar", "vep")
    programs = {}
    for name in names:
        path = tools / (name + ".py")
        path.write_text("#!/usr/bin/env python3\n" + FAKE_TOOL.lstrip("\n"), encoding="utf-8")
        path.chmod(0o755)
        programs[name] = str(path)

    class ExecutableDispatchRunner(CommandRunner):
        def require_executable(self, executable, **kwargs):
            if str(executable) in programs.values():
                # The actual process executable is this verified interpreter;
                # a .py fixture is not falsely declared a Windows PE binary.
                super().require_executable(sys.executable, **kwargs)
                if not Path(executable).is_file():
                    raise FileNotFoundError(executable)
                return Path(executable)
            return super().require_executable(executable, **kwargs)

        def run(self, command, **kwargs):
            args = tuple(map(str, command))
            if args[0] in programs.values():
                args = (sys.executable, "-I", "-B", args[0], args[0], *args[1:])
            safe = {key: value for key, value in os.environ.items() if key.upper() in {"SYSTEMROOT", "WINDIR"}}
            safe.update(PYTHONUTF8="1", PYTHONIOENCODING="utf-8", PYTHONDONTWRITEBYTECODE="1")
            kwargs.update(inherit_env=False, env=safe)
            return super().run(args, **kwargs)

    for module in ("deepvariant", "sawfish", "sniffles2", "pbsv", "cutesv", "sv", "trgt", "hiphase", "hifiasm", "pav", "svim_asm", "jasmine", "bnd_representation", "truvari", "glnexus", "annovar", "vep"):
        monkeypatch.setattr(importlib.import_module("hifivar." + module), "CommandRunner", ExecutableDispatchRunner)
    config = yaml.safe_load((ROOT / "src/hifivar/resources/configs/default.yaml").read_text(encoding="utf-8"))
    config["reference"].update(fasta=str(reference), build="GRCh38")
    config["samples"]["sheet"] = str(sheet)
    config["paths"].update(workdir=str(root / "work"), outdir=str(root / "results"))
    config["small"].update(enabled=True, execution_mode="native", deepvariant_executable=programs["deepvariant"])
    for name in ("sv", "tr", "phasing", "assembly", "assembly_sv", "annotation", "cohort"):
        config[name]["enabled"] = True
    for name in ("sawfish", "sniffles2", "pbsv", "cutesv"):
        config["sv"][name].update(enabled=True, executable=programs[name])
    config["sv"]["finalization"].update(bgzip_executable=programs["bgzip"], tabix_executable=programs["tabix"])
    config["sv"]["harmonization"].update(enabled=True, jasmine_executable=programs["jasmine"], bgzip_executable=programs["bgzip"], tabix_executable=programs["tabix"], truvari_executable=programs["truvari"])
    config["tr"].update(executable=programs["trgt"], bcftools_executable=programs["bcftools"], samtools_executable=programs["samtools"], catalog_reference_build="GRCh38", karyotype="XX")
    config["phasing"].update(executable=programs["hiphase"], tabix_executable=programs["tabix"])
    config["assembly"]["executable"] = programs["hifiasm"]
    for caller in ("pav", "svim_asm"):
        config["assembly_sv"][caller].update(enabled=True, executable=programs[caller], bgzip_executable=programs["bgzip"], tabix_executable=programs["tabix"])
    config["assembly_sv"]["pav"].update(snakefile=str(pav_snakefile), version="2.4.6")
    config["assembly_sv"]["svim_asm"].update(minimap2_executable=programs["minimap2"], samtools_executable=programs["samtools"])
    config["annotation"].update(annovar_enabled=True, vep_enabled=True, overlap_enabled=True,
        annovar_executable=programs["annovar"], annovar_database_root=str(database), annovar_database_version="synthetic-snapshot",
        annovar_version="2020-06-08", annovar_protocols=["refGene"], annovar_operations=["g"],
        vep_executable=programs["vep"], vep_cache_directory=str(cache), vep_cache_version="116")
    config["cohort"].update(cohort_id="TEST_COHORT", overwrite=False)
    config["cohort"]["small_variants"].update(enabled=True, glnexus_executable=programs["glnexus"], bcftools_executable=programs["bcftools"], memory_gb=8)
    config["cohort"]["sv"]["enabled"] = config["cohort"]["tr"]["enabled"] = True
    return SimpleNamespace(root=root, external=external, reference=reference, sheet=sheet, config=config,
                           programs=programs, catalog=catalog, fastq=fastq, monkeypatch=monkeypatch)


def owned_context(bridge):
    config = bridge.config
    effective = bridge.root / "effective.yaml"
    effective.write_text(yaml.safe_dump(config), encoding="utf-8")
    runtime = launcher.runtime_identity(installed_workflow_root())
    identity = {"schema": 2, "run_dir": str(bridge.root), "config": config,
                "sample_sheet": bridge.sheet.read_text(encoding="utf-8"), "runtime": runtime,
                "inputs": [artifact_identity(path) for path in bridge.external.iterdir() if path.is_file()]}
    controls = [effective, bridge.sheet, *(bridge.root / name for name in ("annotation.inputs.tsv", "cohort.inputs.tsv")
                                           if (bridge.root / name).is_file())]
    (bridge.root / "launch_manifest.json").write_text(json.dumps({"identity": identity,
        "control_sha256": {path.name: artifact_identity(path)["sha256"] for path in controls}}), encoding="utf-8")
    attempt = bridge.root / "attempts/bridge_test"
    attempt.mkdir(parents=True)
    (attempt / "status.json").write_text(json.dumps({"status": "DAG_RUNNING", "dry_run": False}), encoding="utf-8")
    context = create_attempt_context(bridge.root, attempt, dry_run=False)
    bridge.monkeypatch.setenv(ATTEMPT_CONTEXT_ENV, str(context))
    return attempt


def job(bridge, **values):
    log = bridge.root / "logs/bridge.log"
    log.parent.mkdir(exist_ok=True)
    base = dict(config=bridge.config, threads=1, resources=SimpleNamespace(mem_mb=1000, runtime_min=1),
                log=[str(log)], wildcards=SimpleNamespace(sample="S1"), params=SimpleNamespace(),
                input=SimpleNamespace(), output=SimpleNamespace())
    base.update(values)
    return SimpleNamespace(**base)


def alignment_input(bridge, **extra):
    return SimpleNamespace(reference=str(bridge.reference), alignment=str(bridge.external / "S1.bam"),
                           alignment_index=str(bridge.external / "S1.bam.bai"), **extra)


def run_bridge(name, snakemake):
    return runpy.run_path(str(SCRIPTS / name), init_globals={"snakemake": snakemake})


def assert_native_receipt(bridge, attempt, track, roles, samples=("S1",), unit_id=None):
    receipt = read_step_receipt(bridge.root, attempt.name, track, samples, unit_id=unit_id)
    assert receipt["status"] == "COMPLETE" and receipt["execution_kind"] == "NATIVE"
    assert set(receipt["outputs"]) == set(roles)
    assert receipt["commands"] and all(row["executed"] and row["returncode"] == 0 for row in receipt["commands"])
    assert any(str(sys.executable) == row["args"][0] for row in receipt["commands"])
    assert all("stdout" not in row and "stderr" not in row for row in receipt["commands"])
    assert receipt["scientific_accuracy"] == "NOT_ASSESSED"
    assert receipt["validation"]["validated"] is True
    return receipt


@pytest.mark.parametrize("name", BRIDGES)
def test_bridge_source_has_preexecution_output_declaration_without_subprocess(name):
    tree = ast.parse((SCRIPTS / name).read_text(encoding="utf-8"))
    steps = [node for node in ast.walk(tree) if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "attempt_step"]
    assert steps and all("outputs" in {key.arg for key in node.keywords} for node in steps)
    assert not any(isinstance(node, ast.Import) and any(alias.name == "subprocess" for alias in node.names) for node in ast.walk(tree))


def small_job(bridge):
    vcf, gvcf = bridge.root / "results/S1.small.vcf.gz", bridge.root / "results/S1.g.vcf.gz"
    return job(bridge, input=alignment_input(bridge), output=SimpleNamespace(vcf=str(vcf), gvcf=str(gvcf)),
               params=SimpleNamespace(region="chr1:1-100", intermediate=str(bridge.root / "work/deepvariant"), tmpdir=str(bridge.root / "tmp")))


def test_deepvariant_bridge_observes_real_processes_and_all_four_outputs(bridge):
    attempt = owned_context(bridge)
    run_bridge("run_deepvariant.py", small_job(bridge))
    receipt = assert_native_receipt(bridge, attempt, "small", ("vcf", "index", "gvcf", "gvcf_index"))
    assert receipt["tool_versions"]["deepvariant"] == "1.10.0"
    assert receipt["metadata"]["region"] == "chr1:1-100"


def test_deepvariant_failure_never_emits_completion(bridge):
    attempt = owned_context(bridge)
    Path(bridge.programs["deepvariant"]).with_suffix(".fail").write_text("fail", encoding="utf-8")
    with pytest.raises(CommandExecutionError):
        run_bridge("run_deepvariant.py", small_job(bridge))
    assert read_step_receipt(bridge.root, attempt.name, "small", ("S1",)) is None


def test_preexisting_deepvariant_artifact_is_not_adopted(bridge):
    attempt = owned_context(bridge)
    snakemake = small_job(bridge)
    path = Path(snakemake.output.vcf)
    path.parent.mkdir()
    path.write_bytes(b"old-output")
    with pytest.raises(OutputValidationError):
        run_bridge("run_deepvariant.py", snakemake)
    assert path.read_bytes() == b"old-output"
    assert read_step_receipt(bridge.root, attempt.name, "small", ("S1",)) is None


@pytest.mark.parametrize("caller", ("sawfish", "sniffles2", "pbsv", "cutesv"))
def test_read_sv_bridge_seals_only_after_caller_and_finalizer_validation(bridge, caller):
    attempt = owned_context(bridge)
    snakemake = job(bridge, wildcards=SimpleNamespace(sample="S1", caller=caller), input=alignment_input(bridge),
                    output=SimpleNamespace(vcf=str(bridge.root / f"results/S1.{caller}.sv.vcf.gz")),
                    params=SimpleNamespace(workdir=str(bridge.root / ("work/" + caller))))
    run_bridge("run_sv_caller.py", snakemake)
    receipt = assert_native_receipt(bridge, attempt, "sv." + caller, ("vcf", "index"))
    assert receipt["metadata"]["evidence_source"] == "READ"
    if caller in {"pbsv", "cutesv"}:
        assert {"bgzip", "tabix"}.issubset(receipt["tool_versions"])


def test_trgt_records_real_catalog_hash_and_independent_spanning_outputs(bridge):
    attempt = owned_context(bridge)
    snakemake = job(bridge, input=alignment_input(bridge, catalog=str(bridge.catalog)),
        output=SimpleNamespace(vcf=str(bridge.root / "results/S1.tr.vcf.gz"), spanning_bam=str(bridge.root / "results/S1.tr.spanning.bam")),
        params=SimpleNamespace(raw_prefix=str(bridge.root / "work/S1.trgt")))
    run_bridge("run_trgt.py", snakemake)
    receipt = assert_native_receipt(bridge, attempt, "tr", ("vcf", "index", "spanning_bam", "spanning_bam_index"))
    assert receipt["metadata"]["catalog_id"] == "sha256:" + hashlib.sha256(bridge.catalog.read_bytes()).hexdigest()


def test_hifiasm_records_three_actual_converted_fastas(bridge):
    attempt = owned_context(bridge)
    snakemake = job(bridge, input=SimpleNamespace(reads=[str(bridge.fastq)]),
                    params=SimpleNamespace(output_prefix=str(bridge.root / "work/S1.asm"), assembly_directory=str(bridge.root / "results/assembly")))
    run_bridge("run_hifiasm.py", snakemake)
    receipt = assert_native_receipt(bridge, attempt, "assembly", ("primary_fasta", "haplotype1_fasta", "haplotype2_fasta"))
    assert receipt["metadata"]["reference_independent"] is True


@pytest.mark.parametrize("caller", ("pav", "svim_asm"))
def test_assembly_sv_receipts_preserve_tool_version_source(bridge, caller):
    attempt = owned_context(bridge)
    run_bridge("run_hifiasm.py", job(bridge, input=SimpleNamespace(reads=[str(bridge.fastq)]),
        params=SimpleNamespace(output_prefix=str(bridge.root / "work/S1.asm"), assembly_directory=str(bridge.root / "results/assembly"))))
    parent = read_step_receipt(bridge.root, attempt.name, "assembly", ("S1",))
    snakemake = job(bridge, input=SimpleNamespace(reference=str(bridge.reference),
        hap1=parent["outputs"]["haplotype1_fasta"]["path"], hap2=parent["outputs"]["haplotype2_fasta"]["path"]),
        params=SimpleNamespace(work=str(bridge.root / ("work/" + caller))),
        output=SimpleNamespace(vcf=str(bridge.root / f"results/S1.{caller}.assembly.sv.vcf.gz")))
    run_bridge("run_" + caller + ".py", snakemake)
    receipt = assert_native_receipt(bridge, attempt, "assembly_sv." + caller, ("vcf", "index"))
    assert receipt["metadata"]["evidence_source"] == "ASSEMBLY"
    binding = receipt["metadata"]["assembly_role_binding"]
    assert [item["role"] for item in binding["assemblies"]] == ["haplotype1", "haplotype2"]
    assert binding["phase_accuracy"] == "NOT_ASSESSED"
    if caller == "pav":
        assert receipt["metadata"]["tool_version_sources"]["pav"] == "USER_DECLARED_NOT_NATIVE_PROBED"


def seed_bgzf(path, payload):
    compressor = zlib.compressobj(wbits=-15)
    compressed = compressor.compress(payload) + compressor.flush()
    size = 18 + len(compressed) + 8
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"\x1f\x8b\x08\x04\0\0\0\0\0\xff\x06\0BC\x02\0" + struct.pack("<H", size - 1)
                     + compressed + struct.pack("<II", zlib.crc32(payload), len(payload)))


def seed_vcf(path, sample, category="sv"):
    extra = ""
    info = "SVTYPE=DEL;END=70;SVLEN=-60"
    if category == "tr":
        extra = '##INFO=<ID=TRID,Number=1,Type=String,Description="ID">\n##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motifs">\n'
        info = "TRID=TR1;MOTIFS=C;END=20"
    elif category == "gvcf":
        extra = '##DeepVariant_version=1.10.0\n##FILTER=<ID=RefCall,Description="Reference">\n##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Depth">\n'
    header = VCF_HEADER + extra + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n"
    text = header + "chr1\t10\tv1\tA\t<DEL>\t.\tPASS\t" + info + "\tGT\t0/1\n"
    if str(path).endswith(".gz"):
        seed_bgzf(path, text.encode())
        seed_bgzf(Path(str(path) + ".tbi"), b"TBI\x01")
    else:
        path.write_text(text, encoding="utf-8")
    return path


def test_hiphase_bridge_keeps_phased_artifacts_separate(bridge):
    small = seed_vcf(bridge.external / "S1.small.vcf.gz", "S1", "small")
    gvcf = seed_vcf(bridge.external / "S1.g.vcf.gz", "S1", "gvcf")
    attempt = owned_context(bridge)
    snakemake = job(bridge, input=alignment_input(bridge, vcf=str(small), gvcf=str(gvcf),
        vcf_index=str(small) + ".tbi", gvcf_index=str(gvcf) + ".tbi"),
        output=SimpleNamespace(vcf=str(bridge.root / "results/S1.phased.vcf.gz")))
    run_bridge("run_hiphase.py", snakemake)
    receipt = assert_native_receipt(bridge, attempt, "phasing", ("vcf", "index"))
    assert receipt["tool_versions"]["hiphase"] == "1.5.0"


@pytest.mark.parametrize("bnd_enabled", [False, True])
def test_harmonization_bridge_keeps_each_actual_comparison_and_evidence(bridge, bnd_enabled):
    if bnd_enabled:
        from hifivar.breakends import BndPolicy, POLICY, POLICY_SHA256
        bridge.config["sv"]["harmonization"]["bnd_representation"] = BndPolicy(POLICY["policy_id"], 1, POLICY_SHA256, POLICY["mode"]).to_dict()
    sources = []
    for caller in ("sawfish", "sniffles2"):
        path = seed_vcf(bridge.external / f"S1.{caller}.sv.vcf.gz", "S1")
        sources.extend((str(path), str(path) + ".tbi"))
    for caller in ("pbsv", "cutesv"):
        bridge.config["sv"][caller]["enabled"] = False
    bridge.config["assembly_sv"]["enabled"] = False
    bridge.config["sv"]["enabled"] = False
    bridge.config["sv"]["harmonization"]["input_vcfs"] = {"S1": {
        caller: str(bridge.external / f"S1.{caller}.sv.vcf.gz") for caller in ("sawfish", "sniffles2")}}
    attempt = owned_context(bridge)
    directory = bridge.root / "results/harmonized/S1"
    snakemake = job(bridge, input=SimpleNamespace(sources=sources),
        output=SimpleNamespace(vcf=str(directory / "S1.harmonized.sv.vcf.gz"),
                               index=str(directory / "S1.harmonized.sv.vcf.gz.tbi"),
                               evidence=str(directory / "S1.sv.evidence.tsv"),
                               provenance=str(directory / "S1.phase9.provenance.json"),
                               bnd_bundle=str(directory / "bndrepr-v1")),
        params=SimpleNamespace(work_root=str(bridge.root / "work/harmonization")))
    run_bridge("run_harmonization.py", snakemake)
    receipt = read_step_receipt(bridge.root, attempt.name, "sv.harmonization", ("S1",))
    assert len(receipt["outputs"]) == (24 if bnd_enabled else 14)
    assert receipt["metadata"]["source_order"] == ["sawfish", "sniffles2"]
    assert receipt["metadata"]["caller_count_is_confidence"] is False
    assert receipt["metadata"]["harmonization_is_truth"] is False
    assert receipt["metadata"]["bnd_representation"]["enabled"] is bnd_enabled
    if bnd_enabled:
        bundle = directory / "bndrepr-v1"
        manifest = json.loads((bundle / "manifest.json").read_text())
        assert manifest["status"] == "completed" and manifest["summary"]["round_trip"]["status"] == "PASS"
        assert receipt["outputs"]["bnd_manifest"]["sha256"] == hashlib.sha256((bundle / "manifest.json").read_bytes()).hexdigest()
        assert receipt["outputs"]["vcf"]["path"] != receipt["outputs"]["bnd_vcf"]["path"]
        assert manifest["standalone_derived_VCF_is_lossless"] is False
    else:
        assert not (directory / "bndrepr-v1").exists()


def annotation_job(bridge, category):
    source = seed_vcf(bridge.external / f"S1.{category}.vcf", "S1", category)
    manifest = bridge.root / "annotation.inputs.tsv"
    previous = manifest.read_text(encoding="utf-8") if manifest.exists() else "sample\tvariant_category\tsource_vcf\tsource_tool\tsource_variant_ids\n"
    manifest.write_text(previous + f"S1\t{category}\t{source}\tsynthetic-caller\tv1\n", encoding="utf-8")
    bridge.config["annotation"]["input_manifest"] = str(manifest)
    return job(bridge, wildcards=SimpleNamespace(sample="S1", category=category),
               input=SimpleNamespace(manifest=str(manifest)),
               output=SimpleNamespace(result=str(bridge.root / ("results/annotation/" + category))))


def test_annotation_backend_receipts_are_separate_and_category_scoped(bridge):
    gene = bridge.external / "genes.bed"
    gene.write_text("chr1\t1\t80\tGENE1\n", encoding="utf-8")
    bridge.config["annotation"].update(gene_bed=str(gene), gene_version="test-declared")
    small, sv = annotation_job(bridge, "small"), annotation_job(bridge, "sv")
    attempt = owned_context(bridge)
    run_bridge("run_annotation.py", small)
    run_bridge("run_annotation.py", sv)
    for category in ("small", "sv"):
        annovar = assert_native_receipt(bridge, attempt, "annotation.annovar", ("tsv", "vcf"), unit_id=category)
        vep = assert_native_receipt(bridge, attempt, "annotation.vep", ("tsv",), unit_id=category)
        assert annovar["metadata"]["tool_version_sources"]["annovar"] == "USER_DECLARED_NOT_NATIVE_PROBED"
        assert vep["metadata"]["tool_version_sources"]["vep"] == "COMMAND_PROBED"
        assert annovar["metadata"]["source_variant_ids_filter_applied"] is False
    overlap = read_step_receipt(bridge.root, attempt.name, "annotation.overlap", ("S1",), unit_id="sv")
    assert overlap["execution_kind"] == "IN_PROCESS" and overlap["commands"] == []
    assert set(overlap["outputs"]) == {"tsv"}


def test_annotation_failed_backend_never_receives_false_receipt(bridge):
    snakemake = annotation_job(bridge, "small")
    attempt = owned_context(bridge)
    Path(bridge.programs["vep"]).with_suffix(".fail").write_text("fail", encoding="utf-8")
    with pytest.raises(CommandExecutionError):
        run_bridge("run_annotation.py", snakemake)
    assert read_step_receipt(bridge.root, attempt.name, "annotation.annovar", ("S1",), unit_id="small") is not None
    assert read_step_receipt(bridge.root, attempt.name, "annotation.vep", ("S1",), unit_id="small") is None


def cohort_job(bridge, track):
    header = "sample\ttrack\tstate\tsource_path\tindex_path\tsource_tool\tsource_version\treference_build\tcatalog_id\n"
    rows = []
    catalog_id = "sha256:" + hashlib.sha256(bridge.catalog.read_bytes()).hexdigest()
    for sample in ("S2", "S1"):
        path = seed_vcf(bridge.external / f"{sample}.{track}.vcf.gz", sample, "gvcf" if track == "small_variants" else track)
        rows.append(f"{sample}\t{track}\tCALLED\t{path}\t{path}.tbi\tsynthetic-caller\t1\tGRCh38\t{catalog_id if track == 'tr' else ''}\n")
    manifest = bridge.root / "cohort.inputs.tsv"
    manifest.write_text(header + "".join(rows), encoding="utf-8")
    bridge.config["cohort"]["input_manifest"] = str(manifest)
    directory = bridge.root / ("results/cohort/" + track)
    outputs = {"result": str(directory / "track_result.json")}
    if track == "small_variants":
        outputs.update(bcf=str(directory / "TEST_COHORT.small.bcf"), vcf=str(directory / "TEST_COHORT.small.vcf.gz"))
    elif track == "sv":
        outputs.update(sites=str(directory / "sites.tsv"), matrix=str(directory / "matrix.tsv"))
    else:
        outputs.update(loci=str(directory / "loci.tsv"), matrix=str(directory / "matrix.tsv"))
    return job(bridge, input=SimpleNamespace(manifest=str(manifest)), output=SimpleNamespace(**outputs))


def test_glnexus_bridge_preserves_declared_order_and_actual_named_output_qc(bridge):
    snakemake = cohort_job(bridge, "small_variants")
    attempt = owned_context(bridge)
    run_bridge("run_cohort_small.py", snakemake)
    receipt = assert_native_receipt(bridge, attempt, "cohort.small_variants",
                                    ("bcf", "bcf_index", "vcf", "index", "provenance"), samples=("S2", "S1"))
    assert receipt["metadata"]["declared_sample_order"] == ["S2", "S1"]
    assert receipt["metadata"]["output_sample_order"] == ["S1", "S2"]
    assert receipt["metadata"]["memory_gb"] == 8


@pytest.mark.parametrize("track,roles", [("sv", {"sites_tsv", "sample_matrix_tsv", "provenance"}),
                                         ("tr", {"loci_tsv", "sample_matrix_tsv", "provenance"})])
def test_pure_cohort_table_bridges_are_in_process_not_native_calling(bridge, track, roles):
    snakemake = cohort_job(bridge, track)
    attempt = owned_context(bridge)
    run_bridge("run_cohort_" + track + ".py", snakemake)
    receipt = read_step_receipt(bridge.root, attempt.name, "cohort." + track, ("S2", "S1"))
    assert receipt["execution_kind"] == "IN_PROCESS" and receipt["commands"] == []
    assert set(receipt["outputs"]) == roles
    assert receipt["metadata"]["input_states"] == {"S2": "CALLED", "S1": "CALLED"}


def test_standalone_snakemake_bridge_remains_unlinked(bridge):
    bridge.monkeypatch.delenv(ATTEMPT_CONTEXT_ENV, raising=False)
    snakemake = small_job(bridge)
    run_bridge("run_deepvariant.py", snakemake)
    assert Path(snakemake.output.vcf).is_file()
    assert not (bridge.root / "attempts").exists()


@pytest.mark.parametrize('caller',('pav','svim_asm'))
def test_assembly_sv_bridge_rejects_unlinked_external_fastas_before_execution(bridge,caller):
    attempt=owned_context(bridge)
    output=bridge.root/f'results/S1.{caller}.assembly.sv.vcf.gz'
    sm=job(bridge,input=SimpleNamespace(reference=str(bridge.reference),hap1=str(bridge.external/'hap1.fa'),hap2=str(bridge.external/'hap2.fa')),
        params=SimpleNamespace(work=str(bridge.root/('work/'+caller))),output=SimpleNamespace(vcf=str(output)))
    with pytest.raises(InputValidationError,match='current-attempt assembly receipt'):
        run_bridge('run_'+caller+'.py',sm)
    assert not output.exists()
    assert read_step_receipt(bridge.root,attempt.name,'assembly_sv.'+caller,('S1',)) is None
