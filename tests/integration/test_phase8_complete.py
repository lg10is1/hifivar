from pathlib import Path
import struct
import zlib
import pytest
from hifivar.assembly import AssemblyArtifact, AssemblyRole, HaplotypeAssemblyArtifact
from hifivar.assembly_sv import create_assembly_sv_artifact
from hifivar.command import format_command
from hifivar.context import AnalysisContext
from hifivar.exceptions import InputValidationError
from hifivar.pav import PavCommandPlan, PavResult, PavResultStatus
from hifivar.phase8 import run_phase8
from hifivar.reference import ReferenceGenome
from hifivar.sample_sheet import SampleRecord
from hifivar.sample import InputDataset, Sample
from hifivar.svim_asm import SvimAsmCommandPlan, SvimAsmResult, SvimAsmResultStatus

def bgzf(path, payload):
    obj=zlib.compressobj(level=6,wbits=-15); data=obj.compress(payload)+obj.flush(); size=18+len(data)+8
    head=b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"+struct.pack("<H",size-1)
    path.parent.mkdir(parents=True,exist_ok=True); path.write_bytes(head+data+struct.pack("<II",zlib.crc32(payload),len(payload)))

def output(request):
    text=("##fileformat=VCFv4.3\n##contig=<ID=chr1,length=4>\n"
          '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
          f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{request.sample_id}\n")
    bgzf(request.output_vcf,text.encode()); bgzf(request.output_index,b"TBI\x01")

class FakePav:
    def run(self, request, **kwargs):
        output(request); cmd=("pav-adapter","--sample",request.sample_id)
        artifact=create_assembly_sv_artifact(request,raw_vcf=request.output_vcf,intermediate_files=(),
            caller_version="2.4.6",backend="fake",commands=(cmd,))
        return PavResult(request,PavResultStatus.COMPLETED,PavCommandPlan(cmd,format_command(cmd)),"2.4.6","8.0",.1,artifact)

class FakeSvim:
    def run(self, request, **kwargs):
        raw=request.work_directory/"native"/"variants.vcf"; raw.parent.mkdir(parents=True)
        raw.write_text("##fileformat=VCFv4.3\n"); output(request)
        cmd=SvimAsmCommandPlan("svim_asm",("svim-asm","diploid"))
        artifact=create_assembly_sv_artifact(request,raw_vcf=raw,intermediate_files=(raw,),
            caller_version="1.0.3",backend="fake",commands=(cmd.args,))
        return SvimAsmResult(request,SvimAsmResultStatus.COMPLETED,(cmd,),{"svim_asm":"1.0.3"},.1,artifact)

def setup(tmp_path):
    ref=tmp_path/"ref.fa"; ref.write_text(">chr1\nACGT\n"); Path(f"{ref}.fai").write_text("chr1\t4\t6\t4\t5\n")
    reads=tmp_path/"reads.fastq"; reads.write_text("@r\nACGT\n+\n!!!!\n")
    sample=Sample("S1",InputDataset.from_files([reads]))
    reference=ReferenceGenome.from_fasta(ref,build="GRCh38")
    context=AnalysisContext(reference,(SampleRecord(sample),),{})
    items=[]
    for role,name in ((AssemblyRole.HAPLOTYPE1,"hap1"),(AssemblyRole.HAPLOTYPE2,"hap2")):
        path=tmp_path/f"S1.{name}.fa"; path.write_text(f">{name}\nACGT\n")
        gfa=tmp_path/f"{name}.gfa"; gfa.write_text(f"S\t{name}\tACGT\n")
        items.append(HaplotypeAssemblyArtifact("S1",role,path,gfa,"0.25",path.stat().st_size))
    assembly=AssemblyArtifact("S1",tuple(x.source_gfa for x in items),tuple(items),(),"0.25",("hifiasm",))
    config={"assembly_sv":{"enabled":True,"overwrite":False,
        "pav":{"enabled":True,"executable":"snakemake","snakefile":str(tmp_path/"Snakefile"),"version":"2.4.6","threads":2,"memory_mb":1000,"runtime_minutes":30},
        "svim_asm":{"enabled":True,"executable":"svim-asm","minimap2_executable":"minimap2","samtools_executable":"samtools","bgzip_executable":"bgzip","tabix_executable":"tabix","threads":2,"memory_mb":1000,"runtime_minutes":30}}}
    return context,assembly,config

def test_phase8_two_independent_artifacts_and_provenance(tmp_path):
    context,assembly,config=setup(tmp_path)
    report=run_phase8(context,{"S1":assembly},output_directory=tmp_path/"out",work_directory=tmp_path/"work",
        config=config,pav_wrapper=FakePav(),svim_wrapper=FakeSvim())
    artifacts=report.sample_results[0].collection.artifacts
    assert len(artifacts)==2 and artifacts[0].raw_vcf != artifacts[1].raw_vcf
    assert all(not item.harmonized for item in artifacts)
    path=report.write_json(tmp_path/"phase8.json")
    assert path.exists() and '"evidence_source": "assembly"' in path.read_text()

def test_phase8_missing_assembly_is_not_no_calls(tmp_path):
    context,_,config=setup(tmp_path)
    with pytest.raises(InputValidationError):
        run_phase8(context,{},output_directory=tmp_path/"out",work_directory=tmp_path/"work",config=config,
                   pav_wrapper=FakePav(),svim_wrapper=FakeSvim())
