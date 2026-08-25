from pathlib import Path
import struct
import zlib
import pytest
import hifivar.jasmine as jasmine_module
from hifivar.assembly_sv import SVEvidenceSource
from hifivar.command import CommandResult
from hifivar.exceptions import CommandExecutionError, InputValidationError, OutputValidationError
from hifivar.harmonization import EvidenceClass, EvidenceRunStatus, SVEvidenceSourceArtifact, SVHarmonizationRequest, iter_sv_evidence, write_evidence_table
from hifivar.jasmine import JasmineResultStatus, JasmineWrapper
from hifivar.reference import ReferenceGenome
from hifivar.truvari import TruvariRequest, TruvariResultStatus, TruvariWrapper

def bgzf(path,payload):
    obj=zlib.compressobj(level=6,wbits=-15); data=obj.compress(payload)+obj.flush(); size=18+len(data)+8
    head=b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"+struct.pack("<H",size-1)
    path.parent.mkdir(parents=True,exist_ok=True); path.write_bytes(head+data+struct.pack("<II",zlib.crc32(payload),len(payload)))

def vcf(sample="S1",records=""):
    return ("##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n"
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"+records).encode()

@pytest.fixture
def reference(tmp_path):
    path=tmp_path/"ref.fa"; path.write_text(">chr1\n"+"A"*1000+"\n")
    Path(f"{path}.fai").write_text("chr1\t1000\t6\t1000\t1001\n")
    return ReferenceGenome.from_fasta(path,build="GRCh38")

def source(tmp_path,reference,caller,kind,records,status=EvidenceRunStatus.COMPLETED,haps=()):
    path=tmp_path/f"{caller}.vcf.gz"; bgzf(path,vcf(records=records)); bgzf(Path(f"{path}.tbi"),b"TBI\x01")
    return SVEvidenceSourceArtifact("S1",kind,caller,path,Path(f"{path}.tbi"),status,haps)

def request(tmp_path,reference,sources):
    return SVHarmonizationRequest("S1",reference,tuple(sources),tmp_path/"work",
        tmp_path/"out"/"S1.harmonized.sv.vcf.gz",tmp_path/"out"/"S1.sv.evidence.tsv",500,"linear")

def test_streaming_normalization_preserves_bnd_ins_and_unknown(tmp_path,reference):
    records=("chr1\t10\tb1\tA\tA]chr1:30]\t.\tPASS\tSVTYPE=BND;CIPOS=-2,2\tGT\t0/1\n"
             "chr1\t20\ti1\tA\tATGC\t.\tPASS\tSVTYPE=INS;SVLEN=3\tGT\t0/1\n"
             "chr1\t30\tc1\tA\t<CPLX>\t.\tPASS\tSVTYPE=CPX\tGT\t0/1\n")
    rows=list(iter_sv_evidence(source(tmp_path,reference,"sawfish",SVEvidenceSource.READ,records)))
    assert rows[0].svtype=="BND" and dict(rows[0].native_info)["CIPOS"]=="-2,2"
    assert rows[1].insertion_sequence_length==3
    assert rows[2].svtype=="UNRESOLVED" and rows[2].native_svtype=="CPX" and rows[2].unresolved

def test_partial_statuses_are_not_no_calls(tmp_path,reference):
    read=source(tmp_path,reference,"sawfish",SVEvidenceSource.READ,"")
    disabled=SVEvidenceSourceArtifact("S1",SVEvidenceSource.READ,"pbsv",None,None,EvidenceRunStatus.DISABLED)
    failed=SVEvidenceSourceArtifact("S1",SVEvidenceSource.ASSEMBLY,"pav",None,None,EvidenceRunStatus.FAILED,error="exit 1")
    item=request(tmp_path,reference,(read,disabled,failed))
    assert [x.caller for x in item.runnable_sources]==["sawfish"]
    assert disabled.status is not EvidenceRunStatus.NO_CALLS and failed.status is not EvidenceRunStatus.NO_CALLS
    with pytest.raises(InputValidationError):
        request(tmp_path,reference,(disabled,failed))

class FakeRunner:
    def __init__(self,request): self.request=request
    def require_executable(self,executable): return Path(executable)
    def run(self,command,**kwargs):
        args=tuple(str(x) for x in command); dry=kwargs.get("dry_run",False)
        if "--version" in args or args[-1:] == ("version",):
            return CommandResult(args,0,f"{args[0]} 1.1.5","",.1,None,True)
        if not dry and args[0]=="jasmine":
            raw=self.request.work_directory/"S1.jasmine.raw.vcf"
            raw.write_bytes(vcf(records="chr1\t10\tJ1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=20;IDLIST=sawfish:r1,pav:a1\tGT\t0/1\n"))
        elif not dry and args[0]=="bgzip":
            bgzf(Path(kwargs["stdout_path"]),Path(args[-1]).read_bytes())
        elif not dry and args[0]=="tabix":
            bgzf(self.request.output_index,b"TBI\x01")
        elif not dry and args[:2]==("truvari","bench"):
            out=Path(args[args.index("-o")+1]); out.mkdir(parents=True)
            (out/"summary.json").write_text('{"comparison_only": true}\n')
        return CommandResult(args,0,"","",.1,None,not dry)

def test_jasmine_evidence_class_is_not_truth(tmp_path,reference):
    read=source(tmp_path,reference,"sawfish",SVEvidenceSource.READ,
                "chr1\t10\tr1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=20\tGT\t0/1\n")
    assembly=source(tmp_path,reference,"pav",SVEvidenceSource.ASSEMBLY,
                    "chr1\t11\ta1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=20\tGT\t0/1\n",haps=("haplotype1",))
    item=request(tmp_path,reference,(read,assembly))
    result=JasmineWrapper(runner=FakeRunner(item)).run(item)
    assert result.status is JasmineResultStatus.COMPLETED and result.artifact is not None
    table=item.evidence_table.read_text()
    assert EvidenceClass.READ_AND_ASSEMBLY.value in table
    assert "sawfish:r1,pav:a1" in table and "HIGH_CONFIDENCE" not in table
    assert result.artifact.truth_label is None
    assert read.vcf_path.exists() and assembly.vcf_path.exists()

def test_jasmine_command_parameters_are_explicit_and_dry_run_has_no_writes(tmp_path,reference):
    read=source(tmp_path,reference,"sawfish",SVEvidenceSource.READ,"")
    item=request(tmp_path,reference,(read,))
    wrapper=JasmineWrapper(runner=FakeRunner(item))
    command=wrapper.plan_commands(item)[0].args
    assert "max_dist=500" in command and not any("confidence" in x.lower() for x in command)
    assert wrapper.run(item,dry_run=True).status is JasmineResultStatus.PLANNED
    assert not item.work_directory.exists()

def test_truvari_is_comparison_only(tmp_path,reference):
    base=source(tmp_path,reference,"base",SVEvidenceSource.READ,"")
    comp=source(tmp_path,reference,"comp",SVEvidenceSource.ASSEMBLY,"")
    item=TruvariRequest("S1",reference,base.vcf_path,comp.vcf_path,tmp_path/"truvari")
    result=TruvariWrapper(runner=FakeRunner(item)).run(item)
    assert result.status is TruvariResultStatus.COMPLETED
    assert result.to_dict()["interpretation"]=="comparison_only_not_truth"
    assert result.command[1]=="bench"
def test_no_calls_empty_vcf_is_distinct_from_failure(tmp_path,reference):
    empty=source(tmp_path,reference,"sawfish",SVEvidenceSource.READ,"",status=EvidenceRunStatus.NO_CALLS)
    assert list(iter_sv_evidence(empty)) == []
    assert request(tmp_path,reference,(empty,)).runnable_sources == (empty,)

def test_unmapped_jasmine_provenance_is_explicitly_unresolved(tmp_path,reference):
    read=source(tmp_path,reference,"sawfish",SVEvidenceSource.READ,"")
    item=request(tmp_path,reference,(read,))
    merged=tmp_path/"merged.vcf"
    merged.write_bytes(vcf(records="chr1\t10\tJ1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=20;IDLIST=native1\tGT\t0/1\n"))
    write_evidence_table(item,merged)
    assert EvidenceClass.UNRESOLVED.value in item.evidence_table.read_text()

def test_harmonization_refuses_silent_overwrite(tmp_path,reference):
    read=source(tmp_path,reference,"sawfish",SVEvidenceSource.READ,"")
    output=tmp_path/"out"/"S1.harmonized.sv.vcf.gz"
    output.parent.mkdir()
    output.write_bytes(b"existing")
    with pytest.raises(OutputValidationError):
        request(tmp_path,reference,(read,))

class FailedRunner(FakeRunner):
    def run(self,command,**kwargs):
        args=tuple(str(x) for x in command)
        if args and args[0]=="jasmine":
            raise CommandExecutionError("jasmine exit code 1")
        return super().run(command,**kwargs)

def test_external_failure_is_not_silently_converted_to_no_calls(tmp_path,reference):
    read=source(tmp_path,reference,"sawfish",SVEvidenceSource.READ,"")
    item=request(tmp_path,reference,(read,))
    with pytest.raises(CommandExecutionError,match="exit code 1"):
        JasmineWrapper(runner=FailedRunner(item)).run(item)


class LauncherRunner:
    def __init__(self, launcher):
        self.launcher = launcher

    def require_executable(self, executable):
        if executable == "jasmine":
            return self.launcher
        return Path(executable)


def test_jasmine_launcher_without_shebang_is_run_through_bash(tmp_path):
    launcher = tmp_path / "jasmine"
    launcher.write_text("JAR_DIR=/opt/jasmine\njava -jar jasmine.jar\n")
    wrapper = JasmineWrapper(runner=LauncherRunner(launcher))
    assert wrapper._resolve_jasmine_prefix() == (
        "bash",
        str(launcher.absolute()),
    )


class UnsortedJasmineRunner(FakeRunner):
    def run(self, command, **kwargs):
        args = tuple(str(item) for item in command)
        if (
            not kwargs.get("dry_run", False)
            and args[0] == "jasmine"
            and "--version" not in args
        ):
            raw = self.request.work_directory / "S1.jasmine.raw.vcf"
            raw.write_text(
                "##fileformat=VCFv4.3\n"
                "##contig=<ID=chr1,length=1000>\n"
                "##contig=<ID=chr2,length=1000>\n"
                '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                "chr2\t5\tJ2\tA\t<DEL>\t.\tPASS\t"
                "SVTYPE=DEL;END=8;IDLIST=sawfish:r2\tAD:GT\t3,2:0/1\n"
                "chr1\t10\tJ1\tA\t<DEL>\t.\tPASS\t"
                "SVTYPE=DEL;END=20;IDLIST=sawfish:r1\tAD:GT\t4,3:0/1\n"
            )
            return CommandResult(args, 0, "", "", .1, None, True)
        return super().run(command, **kwargs)


def test_jasmine_unzips_inputs_restores_format_and_sorts_by_contig(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(jasmine_module, "_SORT_CHUNK_RECORDS", 1)
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + "A" * 1000 + "\n>chr2\n" + "A" * 1000 + "\n")
    Path(f"{fasta}.fai").write_text(
        "chr1\t1000\t6\t1000\t1001\nchr2\t1000\t1013\t1000\t1001\n"
    )
    ref = ReferenceGenome.from_fasta(fasta, build="GRCh38")
    source_vcf = tmp_path / "sawfish.vcf.gz"
    payload = (
        "##fileformat=VCFv4.3\n"
        "##contig=<ID=chr1,length=1000>\n"
        "##contig=<ID=chr2,length=1000>\n"
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t10\tr1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=20\tAD:GT\t4,3:0/1\n"
    ).encode()
    bgzf(source_vcf, payload)
    bgzf(Path(f"{source_vcf}.tbi"), b"TBI\x01")
    evidence = SVEvidenceSourceArtifact(
        "S1", SVEvidenceSource.READ, "sawfish", source_vcf,
        Path(f"{source_vcf}.tbi"), EvidenceRunStatus.COMPLETED,
    )
    item = request(tmp_path, ref, (evidence,))
    result = JasmineWrapper(runner=UnsortedJasmineRunner(item)).run(item)
    assert result.artifact is not None
    raw = item.work_directory / "S1.jasmine.raw.vcf"
    sorted_vcf = item.work_directory / "S1.jasmine.sorted.vcf"
    assert raw.exists() and sorted_vcf.exists()
    sorted_text = sorted_vcf.read_text()
    assert "##FORMAT=<ID=AD" in sorted_text
    records = [line for line in sorted_text.splitlines() if not line.startswith("#")]
    assert [line.split("\t", 1)[0] for line in records] == ["chr1", "chr2"]
    listing = (item.work_directory / "jasmine.inputs.txt").read_text()
    assert ".unzipped.vcf" in listing and ".vcf.gz" not in listing


class TruvariVersionRunner:
    def __init__(self, *, fail_primary=False):
        self.fail_primary = fail_primary
        self.calls = []

    def require_executable(self, executable):
        return Path(executable)

    def run(self, command, **kwargs):
        args = tuple(str(item) for item in command)
        self.calls.append(args)
        if args[-1] == "version" and self.fail_primary:
            raise CommandExecutionError("unsupported version command")
        return CommandResult(args, 0, "Truvari v5.4.0", "", .1, None, True)


def test_truvari_uses_version_subcommand_and_legacy_fallback():
    primary = TruvariVersionRunner()
    assert TruvariWrapper(runner=primary).detect_version() == "5.4.0"
    assert primary.calls == [("truvari", "version")]

    fallback = TruvariVersionRunner(fail_primary=True)
    assert TruvariWrapper(runner=fallback).detect_version() == "5.4.0"
    assert fallback.calls == [
        ("truvari", "version"),
        ("truvari", "--version"),
    ]
