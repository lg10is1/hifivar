import json
import gzip
from pathlib import Path
import struct
import zlib
import pytest
from hifivar.assembly import AssemblyRole, HaplotypeAssemblyArtifact
from hifivar.assembly_sv import AssemblySvCaller, AssemblySvRequest, AssemblySvResources
from hifivar.command import CommandResult
from hifivar.exceptions import CommandExecutionError, InputValidationError, OutputValidationError, ToolVersionError
from hifivar.pav import PavResultStatus, PavWrapper
from hifivar.reference import ReferenceGenome
from hifivar.svim_asm import SvimAsmResultStatus, SvimAsmWrapper

def write_bgzf(path: Path, payload: bytes) -> None:
    obj = zlib.compressobj(level=6, wbits=-15)
    data = obj.compress(payload) + obj.flush()
    size = 18 + len(data) + 8
    head = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", size - 1)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(head + data + struct.pack("<II", zlib.crc32(payload), len(payload)))

def vcf(sample="S1"):
    return ("##fileformat=VCFv4.3\n##contig=<ID=chr1,length=4>\n"
            "##source=PAV 2.4.6.0\n"
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
            '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length">\n'
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
            "chr1\t1\tsnv1\tA\tC\t.\tPASS\tSVTYPE=SNV\tGT\t0/1\n"
            "chr1\t2\tins49\tA\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=49\tGT\t0/1\n"
            "chr1\t2\tins50\tA\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=50\tGT\t0/1\n"
            "chr1\t2\tdel50\tA\t<DEL>\t.\tCOMPOUND\tSVTYPE=DEL;SVLEN=-50\tGT\t0/1\n"
            "chr1\t2\tinv1\tA\t<INV>\t.\tPASS\tSVTYPE=INV;SVLEN=100\tGT\t0/1\n").encode()

@pytest.fixture
def reference(tmp_path):
    fasta = tmp_path / ("ref_" + chr(0x53C2) + ".fa")
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    return ReferenceGenome.from_fasta(fasta, build="GRCh38")

@pytest.fixture
def assemblies(tmp_path):
    items = []
    for role, name in ((AssemblyRole.HAPLOTYPE1, "hap1"), (AssemblyRole.HAPLOTYPE2, "hap2")):
        fasta = tmp_path / ("asm_" + chr(0x7EC4)) / f"S1.{name}.fa"
        fasta.parent.mkdir(parents=True, exist_ok=True)
        fasta.write_text(f">{name}\nACGT\n", encoding="utf-8")
        gfa = tmp_path / f"{name}.gfa"; gfa.write_text(f"S\t{name}\tACGT\n", encoding="utf-8")
        items.append(HaplotypeAssemblyArtifact("S1", role, fasta, gfa, "0.25.0-r726", fasta.stat().st_size))
    return tuple(items)

def request(tmp_path, reference, assemblies, caller, *, overwrite=False):
    return AssemblySvRequest("S1", caller, reference, assemblies, tmp_path / caller.value,
        tmp_path / ("result_" + chr(0x7ED3)) / f"S1.{caller.value}.assembly.sv.vcf.gz",
        AssemblySvResources(2, 1000, 30), overwrite)

class FakeRunner:
    def __init__(self, request): self.request, self.calls = request, []
    def require_executable(self, executable): return Path(executable)
    def run(self, command, **kwargs):
        args = tuple(str(x) for x in command); self.calls.append((args, kwargs))
        dry = kwargs.get("dry_run", False)
        if "--version" in args: return CommandResult(args, 0, f"{args[0]} 1.21.0", "", .1, None, True)
        if not dry and args[0] == "snakemake":
            raw = self.request.work_directory / "S1.vcf.gz"
            write_bgzf(raw, vcf()); write_bgzf(Path(f"{raw}.tbi"), b"TBI\x01")
        elif not dry and args[0] == "minimap2": Path(kwargs["stdout_path"]).write_text("@HD\tVN:1.6\n")
        elif not dry and args[:2] == ("samtools", "sort"): Path(args[args.index("-o") + 1]).write_bytes(b"BAM")
        elif not dry and args[:2] == ("samtools", "index"): Path(f"{args[-1]}.bai").write_bytes(b"BAI")
        elif not dry and args[0] == "svim-asm": (self.request.work_directory / "native" / "variants.vcf").write_bytes(vcf())
        elif not dry and args[0] == "bgzip": write_bgzf(Path(kwargs["stdout_path"]), Path(args[-1]).read_bytes())
        elif not dry and args[0] == "tabix": write_bgzf(self.request.output_index, b"TBI\x01")
        return CommandResult(args, 0, "", "", .1, None, not dry)

def test_pav_dry_run_and_complete(tmp_path, reference, assemblies):
    snakefile = tmp_path / "pav_site" / "Snakefile"
    snakefile.parent.mkdir(); snakefile.write_text("# PAV\n")
    item = request(tmp_path, reference, assemblies, AssemblySvCaller.PAV)
    assert PavWrapper(snakefile=snakefile, runner=FakeRunner(item)).run(item, dry_run=True).status is PavResultStatus.PLANNED
    assert not item.work_directory.exists()
    result = PavWrapper(snakefile=snakefile, pav_version="2.4.6", runner=FakeRunner(item)).run(item)
    assert result.status is PavResultStatus.COMPLETED and result.artifact is not None
    assert result.artifact.raw_vcf == item.work_directory / "S1.vcf.gz"
    assert result.artifact.raw_vcf != result.artifact.vcf_path and not result.artifact.harmonized
    assert result.command.args == (
        "snakemake",
        "--snakefile",
        str(snakefile.absolute()),
        "--directory",
        str(item.work_directory.absolute()),
        "--cores",
        "2",
    )
    assert result.pav_version == "2.4.6"
    assert result.pav_version_source == "config"
    assert result.to_dict()["pav_version_source"] == "config"
    assert result.selection is not None
    assert result.selection.total_records == 5
    assert result.selection.selected_records == 3
    assert result.selection.policy.startswith("PAV_2.4.6_VARTYPE")
    assert result.finalizer_versions == {"bgzip": "1.21.0", "tabix": "1.21.0"}
    with gzip.open(result.artifact.raw_vcf, "rt", encoding="utf-8") as handle:
        raw_text = handle.read()
    with gzip.open(result.artifact.vcf_path, "rt", encoding="utf-8") as handle:
        selected_text = handle.read()
    assert "snv1" in raw_text and "ins49" in raw_text
    assert "snv1" not in selected_text and "ins49" not in selected_text
    for record_id in ("ins50", "del50", "inv1"):
        assert record_id in selected_text
    assert "##hifivar_pav_sv_selection=PAV_2.4.6_VARTYPE" in selected_text
    assert "COMPOUND" in selected_text
    table = (item.work_directory / "assemblies.tsv").read_text()
    assert "HAP_h1" in table and "HAP_h2" in table

def test_pav_overwrite_replaces_owned_inputs_atomically(tmp_path, reference, assemblies):
    snakefile = tmp_path / "pav_site" / "Snakefile"
    snakefile.parent.mkdir(); snakefile.write_text("# PAV\n")
    item = request(tmp_path, reference, assemblies, AssemblySvCaller.PAV, overwrite=True)
    item.work_directory.mkdir(parents=True)
    config = item.work_directory / "config.json"
    table = item.work_directory / "assemblies.tsv"
    config.write_text("stale config\n", encoding="utf-8")
    table.write_text("stale table\n", encoding="utf-8")

    PavWrapper(snakefile=snakefile, pav_version="2.4.6", runner=FakeRunner(item))._prepare(item)

    assert "stale" not in config.read_text(encoding="utf-8")
    assert json.loads(config.read_text(encoding="utf-8"))["reference"] == str(reference.fasta.absolute())
    assert table.read_text(encoding="utf-8").startswith("NAME\tHAP_h1\tHAP_h2\nS1\t")

def test_pav_real_execution_requires_explicit_version(tmp_path, reference, assemblies):
    snakefile = tmp_path / "pav_site" / "Snakefile"
    snakefile.parent.mkdir(); snakefile.write_text("# PAV\n")
    item = request(tmp_path, reference, assemblies, AssemblySvCaller.PAV)
    wrapper = PavWrapper(snakefile=snakefile, runner=FakeRunner(item))

    planned = wrapper.run(item, dry_run=True)
    assert planned.pav_version_source == "unresolved"
    with pytest.raises(ToolVersionError, match="explicit numeric release"):
        wrapper.run(item)
    assert not item.work_directory.exists()

def test_svim_diploid_chain_and_complete(tmp_path, reference, assemblies):
    item = request(tmp_path, reference, assemblies, AssemblySvCaller.SVIM_ASM)
    runner = FakeRunner(item); wrapper = SvimAsmWrapper(runner=runner)
    caller = next(x for x in wrapper.plan_commands(item) if x.step == "svim_asm")
    assert caller.args[1] == "diploid"
    assert "haplotype1.sorted.bam" in caller.args[-3] and "haplotype2.sorted.bam" in caller.args[-2]
    result = wrapper.run(item)
    assert result.status is SvimAsmResultStatus.COMPLETED
    assert result.artifact is not None and result.artifact.raw_vcf.name == "variants.vcf"

def test_haploid_and_invalid_contracts(tmp_path, reference, assemblies):
    item = request(tmp_path, reference, assemblies[:1], AssemblySvCaller.SVIM_ASM)
    caller = next(x for x in SvimAsmWrapper(runner=FakeRunner(item)).plan_commands(item) if x.step == "svim_asm")
    assert caller.args[1] == "haploid"
    primary = HaplotypeAssemblyArtifact("S1", AssemblyRole.PRIMARY, assemblies[0].path, assemblies[0].source_gfa, "1", 1)
    with pytest.raises(InputValidationError): request(tmp_path, reference, (primary,), AssemblySvCaller.PAV)
    wrong = HaplotypeAssemblyArtifact("OTHER", assemblies[0].role, assemblies[0].path, assemblies[0].source_gfa, "1", 1)
    with pytest.raises(InputValidationError): request(tmp_path, reference, (wrong,), AssemblySvCaller.PAV)

def test_missing_and_overwrite_are_rejected(tmp_path, reference, assemblies):
    missing = HaplotypeAssemblyArtifact("S1", AssemblyRole.HAPLOTYPE1, tmp_path / "missing.fa", tmp_path / "missing.gfa", "1", 1)
    item = request(tmp_path, reference, (missing,), AssemblySvCaller.SVIM_ASM)
    with pytest.raises(OutputValidationError): SvimAsmWrapper(runner=FakeRunner(item)).plan_commands(item)
    output = tmp_path / ("result_" + chr(0x7ED3)) / "S1.pav.assembly.sv.vcf.gz"
    output.parent.mkdir(parents=True); output.write_bytes(b"old")
    with pytest.raises(OutputValidationError): request(tmp_path, reference, assemblies, AssemblySvCaller.PAV)
class FailedRunner(FakeRunner):
    def run(self,command,**kwargs):
        args=tuple(str(x) for x in command)
        if args and args[0]=="svim-asm":
            raise CommandExecutionError("svim-asm exit code 1")
        return super().run(command,**kwargs)

class MissingOutputRunner(FakeRunner):
    def run(self,command,**kwargs):
        args=tuple(str(x) for x in command)
        if args and args[0]=="snakemake" and "--version" not in args:
            return CommandResult(args,0,"","",.1,None,True)
        return super().run(command,**kwargs)

class InvalidVcfRunner(FakeRunner):
    def run(self,command,**kwargs):
        args=tuple(str(x) for x in command)
        if args and args[0]=="snakemake" and "--version" not in args:
            raw = self.request.work_directory / "S1.vcf.gz"
            write_bgzf(raw,b"not a VCF\n")
            write_bgzf(Path(f"{raw}.tbi"),b"TBI\x01")
            return CommandResult(args,0,"","",.1,None,True)
        return super().run(command,**kwargs)

class MissingSvlenRunner(FakeRunner):
    def run(self, command, **kwargs):
        args = tuple(str(x) for x in command)
        if args and args[0] == "snakemake" and "--version" not in args:
            raw = self.request.work_directory / "S1.vcf.gz"
            payload = vcf().replace(b"SVTYPE=INS;SVLEN=50", b"SVTYPE=INS")
            write_bgzf(raw, payload)
            write_bgzf(Path(f"{raw}.tbi"), b"TBI\x01")
            return CommandResult(args, 0, "", "", .1, None, True)
        return super().run(command, **kwargs)

class WrongPavVersionRunner(FakeRunner):
    def run(self, command, **kwargs):
        args = tuple(str(x) for x in command)
        if args and args[0] == "snakemake" and "--version" not in args:
            raw = self.request.work_directory / "S1.vcf.gz"
            write_bgzf(raw, vcf().replace(b"PAV 2.4.6.0", b"PAV 2.4.5"))
            write_bgzf(Path(f"{raw}.tbi"), b"TBI\x01")
            return CommandResult(args, 0, "", "", .1, None, True)
        return super().run(command, **kwargs)

def test_external_failure_is_propagated(tmp_path,reference,assemblies):
    item=request(tmp_path,reference,assemblies,AssemblySvCaller.SVIM_ASM)
    with pytest.raises(CommandExecutionError,match="exit code 1"):
        SvimAsmWrapper(runner=FailedRunner(item)).run(item)

@pytest.mark.parametrize("runner_type",(MissingOutputRunner,InvalidVcfRunner))
def test_pav_missing_or_invalid_output_is_rejected(tmp_path,reference,assemblies,runner_type):
    snakefile=tmp_path/"pav_site"/"Snakefile"
    snakefile.parent.mkdir(); snakefile.write_text("# PAV\n")
    item=request(tmp_path,reference,assemblies,AssemblySvCaller.PAV)
    with pytest.raises(OutputValidationError):
        PavWrapper(snakefile=snakefile,pav_version="2.4.6",runner=runner_type(item)).run(item)

def test_pav_sv_only_selection_rejects_missing_svlen(tmp_path, reference, assemblies):
    snakefile = tmp_path / "pav_site" / "Snakefile"
    snakefile.parent.mkdir(); snakefile.write_text("# PAV\n")
    item = request(tmp_path, reference, assemblies, AssemblySvCaller.PAV)
    with pytest.raises(OutputValidationError, match="lacks one scalar SVLEN"):
        PavWrapper(
            snakefile=snakefile,
            pav_version="2.4.6",
            runner=MissingSvlenRunner(item),
        ).run(item)
    assert (item.work_directory / "S1.vcf.gz").exists()
    assert not item.output_vcf.exists()

def test_pav_vcf_version_must_match_config(tmp_path, reference, assemblies):
    snakefile = tmp_path / "pav_site" / "Snakefile"
    snakefile.parent.mkdir(); snakefile.write_text("# PAV\n")
    item = request(tmp_path, reference, assemblies, AssemblySvCaller.PAV)
    with pytest.raises(OutputValidationError, match="differs from configured"):
        PavWrapper(
            snakefile=snakefile,
            pav_version="2.4.6",
            runner=WrongPavVersionRunner(item),
        ).run(item)

def test_pav_sv_selection_rejects_unvalidated_version(tmp_path, reference, assemblies):
    snakefile = tmp_path / "pav_site" / "Snakefile"
    snakefile.parent.mkdir(); snakefile.write_text("# PAV\n")
    item = request(tmp_path, reference, assemblies, AssemblySvCaller.PAV)
    with pytest.raises(ToolVersionError, match="validated only for PAV 2.4.6"):
        PavWrapper(
            snakefile=snakefile,
            pav_version="2.5.0",
            runner=FakeRunner(item),
        ).run(item)
