import json
from pathlib import Path
import pytest
from hifivar.benchmark import BenchmarkVariantClass
from hifivar.reference import ReferenceGenome
from hifivar.truvari import TruvariRequest, TruvariThresholds, TruvariWrapper, parse_truvari_summary, stratify_truvari_outputs
import gzip
from hifivar.command import CommandResult


def test_truvari_benchmark_policy_and_official_summary_fields(tmp_path: Path) -> None:
    fasta=tmp_path/"ref.fa"; fasta.write_text(">chr1\nACGT\n"); Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n")
    ref=ReferenceGenome.from_fasta(fasta,build="GRCh38"); base=tmp_path/"truth.vcf.gz"; comp=tmp_path/"query.vcf.gz"; bed=tmp_path/"highconf.bed"
    policy=TruvariThresholds(refdist=500,pctseq=.7,pctsize=.7,sizemin=50,pass_only=True)
    request=TruvariRequest("S1",ref,base,comp,tmp_path/"out",confident_regions=bed,thresholds=policy)
    command=TruvariWrapper().plan_command(request)
    assert command[1]=="bench" and "--includebed" in command and "--refdist" in command and "--passonly" in command
    summary=tmp_path/"summary.json"; summary.write_text(json.dumps({"TP-base":8,"TP-comp":7,"FP":2,"FN":3,"precision":.777,"recall":.727,"f1":.75}))
    metrics={m.name:(m.value,m.source_field) for m in parse_truvari_summary(summary,variant_class=BenchmarkVariantClass.ASSEMBLY_SV)}
    assert metrics["tp_call"] == (7,"TP-comp") and metrics["tp_base"] == (8,"TP-base")


def test_truvari_thresholds_reject_invalid_fraction() -> None:
    with pytest.raises(Exception,match="between 0 and 1"): TruvariThresholds(pctseq=1.1)

def test_truvari_stratification_reuses_assignments_and_excludes_bnd_from_size(tmp_path: Path) -> None:
    root=tmp_path/"truvari"; root.mkdir()
    rows={"tp-base.vcf.gz":"chr1\t10\td1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-80\nchr1\t30\tb1\tA\tA]chr2:40]\t.\tPASS\tSVTYPE=BND\n",
          "tp-comp.vcf.gz":"chr1\t10\td1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-80\nchr1\t30\tb1\tA\tA]chr2:40]\t.\tPASS\tSVTYPE=BND\n",
          "fp.vcf.gz":"chr1\t50\ti1\tA\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=200\n",
          "fn.vcf.gz":""}
    for name,records in rows.items():
        with gzip.open(root/name,"wt",encoding="utf-8") as handle: handle.write("##fileformat=VCFv4.3\n"+records)
    metrics,unsupported=stratify_truvari_outputs(root,size_bins=(50,100,500))
    assert "BND" in unsupported
    strata={metric.stratum for metric in metrics}
    assert "SVTYPE:DEL" in strata and "SIZE:50-99bp" in strata and not any("BND" in value and value.startswith("SIZE") for value in strata)

class _Runner:
    def require_executable(self, executable): return Path(executable)
    def run(self, command, **kwargs):
        args=tuple(str(item) for item in command)
        if args[-1:] == ("version",): return CommandResult(args,0,"Truvari v5.4.0","",.1,None,True)
        out=Path(args[args.index("-o")+1]); out.mkdir(parents=True)
        (out/"summary.json").write_text(json.dumps({"TP-base":1,"TP-comp":1,"FP":0,"FN":0,"precision":1.0,"recall":1.0,"f1":1.0}))
        for name in ("tp-base.vcf.gz","tp-comp.vcf.gz","fp.vcf.gz","fn.vcf.gz"):
            with gzip.open(out/name,"wt",encoding="utf-8") as handle: handle.write("##fileformat=VCFv4.3\n" + ("chr1\t10\td1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-80\n" if name.startswith("tp-") else ""))
        return CommandResult(args,0,"","",.2,None,True)

def test_truvari_fake_end_to_end_uses_existing_wrapper(tmp_path: Path) -> None:
    fasta=tmp_path/"ref.fa"; fasta.write_text(">chr1\nACGT\n"); Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n")
    ref=ReferenceGenome.from_fasta(fasta,build="GRCh38"); truth=tmp_path/"truth.vcf.gz"; query=tmp_path/"query.vcf.gz"
    for path in (truth,query): path.write_bytes(b"vcf"); Path(f"{path}.tbi").write_bytes(b"index")
    item=TruvariRequest("S1",ref,truth,query,tmp_path/"out")
    result=TruvariWrapper(runner=_Runner()).run(item)
    assert result.version=="5.4.0" and {m.name:m.value for m in parse_truvari_summary(result.summary_path)}["f1"]==1.0
    stratified,unsupported=stratify_truvari_outputs(item.output_directory,size_bins=(50,100))
    assert any(m.stratum=="SIZE:50-99bp" for m in stratified) and unsupported==()
