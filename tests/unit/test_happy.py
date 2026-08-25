from __future__ import annotations
import gzip
import json
from pathlib import Path
import pytest
from hifivar.benchmark import BenchmarkRegion, TruthSet
from hifivar.command import CommandResult
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar.happy import (
    HappyRequest,
    HappyResultStatus,
    HappyWrapper,
    discover_happy_metrics,
    parse_happy_metrics_json,
    parse_happy_summary,
)
from hifivar.reference import ReferenceGenome


def reference(tmp_path: Path) -> ReferenceGenome:
    fasta=tmp_path/"参考.fa"; fasta.write_text(">chr1\nACGT\n",encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n",encoding="utf-8")
    return ReferenceGenome.from_fasta(fasta,build="GRCh38")


def request(tmp_path: Path) -> HappyRequest:
    ref=reference(tmp_path); query=tmp_path/"query.vcf.gz"; expected=tmp_path/"truth.vcf.gz"; bed=tmp_path/"confident.bed"
    for path in (query,expected): path.write_bytes(b"vcf"); Path(f"{path}.tbi").write_bytes(b"tbi")
    bed.write_text("chr1\t0\t4\n",encoding="utf-8")
    return HappyRequest("B1","S1",ref,query,TruthSet(expected,"GIAB-v4.2.1","GRCh38","GIAB"),
        BenchmarkRegion("confident",bed,"GIAB-v4.2.1"),tmp_path/"输出"/"happy",threads=2)


class Runner:
    def __init__(self,item, *, version_output="hap.py 0.3.15", metrics_compression="plain"):
        self.item=item; self.calls=[]; self.version_output=version_output; self.metrics_compression=metrics_compression
    def require_executable(self,executable): return Path(executable)
    def run(self,command,**kwargs):
        args=tuple(str(x) for x in command); self.calls.append((args,kwargs))
        if "--version" in args: return CommandResult(args,0,self.version_output,"",.1,None,True)
        if not kwargs.get("dry_run"):
            self.item.summary_csv.write_text("Type,Filter,METRIC.Recall,METRIC.Precision\nSNP,PASS,0.8,1.0\nINDEL,PASS,0.5,0.5\n",encoding="utf-8")
            if self.metrics_compression == "gzip":
                with gzip.open(self.item.metrics_json_gz, "wt", encoding="utf-8") as handle: json.dump({"format":"gzip"}, handle)
            elif self.metrics_compression == "plain":
                self.item.metrics_json.write_text("{}\n",encoding="utf-8")
            elif self.metrics_compression == "both":
                self.item.metrics_json.write_text("{}\n",encoding="utf-8")
                with gzip.open(self.item.metrics_json_gz, "wt", encoding="utf-8") as handle: json.dump({}, handle)
        return CommandResult(args,0,"","",.2,None,not kwargs.get("dry_run",False))


def test_command_is_deterministic_and_uses_command_runner(tmp_path: Path) -> None:
    item=request(tmp_path); runner=Runner(item); wrapper=HappyWrapper(runner=runner)
    command=wrapper.plan_command(item)
    assert command[:3] == ("hap.py",str(item.truth_set.path.absolute()),str(item.query_vcf.absolute()))
    assert "--threads" in command and "--engine=xcmp" in command and command[command.index("-f")+1].endswith("confident.bed")
    result=wrapper.run(item)
    assert result.status is HappyResultStatus.COMPLETED and result.version=="0.3.15"
    assert result.version_source == "command" and result.metrics_artifact == item.metrics_json
    assert result.to_dict()["version_source"] == "command"
    metrics={(m.stratum,m.name):m.value for m in result.metrics}
    assert metrics[("SNP","f1")] == pytest.approx(8/9)


def test_dry_run_creates_no_output_or_stratification(tmp_path: Path) -> None:
    item=request(tmp_path); result=HappyWrapper(runner=Runner(item)).run(item,dry_run=True)
    assert result.status is HappyResultStatus.PLANNED and not item.output_prefix.parent.exists()


def test_parser_uses_columns_and_named_rows_not_positions(tmp_path: Path) -> None:
    summary=tmp_path/"summary.csv"
    summary.write_text("Filter,METRIC.F1_Score,Type,METRIC.Precision,METRIC.Recall\nPASS,0.2,INDEL,0.3,0.4\nPASS,0.9,SNP,0.8,0.7\n",encoding="utf-8")
    parsed={(m.stratum,m.name):m.value for m in parse_happy_summary(summary)}
    assert parsed[("SNP","recall")] == 0.7 and parsed[("INDEL","f1")] == 0.2


def test_parser_rejects_missing_metric_column(tmp_path: Path) -> None:
    summary=tmp_path/"summary.csv"; summary.write_text("Type,Filter,METRIC.Recall\nSNP,PASS,1\nINDEL,PASS,1\n",encoding="utf-8")
    with pytest.raises(OutputValidationError,match="Precision"): parse_happy_summary(summary)

def test_missing_truth_and_confident_bed_fail_before_external_execution(tmp_path: Path) -> None:
    item=request(tmp_path); item.truth_set.path.unlink()
    with pytest.raises(InputValidationError): HappyWrapper(runner=Runner(item)).run(item)
    item.truth_set.path.write_bytes(b"vcf"); item.confident_regions.path.unlink()
    with pytest.raises(InputValidationError): HappyWrapper(runner=Runner(item)).run(item)

def test_truth_reference_build_mismatch_is_rejected(tmp_path: Path) -> None:
    item=request(tmp_path)
    with pytest.raises(InputValidationError,match="build mismatch"):
        HappyRequest("B","S1",item.reference,item.query_vcf,
            TruthSet(item.truth_set.path,"v1","GRCh37","truth"),item.confident_regions,tmp_path/"other")


def test_version_detection_rejects_empty_output_without_trustworthy_fallback(tmp_path: Path) -> None:
    item=request(tmp_path)
    wrapper=HappyWrapper(runner=Runner(item,version_output="Hap.py "))
    with pytest.raises(ToolVersionError,match="happy_version"):
        wrapper.detect_version()


@pytest.mark.parametrize("empty_version_output", ["", "Hap.py "])
def test_explicit_version_fallback_is_recorded_in_provenance(
    tmp_path: Path, empty_version_output: str
) -> None:
    item=request(tmp_path)
    wrapper=HappyWrapper(
        configured_version="0.3.15",
        runner=Runner(item,version_output=empty_version_output,metrics_compression="gzip"),
    )
    result=wrapper.run(item)
    assert result.version == "0.3.15" and result.version_source == "config"
    assert result.metrics_artifact == item.metrics_json_gz
    assert result.to_dict()["version_source"] == "config"


def test_configured_version_must_be_an_explicit_release() -> None:
    with pytest.raises(ToolVersionError,match="explicit numeric release"):
        HappyWrapper(configured_version="unknown")


def test_metrics_discovery_accepts_each_format_and_rejects_ambiguity_or_missing(tmp_path: Path) -> None:
    prefix=tmp_path/"happy"
    plain=Path(f"{prefix}.metrics.json"); compressed=Path(f"{prefix}.metrics.json.gz")
    with pytest.raises(OutputValidationError,match="no metrics artifact"):
        discover_happy_metrics(prefix)
    plain.write_text("{}",encoding="utf-8")
    assert discover_happy_metrics(prefix) == plain
    plain.unlink()
    with gzip.open(compressed,"wt",encoding="utf-8") as handle: json.dump({"ok":True},handle)
    assert discover_happy_metrics(prefix) == compressed
    plain.write_text("{}",encoding="utf-8")
    with pytest.raises(OutputValidationError,match="ambiguous"):
        discover_happy_metrics(prefix)


def test_metrics_parser_transparently_reads_plain_and_gzip_json(tmp_path: Path) -> None:
    plain=tmp_path/"metrics.json"; compressed=tmp_path/"metrics.json.gz"
    plain.write_text('{"kind":"plain"}',encoding="utf-8")
    with gzip.open(compressed,"wt",encoding="utf-8") as handle: json.dump({"kind":"gzip"},handle)
    assert parse_happy_metrics_json(plain) == {"kind":"plain"}
    assert parse_happy_metrics_json(compressed) == {"kind":"gzip"}
