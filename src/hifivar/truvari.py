"""Minimal Truvari concordance wrapper; never a truth generator."""

from __future__ import annotations
import json
import gzip
from collections import defaultdict
import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from hifivar.command import CommandRunner, format_command
from hifivar.exceptions import (
    CommandExecutionError,
    InputValidationError,
    OutputValidationError,
    ToolVersionError,
)
from hifivar.reference import ReferenceGenome
from hifivar.validation import validate_output_file
from hifivar.benchmark import BenchmarkMetric, BenchmarkVariantClass

_VERSION = re.compile(r"(\d+(?:\.\d+)+)")

@dataclass(frozen=True, slots=True)
class TruvariThresholds:
    """Explicit Truvari bench policy; ``None`` preserves tool defaults."""
    refdist: int | None = None
    pctseq: float | None = None
    pctsize: float | None = None
    pctovl: float | None = None
    sizemin: int | None = None
    sizemax: int | None = None
    bnddist: int | None = None
    pass_only: bool = False

    def __post_init__(self):
        for name in ("refdist", "sizemin", "sizemax", "bnddist"):
            value=getattr(self,name)
            if value is not None and (isinstance(value,bool) or not isinstance(value,int) or value < 0):
                raise InputValidationError(f"Truvari {name} must be a non-negative integer or null.")
        for name in ("pctseq", "pctsize", "pctovl"):
            value=getattr(self,name)
            if value is not None and (isinstance(value,bool) or not isinstance(value,(int,float)) or not 0 <= value <= 1):
                raise InputValidationError(f"Truvari {name} must be between 0 and 1 or null.")

    def to_dict(self):
        return {name:getattr(self,name) for name in ("refdist","pctseq","pctsize","pctovl","sizemin","sizemax","bnddist","pass_only")}

@dataclass(frozen=True, slots=True)
class TruvariRequest:
    sample_id: str
    reference: ReferenceGenome
    base_vcf: Path
    comparison_vcf: Path
    output_directory: Path
    overwrite: bool = False
    confident_regions: Path | None = None
    thresholds: TruvariThresholds = TruvariThresholds()

    def __post_init__(self):
        if not self.sample_id: raise InputValidationError("Truvari sample_id must be non-empty.")
        object.__setattr__(self, "base_vcf", Path(self.base_vcf))
        object.__setattr__(self, "comparison_vcf", Path(self.comparison_vcf))
        object.__setattr__(self, "output_directory", Path(self.output_directory))
        if self.confident_regions is not None:
            object.__setattr__(self, "confident_regions", Path(self.confident_regions))
        if self.output_directory.exists() and not self.overwrite:
            raise OutputValidationError(f"Truvari output directory exists: '{self.output_directory}'.")

    def to_dict(self):
        return {"sample_id": self.sample_id, "reference": self.reference.to_dict(),
                "base_vcf": str(self.base_vcf), "comparison_vcf": str(self.comparison_vcf),
                "output_directory": str(self.output_directory), "overwrite": self.overwrite,
                "confident_regions": str(self.confident_regions) if self.confident_regions else None,
                "thresholds": self.thresholds.to_dict()}

class TruvariResultStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"

@dataclass(frozen=True, slots=True)
class TruvariResult:
    request: TruvariRequest
    status: TruvariResultStatus
    command: tuple[str, ...]
    version: str | None = None
    runtime_seconds: float = 0.0
    summary_path: Path | None = None

    def to_dict(self):
        return {"request": self.request.to_dict(), "status": self.status.value,
                "command": list(self.command), "display_command": format_command(self.command),
                "version": self.version, "runtime_seconds": self.runtime_seconds,
                "summary_path": str(self.summary_path) if self.summary_path else None,
                "interpretation": "comparison_only_not_truth"}

class TruvariWrapper:
    def __init__(self, *, executable="truvari", runner: CommandRunner | None = None):
        self.executable, self.runner = executable, runner or CommandRunner()

    def plan_command(self, request: TruvariRequest):
        command = [self.executable, "bench", "-b", str(request.base_vcf.absolute()),
                "-c", str(request.comparison_vcf.absolute()), "-f",
                str(request.reference.fasta.absolute()), "-o",
                str(request.output_directory.absolute())]
        if request.confident_regions is not None:
            command.extend(("--includebed", str(request.confident_regions.absolute())))
        flags=(("refdist","--refdist"),("pctseq","--pctseq"),("pctsize","--pctsize"),
               ("pctovl","--pctovl"),("sizemin","--sizemin"),("sizemax","--sizemax"),
               ("bnddist","--bnddist"))
        for field,flag in flags:
            value=getattr(request.thresholds,field)
            if value is not None: command.extend((flag,str(value)))
        if request.thresholds.pass_only: command.append("--passonly")
        return tuple(command)

    def detect_version(self):
        self.runner.require_executable(self.executable)
        try:
            result = self.runner.run([self.executable, "version"])
        except CommandExecutionError:
            result = self.runner.run([self.executable, "--version"])
        match=_VERSION.search("\n".join(x for x in (result.stdout,result.stderr) if x))
        if match is None: raise ToolVersionError("Unable to parse Truvari version.")
        return match.group(1)

    def run(self, request: TruvariRequest, *, dry_run=False, stderr_path: Path | None = None):
        command=self.plan_command(request)
        if dry_run:
            self.runner.run(command,dry_run=True)
            return TruvariResult(request,TruvariResultStatus.PLANNED,command)
        for path in (request.reference.fasta,request.base_vcf,Path(f"{request.base_vcf}.tbi"),
                     request.comparison_vcf,Path(f"{request.comparison_vcf}.tbi")):
            validate_output_file(path)
        if request.confident_regions is not None: validate_output_file(request.confident_regions)
        version=self.detect_version()
        request.output_directory.parent.mkdir(parents=True,exist_ok=True)
        result=self.runner.run(command,stderr_path=stderr_path)
        summary=request.output_directory/"summary.json"
        validate_output_file(summary)
        return TruvariResult(request,TruvariResultStatus.COMPLETED,command,version,result.duration_seconds,summary)

def parse_truvari_summary(path: Path, *, variant_class: BenchmarkVariantClass = BenchmarkVariantClass.SV):
    """Parse official Truvari summary keys and expose TP-call from TP-comp."""
    try:
        payload=json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError,json.JSONDecodeError) as error:
        raise OutputValidationError(f"Unable to read Truvari summary '{path}': {error}") from error
    fields=(("tp_base","TP-base"),("tp_call","TP-comp"),("fp","FP"),("fn","FN"),
            ("precision","precision"),("recall","recall"),("f1","f1"))
    metrics=[]
    for name,source in fields:
        value=payload.get(source)
        if isinstance(value,bool) or not isinstance(value,(int,float)):
            raise OutputValidationError(f"Truvari summary lacks numeric field '{source}'.")
        metrics.append(BenchmarkMetric(name,value,variant_class,source_field=source))
    return tuple(metrics)

def stratify_truvari_outputs(output_directory: Path, *, size_bins: tuple[int, ...] = (),
                             variant_class: BenchmarkVariantClass = BenchmarkVariantClass.SV):
    """Derive descriptive SVTYPE/size-bin metrics from Truvari-assigned VCFs.

    This does not perform matching. BND and unresolved/complex records are
    counted by type but explicitly excluded from length bins.
    """
    bins=tuple(size_bins)
    if any(isinstance(value,bool) or not isinstance(value,int) or value <= 0 for value in bins) or tuple(sorted(set(bins))) != bins:
        raise InputValidationError("Truvari size bins must be unique ascending positive integers.")
    roles={"TP-base":"tp-base.vcf.gz","TP-comp":"tp-comp.vcf.gz","FP":"fp.vcf.gz","FN":"fn.vcf.gz"}
    counts: dict[str,dict[str,int]]=defaultdict(lambda:defaultdict(int))
    unsupported:set[str]=set()
    for role,name in roles.items():
        path=Path(output_directory)/name
        validate_output_file(path)
        opener=gzip.open if str(path).endswith(".gz") else open
        try:
            with opener(path,"rt",encoding="utf-8",newline="") as handle:
                for line in handle:
                    if line.startswith("#") or not line.strip(): continue
                    fields=line.rstrip().split("\t")
                    if len(fields) < 8: raise OutputValidationError(f"Malformed Truvari VCF record in '{path}'.")
                    info={part.partition("=")[0]:part.partition("=")[2] for part in fields[7].split(";")}
                    kind=(info.get("SVTYPE") or "UNRESOLVED").upper()
                    counts[f"SVTYPE:{kind}"][role]+=1
                    if kind in {"BND","TRA","CPX","CTX","UNRESOLVED"}:
                        unsupported.add(kind); continue
                    if not bins: continue
                    length=_sv_length(fields,info)
                    if length is None: unsupported.add(kind); continue
                    counts[f"SIZE:{_size_label(length,bins)}"][role]+=1
        except (OSError,UnicodeError) as error:
            raise OutputValidationError(f"Unable to stream Truvari VCF '{path}': {error}") from error
    metrics=[]
    for stratum,values in sorted(counts.items()):
        tp_base, tp_call, fp, fn=(values.get(key,0) for key in ("TP-base","TP-comp","FP","FN"))
        precision=tp_call/(tp_call+fp) if tp_call+fp else None
        recall=tp_base/(tp_base+fn) if tp_base+fn else None
        f1=2*precision*recall/(precision+recall) if precision is not None and recall is not None and precision+recall else None
        for name,value in (("tp_base",tp_base),("tp_call",tp_call),("fp",fp),("fn",fn),("precision",precision),("recall",recall),("f1",f1)):
            metrics.append(BenchmarkMetric(name,value,variant_class,stratum))
    return tuple(metrics), tuple(sorted(unsupported))

def _sv_length(fields,info):
    try:
        if info.get("SVLEN") not in {None,"","."}: return abs(int(str(info["SVLEN"]).split(",")[0]))
        if info.get("END") not in {None,"","."}: return abs(int(info["END"])-int(fields[1]))
        if not fields[4].startswith("<") and "[" not in fields[4] and "]" not in fields[4]: return abs(len(fields[4])-len(fields[3]))
    except ValueError: return None
    return None

def _size_label(length,bins):
    lower=0
    for upper in bins:
        if length < upper: return f"{lower}-{upper-1}bp"
        lower=upper
    return f">={lower}bp"

__all__=["TruvariRequest","TruvariResult","TruvariResultStatus","TruvariThresholds","TruvariWrapper","parse_truvari_summary","stratify_truvari_outputs"]
