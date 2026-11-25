"""Minimal Truvari concordance wrapper; never a truth generator."""

from __future__ import annotations
import json
import gzip
import math
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
from hifivar.validation import compute_sha256, validate_output_file
from hifivar.benchmark import BenchmarkMetric, BenchmarkVariantClass

_VERSION = re.compile(r"(\d+(?:\.\d+)+)")

@dataclass(frozen=True, slots=True)
class TruvariThresholds:
    """Explicit Truvari bench policy; ``None`` preserves tool defaults.

    ``sizemax=-1`` explicitly disables the native upper size limit. This does
    not change the native default (50,000 bp) or select an unlimited policy.
    """
    refdist: int | None = None
    pctseq: float | None = None
    pctsize: float | None = None
    pctovl: float | None = None
    sizemin: int | None = None
    sizemax: int | None = None
    bnddist: int | None = None
    pass_only: bool = False
    pick: str | None = None
    chunksize: int | None = None
    refine: bool = False

    def __post_init__(self):
        for name in ("refdist", "sizemin", "bnddist"):
            value=getattr(self,name)
            if value is not None and (isinstance(value,bool) or not isinstance(value,int) or value < 0):
                raise InputValidationError(f"Truvari {name} must be a non-negative integer or null.")
        if self.sizemax is not None and (
            isinstance(self.sizemax, bool) or not isinstance(self.sizemax, int)
            or self.sizemax < -1
        ):
            raise InputValidationError("Truvari sizemax must be -1 (no upper limit), a non-negative integer or null.")
        for name in ("pctseq", "pctsize", "pctovl"):
            value=getattr(self,name)
            if value is not None and (isinstance(value,bool) or not isinstance(value,(int,float)) or not 0 <= value <= 1):
                raise InputValidationError(f"Truvari {name} must be between 0 and 1 or null.")
        for name in ("pass_only", "refine"):
            if not isinstance(getattr(self, name), bool):
                raise InputValidationError(f"Truvari {name} must be a boolean.")
        if self.pick is not None and (
            not isinstance(self.pick, str)
            or self.pick not in {"single", "multi", "ac"}
        ):
            raise InputValidationError("Truvari pick must be single, multi, ac or null.")
        if self.chunksize is not None and (
            isinstance(self.chunksize, bool)
            or not isinstance(self.chunksize, int)
            or self.chunksize <= 0
        ):
            raise InputValidationError("Truvari chunksize must be a positive integer or null.")
        if self.sizemin is not None and self.sizemax not in {None, -1} and self.sizemin > self.sizemax:
            raise InputValidationError("Truvari sizemin must not exceed sizemax.")
        if self.chunksize is not None and self.refdist is not None and self.chunksize < self.refdist:
            raise InputValidationError("Truvari chunksize must be at least refdist when both are explicit.")

    def to_dict(self):
        return {name:getattr(self,name) for name in ("refdist","pctseq","pctsize","pctovl","sizemin","sizemax","bnddist","pass_only","pick","chunksize","refine")}

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
                "thresholds": self.thresholds.to_dict(),
                "bench_summary_candidate": str(self.output_directory / "summary.json"),
                "refinement_output_candidates": {
                    name: str(path) for name, path in self.refinement_output_candidates.items()
                }}

    @property
    def refinement_output_candidates(self) -> dict[str, Path]:
        """Official Truvari 5.4 output names, not evidence of successful refinement."""
        if not self.thresholds.refine:
            return {}
        return {
            "variant_summary": self.output_directory / "refine.variant_summary.json",
            "region_summary": self.output_directory / "refine.region_summary.json",
            "base_vcf": self.output_directory / "refine.base.vcf.gz",
            "base_index": self.output_directory / "refine.base.vcf.gz.tbi",
            "comparison_vcf": self.output_directory / "refine.comp.vcf.gz",
            "comparison_index": self.output_directory / "refine.comp.vcf.gz.tbi",
        }

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
    refined_summary_path: Path | None = None
    refined_region_summary_path: Path | None = None
    refined_base_vcf_path: Path | None = None
    refined_comparison_vcf_path: Path | None = None
    refinement_status: str | None = None
    refinement_proof: dict | None = None
    effective_summary_path: Path | None = None
    effective_summary_sha256: str | None = None

    def to_dict(self):
        refinement_status = self.refinement_status or ("PLANNED" if self.status is TruvariResultStatus.PLANNED else "COMPLETE" if self.refined_summary_path else "NOT_REQUESTED")
        effective = self.effective_summary_path or self.refined_summary_path or self.summary_path
        return {"request": self.request.to_dict(), "status": self.status.value,
                "command": list(self.command), "display_command": format_command(self.command),
                "version": self.version, "runtime_seconds": self.runtime_seconds,
                "summary_path": str(self.summary_path) if self.summary_path else None,
                "summary_kind": "bench",
                "refinement_requested": self.request.thresholds.refine,
                "refinement_status": refinement_status,
                "refinement_proof": self.refinement_proof,
                "effective_summary_kind": ("refined" if self.refined_summary_path else "bench") if effective else None,
                "effective_summary_path": str(effective) if effective else None,
                "effective_summary_sha256": self.effective_summary_sha256,
                "refined_summary_path": str(self.refined_summary_path) if self.refined_summary_path else None,
                "refined_summary_kind": "variant" if self.refined_summary_path else None,
                "refined_region_summary_path": str(self.refined_region_summary_path) if self.refined_region_summary_path else None,
                "refined_base_vcf_path": str(self.refined_base_vcf_path) if self.refined_base_vcf_path else None,
                "refined_comparison_vcf_path": str(self.refined_comparison_vcf_path) if self.refined_comparison_vcf_path else None,
                "refined_stratification_status": ("NOT_NEEDED_NATIVE_ZERO_REGION_PERFECT_BENCH" if refinement_status == "NOT_NEEDED" else "NOT_RUN_NEEDS_REFINED_ASSIGNMENTS") if self.request.thresholds.refine else "NOT_REQUESTED",
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
        if request.thresholds.pick is not None:
            command.extend(("--pick", request.thresholds.pick))
        if request.thresholds.chunksize is not None:
            command.extend(("--chunksize", str(request.thresholds.chunksize)))
        if request.thresholds.refine:
            command.append("--refine")
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
        output_absent_before = not request.output_directory.exists()
        result=self.runner.run(command,stderr_path=stderr_path)
        summary=request.output_directory/"summary.json"
        validate_output_file(summary)
        refined = request.refinement_output_candidates
        refinement_status, proof = "NOT_REQUESTED", None
        if refined and all(path.is_file() for path in refined.values()):
            for path in refined.values():
                validate_output_file(path)
            parse_truvari_refined_summary(refined["variant_summary"])
            try:
                region_summary = json.loads(refined["region_summary"].read_text(encoding="utf-8"))
            except (OSError, UnicodeError, json.JSONDecodeError) as error:
                raise OutputValidationError(f"Unable to read Truvari region summary: {error}") from error
            if not isinstance(region_summary, dict) or not region_summary:
                raise OutputValidationError("Truvari refined region summary must be a non-empty JSON object.")
            refinement_status = "COMPLETE"
        elif refined:
            from hifivar.truvari_refinement import verify_noop_refinement
            if getattr(result, "args", None) != command:
                raise OutputValidationError("Native no-op command result is not bound to this invocation.")
            proof = verify_noop_refinement(
                request.to_dict(), command, version=version,
                returncode=getattr(result, "returncode", None), executed=getattr(result, "executed", None),
                stderr_path=getattr(result, "stderr_path", None) or stderr_path,
                output_absent_before=output_absent_before,
            )
            refinement_status, refined = "NOT_NEEDED", {}
        effective = refined.get("variant_summary", summary)
        return TruvariResult(
            request, TruvariResultStatus.COMPLETED, command, version, result.duration_seconds, summary,
            refined.get("variant_summary"), refined.get("region_summary"),
            refined.get("base_vcf"), refined.get("comparison_vcf"),
            refinement_status=refinement_status, refinement_proof=proof,
            effective_summary_path=effective, effective_summary_sha256=compute_sha256(effective),
        )

def parse_truvari_summary(path: Path, *, variant_class: BenchmarkVariantClass = BenchmarkVariantClass.SV):
    """Parse bench metrics, keeping distinct TP sides and native null rates.

    Native summaries carry both totals and receive denominator/rate consistency
    checks. Legacy abbreviated summaries without either total remain readable,
    including historically rounded rates; their rate formulas are not certified.
    A refined or regional summary cannot be substituted for the bench summary.
    """
    source = Path(path)
    if source.name in {"refine.variant_summary.json", "refine.region_summary.json"}:
        raise OutputValidationError("Bench metrics require the initial bench summary, not a refined summary.")
    return _parse_summary_metrics(source, variant_class=variant_class, require_totals=False)


def parse_truvari_refined_summary(
    path: Path, *, variant_class: BenchmarkVariantClass = BenchmarkVariantClass.SV,
) -> tuple[BenchmarkMetric, ...]:
    """Read native refined variant metrics, never the regional confusion matrix.

    Truvari 5.4 ``make_ga4gh`` preserves distinct truth/query TP counts and
    allows null rates when their denominators are empty. The initial bench
    summary is intentionally parsed by a different entry point.
    """
    source = Path(path)
    if source.name != "refine.variant_summary.json":
        raise OutputValidationError("Refined variant metrics require refine.variant_summary.json, not bench or region summary.")
    return _parse_summary_metrics(source, variant_class=variant_class, require_totals=True)


def _parse_summary_metrics(
    source: Path, *, variant_class: BenchmarkVariantClass, require_totals: bool,
) -> tuple[BenchmarkMetric, ...]:
    """Validate native metric semantics without replacing or inventing rates."""
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise OutputValidationError(f"Unable to read Truvari summary '{source}': {error}") from error
    if not isinstance(payload, dict):
        raise OutputValidationError("Truvari summary must be a JSON object.")
    totals = ("base cnt", "comp cnt")
    has_totals = all(name in payload for name in totals)
    if (require_totals or any(name in payload for name in totals)) and not has_totals:
        raise OutputValidationError("Native Truvari summary requires both truth/query totals.")
    count_fields = ("TP-base", "TP-comp", "FP", "FN") + (totals if has_totals else ())
    for field in count_fields:
        value = payload.get(field)
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise OutputValidationError(f"Truvari summary requires a non-negative integer '{field}'.")
    base_count = payload["TP-base"] + payload["FN"]
    comp_count = payload["TP-comp"] + payload["FP"]
    if has_totals and (payload["base cnt"] != base_count or payload["comp cnt"] != comp_count):
        raise OutputValidationError("Truvari summary truth/query denominators do not reconcile.")
    # Match Truvari 5.4 utils.performance_metrics: either empty side makes ALL
    # three rates null, and zero P+R has a null F1 even with nonempty sides.
    expected = {"precision": None, "recall": None, "f1": None}
    if base_count and comp_count:
        precision = payload["TP-comp"] / comp_count
        recall = payload["TP-base"] / base_count
        expected.update(precision=precision, recall=recall)
        if precision + recall:
            expected["f1"] = 2 * precision * recall / (precision + recall)
    for field in ("precision", "recall", "f1"):
        if field not in payload:
            raise OutputValidationError(f"Truvari summary lacks '{field}'.")
        value = payload[field]
        if value is not None and (
            isinstance(value, bool) or not isinstance(value, (int, float))
            or not math.isfinite(value) or not 0 <= value <= 1
        ):
            raise OutputValidationError(f"Truvari summary '{field}' must be a finite rate or null.")
        if (value is None) != (expected[field] is None):
            raise OutputValidationError(f"Truvari summary '{field}' null state conflicts with native denominator semantics.")
        if has_totals and value is not None and not math.isclose(value, expected[field], rel_tol=1e-6, abs_tol=1e-6):
            raise OutputValidationError(f"Truvari summary '{field}' does not reconcile with native truth/query counts.")
    fields = (("tp_base", "TP-base"), ("tp_call", "TP-comp"), ("fp", "FP"), ("fn", "FN"),
              ("precision", "precision"), ("recall", "recall"), ("f1", "f1"))
    return tuple(
        BenchmarkMetric(name, payload[field], variant_class, source_field=field)
        for name, field in fields
    )

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

__all__=["TruvariRequest","TruvariResult","TruvariResultStatus","TruvariThresholds","TruvariWrapper","parse_truvari_summary","parse_truvari_refined_summary","stratify_truvari_outputs"]
