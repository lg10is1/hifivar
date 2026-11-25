"""Phase 9 per-sample SV harmonization and concordance provenance."""

from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Final
from hifivar import __version__
from hifivar.bnd_representation import BndRepresentationRequest, BndRepresentationResult, BndRepresentationWrapper
from hifivar.breakends import BndPolicy
from hifivar.context import AnalysisContext
from hifivar.exceptions import ConfigurationError, InputValidationError, OutputValidationError
from hifivar.harmonization import EvidenceRunStatus, SVEvidenceSourceArtifact, SVHarmonizationRequest
from hifivar.jasmine import JasmineResult, JasmineWrapper
from hifivar.sample import validate_sample_id
from hifivar.serialization import standardize_data, utc_now_iso8601, write_json_atomic, write_yaml_atomic
from hifivar.truvari import TruvariRequest, TruvariResult, TruvariWrapper
from hifivar.sv_input_identity import SourceIdentityPolicy
from hifivar.sv_support_selection import (
    SupportSelectionRequest, SupportSelectionResult, SupportSelectionWrapper,
    support_selection_paths,
)

PHASE9_REPORT_SCHEMA_VERSION: Final[str] = "1.0"

@dataclass(frozen=True, slots=True)
class Phase9Settings:
    jasmine_executable: str
    truvari_executable: str
    bgzip_executable: str
    tabix_executable: str
    threads: int
    memory_mb: int
    runtime_minutes: int
    max_dist: int
    distance_type: str
    overwrite: bool
    bnd_policy: BndPolicy | None = None
    input_identity_policy: SourceIdentityPolicy | None = None
    minimum_caller_support: int = 1

    def __post_init__(self):
        if type(self.minimum_caller_support) is not int or not 1 <= self.minimum_caller_support <= 6:
            raise ConfigurationError("sv.harmonization.minimum_caller_support must be an integer between 1 and 6.")

    @classmethod
    def from_config(cls, config: Mapping[str, object]):
        sv=config.get("sv")
        if not isinstance(sv,Mapping): raise ConfigurationError("Phase 9 requires sv mapping.")
        section=sv.get("harmonization")
        if not isinstance(section,Mapping): raise ConfigurationError("Phase 9 requires sv.harmonization mapping.")
        try:
            if section["backend"] != "jasmine": raise ConfigurationError("Primary harmonizer must be Jasmine.")
            return cls(str(section["jasmine_executable"]),str(section["truvari_executable"]),
                str(section["bgzip_executable"]),str(section["tabix_executable"]),
                int(section["threads"]),int(section["memory_mb"]),int(section["runtime_minutes"]),
                int(section["max_dist"]),str(section["distance_type"]),bool(section["overwrite"]),
                BndPolicy.from_config(section.get("bnd_representation", {"enabled": False})),
                SourceIdentityPolicy.from_config(section.get("input_identity", {"enabled": False})),
                section.get("minimum_caller_support", 1))
        except (KeyError,TypeError,ValueError) as error:
            raise ConfigurationError(f"Invalid Phase 9 settings: {error}") from error

    def to_dict(self):
        return {"backend":"jasmine","jasmine_executable":self.jasmine_executable,
                "truvari_executable":self.truvari_executable,
                "bgzip_executable":self.bgzip_executable,"tabix_executable":self.tabix_executable,
                "threads":self.threads,"memory_mb":self.memory_mb,
                "runtime_minutes":self.runtime_minutes,"max_dist":self.max_dist,
                "distance_type":self.distance_type,"overwrite":self.overwrite,
                "minimum_caller_support":self.minimum_caller_support,
                "bnd_representation":self.bnd_policy.to_dict() if self.bnd_policy else {"enabled": False},
                "input_identity":self.input_identity_policy.to_dict() if self.input_identity_policy else {"enabled": False}}


def phase9_identity_paths(output_directory: str | Path, sample_id: str,
                          callers: tuple[str, ...]) -> dict[str, Path]:
    """Deterministic retained identity artifacts, separate from canonical VCFs."""
    validate_sample_id(sample_id)
    import re
    if not callers or len(set(callers)) != len(callers) or any(re.fullmatch(r"[A-Za-z0-9_]+", c) is None for c in callers):
        raise InputValidationError("Source identity paths require unique safe caller names.")
    directory = Path(output_directory).absolute() / sample_id / "source-identity-v1"
    paths = {"directory": directory, "manifest": directory / "manifest.json", "database": directory / "records.sqlite3"}
    for i, caller in enumerate(callers):
        stem = f"{i:02d}.{caller}.analysis"
        paths.update({f"{caller}_vcf": directory / f"{stem}.vcf.gz",
                      f"{caller}_index": directory / f"{stem}.vcf.gz.tbi",
                      f"{caller}_ledger": directory / f"{stem}.records.jsonl"})
    return paths


def phase9_bnd_paths(output_directory: str | Path, sample_id: str) -> dict[str, Path]:
    """Named files for the optional, separate mandatory native-plus-ledger bundle."""
    validate_sample_id(sample_id)
    directory = Path(output_directory).absolute() / sample_id / "bndrepr-v1"
    return {
        "directory": directory, "manifest": directory / "manifest.json",
        "vcf": directory / f"{sample_id}.bndrepr-v1.sv.vcf.gz",
        "index": directory / f"{sample_id}.bndrepr-v1.sv.vcf.gz.tbi",
        "ledger": directory / "records.jsonl", "junctions": directory / "junctions.tsv",
        "native_result": directory / "native_result.json",
        "native_vcf": directory / "native" / f"{sample_id}.harmonized.sv.vcf.gz",
        "native_index": directory / "native" / f"{sample_id}.harmonized.sv.vcf.gz.tbi",
        "native_raw": directory / "native" / f"{sample_id}.jasmine.raw.vcf",
        "native_evidence": directory / "native" / f"{sample_id}.sv.evidence.tsv",
    }

@dataclass(frozen=True, slots=True)
class Phase9SampleResult:
    sample_id: str
    sources: tuple[SVEvidenceSourceArtifact, ...]
    jasmine: JasmineResult
    truvari: tuple[TruvariResult, ...]
    bnd_representation: BndRepresentationResult | None = None
    bnd_plan: dict[str, object] | None = None
    support_selection: SupportSelectionResult | dict[str, object] | None = None

    def to_dict(self):
        return {"sample_id":self.sample_id,"sources":[x.to_dict() for x in self.sources],
                "jasmine":self.jasmine.to_dict(),"truvari":[x.to_dict() for x in self.truvari],
                "bnd_representation":self.bnd_representation.to_dict() if self.bnd_representation else self.bnd_plan,
                "support_selection":self.support_selection.to_dict() if isinstance(self.support_selection, SupportSelectionResult) else self.support_selection}

@dataclass(frozen=True, slots=True)
class Phase9RunReport:
    context: AnalysisContext
    settings: Phase9Settings
    sample_results: tuple[Phase9SampleResult, ...]
    dry_run: bool
    created_at: str
    hifivar_version: str = __version__
    schema_version: str = PHASE9_REPORT_SCHEMA_VERSION

    def to_dict(self):
        payload={"schema_version":self.schema_version,"hifivar_version":self.hifivar_version,
                 "created_at":self.created_at,"status":"planned" if self.dry_run else "completed",
                 "dry_run":self.dry_run,"context":self.context.to_dict(),
                 "settings":self.settings.to_dict(),
                 "sample_results":[x.to_dict() for x in self.sample_results],
                 "scientific_policy":{"harmonization_is_truth":False,
                    "caller_count_is_confidence":False,
                    "read_and_assembly_is_evidence_class_only":True}}
        value=standardize_data(payload,context="Phase 9 report value")
        if not isinstance(value,dict): raise InputValidationError("Phase 9 serialization failed.")
        return value

    def write_json(self,path:Path,*,overwrite=False):
        return write_json_atomic(self.to_dict(),path,overwrite=overwrite,artifact_name="Phase 9 report")

    def write_yaml(self,path:Path,*,overwrite=False):
        return write_yaml_atomic(self.to_dict(),path,overwrite=overwrite,artifact_name="Phase 9 report")


def run_phase9(
    context: AnalysisContext, evidence: Mapping[str, tuple[SVEvidenceSourceArtifact, ...]], *,
    output_directory: str | Path, work_directory: str | Path,
    config: Mapping[str, object], dry_run=False,
    jasmine_wrapper: JasmineWrapper | None = None,
    truvari_wrapper: TruvariWrapper | None = None,
    bnd_wrapper: BndRepresentationWrapper | None = None,
    support_wrapper: SupportSelectionWrapper | None = None,
) -> Phase9RunReport:
    sv=config.get("sv")
    if not isinstance(sv,Mapping) or not isinstance(sv.get("harmonization"),Mapping) or sv["harmonization"].get("enabled") is not True:
        raise ConfigurationError("Phase 9 requires sv.harmonization.enabled: true.")
    settings=Phase9Settings.from_config(config)
    jasmine=jasmine_wrapper or JasmineWrapper(executable=settings.jasmine_executable,
        bgzip_executable=settings.bgzip_executable,tabix_executable=settings.tabix_executable)
    truvari=truvari_wrapper or TruvariWrapper(executable=settings.truvari_executable)
    output_root,work_root=Path(output_directory),Path(work_directory)
    unknown=sorted(set(evidence).difference(context.sample_ids))
    if unknown: raise InputValidationError(f"Phase 9 received evidence for unknown samples: {unknown!r}.")
    if settings.bnd_policy is not None:
        # Fail before calling Jasmine for any sample if a protected bundle would
        # be reused. Native overwrite never authorizes bundle replacement.
        for sample in context.sample_ids:
            directory = phase9_bnd_paths(output_root, sample)["directory"]
            if directory.exists() or directory.is_symlink():
                raise OutputValidationError("BND representation requires a fresh per-sample bundle; overwrite is unsupported.")
    if settings.input_identity_policy is not None:
        for sample in context.sample_ids:
            directory = output_root.absolute() / sample / "source-identity-v1"
            if directory.exists() or directory.is_symlink():
                raise OutputValidationError("SV source identity requires a fresh per-sample bundle; overwrite is unsupported.")
    # Validate every sample before any producer or optional bundle is executed.
    for sample in context.sample_ids:
        available = evidence.get(sample)
        if not available:
            raise InputValidationError(f"Phase 9 evidence is missing for sample '{sample}'.")
        runnable = tuple(item for item in available if item.status in {EvidenceRunStatus.COMPLETED, EvidenceRunStatus.NO_CALLS})
        if settings.minimum_caller_support > len(runnable):
            raise ConfigurationError(f"minimum_caller_support exceeds the {len(runnable)} actual input callers for sample '{sample}'.")
        if settings.minimum_caller_support > 1:
            if len({item.caller for item in runnable}) != len(runnable):
                raise InputValidationError("Caller-support selection requires unique actual caller sources.")
            directory = support_selection_paths(output_root, sample, settings.minimum_caller_support)["directory"]
            if directory.exists() or directory.is_symlink():
                raise OutputValidationError("Caller-support selection requires a fresh per-sample bundle; overwrite is unsupported.")
    results=[]
    for record in context.samples:
        sample=record.sample.sample_id
        sources=evidence.get(sample)
        if not sources: raise InputValidationError(f"Phase 9 evidence is missing for sample '{sample}'.")
        request=SVHarmonizationRequest(sample,context.reference,sources,work_root/"jasmine"/sample,
            output_root/sample/f"{sample}.harmonized.sv.vcf.gz",
            output_root/sample/f"{sample}.sv.evidence.tsv",
            settings.max_dist,settings.distance_type,settings.overwrite,settings.input_identity_policy)
        jasmine_result=jasmine.run(request,dry_run=dry_run,
            stderr_path=work_root/"jasmine"/sample/f"{sample}.jasmine.log")
        support_selection = None
        if settings.minimum_caller_support > 1:
            support_paths = support_selection_paths(output_root, sample, settings.minimum_caller_support)
            if dry_run:
                support_selection = {"status": "planned", "executed": False,
                    "minimum_caller_support": settings.minimum_caller_support, "operator": ">=",
                    "source_order": [item.caller for item in request.runnable_sources],
                    "output_paths": {key: str(path) for key, path in support_paths.items()},
                    "native_output_role_unchanged": True, "caller_count_is_confidence": False,
                    "bnd_mate_completeness": "NOT_ASSESSED"}
            else:
                if jasmine_result.artifact is None:
                    raise OutputValidationError("Support selection requires a completed native Jasmine artifact.")
                selector = support_wrapper or SupportSelectionWrapper(
                    bgzip_executable=settings.bgzip_executable, tabix_executable=settings.tabix_executable)
                support_selection = selector.run(SupportSelectionRequest(
                    jasmine_result.artifact, context.reference, settings.minimum_caller_support,
                    support_paths["directory"]), stderr_path=work_root/"jasmine"/sample/f"{sample}.support-selection.log")
        bnd_result, bnd_plan = None, None
        if settings.bnd_policy is not None:
            paths = phase9_bnd_paths(output_root, sample)
            if dry_run:
                bnd_plan = {"status": "planned", "executed": False,
                    "policy": settings.bnd_policy.to_dict(),
                    "output_paths": {key: str(path) for key, path in paths.items()},
                    "native_output_role_unchanged": True, "consumer_opt_in_required": True}
            else:
                wrapper = bnd_wrapper or BndRepresentationWrapper(
                    bgzip_executable=settings.bgzip_executable, tabix_executable=settings.tabix_executable)
                bnd_result = wrapper.run(BndRepresentationRequest(jasmine_result, paths["directory"], settings.bnd_policy))
        comparisons=[]
        for source in request.runnable_sources:
            if source.vcf_path is None: continue
            tr_request=TruvariRequest(sample,context.reference,request.output_vcf,source.vcf_path,
                output_root/sample/"truvari"/source.caller,settings.overwrite)
            comparisons.append(truvari.run(tr_request,dry_run=dry_run,
                stderr_path=work_root/"truvari"/sample/f"{source.caller}.log"))
        if jasmine_result.source_identity is not None:
            jasmine_result.source_identity.validate()
        if isinstance(support_selection, SupportSelectionResult):
            support_selection.validate()
        results.append(Phase9SampleResult(sample,sources,jasmine_result,tuple(comparisons),bnd_result,bnd_plan,support_selection))
    return Phase9RunReport(context,settings,tuple(results),dry_run,utc_now_iso8601())


__all__=["PHASE9_REPORT_SCHEMA_VERSION","Phase9RunReport","Phase9SampleResult","Phase9Settings","phase9_bnd_paths","phase9_identity_paths","run_phase9"]
