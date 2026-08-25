"""Phase 9 per-sample SV harmonization and concordance provenance."""

from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Final
from hifivar import __version__
from hifivar.context import AnalysisContext
from hifivar.exceptions import ConfigurationError, InputValidationError
from hifivar.harmonization import SVEvidenceSourceArtifact, SVHarmonizationRequest
from hifivar.jasmine import JasmineResult, JasmineWrapper
from hifivar.serialization import standardize_data, utc_now_iso8601, write_json_atomic, write_yaml_atomic
from hifivar.truvari import TruvariRequest, TruvariResult, TruvariWrapper

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
                int(section["max_dist"]),str(section["distance_type"]),bool(section["overwrite"]))
        except (KeyError,TypeError,ValueError) as error:
            raise ConfigurationError(f"Invalid Phase 9 settings: {error}") from error

    def to_dict(self):
        return {"backend":"jasmine","jasmine_executable":self.jasmine_executable,
                "truvari_executable":self.truvari_executable,
                "bgzip_executable":self.bgzip_executable,"tabix_executable":self.tabix_executable,
                "threads":self.threads,"memory_mb":self.memory_mb,
                "runtime_minutes":self.runtime_minutes,"max_dist":self.max_dist,
                "distance_type":self.distance_type,"overwrite":self.overwrite}

@dataclass(frozen=True, slots=True)
class Phase9SampleResult:
    sample_id: str
    sources: tuple[SVEvidenceSourceArtifact, ...]
    jasmine: JasmineResult
    truvari: tuple[TruvariResult, ...]

    def to_dict(self):
        return {"sample_id":self.sample_id,"sources":[x.to_dict() for x in self.sources],
                "jasmine":self.jasmine.to_dict(),"truvari":[x.to_dict() for x in self.truvari]}

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
    results=[]
    for record in context.samples:
        sample=record.sample.sample_id
        sources=evidence.get(sample)
        if not sources: raise InputValidationError(f"Phase 9 evidence is missing for sample '{sample}'.")
        request=SVHarmonizationRequest(sample,context.reference,sources,work_root/"jasmine"/sample,
            output_root/sample/f"{sample}.harmonized.sv.vcf.gz",
            output_root/sample/f"{sample}.sv.evidence.tsv",
            settings.max_dist,settings.distance_type,settings.overwrite)
        jasmine_result=jasmine.run(request,dry_run=dry_run,
            stderr_path=work_root/"jasmine"/sample/f"{sample}.jasmine.log")
        comparisons=[]
        for source in request.runnable_sources:
            if source.vcf_path is None: continue
            tr_request=TruvariRequest(sample,context.reference,request.output_vcf,source.vcf_path,
                output_root/sample/"truvari"/source.caller,settings.overwrite)
            comparisons.append(truvari.run(tr_request,dry_run=dry_run,
                stderr_path=work_root/"truvari"/sample/f"{source.caller}.log"))
        results.append(Phase9SampleResult(sample,sources,jasmine_result,tuple(comparisons)))
    return Phase9RunReport(context,settings,tuple(results),dry_run,utc_now_iso8601())


__all__=["PHASE9_REPORT_SCHEMA_VERSION","Phase9RunReport","Phase9SampleResult","Phase9Settings","run_phase9"]
