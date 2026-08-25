"""Phase 4 read-based structural-variant orchestration and provenance."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Final

from hifivar import __version__
from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import AlignmentArtifact, AlignmentSortOrder, AlignmentSource, find_alignment_index, validate_alignment_artifact
from hifivar.context import AnalysisContext
from hifivar.cutesv import CuteSvRequest, CuteSvResources, CuteSvResult, CuteSvWrapper
from hifivar.exceptions import ConfigurationError, InputValidationError, ReferenceError
from hifivar.pbsv import PbsvRequest, PbsvResources, PbsvResult, PbsvWrapper
from hifivar.sample import InputType
from hifivar.sawfish import SawfishRequest, SawfishResources, SawfishResult, SawfishWrapper
from hifivar.serialization import standardize_data, utc_now_iso8601, write_json_atomic, write_yaml_atomic
from hifivar.sniffles2 import Sniffles2Request, Sniffles2Resources, Sniffles2Result, Sniffles2Wrapper
from hifivar.sv import BgzipTabixWrapper, StructuralVariantArtifact, SvCaller, VcfFinalizeRequest, VcfFinalizeResult, create_structural_variant_artifact


PHASE4_REPORT_SCHEMA_VERSION: Final[str] = "1.0"


class Phase4RunStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class Phase4Settings:
    enabled: bool
    overwrite: bool
    sawfish_enabled: bool
    sawfish_executable: str
    sawfish_resources: SawfishResources
    sawfish_disable_cnv: bool
    sniffles2_enabled: bool
    sniffles2_executable: str
    sniffles2_resources: Sniffles2Resources
    sniffles2_minimum_support: int | None
    sniffles2_minimum_sv_length: int
    pbsv_enabled: bool
    pbsv_executable: str
    pbsv_resources: PbsvResources
    cutesv_enabled: bool
    cutesv_executable: str
    cutesv_resources: CuteSvResources
    cutesv_minimum_support: int
    cutesv_minimum_sv_size: int
    cutesv_max_cluster_bias_ins: int
    cutesv_diff_ratio_merging_ins: float
    cutesv_max_cluster_bias_del: int
    cutesv_diff_ratio_merging_del: float
    cutesv_genotype: bool
    bgzip_executable: str
    tabix_executable: str

    @classmethod
    def from_config(cls, config: Mapping[str, object]) -> "Phase4Settings":
        sv = _mapping(config, "sv")
        sawfish = _mapping(sv, "sawfish")
        sniffles2 = _mapping(sv, "sniffles2")
        pbsv = _mapping(sv, "pbsv")
        cutesv = _mapping(sv, "cutesv")
        finalization = _mapping(sv, "finalization")
        try:
            return cls(
                enabled=bool(sv["enabled"]),
                overwrite=bool(sv["overwrite"]),
                sawfish_enabled=bool(sawfish["enabled"]),
                sawfish_executable=str(sawfish["executable"]),
                sawfish_resources=SawfishResources(int(sawfish["threads"]), int(sawfish["memory_mb"]), int(sawfish["runtime_minutes"])),
                sawfish_disable_cnv=bool(sawfish["disable_cnv"]),
                sniffles2_enabled=bool(sniffles2["enabled"]),
                sniffles2_executable=str(sniffles2["executable"]),
                sniffles2_resources=Sniffles2Resources(int(sniffles2["threads"]), int(sniffles2["memory_mb"]), int(sniffles2["runtime_minutes"])),
                sniffles2_minimum_support=(None if sniffles2["minimum_support"] is None else int(sniffles2["minimum_support"])),
                sniffles2_minimum_sv_length=int(sniffles2["minimum_sv_length"]),
                pbsv_enabled=bool(pbsv["enabled"]),
                pbsv_executable=str(pbsv["executable"]),
                pbsv_resources=PbsvResources(int(pbsv["threads"]), int(pbsv["memory_mb"]), int(pbsv["runtime_minutes"])),
                cutesv_enabled=bool(cutesv["enabled"]),
                cutesv_executable=str(cutesv["executable"]),
                cutesv_resources=CuteSvResources(int(cutesv["threads"]), int(cutesv["memory_mb"]), int(cutesv["runtime_minutes"])),
                cutesv_minimum_support=int(cutesv["minimum_support"]),
                cutesv_minimum_sv_size=int(cutesv["minimum_sv_size"]),
                cutesv_max_cluster_bias_ins=int(cutesv["max_cluster_bias_ins"]),
                cutesv_diff_ratio_merging_ins=float(cutesv["diff_ratio_merging_ins"]),
                cutesv_max_cluster_bias_del=int(cutesv["max_cluster_bias_del"]),
                cutesv_diff_ratio_merging_del=float(cutesv["diff_ratio_merging_del"]),
                cutesv_genotype=bool(cutesv["genotype"]),
                bgzip_executable=str(finalization["bgzip_executable"]),
                tabix_executable=str(finalization["tabix_executable"]),
            )
        except (KeyError, TypeError, ValueError) as error:
            raise ConfigurationError(f"Invalid or incomplete Phase 4 sv configuration: {error}") from error

    @property
    def enabled_callers(self) -> tuple[SvCaller, ...]:
        pairs = (
            (SvCaller.SAWFISH, self.sawfish_enabled),
            (SvCaller.SNIFFLES2, self.sniffles2_enabled),
            (SvCaller.PBSV, self.pbsv_enabled),
            (SvCaller.CUTESV, self.cutesv_enabled),
        )
        return tuple(caller for caller, enabled in pairs if enabled)

    def to_dict(self) -> dict[str, object]:
        return {
            "enabled": self.enabled,
            "overwrite": self.overwrite,
            "enabled_callers": [caller.value for caller in self.enabled_callers],
            "sawfish": {"executable": self.sawfish_executable, "resources": self.sawfish_resources.to_dict(), "disable_cnv": self.sawfish_disable_cnv},
            "sniffles2": {"executable": self.sniffles2_executable, "resources": self.sniffles2_resources.to_dict(), "minimum_support": self.sniffles2_minimum_support, "minimum_sv_length": self.sniffles2_minimum_sv_length},
            "pbsv": {"executable": self.pbsv_executable, "resources": self.pbsv_resources.to_dict()},
            "cutesv": {
                "executable": self.cutesv_executable, "resources": self.cutesv_resources.to_dict(),
                "minimum_support": self.cutesv_minimum_support, "minimum_sv_size": self.cutesv_minimum_sv_size,
                "max_cluster_bias_ins": self.cutesv_max_cluster_bias_ins, "diff_ratio_merging_ins": self.cutesv_diff_ratio_merging_ins,
                "max_cluster_bias_del": self.cutesv_max_cluster_bias_del, "diff_ratio_merging_del": self.cutesv_diff_ratio_merging_del,
                "genotype": self.cutesv_genotype,
            },
            "finalization": {"bgzip_executable": self.bgzip_executable, "tabix_executable": self.tabix_executable},
        }


CallerResult = SawfishResult | Sniffles2Result | PbsvResult | CuteSvResult


@dataclass(frozen=True, slots=True)
class Phase4SampleResult:
    alignment: AlignmentArtifact
    caller_results: tuple[CallerResult, ...]
    finalization_results: tuple[VcfFinalizeResult, ...]
    artifacts: tuple[StructuralVariantArtifact, ...]

    @property
    def sample_id(self) -> str:
        return self.alignment.sample_id

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "alignment": self.alignment.to_dict(),
            "caller_results": [result.to_dict() for result in self.caller_results],
            "finalization_results": [result.to_dict() for result in self.finalization_results],
            "artifacts": [artifact.to_dict() for artifact in self.artifacts],
        }


@dataclass(frozen=True, slots=True)
class Phase4RunReport:
    context: AnalysisContext
    settings: Phase4Settings
    sample_results: tuple[Phase4SampleResult, ...]
    dry_run: bool
    schema_version: str = PHASE4_REPORT_SCHEMA_VERSION
    hifivar_version: str = __version__
    created_at: str = field(default_factory=utc_now_iso8601)
    status: Phase4RunStatus = field(init=False)

    def __post_init__(self) -> None:
        if tuple(result.sample_id for result in self.sample_results) != self.context.sample_ids:
            raise InputValidationError("Phase 4 result order differs from AnalysisContext sample order.")
        expected_artifacts = 0 if self.dry_run else len(self.settings.enabled_callers)
        if any(len(result.artifacts) != expected_artifacts for result in self.sample_results):
            raise InputValidationError("Phase 4 artifact count conflicts with enabled callers/dry-run state.")
        object.__setattr__(self, "status", Phase4RunStatus.PLANNED if self.dry_run else Phase4RunStatus.COMPLETED)

    @property
    def tool_versions(self) -> dict[str, str | None]:
        versions: dict[str, str | None] = {caller.value: None for caller in self.settings.enabled_callers}
        for sample in self.sample_results:
            for artifact in sample.artifacts:
                versions[artifact.caller.value] = artifact.caller_version
            for finalization in sample.finalization_results:
                if finalization.tool_versions:
                    versions.update(finalization.tool_versions)
        return versions

    def to_dict(self) -> dict[str, object]:
        payload = {
            "schema_version": self.schema_version,
            "hifivar_version": self.hifivar_version,
            "created_at": self.created_at,
            "status": self.status.value,
            "dry_run": self.dry_run,
            "settings": self.settings.to_dict(),
            "tool_versions": self.tool_versions,
            "analysis_context": self.context.to_dict(),
            "sample_results": [result.to_dict() for result in self.sample_results],
        }
        standardized = standardize_data(payload, context="Phase 4 report value")
        if not isinstance(standardized, dict):
            raise InputValidationError("Phase 4 report serialization failed.")
        return standardized

    def write_json(self, path: str | Path, *, overwrite: bool = False) -> Path:
        return write_json_atomic(self.to_dict(), path, overwrite=overwrite, artifact_name="Phase 4 report")

    def write_yaml(self, path: str | Path, *, overwrite: bool = False) -> Path:
        return write_yaml_atomic(self.to_dict(), path, overwrite=overwrite, artifact_name="Phase 4 report")


def run_phase4(
    context: AnalysisContext,
    output_directory: str | Path,
    *,
    alignment_artifacts: Mapping[str, AlignmentArtifact] | None = None,
    settings: Phase4Settings | None = None,
    sawfish_wrapper: SawfishWrapper | None = None,
    sniffles2_wrapper: Sniffles2Wrapper | None = None,
    pbsv_wrapper: PbsvWrapper | None = None,
    cutesv_wrapper: CuteSvWrapper | None = None,
    finalizer: BgzipTabixWrapper | None = None,
    dry_run: bool = False,
) -> Phase4RunReport:
    if not isinstance(context, AnalysisContext):
        raise InputValidationError("run_phase4 requires AnalysisContext.")
    selected = settings or Phase4Settings.from_config(context.config)
    if not selected.enabled:
        raise ConfigurationError("Phase 4 requires sv.enabled: true.")
    if not selected.enabled_callers:
        raise ConfigurationError("Phase 4 requires at least one enabled SV caller.")
    alignments = _resolve_alignments(context, dict(alignment_artifacts or {}))
    if (selected.pbsv_enabled or selected.cutesv_enabled) and any(
        artifact.output_format is AlignmentOutputFormat.CRAM for artifact in alignments
    ):
        raise InputValidationError(
            "Phase 4 CRAM input is incompatible with enabled pbsv/cuteSV; "
            "disable those BAM-only callers or provide indexed BAM."
        )
    output_root = Path(output_directory).expanduser()
    sv_root = output_root / "sv"
    work_root = output_root / "work" / "sv"
    sawfish = sawfish_wrapper or SawfishWrapper(executable=selected.sawfish_executable)
    sniffles2 = sniffles2_wrapper or Sniffles2Wrapper(executable=selected.sniffles2_executable)
    pbsv = pbsv_wrapper or PbsvWrapper(executable=selected.pbsv_executable)
    cutesv = cutesv_wrapper or CuteSvWrapper(executable=selected.cutesv_executable)
    finalize = finalizer or BgzipTabixWrapper(bgzip_executable=selected.bgzip_executable, tabix_executable=selected.tabix_executable)
    sample_results: list[Phase4SampleResult] = []
    for alignment in alignments:
        caller_results: list[CallerResult] = []
        finalization_results: list[VcfFinalizeResult] = []
        artifacts: list[StructuralVariantArtifact] = []
        log_root = output_root / "logs" / "sv"
        if selected.sawfish_enabled:
            request = SawfishRequest.create(alignment, sv_root, work_root / "sawfish", resources=selected.sawfish_resources, overwrite=selected.overwrite, disable_cnv=selected.sawfish_disable_cnv)
            result = sawfish.run(request, dry_run=dry_run, stderr_path=log_root / f"{alignment.sample_id}.sawfish.log")
            caller_results.append(result)
            if not dry_run:
                artifacts.append(create_structural_variant_artifact(caller=SvCaller.SAWFISH, sample_id=alignment.sample_id, reference=alignment.reference, vcf_path=request.output_vcf, caller_version=result.tool_version, commands=tuple(command.args for command in result.commands)))
        if selected.sniffles2_enabled:
            request = Sniffles2Request.create(alignment, sv_root, resources=selected.sniffles2_resources, minimum_support=selected.sniffles2_minimum_support, minimum_sv_length=selected.sniffles2_minimum_sv_length, overwrite=selected.overwrite)
            result = sniffles2.run(request, dry_run=dry_run, stderr_path=log_root / f"{alignment.sample_id}.sniffles2.log")
            caller_results.append(result)
            if not dry_run:
                artifacts.append(create_structural_variant_artifact(caller=SvCaller.SNIFFLES2, sample_id=alignment.sample_id, reference=alignment.reference, vcf_path=request.output_vcf, caller_version=result.tool_version, commands=(result.command.args,)))
        if selected.pbsv_enabled:
            request = PbsvRequest.create(alignment, work_root / "pbsv" / "native", work_root / "pbsv" / "signatures", resources=selected.pbsv_resources, overwrite=selected.overwrite)
            result = pbsv.run(request, dry_run=dry_run, stderr_path=log_root / f"{alignment.sample_id}.pbsv.log")
            caller_results.append(result)
            final = finalize.run(VcfFinalizeRequest(SvCaller.PBSV, alignment.sample_id, alignment.reference, request.raw_vcf, sv_root / f"{alignment.sample_id}.pbsv.sv.vcf.gz", result.tool_version, tuple(command.args for command in result.commands), selected.overwrite), dry_run=dry_run, stderr_path=log_root / f"{alignment.sample_id}.pbsv.finalize.log")
            finalization_results.append(final)
            if final.artifact is not None:
                artifacts.append(final.artifact)
        if selected.cutesv_enabled:
            request = CuteSvRequest.create(
                alignment, work_root / "cutesv" / "native", work_root / "cutesv" / "tmp",
                resources=selected.cutesv_resources, minimum_support=selected.cutesv_minimum_support,
                minimum_sv_size=selected.cutesv_minimum_sv_size, max_cluster_bias_ins=selected.cutesv_max_cluster_bias_ins,
                diff_ratio_merging_ins=selected.cutesv_diff_ratio_merging_ins, max_cluster_bias_del=selected.cutesv_max_cluster_bias_del,
                diff_ratio_merging_del=selected.cutesv_diff_ratio_merging_del, genotype=selected.cutesv_genotype,
                overwrite=selected.overwrite,
            )
            result = cutesv.run(request, dry_run=dry_run, stderr_path=log_root / f"{alignment.sample_id}.cutesv.log")
            caller_results.append(result)
            final = finalize.run(VcfFinalizeRequest(SvCaller.CUTESV, alignment.sample_id, alignment.reference, request.raw_vcf, sv_root / f"{alignment.sample_id}.cutesv.sv.vcf.gz", result.tool_version, (result.command.args,), selected.overwrite), dry_run=dry_run, stderr_path=log_root / f"{alignment.sample_id}.cutesv.finalize.log")
            finalization_results.append(final)
            if final.artifact is not None:
                artifacts.append(final.artifact)
        sample_results.append(Phase4SampleResult(alignment, tuple(caller_results), tuple(finalization_results), tuple(artifacts)))
    return Phase4RunReport(context, selected, tuple(sample_results), dry_run)


def _resolve_alignments(context: AnalysisContext, supplied: Mapping[str, AlignmentArtifact]) -> tuple[AlignmentArtifact, ...]:
    unknown = sorted(set(supplied).difference(context.sample_ids))
    if unknown:
        raise InputValidationError(f"Phase 4 received artifacts for unknown samples: {unknown!r}.")
    resolved: list[AlignmentArtifact] = []
    for record in context.samples:
        sample = record.sample
        artifact = supplied.get(sample.sample_id)
        if artifact is None:
            if sample.input.input_type is InputType.FASTQ:
                raise InputValidationError(f"Phase 4 sample '{sample.sample_id}' has raw FASTQ; provide a completed Phase 2 AlignmentArtifact.")
            path = sample.input.files[0]
            artifact = AlignmentArtifact(
                sample.sample_id, path,
                AlignmentOutputFormat.BAM if sample.input.input_type is InputType.BAM else AlignmentOutputFormat.CRAM,
                context.reference, AlignmentSource.EXISTING, AlignmentSortOrder.UNKNOWN, find_alignment_index(path),
            )
        if artifact.sample_id != sample.sample_id:
            raise InputValidationError(f"Phase 4 alignment sample mismatch for '{sample.sample_id}'.")
        if artifact.reference.fasta.absolute() != context.reference.fasta.absolute() or artifact.reference.build != context.reference.build or artifact.reference.contigs != context.reference.contigs:
            raise ReferenceError(f"Phase 4 alignment reference conflicts for sample '{sample.sample_id}'.")
        resolved.append(validate_alignment_artifact(artifact, require_index=True))
    return tuple(resolved)


def _mapping(mapping: Mapping[str, object], key: str) -> Mapping[str, object]:
    value = mapping.get(key)
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"Phase 4 requires configuration mapping {key}.")
    return value


__all__ = ["PHASE4_REPORT_SCHEMA_VERSION", "Phase4RunReport", "Phase4RunStatus", "Phase4SampleResult", "Phase4Settings", "run_phase4"]
