"""Phase 2 QC, alignment, indexing, and provenance orchestration."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

from hifivar import __version__
from hifivar.alignment import (
    AlignmentAction,
    AlignmentOutputFormat,
    AlignmentPlan,
    AlignmentResources,
    AlignmentResult,
    AlignmentResultStatus,
    AlignmentTool,
    build_alignment_plans,
)
from hifivar.alignment_postprocess import (
    AlignmentArtifact,
    AlignmentIndexFormat,
    AlignmentIndexRequest,
    AlignmentSortOrder,
    AlignmentSource,
    run_alignment_qc,
    validate_alignment_artifact,
)
from hifivar.context import AnalysisContext
from hifivar.exceptions import ConfigurationError, InputValidationError
from hifivar.pbmm2 import Pbmm2Options, Pbmm2Wrapper
from hifivar.qc import (
    QCResult,
    QCStatus,
    RunQCReport,
    aggregate_qc_status,
    run_input_qc,
)
from hifivar.samtools import SamtoolsIndexResult, SamtoolsWrapper
from hifivar.serialization import (
    standardize_data,
    utc_now_iso8601,
    write_json_atomic,
    write_yaml_atomic,
)


PHASE2_REPORT_SCHEMA_VERSION: Final[str] = "1.0"


@dataclass(frozen=True, slots=True)
class Phase2Settings:
    """Effective settings consumed by the implemented Phase 2 pipeline."""

    tool: AlignmentTool
    output_format: AlignmentOutputFormat
    resources: AlignmentResources
    overwrite: bool
    index_threads: int
    bam_index_format: AlignmentIndexFormat | None = None

    def __post_init__(self) -> None:
        """Restrict orchestration to the concrete Phase 2 implementation."""
        if self.tool is not AlignmentTool.PBMM2:
            raise ConfigurationError(
                "Phase 2 execution currently supports alignment.tool=pbmm2 only."
            )
        if self.output_format is not AlignmentOutputFormat.BAM:
            raise ConfigurationError(
                "Phase 2 pbmm2 execution currently requires output_format=bam."
            )
        if not isinstance(self.resources, AlignmentResources):
            raise ConfigurationError(
                "Phase 2 resources must be AlignmentResources."
            )
        if not isinstance(self.overwrite, bool):
            raise ConfigurationError("Phase 2 overwrite must be a boolean.")
        if (
            not isinstance(self.index_threads, int)
            or isinstance(self.index_threads, bool)
            or self.index_threads <= 0
        ):
            raise ConfigurationError(
                "Phase 2 index_threads must be a positive integer."
            )
        if self.bam_index_format is not None and self.bam_index_format not in {
            AlignmentIndexFormat.BAI,
            AlignmentIndexFormat.CSI,
        }:
            raise ConfigurationError(
                "Phase 2 BAM index format must be BAI, CSI, or automatic."
            )

    @classmethod
    def from_config(cls, config: Mapping[str, object]) -> Phase2Settings:
        """Read validated Phase 2 fields from an effective config mapping."""
        if not isinstance(config, Mapping):
            raise ConfigurationError("Phase 2 settings require a config mapping.")
        section = config.get("alignment")
        if not isinstance(section, Mapping):
            raise ConfigurationError(
                "Phase 2 settings require configuration section alignment."
            )
        try:
            tool = AlignmentTool(section["tool"])
            output_format = AlignmentOutputFormat(section["output_format"])
            resources = AlignmentResources(
                threads=section["threads"],  # type: ignore[arg-type]
                memory_mb=section["memory_mb"],  # type: ignore[arg-type]
                runtime_minutes=section["runtime_minutes"],  # type: ignore[arg-type]
            )
            overwrite = section["overwrite"]
            index_threads = section["index_threads"]
            configured_index = section["bam_index_format"]
        except (KeyError, TypeError, ValueError) as error:
            raise ConfigurationError(
                f"Invalid or incomplete Phase 2 alignment configuration: {error}"
            ) from error
        if not isinstance(overwrite, bool):
            raise ConfigurationError("alignment.overwrite must be a boolean.")
        if not isinstance(configured_index, str):
            raise ConfigurationError(
                "alignment.bam_index_format must be a string."
            )
        try:
            index_format = (
                None
                if configured_index.casefold() == "auto"
                else AlignmentIndexFormat(configured_index.casefold())
            )
        except ValueError as error:
            raise ConfigurationError(
                "alignment.bam_index_format must be auto, bai, or csi."
            ) from error
        return cls(
            tool=tool,
            output_format=output_format,
            resources=resources,
            overwrite=overwrite,
            index_threads=index_threads,  # type: ignore[arg-type]
            bam_index_format=index_format,
        )

    def to_dict(self) -> dict[str, object]:
        """Return standard effective Phase 2 settings."""
        return {
            "tool": self.tool.value,
            "output_format": self.output_format.value,
            "resources": self.resources.to_dict(),
            "overwrite": self.overwrite,
            "index_threads": self.index_threads,
            "bam_index_format": (
                self.bam_index_format.value
                if self.bam_index_format is not None
                else "auto"
            ),
        }


@dataclass(frozen=True, slots=True)
class Phase2SampleResult:
    """Ordered per-sample Phase 2 handoff and QC provenance."""

    input_qc: QCResult
    plan: AlignmentPlan
    alignment_result: AlignmentResult
    alignment_qc: QCResult
    artifact: AlignmentArtifact | None = None
    index_result: SamtoolsIndexResult | None = None

    def __post_init__(self) -> None:
        """Require all per-sample records to refer to one sample."""
        if not isinstance(self.input_qc, QCResult):
            raise InputValidationError("Phase 2 input_qc must be QCResult.")
        if not isinstance(self.plan, AlignmentPlan):
            raise InputValidationError("Phase 2 plan must be AlignmentPlan.")
        if not isinstance(self.alignment_result, AlignmentResult):
            raise InputValidationError(
                "Phase 2 alignment_result must be AlignmentResult."
            )
        if not isinstance(self.alignment_qc, QCResult):
            raise InputValidationError("Phase 2 alignment_qc must be QCResult.")
        sample_id = self.plan.sample_id
        if self.input_qc.sample_id != sample_id:
            raise InputValidationError(
                "Phase 2 input QC sample does not match its alignment plan."
            )
        if self.alignment_result.plan != self.plan:
            raise InputValidationError(
                "Phase 2 alignment result does not match its plan."
            )
        if self.alignment_qc.sample_id != sample_id:
            raise InputValidationError(
                "Phase 2 alignment QC sample does not match its plan."
            )
        if self.artifact is not None and self.artifact.sample_id != sample_id:
            raise InputValidationError(
                "Phase 2 artifact sample does not match its plan."
            )

    @property
    def sample_id(self) -> str:
        """Return the shared sample identifier."""
        return self.plan.sample_id

    def to_dict(self) -> dict[str, object]:
        """Return standard per-sample pipeline provenance."""
        return {
            "sample_id": self.sample_id,
            "input_qc": self.input_qc.to_dict(),
            "alignment_plan": self.plan.to_dict(),
            "alignment_result": self.alignment_result.to_dict(),
            "index_result": (
                self.index_result.to_dict()
                if self.index_result is not None
                else None
            ),
            "artifact": self.artifact.to_dict() if self.artifact else None,
            "alignment_qc": self.alignment_qc.to_dict(),
        }


@dataclass(frozen=True, slots=True)
class Phase2RunReport:
    """Serializable run-level Phase 2 integration result."""

    input_qc: RunQCReport
    sample_results: tuple[Phase2SampleResult, ...]
    settings: Phase2Settings
    dry_run: bool
    schema_version: str = PHASE2_REPORT_SCHEMA_VERSION
    hifivar_version: str = __version__
    created_at: str = field(default_factory=utc_now_iso8601)
    overall_qc_status: QCStatus = field(init=False)

    def __post_init__(self) -> None:
        """Validate ordered samples and deterministically aggregate QC."""
        if not isinstance(self.input_qc, RunQCReport):
            raise InputValidationError(
                "Phase 2 report input_qc must be RunQCReport."
            )
        results = tuple(self.sample_results)
        if any(not isinstance(item, Phase2SampleResult) for item in results):
            raise InputValidationError(
                "Phase 2 report sample_results must contain Phase2SampleResult."
            )
        expected = tuple(
            result.sample_id for result in self.input_qc.sample_results
        )
        observed = tuple(result.sample_id for result in results)
        if observed != expected:
            raise InputValidationError(
                "Phase 2 report sample order conflicts with input QC."
            )
        if not isinstance(self.settings, Phase2Settings):
            raise InputValidationError(
                "Phase 2 report settings must be Phase2Settings."
            )
        if not isinstance(self.dry_run, bool):
            raise InputValidationError("Phase 2 report dry_run must be boolean.")
        object.__setattr__(self, "sample_results", results)
        object.__setattr__(
            self,
            "overall_qc_status",
            aggregate_qc_status(
                (
                    self.input_qc.overall_status,
                    *(result.alignment_qc.status for result in results),
                )
            ),
        )

    @property
    def tool_versions(self) -> dict[str, str | None]:
        """Return run-level external tool versions when actually executed."""
        pbmm2_version = next(
            (
                item.alignment_result.tool_version
                for item in self.sample_results
                if item.alignment_result.tool_version is not None
            ),
            None,
        )
        samtools_version = next(
            (
                item.index_result.tool_version
                for item in self.sample_results
                if item.index_result is not None
                and item.index_result.tool_version is not None
            ),
            None,
        )
        return {"pbmm2": pbmm2_version, "samtools": samtools_version}

    def to_dict(self) -> dict[str, object]:
        """Return an independent JSON/YAML-friendly report payload."""
        payload = {
            "schema_version": self.schema_version,
            "hifivar_version": self.hifivar_version,
            "created_at": self.created_at,
            "dry_run": self.dry_run,
            "overall_qc_status": self.overall_qc_status.value,
            "settings": self.settings.to_dict(),
            "tool_versions": self.tool_versions,
            "input_qc": self.input_qc.to_dict(),
            "sample_results": [item.to_dict() for item in self.sample_results],
        }
        standardized = standardize_data(payload, context="Phase 2 report value")
        if not isinstance(standardized, dict):  # pragma: no cover
            raise InputValidationError("Phase 2 report serialization failed.")
        return standardized

    def write_json(self, path: str | Path, *, overwrite: bool = False) -> Path:
        """Atomically write a UTF-8 Phase 2 JSON report."""
        return write_json_atomic(
            self.to_dict(),
            path,
            overwrite=overwrite,
            artifact_name="Phase 2 report",
        )

    def write_yaml(self, path: str | Path, *, overwrite: bool = False) -> Path:
        """Atomically write a UTF-8 Phase 2 YAML report."""
        return write_yaml_atomic(
            self.to_dict(),
            path,
            overwrite=overwrite,
            artifact_name="Phase 2 report",
        )


def run_phase2(
    context: AnalysisContext,
    output_directory: str | Path,
    *,
    settings: Phase2Settings | None = None,
    pbmm2_wrapper: Pbmm2Wrapper | None = None,
    samtools_wrapper: SamtoolsWrapper | None = None,
    dry_run: bool = False,
) -> Phase2RunReport:
    """Run the complete Phase 2 path without invoking variant callers."""
    if not isinstance(context, AnalysisContext):
        raise InputValidationError("run_phase2 requires AnalysisContext.")
    if not isinstance(output_directory, (str, Path)) or (
        isinstance(output_directory, str) and not output_directory.strip()
    ):
        raise InputValidationError(
            "Phase 2 output_directory must be a non-empty string or Path."
        )
    if not isinstance(dry_run, bool):
        raise InputValidationError("Phase 2 dry_run must be boolean.")
    selected_settings = settings or Phase2Settings.from_config(context.config)
    if not isinstance(selected_settings, Phase2Settings):
        raise InputValidationError("settings must be Phase2Settings or None.")
    output_root = Path(output_directory).expanduser()
    input_qc = run_input_qc(context)
    plans = build_alignment_plans(
        context,
        output_root,
        tool=selected_settings.tool,
        output_format=selected_settings.output_format,
        resources=selected_settings.resources,
        overwrite=selected_settings.overwrite,
    )
    pbmm2 = pbmm2_wrapper or Pbmm2Wrapper(
        options=Pbmm2Options.from_config(context.config)
    )
    samtools = samtools_wrapper or SamtoolsWrapper()
    input_results = {
        result.sample_id: result for result in input_qc.sample_results
    }

    sample_results: list[Phase2SampleResult] = []
    for plan in plans:
        if plan.action is AlignmentAction.REUSE:
            alignment_result = AlignmentResult(
                plan=plan,
                status=AlignmentResultStatus.REUSED,
            )
            artifact = validate_alignment_artifact(
                AlignmentArtifact.from_result(alignment_result)
            )
            alignment_qc = run_alignment_qc(artifact)
            sample_results.append(
                Phase2SampleResult(
                    input_qc=input_results[plan.sample_id],
                    plan=plan,
                    alignment_result=alignment_result,
                    artifact=artifact,
                    index_result=None,
                    alignment_qc=alignment_qc,
                )
            )
            continue

        if plan.request is None:  # pragma: no cover - model invariant
            raise InputValidationError("ALIGN plan is missing its request.")
        alignment_result = pbmm2.run(
            plan.request,
            dry_run=dry_run,
            stderr_path=output_root / "logs" / f"{plan.sample_id}.pbmm2.log",
        )
        if dry_run:
            provisional = AlignmentArtifact(
                sample_id=plan.sample_id,
                path=plan.alignment_path,
                output_format=plan.output_format,
                reference=plan.reference,
                source=AlignmentSource.GENERATED,
                sort_order=AlignmentSortOrder.COORDINATE,
                tool=AlignmentTool.PBMM2,
            )
            index_request = AlignmentIndexRequest.create(
                provisional,
                index_format=selected_settings.bam_index_format,
                threads=selected_settings.index_threads,
                overwrite=selected_settings.overwrite,
            )
            index_result = samtools.run_index(
                index_request,
                dry_run=True,
                stderr_path=(
                    output_root / "logs" / f"{plan.sample_id}.samtools-index.log"
                ),
            )
            alignment_qc = QCResult(
                status=QCStatus.NOT_CHECKED,
                module="alignment_qc",
                sample_id=plan.sample_id,
            )
            artifact = None
        else:
            artifact = validate_alignment_artifact(
                AlignmentArtifact.from_result(alignment_result)
            )
            index_request = AlignmentIndexRequest.create(
                artifact,
                index_format=selected_settings.bam_index_format,
                threads=selected_settings.index_threads,
                overwrite=selected_settings.overwrite,
            )
            index_result = samtools.run_index(
                index_request,
                stderr_path=(
                    output_root / "logs" / f"{plan.sample_id}.samtools-index.log"
                ),
            )
            artifact = validate_alignment_artifact(
                index_result.artifact,
                require_index=True,
            )
            alignment_qc = run_alignment_qc(artifact)
        sample_results.append(
            Phase2SampleResult(
                input_qc=input_results[plan.sample_id],
                plan=plan,
                alignment_result=alignment_result,
                artifact=artifact,
                index_result=index_result,
                alignment_qc=alignment_qc,
            )
        )

    return Phase2RunReport(
        input_qc=input_qc,
        sample_results=tuple(sample_results),
        settings=selected_settings,
        dry_run=dry_run,
    )


__all__ = [
    "PHASE2_REPORT_SCHEMA_VERSION",
    "Phase2RunReport",
    "Phase2SampleResult",
    "Phase2Settings",
    "run_phase2",
]
