"""Integrated Phase 3 DeepVariant orchestration and run provenance."""

from __future__ import annotations

import os
from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Final

from hifivar import __version__
from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import (
    AlignmentArtifact,
    AlignmentSortOrder,
    AlignmentSource,
    find_alignment_index,
    validate_alignment_artifact,
)
from hifivar.context import AnalysisContext
from hifivar.deepvariant import DeepVariantRuntime, DeepVariantWrapper
from hifivar.exceptions import ConfigurationError, InputValidationError, ReferenceError
from hifivar.phase2 import Phase2RunReport
from hifivar.sample import InputType
from hifivar.serialization import (
    standardize_data,
    utc_now_iso8601,
    write_json_atomic,
    write_yaml_atomic,
)
from hifivar.small import (
    DeepVariantRequest,
    SmallVariantModelType,
    SmallVariantResources,
    SmallVariantResult,
    SmallVariantResultStatus,
)


PHASE3_REPORT_SCHEMA_VERSION: Final[str] = "1.0"


class Phase3RunStatus(str, Enum):
    """Run-level state for planned or completed Phase 3 calls."""

    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class Phase3Settings:
    """Validated effective settings for DeepVariant Phase 3 execution."""

    runtime: DeepVariantRuntime
    resources: SmallVariantResources
    model_type: SmallVariantModelType = SmallVariantModelType.PACBIO
    overwrite: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.runtime, DeepVariantRuntime):
            raise ConfigurationError("Phase 3 runtime must be DeepVariantRuntime.")
        if not isinstance(self.resources, SmallVariantResources):
            raise ConfigurationError(
                "Phase 3 resources must be SmallVariantResources."
            )
        if self.model_type is not SmallVariantModelType.PACBIO:
            raise ConfigurationError("Phase 3 model_type must be PACBIO.")
        if not isinstance(self.overwrite, bool):
            raise ConfigurationError("Phase 3 overwrite must be boolean.")

    @classmethod
    def from_config(cls, config: Mapping[str, object]) -> Phase3Settings:
        if not isinstance(config, Mapping):
            raise ConfigurationError("Phase 3 settings require a config mapping.")
        section = config.get("small")
        if not isinstance(section, Mapping):
            raise ConfigurationError(
                "Phase 3 settings require configuration section small."
            )
        try:
            model_type = SmallVariantModelType(str(section["model_type"]).upper())
            resources = SmallVariantResources(
                threads=section["threads"],  # type: ignore[arg-type]
                memory_mb=section["memory_mb"],  # type: ignore[arg-type]
                runtime_minutes=section["runtime_minutes"],  # type: ignore[arg-type]
            )
            overwrite = section["overwrite"]
        except (KeyError, TypeError, ValueError) as error:
            raise ConfigurationError(
                f"Invalid or incomplete Phase 3 small configuration: {error}"
            ) from error
        if not isinstance(overwrite, bool):
            raise ConfigurationError("small.overwrite must be boolean.")
        return cls(
            runtime=DeepVariantRuntime.from_config(config),
            resources=resources,
            model_type=model_type,
            overwrite=overwrite,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "runtime": self.runtime.to_dict(),
            "resources": self.resources.to_dict(),
            "model_type": self.model_type.value,
            "overwrite": self.overwrite,
        }


@dataclass(frozen=True, slots=True)
class Phase3SampleResult:
    """One ordered alignment-to-small-variant handoff."""

    alignment: AlignmentArtifact
    call: SmallVariantResult

    def __post_init__(self) -> None:
        if not isinstance(self.alignment, AlignmentArtifact):
            raise InputValidationError(
                "Phase 3 sample alignment must be AlignmentArtifact."
            )
        if not isinstance(self.call, SmallVariantResult):
            raise InputValidationError(
                "Phase 3 sample call must be SmallVariantResult."
            )
        if self.alignment.sample_id != self.call.request.sample_id:
            raise InputValidationError(
                "Phase 3 alignment and DeepVariant result samples differ."
            )

    @property
    def sample_id(self) -> str:
        return self.alignment.sample_id

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "alignment": self.alignment.to_dict(),
            "call": self.call.to_dict(),
        }


@dataclass(frozen=True, slots=True)
class Phase3RunReport:
    """Versioned run-level DeepVariant command and artifact provenance."""

    context: AnalysisContext
    settings: Phase3Settings
    sample_results: tuple[Phase3SampleResult, ...]
    dry_run: bool
    schema_version: str = PHASE3_REPORT_SCHEMA_VERSION
    hifivar_version: str = __version__
    created_at: str = field(default_factory=utc_now_iso8601)
    status: Phase3RunStatus = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.context, AnalysisContext):
            raise InputValidationError("Phase 3 report context must be AnalysisContext.")
        if not isinstance(self.settings, Phase3Settings):
            raise InputValidationError("Phase 3 report settings must be Phase3Settings.")
        results = tuple(self.sample_results)
        if any(not isinstance(item, Phase3SampleResult) for item in results):
            raise InputValidationError(
                "Phase 3 report results must contain Phase3SampleResult."
            )
        if tuple(item.sample_id for item in results) != self.context.sample_ids:
            raise InputValidationError(
                "Phase 3 result order differs from AnalysisContext sample order."
            )
        if not isinstance(self.dry_run, bool):
            raise InputValidationError("Phase 3 report dry_run must be boolean.")
        expected = (
            SmallVariantResultStatus.PLANNED
            if self.dry_run
            else SmallVariantResultStatus.COMPLETED
        )
        if any(item.call.status is not expected for item in results):
            raise InputValidationError(
                "Phase 3 result status conflicts with report dry_run state."
            )
        object.__setattr__(self, "sample_results", results)
        object.__setattr__(
            self,
            "status",
            Phase3RunStatus.PLANNED if self.dry_run else Phase3RunStatus.COMPLETED,
        )

    @property
    def tool_versions(self) -> dict[str, str | None]:
        version = next(
            (
                item.call.tool_version
                for item in self.sample_results
                if item.call.tool_version is not None
            ),
            None,
        )
        return {"deepvariant": version}

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
            "sample_results": [item.to_dict() for item in self.sample_results],
        }
        standardized = standardize_data(payload, context="Phase 3 report value")
        if not isinstance(standardized, dict):  # pragma: no cover
            raise InputValidationError("Phase 3 report serialization failed.")
        return standardized

    def write_json(self, path: str | Path, *, overwrite: bool = False) -> Path:
        return write_json_atomic(
            self.to_dict(),
            path,
            overwrite=overwrite,
            artifact_name="Phase 3 report",
        )

    def write_yaml(self, path: str | Path, *, overwrite: bool = False) -> Path:
        return write_yaml_atomic(
            self.to_dict(),
            path,
            overwrite=overwrite,
            artifact_name="Phase 3 report",
        )


def collect_phase2_alignment_artifacts(
    report: Phase2RunReport,
) -> dict[str, AlignmentArtifact]:
    """Collect completed Phase 2 artifacts for an explicit Phase 3 handoff."""
    if not isinstance(report, Phase2RunReport):
        raise InputValidationError(
            "Phase 2 handoff requires a Phase2RunReport."
        )
    artifacts: dict[str, AlignmentArtifact] = {}
    for result in report.sample_results:
        if result.artifact is None:
            raise InputValidationError(
                f"Phase 2 has no completed alignment artifact for sample "
                f"'{result.sample_id}'; dry-run results cannot enter Phase 3."
            )
        artifacts[result.sample_id] = result.artifact
    return artifacts


def run_phase3(
    context: AnalysisContext,
    output_directory: str | Path,
    *,
    alignment_artifacts: Mapping[str, AlignmentArtifact] | None = None,
    settings: Phase3Settings | None = None,
    deepvariant_wrapper: DeepVariantWrapper | None = None,
    dry_run: bool = False,
) -> Phase3RunReport:
    """Call DeepVariant for all context samples in deterministic order."""
    if not isinstance(context, AnalysisContext):
        raise InputValidationError("run_phase3 requires AnalysisContext.")
    if not isinstance(output_directory, (str, Path)) or (
        isinstance(output_directory, str) and not output_directory.strip()
    ):
        raise InputValidationError(
            "Phase 3 output_directory must be a non-empty string or Path."
        )
    if not isinstance(dry_run, bool):
        raise InputValidationError("Phase 3 dry_run must be boolean.")
    selected_settings = settings or Phase3Settings.from_config(context.config)
    if not isinstance(selected_settings, Phase3Settings):
        raise InputValidationError("settings must be Phase3Settings or None.")
    supplied = dict(alignment_artifacts or {})
    unknown = sorted(set(supplied).difference(context.sample_ids))
    if unknown:
        raise InputValidationError(
            f"Phase 3 received alignment artifacts for unknown samples: {unknown!r}."
        )
    alignments = _resolve_alignments(context, supplied)
    output_root = Path(output_directory).expanduser()
    wrapper = deepvariant_wrapper or DeepVariantWrapper(
        runtime=selected_settings.runtime
    )
    results: list[Phase3SampleResult] = []
    for alignment in alignments:
        request = DeepVariantRequest.create(
            alignment,
            output_root / "small",
            resources=selected_settings.resources,
            model_type=selected_settings.model_type,
            overwrite=selected_settings.overwrite,
        )
        call = wrapper.run(
            request,
            dry_run=dry_run,
            stderr_path=output_root / "logs" / "small" / f"{alignment.sample_id}.deepvariant.log",
        )
        results.append(Phase3SampleResult(alignment=alignment, call=call))
    return Phase3RunReport(
        context=context,
        settings=selected_settings,
        sample_results=tuple(results),
        dry_run=dry_run,
    )


def _resolve_alignments(
    context: AnalysisContext,
    supplied: Mapping[str, AlignmentArtifact],
) -> tuple[AlignmentArtifact, ...]:
    resolved: list[AlignmentArtifact] = []
    for record in context.samples:
        sample = record.sample
        artifact = supplied.get(sample.sample_id)
        if artifact is None:
            if sample.input.input_type is InputType.FASTQ:
                raise InputValidationError(
                    f"Phase 3 sample '{sample.sample_id}' has raw FASTQ input; "
                    "provide its completed Phase 2 AlignmentArtifact."
                )
            path = sample.input.files[0]
            artifact = AlignmentArtifact(
                sample_id=sample.sample_id,
                path=path,
                output_format=(
                    AlignmentOutputFormat.BAM
                    if sample.input.input_type is InputType.BAM
                    else AlignmentOutputFormat.CRAM
                ),
                reference=context.reference,
                source=AlignmentSource.EXISTING,
                sort_order=AlignmentSortOrder.UNKNOWN,
                index_path=find_alignment_index(path),
            )
        _validate_handoff(context, sample.sample_id, artifact)
        validated = validate_alignment_artifact(artifact, require_index=True)
        if validated.index_path is None:  # pragma: no cover - validation invariant
            raise InputValidationError("Validated Phase 3 alignment lacks an index.")
        resolved.append(validated)
    return tuple(resolved)


def _validate_handoff(
    context: AnalysisContext,
    sample_id: str,
    artifact: AlignmentArtifact,
) -> None:
    if not isinstance(artifact, AlignmentArtifact):
        raise InputValidationError(
            f"Phase 3 alignment handoff for '{sample_id}' is not AlignmentArtifact."
        )
    if artifact.sample_id != sample_id:
        raise InputValidationError(
            f"Phase 3 alignment handoff sample mismatch: '{artifact.sample_id}' "
            f"!= '{sample_id}'."
        )
    expected = context.reference
    observed = artifact.reference
    if (
        os.path.normcase(os.path.normpath(str(observed.fasta.absolute())))
        != os.path.normcase(os.path.normpath(str(expected.fasta.absolute())))
        or observed.build != expected.build
        or observed.contigs != expected.contigs
    ):
        raise ReferenceError(
            f"Phase 3 alignment reference conflicts with AnalysisContext for "
            f"sample '{sample_id}'."
        )


__all__ = [
    "PHASE3_REPORT_SCHEMA_VERSION",
    "Phase3RunReport",
    "Phase3RunStatus",
    "Phase3SampleResult",
    "Phase3Settings",
    "collect_phase2_alignment_artifacts",
    "run_phase3",
]
