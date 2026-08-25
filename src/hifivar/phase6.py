"""Phase 6 HiPhase orchestration and provenance."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Final

from hifivar import __version__
from hifivar.alignment_postprocess import AlignmentArtifact
from hifivar.context import AnalysisContext
from hifivar.exceptions import ConfigurationError, InputValidationError
from hifivar.hiphase import HiPhaseWrapper, PhasingResult, PhasingResultStatus
from hifivar.phasing import PhasingRequest, PhasingResources
from hifivar.serialization import (
    standardize_data,
    utc_now_iso8601,
    write_json_atomic,
    write_yaml_atomic,
)
from hifivar.small import SmallVariantArtifact


PHASE6_REPORT_SCHEMA_VERSION: Final[str] = "1.0"


class Phase6RunStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class Phase6Settings:
    backend: str
    executable: str
    tabix_executable: str
    resources: PhasingResources
    overwrite: bool = False

    def __post_init__(self) -> None:
        if self.backend != "hiphase":
            raise ConfigurationError("Phase 6 backend must be hiphase.")
        for name in ("executable", "tabix_executable"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ConfigurationError(f"Phase 6 {name} must be non-empty.")
        if not isinstance(self.resources, PhasingResources):
            raise ConfigurationError("Phase 6 resources must be PhasingResources.")
        if not isinstance(self.overwrite, bool):
            raise ConfigurationError("Phase 6 overwrite must be boolean.")

    @classmethod
    def from_config(cls, config: Mapping[str, object]) -> "Phase6Settings":
        section = config.get("phasing")
        if not isinstance(section, Mapping):
            raise ConfigurationError("Phase 6 requires configuration mapping phasing.")
        try:
            return cls(
                backend=str(section["backend"]),
                executable=str(section["executable"]),
                tabix_executable=str(section["tabix_executable"]),
                resources=PhasingResources(
                    threads=section["threads"],  # type: ignore[arg-type]
                    memory_mb=section["memory_mb"],  # type: ignore[arg-type]
                    runtime_minutes=section["runtime_minutes"],  # type: ignore[arg-type]
                ),
                overwrite=section["overwrite"],  # type: ignore[arg-type]
            )
        except (KeyError, TypeError, ValueError) as error:
            raise ConfigurationError(f"Invalid Phase 6 settings: {error}") from error

    def to_dict(self) -> dict[str, object]:
        return {
            "backend": self.backend,
            "executable": self.executable,
            "tabix_executable": self.tabix_executable,
            "resources": self.resources.to_dict(),
            "overwrite": self.overwrite,
        }


@dataclass(frozen=True, slots=True)
class Phase6SampleResult:
    sample_id: str
    phasing: PhasingResult

    def __post_init__(self) -> None:
        if self.sample_id != self.phasing.request.sample_id:
            raise InputValidationError("Phase 6 sample/result IDs differ.")

    def to_dict(self) -> dict[str, object]:
        return {"sample_id": self.sample_id, "phasing": self.phasing.to_dict()}


@dataclass(frozen=True, slots=True)
class Phase6RunReport:
    context: AnalysisContext
    settings: Phase6Settings
    sample_results: tuple[Phase6SampleResult, ...]
    dry_run: bool
    created_at: str
    hifivar_version: str = __version__
    schema_version: str = PHASE6_REPORT_SCHEMA_VERSION

    def __post_init__(self) -> None:
        expected = tuple(record.sample.sample_id for record in self.context.samples)
        observed = tuple(item.sample_id for item in self.sample_results)
        if expected != observed:
            raise InputValidationError("Phase 6 result order differs from AnalysisContext.")
        planned = all(
            item.phasing.status is PhasingResultStatus.PLANNED
            for item in self.sample_results
        )
        if planned != self.dry_run:
            raise InputValidationError("Phase 6 report status conflicts with dry_run.")

    @property
    def status(self) -> Phase6RunStatus:
        return Phase6RunStatus.PLANNED if self.dry_run else Phase6RunStatus.COMPLETED

    def to_dict(self) -> dict[str, object]:
        payload = {
            "schema_version": self.schema_version,
            "hifivar_version": self.hifivar_version,
            "created_at": self.created_at,
            "status": self.status.value,
            "dry_run": self.dry_run,
            "context": self.context.to_dict(),
            "settings": self.settings.to_dict(),
            "sample_results": [item.to_dict() for item in self.sample_results],
        }
        standardized = standardize_data(payload, context="Phase 6 report value")
        if not isinstance(standardized, dict):
            raise InputValidationError("Phase 6 report serialization failed.")
        return standardized

    def write_json(self, path: Path, *, overwrite: bool = False) -> Path:
        return write_json_atomic(
            self.to_dict(), path, overwrite=overwrite, artifact_name="Phase 6 report"
        )

    def write_yaml(self, path: Path, *, overwrite: bool = False) -> Path:
        return write_yaml_atomic(
            self.to_dict(), path, overwrite=overwrite, artifact_name="Phase 6 report"
        )


def run_phase6(
    context: AnalysisContext,
    *,
    alignment_artifacts: Mapping[str, AlignmentArtifact],
    small_variant_artifacts: Mapping[str, SmallVariantArtifact],
    output_directory: str | Path,
    config: Mapping[str, object],
    dry_run: bool = False,
    wrapper: HiPhaseWrapper | None = None,
) -> Phase6RunReport:
    """Phase every context sample in deterministic order."""
    section = config.get("phasing")
    if not isinstance(section, Mapping) or section.get("enabled") is not True:
        raise ConfigurationError("Phase 6 requires phasing.enabled: true.")
    settings = Phase6Settings.from_config(config)
    output_root = Path(output_directory).expanduser()
    engine = wrapper or HiPhaseWrapper(
        executable=settings.executable,
        tabix_executable=settings.tabix_executable,
    )
    expected = {record.sample.sample_id for record in context.samples}
    unknown = sorted(
        (set(alignment_artifacts) | set(small_variant_artifacts)).difference(expected)
    )
    if unknown:
        raise InputValidationError(f"Phase 6 received artifacts for unknown samples: {unknown!r}.")

    results: list[Phase6SampleResult] = []
    for record in context.samples:
        sample_id = record.sample.sample_id
        try:
            alignment = alignment_artifacts[sample_id]
            small = small_variant_artifacts[sample_id]
        except KeyError as error:
            raise InputValidationError(
                f"Phase 6 requires BAM and small-variant artifacts for '{sample_id}'."
            ) from error
        request = PhasingRequest(
            alignment,
            small,
            output_root / f"{sample_id}.phased.vcf.gz",
            settings.resources,
            settings.overwrite,
        )
        result = engine.run(
            request,
            dry_run=dry_run,
            stderr_path=output_root / "logs" / f"{sample_id}.hiphase.log",
        )
        results.append(Phase6SampleResult(sample_id, result))
    return Phase6RunReport(
        context,
        settings,
        tuple(results),
        dry_run,
        utc_now_iso8601(),
    )


__all__ = [
    "PHASE6_REPORT_SCHEMA_VERSION",
    "Phase6RunReport",
    "Phase6RunStatus",
    "Phase6SampleResult",
    "Phase6Settings",
    "run_phase6",
]
