"""Phase 7 hifiasm assembly orchestration and provenance."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Final

from hifivar import __version__
from hifivar.assembly import AssemblyRequest, AssemblyResources
from hifivar.context import AnalysisContext
from hifivar.exceptions import ConfigurationError, InputValidationError
from hifivar.hifiasm import AssemblyResult, AssemblyResultStatus, HifiasmWrapper
from hifivar.sample import InputType
from hifivar.serialization import (
    standardize_data,
    utc_now_iso8601,
    write_json_atomic,
    write_yaml_atomic,
)


PHASE7_REPORT_SCHEMA_VERSION: Final[str] = "1.0"


class Phase7RunStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class Phase7Settings:
    backend: str
    executable: str
    resources: AssemblyResources
    overwrite: bool = False

    def __post_init__(self) -> None:
        if self.backend != "hifiasm":
            raise ConfigurationError("Phase 7 backend must be hifiasm.")
        if not isinstance(self.executable, str) or not self.executable.strip():
            raise ConfigurationError("Phase 7 executable must be non-empty.")
        if not isinstance(self.resources, AssemblyResources):
            raise ConfigurationError("Phase 7 resources must be AssemblyResources.")
        if not isinstance(self.overwrite, bool):
            raise ConfigurationError("Phase 7 overwrite must be boolean.")

    @classmethod
    def from_config(cls, config: Mapping[str, object]) -> "Phase7Settings":
        section = config.get("assembly")
        if not isinstance(section, Mapping):
            raise ConfigurationError("Phase 7 requires configuration mapping assembly.")
        try:
            return cls(
                str(section["backend"]),
                str(section["executable"]),
                AssemblyResources(
                    threads=section["threads"],  # type: ignore[arg-type]
                    memory_mb=section["memory_mb"],  # type: ignore[arg-type]
                    runtime_minutes=section["runtime_minutes"],  # type: ignore[arg-type]
                ),
                section["overwrite"],  # type: ignore[arg-type]
            )
        except (KeyError, TypeError, ValueError) as error:
            raise ConfigurationError(f"Invalid Phase 7 settings: {error}") from error

    def to_dict(self) -> dict[str, object]:
        return {
            "backend": self.backend,
            "executable": self.executable,
            "resources": self.resources.to_dict(),
            "overwrite": self.overwrite,
        }


@dataclass(frozen=True, slots=True)
class Phase7SampleResult:
    sample_id: str
    assembly: AssemblyResult

    def __post_init__(self) -> None:
        if self.sample_id != self.assembly.request.sample_id:
            raise InputValidationError("Phase 7 sample/result IDs differ.")

    def to_dict(self) -> dict[str, object]:
        return {"sample_id": self.sample_id, "assembly": self.assembly.to_dict()}


@dataclass(frozen=True, slots=True)
class Phase7RunReport:
    context: AnalysisContext
    settings: Phase7Settings
    sample_results: tuple[Phase7SampleResult, ...]
    dry_run: bool
    created_at: str
    hifivar_version: str = __version__
    schema_version: str = PHASE7_REPORT_SCHEMA_VERSION

    def __post_init__(self) -> None:
        expected = tuple(record.sample.sample_id for record in self.context.samples)
        observed = tuple(item.sample_id for item in self.sample_results)
        if expected != observed:
            raise InputValidationError("Phase 7 result order differs from AnalysisContext.")
        planned = all(
            item.assembly.status is AssemblyResultStatus.PLANNED
            for item in self.sample_results
        )
        if planned != self.dry_run:
            raise InputValidationError("Phase 7 report status conflicts with dry_run.")

    @property
    def status(self) -> Phase7RunStatus:
        return Phase7RunStatus.PLANNED if self.dry_run else Phase7RunStatus.COMPLETED

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
        standardized = standardize_data(payload, context="Phase 7 report value")
        if not isinstance(standardized, dict):
            raise InputValidationError("Phase 7 report serialization failed.")
        return standardized

    def write_json(self, path: Path, *, overwrite: bool = False) -> Path:
        return write_json_atomic(
            self.to_dict(), path, overwrite=overwrite, artifact_name="Phase 7 report"
        )

    def write_yaml(self, path: Path, *, overwrite: bool = False) -> Path:
        return write_yaml_atomic(
            self.to_dict(), path, overwrite=overwrite, artifact_name="Phase 7 report"
        )


def run_phase7(
    context: AnalysisContext,
    *,
    output_directory: str | Path,
    work_directory: str | Path,
    config: Mapping[str, object],
    dry_run: bool = False,
    wrapper: HifiasmWrapper | None = None,
) -> Phase7RunReport:
    """Assemble every FASTQ sample in context order without reference coupling."""
    section = config.get("assembly")
    if not isinstance(section, Mapping) or section.get("enabled") is not True:
        raise ConfigurationError("Phase 7 requires assembly.enabled: true.")
    settings = Phase7Settings.from_config(config)
    output_root = Path(output_directory).expanduser()
    work_root = Path(work_directory).expanduser()
    engine = wrapper or HifiasmWrapper(executable=settings.executable)
    results: list[Phase7SampleResult] = []
    for record in context.samples:
        sample = record.sample
        if sample.input.input_type is not InputType.FASTQ:
            raise InputValidationError(
                f"Phase 7 sample '{sample.sample_id}' is not applicable: "
                "primary input must be HiFi FASTQ; BAM/CRAM extraction is disabled."
            )
        request = AssemblyRequest(
            sample,
            work_root / sample.sample_id / f"{sample.sample_id}.asm",
            output_root / sample.sample_id,
            settings.resources,
            settings.overwrite,
        )
        result = engine.run(
            request,
            dry_run=dry_run,
            stderr_path=work_root / sample.sample_id / f"{sample.sample_id}.hifiasm.log",
        )
        results.append(Phase7SampleResult(sample.sample_id, result))
    return Phase7RunReport(
        context,
        settings,
        tuple(results),
        dry_run,
        utc_now_iso8601(),
    )


__all__ = [
    "PHASE7_REPORT_SCHEMA_VERSION",
    "Phase7RunReport",
    "Phase7RunStatus",
    "Phase7SampleResult",
    "Phase7Settings",
    "run_phase7",
]
