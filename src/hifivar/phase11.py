"""Phase 11 annotation and optional functional-prioritization orchestration."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from hifivar import __version__
from hifivar.annotation import AnnotationInput, AnnotationResult, RegionOverlapResult
from hifivar.annovar import AnnovarRequest, AnnovarWrapper
from hifivar.exceptions import InputValidationError
from hifivar.functional import (
    FunctionalBackend,
    FunctionalPrioritizationRequest,
    FunctionalPrioritizationResult,
    run_functional_prioritization,
)
from hifivar.serialization import standardize_data, utc_now_iso8601, write_json_atomic, write_yaml_atomic
from hifivar.vep import VepRequest, VepWrapper


PHASE11_REPORT_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True, slots=True)
class AnnotationJob:
    """Independent adapters configured for one immutable source VCF."""

    input: AnnotationInput
    annovar: AnnovarRequest | None = None
    vep: VepRequest | None = None

    def __post_init__(self) -> None:
        if self.annovar is None and self.vep is None:
            raise InputValidationError("Annotation job requires ANNOVAR and/or VEP.")
        if self.annovar is not None and self.annovar.input != self.input:
            raise InputValidationError("ANNOVAR request does not match annotation job input.")
        if self.vep is not None and self.vep.input != self.input:
            raise InputValidationError("VEP request does not match annotation job input.")


@dataclass(frozen=True, slots=True)
class Phase11RunReport:
    annotation_results: tuple[AnnotationResult, ...]
    region_overlaps: tuple[RegionOverlapResult, ...]
    functional_result: FunctionalPrioritizationResult | None
    dry_run: bool
    created_at: str
    hifivar_version: str = __version__
    schema_version: str = PHASE11_REPORT_SCHEMA_VERSION

    def to_dict(self) -> dict[str, object]:
        payload = {
            "schema_version": self.schema_version,
            "hifivar_version": self.hifivar_version,
            "created_at": self.created_at,
            "status": "planned" if self.dry_run else "completed",
            "dry_run": self.dry_run,
            "annotation_results": [item.to_dict() for item in self.annotation_results],
            "region_overlaps": [item.to_dict() for item in self.region_overlaps],
            "functional_result": self.functional_result.to_dict() if self.functional_result else None,
            "scientific_policy": {
                "variant_call_confidence_is_functional_impact": False,
                "functional_impact_is_variant_call_confidence": False,
                "raw_caller_outputs_modified": False,
                "variant_classes_remain_separate": True,
                "functional_selection_is_explicit": True,
            },
        }
        standardized = standardize_data(payload, context="Phase 11 report value")
        if not isinstance(standardized, dict):
            raise InputValidationError("Phase 11 report serialization failed.")
        return standardized

    def write_json(self, path: str | Path, *, overwrite: bool = False) -> Path:
        return write_json_atomic(
            self.to_dict(), path, overwrite=overwrite,
            artifact_name="Phase 11 report",
        )

    def write_yaml(self, path: str | Path, *, overwrite: bool = False) -> Path:
        return write_yaml_atomic(
            self.to_dict(), path, overwrite=overwrite,
            artifact_name="Phase 11 report",
        )


def run_phase11(
    jobs: Sequence[AnnotationJob],
    *,
    annovar_wrapper: AnnovarWrapper | None = None,
    vep_wrapper: VepWrapper | None = None,
    region_overlaps: Sequence[RegionOverlapResult] = (),
    functional_request: FunctionalPrioritizationRequest | None = None,
    functional_backend: FunctionalBackend | None = None,
    log_directory: str | Path | None = None,
    dry_run: bool = False,
) -> Phase11RunReport:
    """Run independent annotations, then only an explicit functional selection."""
    ordered_jobs = tuple(jobs)
    if not ordered_jobs and not region_overlaps:
        raise InputValidationError("Phase 11 requires a tool annotation or region-overlap result.")
    keys = [(job.input.sample_id, job.input.variant_category.value) for job in ordered_jobs]
    if len(keys) != len(set(keys)):
        raise InputValidationError("Phase 11 jobs duplicate a sample/category input.")
    if (functional_request is None) != (functional_backend is None):
        raise InputValidationError(
            "Functional request and backend must either both be supplied or both be absent."
        )

    annovar = annovar_wrapper or AnnovarWrapper()
    vep = vep_wrapper or VepWrapper()
    logs = Path(log_directory).expanduser() if log_directory is not None else None
    results: list[AnnotationResult] = []
    for job in ordered_jobs:
        stem = f"{job.input.sample_id}.{job.input.variant_category.value}"
        if job.annovar is not None:
            results.append(
                annovar.run(
                    job.annovar,
                    dry_run=dry_run,
                    stdout_path=logs / f"{stem}.annovar.stdout.log" if logs and not dry_run else None,
                    stderr_path=logs / f"{stem}.annovar.stderr.log" if logs and not dry_run else None,
                )
            )
        if job.vep is not None:
            results.append(
                vep.run(
                    job.vep,
                    dry_run=dry_run,
                    stdout_path=logs / f"{stem}.vep.stdout.log" if logs and not dry_run else None,
                    stderr_path=logs / f"{stem}.vep.stderr.log" if logs and not dry_run else None,
                )
            )

    functional = None
    if functional_request is not None and functional_backend is not None:
        if dry_run:
            functional = None
        else:
            functional = run_functional_prioritization(
                functional_request,
                backend=functional_backend,
            )
    return Phase11RunReport(
        tuple(results), tuple(region_overlaps), functional, dry_run,
        utc_now_iso8601(),
    )


__all__ = [
    "AnnotationJob", "PHASE11_REPORT_SCHEMA_VERSION", "Phase11RunReport",
    "run_phase11",
]
