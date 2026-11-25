"""Phase 11 annotation and optional functional-prioritization orchestration."""

from __future__ import annotations

from dataclasses import dataclass
import logging
from pathlib import Path
from typing import Sequence

from hifivar import __version__
from hifivar.annotation import AnnotationInput, AnnotationResult, RegionOverlapResult
from hifivar.annovar import AnnovarRequest, AnnovarWrapper
from hifivar.exceptions import HiFiVarError, InputValidationError
from hifivar.functional import (
    FunctionalBackend,
    FunctionalPrioritizationRequest,
    FunctionalPrioritizationResult,
    run_functional_prioritization,
)
from hifivar.serialization import standardize_data, utc_now_iso8601, write_json_atomic, write_yaml_atomic
from hifivar.vep import VepRequest, VepWrapper


PHASE11_REPORT_SCHEMA_VERSION = "1.1"
_LOGGER = logging.getLogger(__name__)


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
class AnnotationBackendFailure:
    """An expected failure retained separately from successful annotations."""

    sample_id: str
    variant_category: str
    backend: str
    error_type: str
    message: str

    def to_dict(self) -> dict[str, object]:
        return {
            "sample": self.sample_id,
            "variant_category": self.variant_category,
            "backend": self.backend,
            "status": "failed",
            "error_type": self.error_type,
            "message": self.message,
        }


@dataclass(frozen=True, slots=True)
class Phase11RunReport:
    annotation_results: tuple[AnnotationResult, ...]
    region_overlaps: tuple[RegionOverlapResult, ...]
    functional_result: FunctionalPrioritizationResult | None
    dry_run: bool
    created_at: str
    hifivar_version: str = __version__
    schema_version: str = PHASE11_REPORT_SCHEMA_VERSION
    backend_failures: tuple[AnnotationBackendFailure, ...] = ()
    isolation_enabled: bool = False
    functional_requested: bool = False

    def to_dict(self) -> dict[str, object]:
        if self.backend_failures:
            status = "partial" if self.annotation_results or self.region_overlaps else "failed"
        else:
            status = "planned" if self.dry_run else "completed"
        if self.functional_result is not None:
            functional_status = "completed"
        elif not self.functional_requested:
            functional_status = "disabled"
        elif self.backend_failures:
            functional_status = "not_run_annotation_failure"
        else:
            functional_status = "planned" if self.dry_run else "not_run"
        payload = {
            "schema_version": self.schema_version,
            "hifivar_version": self.hifivar_version,
            "created_at": self.created_at,
            "status": status,
            "dry_run": self.dry_run,
            "annotation_results": [item.to_dict() for item in self.annotation_results],
            "region_overlaps": [item.to_dict() for item in self.region_overlaps],
            "functional_result": self.functional_result.to_dict() if self.functional_result else None,
            "functional_status": functional_status,
            "failure_isolation_enabled": self.isolation_enabled,
            "backend_failures": [item.to_dict() for item in self.backend_failures],
            "tool_annotation_scope": "whole_source_vcf",
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
    isolate_failures: bool = False,
) -> Phase11RunReport:
    """Run independent annotations, then only an explicit functional selection.

    Default execution remains fail-fast. With explicit ``isolate_failures``,
    expected HiFiVar errors are recorded per backend and later independent
    annotations continue. Programming errors still propagate. An incomplete
    annotation attempt never triggers functional prediction automatically.
    This API does not change the Snakemake scheduler's failure policy.
    """
    if not isinstance(isolate_failures, bool):
        raise InputValidationError("isolate_failures must be boolean.")
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
    failures: list[AnnotationBackendFailure] = []
    for job in ordered_jobs:
        stem = f"{job.input.sample_id}.{job.input.variant_category.value}"
        for backend, wrapper, request in (
            ("annovar", annovar, job.annovar), ("vep", vep, job.vep),
        ):
            if request is None:
                continue
            try:
                result = wrapper.run(
                    request, dry_run=dry_run,
                    stdout_path=logs / f"{stem}.{backend}.stdout.log" if logs and not dry_run else None,
                    stderr_path=logs / f"{stem}.{backend}.stderr.log" if logs and not dry_run else None,
                )
            except HiFiVarError as exc:
                if not isolate_failures:
                    raise
                failures.append(AnnotationBackendFailure(
                    job.input.sample_id, job.input.variant_category.value,
                    backend, type(exc).__name__, str(exc),
                ))
                _LOGGER.error("Annotation backend failed: sample=%s category=%s backend=%s error=%s",
                              job.input.sample_id, job.input.variant_category.value,
                              backend, exc)
            else:
                results.append(result)

    functional = None
    if functional_request is not None and functional_backend is not None and not failures:
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
        backend_failures=tuple(failures), isolation_enabled=isolate_failures,
        functional_requested=functional_request is not None,
    )


__all__ = [
    "AnnotationBackendFailure", "AnnotationJob", "PHASE11_REPORT_SCHEMA_VERSION", "Phase11RunReport",
    "run_phase11",
]
