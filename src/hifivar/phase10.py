"""Phase 10 orchestration for optional IGV/manual variant review."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Sequence

from hifivar.exceptions import InputValidationError
from hifivar.igv import IgvRunResult, IgvRunStatus, IgvWrapper
from hifivar.review import ReviewManifest, ReviewResult, ReviewTarget


class Phase10Status(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class Phase10RunReport:
    status: Phase10Status
    igv: IgvRunResult
    manifest: ReviewManifest
    output_directory: Path

    def to_dict(self) -> dict[str, object]:
        return {
            "phase": "10",
            "status": self.status.value,
            "output_directory": str(self.output_directory),
            "scientific_policy": {
                "manual_status_is_truth": False,
                "manual_status_is_pathogenicity": False,
                "selection_is_explicit": True,
                "raw_variant_artifacts_modified": False,
            },
            "igv": self.igv.to_dict(),
            "manifest": self.manifest.to_dict(),
        }


def run_phase10(
    targets: Sequence[ReviewTarget],
    *,
    output_directory: str | Path,
    igv_wrapper: IgvWrapper | None = None,
    dry_run: bool = False,
    overwrite: bool = False,
) -> Phase10RunReport:
    """Generate evidence for an explicit target list without changing callers."""
    output = Path(output_directory).expanduser()
    ordered = tuple(targets)
    if any(target.output_directory != output for target in ordered):
        raise InputValidationError(
            "Every Phase 10 target output_directory must match the run output directory."
        )
    wrapper = igv_wrapper or IgvWrapper()
    plan = wrapper.plan(
        ordered,
        batch_path=output / "review.igv.batch",
        snapshot_directory=output / "screenshots",
    )
    igv_result = wrapper.run(
        plan,
        dry_run=dry_run,
        overwrite=overwrite,
        stdout_path=output / "logs" / "igv.stdout.log" if not dry_run else None,
        stderr_path=output / "logs" / "igv.stderr.log" if not dry_run else None,
    )
    results = tuple(
        ReviewResult(target=evidence.target, evidence=evidence)
        for evidence in igv_result.evidence
    )
    manifest = ReviewManifest(results)
    status = (
        Phase10Status.PLANNED
        if igv_result.status is IgvRunStatus.PLANNED
        else Phase10Status.COMPLETED
    )
    report = Phase10RunReport(status, igv_result, manifest, output)
    if not dry_run:
        manifest.write_json(output / "review_manifest.json", overwrite=overwrite)
        manifest.write_yaml(output / "review_manifest.yaml", overwrite=overwrite)
        manifest.write_tsv(output / "review_manifest.tsv", overwrite=overwrite)
    return report


__all__ = ["Phase10RunReport", "Phase10Status", "run_phase10"]
