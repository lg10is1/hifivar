"""Phase 14 final run reports with machine and offline human renderings."""

from __future__ import annotations

import html
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Mapping

from hifivar import __version__
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.serialization import (
    redact_sensitive_data,
    standardize_data,
    utc_now_iso8601,
    write_json_atomic,
    write_yaml_atomic,
)

FINAL_REPORT_SCHEMA_VERSION = "1.0"


class FinalStatus(str, Enum):
    COMPLETE = "COMPLETE"
    PARTIAL = "PARTIAL"
    FAILED = "FAILED"
    NOT_RUN = "NOT_RUN"
    DISABLED = "DISABLED"


@dataclass(frozen=True, slots=True)
class ReportArtifact:
    path: Path
    role: str
    artifact_type: str
    sample_id: str | None = None
    sha256: str | None = None
    selected_for_bundle: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "path", Path(self.path))
        for name in ("role", "artifact_type"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise InputValidationError(f"Report artifact {name} must be non-empty.")

    def to_dict(self) -> dict[str, object]:
        return {
            "path": str(self.path), "role": self.role,
            "artifact_type": self.artifact_type, "sample_id": self.sample_id,
            "sha256": self.sha256, "selected_for_bundle": self.selected_for_bundle,
        }


@dataclass(frozen=True, slots=True)
class ToolRecord:
    name: str
    version: str | None
    executable_or_image: str | None = None
    backend: str | None = None
    validation_status: str | None = None

    def to_dict(self) -> dict[str, object]:
        return {"name": self.name, "version": self.version,
                "executable_or_image": self.executable_or_image,
                "backend": self.backend, "validation_status": self.validation_status}


@dataclass(frozen=True, slots=True)
class TrackReport:
    name: str
    phase: str
    category: str
    status: FinalStatus
    enabled: bool
    artifacts: tuple[ReportArtifact, ...] = ()
    samples: tuple[str, ...] = ()
    qc: Mapping[str, object] = field(default_factory=dict)
    benchmark: Mapping[str, object] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()
    message: str | None = None

    def __post_init__(self) -> None:
        if not self.name.strip() or not self.phase.strip() or not self.category.strip():
            raise InputValidationError("Report track name, phase, and category must be non-empty.")
        if self.enabled and self.status is FinalStatus.DISABLED:
            raise InputValidationError("An enabled report track cannot have DISABLED status.")
        if not self.enabled and self.status is not FinalStatus.DISABLED:
            raise InputValidationError("A disabled report track must have DISABLED status.")

    def to_dict(self) -> dict[str, object]:
        return {"name": self.name, "phase": self.phase, "category": self.category,
                "status": self.status.value, "enabled": self.enabled,
                "samples": list(self.samples),
                "artifacts": [item.to_dict() for item in self.artifacts],
                "qc": dict(self.qc), "benchmark": dict(self.benchmark),
                "warnings": list(self.warnings), "message": self.message}


@dataclass(frozen=True, slots=True)
class FinalRunReport:
    run_id: str
    git_sha: str | None
    reference: Mapping[str, object]
    samples: tuple[Mapping[str, object], ...]
    tracks: tuple[TrackReport, ...]
    config: Mapping[str, object]
    tools: tuple[ToolRecord, ...] = ()
    cohort: Mapping[str, object] | None = None
    provenance: Mapping[str, object] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()
    created_at: str = field(default_factory=utc_now_iso8601)
    schema_version: str = FINAL_REPORT_SCHEMA_VERSION
    project_version: str = __version__

    def __post_init__(self) -> None:
        if not isinstance(self.run_id, str) or not self.run_id.strip():
            raise InputValidationError("Final report run_id must be non-empty.")
        names = [track.name for track in self.tracks]
        if len(names) != len(set(names)):
            raise InputValidationError("Final report track names must be unique.")

    @property
    def status(self) -> FinalStatus:
        states = {track.status for track in self.tracks}
        if FinalStatus.FAILED in states:
            return FinalStatus.FAILED
        if FinalStatus.PARTIAL in states:
            return FinalStatus.PARTIAL
        if FinalStatus.COMPLETE in states and FinalStatus.NOT_RUN in states:
            return FinalStatus.PARTIAL
        if FinalStatus.COMPLETE in states:
            return FinalStatus.COMPLETE
        if FinalStatus.NOT_RUN in states:
            return FinalStatus.NOT_RUN
        return FinalStatus.DISABLED

    @property
    def enabled_phases(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys(track.phase for track in self.tracks if track.enabled))

    def to_dict(self) -> dict[str, object]:
        payload = {
            "schema_version": self.schema_version, "project_version": self.project_version,
            "run_id": self.run_id, "created_at": self.created_at,
            "git_sha": self.git_sha, "status": self.status.value,
            "enabled_phases": list(self.enabled_phases),
            "reference": dict(self.reference), "samples": [dict(item) for item in self.samples],
            "cohort": dict(self.cohort) if self.cohort is not None else None,
            "config": redact_sensitive_data(self.config),
            "tools": [tool.to_dict() for tool in self.tools],
            "tracks": [track.to_dict() for track in self.tracks],
            "warnings": list(self.warnings), "provenance": dict(self.provenance),
            "interpretation_policy": {
                "clinical_interpretation_performed": False,
                "benchmark_is_not_pathogenicity": True,
                "manual_review_is_not_biological_truth": True,
            },
        }
        standardized = standardize_data(redact_sensitive_data(payload), context="Final report value")
        if not isinstance(standardized, dict):
            raise InputValidationError("Final report serialization failed.")
        return standardized

    def write_json(self, path: Path, *, overwrite: bool = False) -> Path:
        return write_json_atomic(self.to_dict(), path, overwrite=overwrite, artifact_name="Final report")

    def write_yaml(self, path: Path, *, overwrite: bool = False) -> Path:
        return write_yaml_atomic(self.to_dict(), path, overwrite=overwrite, artifact_name="Final report")

    def write_markdown(self, path: Path, *, overwrite: bool = False) -> Path:
        return _write_text(path, self._markdown(), overwrite=overwrite, name="Final report")

    def write_html(self, path: Path, *, overwrite: bool = False) -> Path:
        return _write_text(path, self._html(), overwrite=overwrite, name="Final report")

    def _markdown(self) -> str:
        lines = [f"# HiFiVar run report: {self.run_id}", "",
                 f"- Status: **{self.status.value}**", f"- HiFiVar: `{self.project_version}`",
                 f"- Git SHA: `{self.git_sha or 'NOT_RECORDED'}`", f"- Created: `{self.created_at}`", "",
                 "> Research-use report. Status, benchmark, annotation, and manual review do not constitute clinical interpretation.", "",
                 "## Summary", "", f"Samples: {len(self.samples)}; enabled phases: {', '.join(self.enabled_phases) or 'none'}.", ""]
        for category, title in _REPORT_SECTIONS:
            items = [track for track in self.tracks if track.category == category]
            lines.extend((f"## {title}", ""))
            if not items:
                lines.extend(("NOT_RUN", "")); continue
            lines.extend(("| Track | Phase | Status | Samples | Artifacts |", "|---|---|---|---:|---:|"))
            lines.extend(f"| {item.name} | {item.phase} | {item.status.value} | {len(item.samples)} | {len(item.artifacts)} |" for item in items)
            for item in items:
                if item.message: lines.append(f"- {item.name}: {item.message}")
                lines.extend(f"- WARNING ({item.name}): {warning}" for warning in item.warnings)
            lines.append("")
        lines.extend(("## Warnings and limitations", ""))
        lines.extend(f"- {warning}" for warning in self.warnings)
        if not self.warnings: lines.append("- None recorded.")
        return "\n".join(lines) + "\n"

    def _html(self) -> str:
        sections=[]
        for category,title in _REPORT_SECTIONS:
            items=[track for track in self.tracks if track.category==category]
            rows="".join(f"<tr><td>{html.escape(item.name)}</td><td>{html.escape(item.phase)}</td><td class='status {item.status.value.lower()}'>{item.status.value}</td><td>{len(item.samples)}</td><td>{len(item.artifacts)}</td></tr>" for item in items)
            body=f"<table><thead><tr><th>Track</th><th>Phase</th><th>Status</th><th>Samples</th><th>Artifacts</th></tr></thead><tbody>{rows}</tbody></table>" if rows else "<p>NOT_RUN</p>"
            sections.append(f"<section><h2>{html.escape(title)}</h2>{body}</section>")
        warning_items="".join(f"<li>{html.escape(item)}</li>" for item in self.warnings) or "<li>None recorded.</li>"
        return "<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>HiFiVar report</title><style>body{font-family:system-ui,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#17202a}table{border-collapse:collapse;width:100%}th,td{border:1px solid #ccd1d1;padding:.45rem;text-align:left}.status{font-weight:700}.complete{color:#18743b}.partial,.not_run{color:#936600}.failed{color:#b42318}.disabled{color:#667085}.notice{padding:1rem;background:#f2f4f7;border-left:4px solid #667085}</style></head><body>" + f"<h1>HiFiVar run report: {html.escape(self.run_id)}</h1><p><strong>Status:</strong> {self.status.value}<br><strong>Version:</strong> {html.escape(self.project_version)}<br><strong>Git SHA:</strong> {html.escape(self.git_sha or 'NOT_RECORDED')}<br><strong>Created:</strong> {html.escape(self.created_at)}</p><p class='notice'>Research use only. No clinical interpretation is performed.</p><h2>Summary</h2><p>Samples: {len(self.samples)}; enabled phases: {html.escape(', '.join(self.enabled_phases) or 'none')}.</p>" + "".join(sections) + f"<section><h2>Warnings and limitations</h2><ul>{warning_items}</ul></section></body></html>"


_REPORT_SECTIONS = (
    ("sample", "Sample/cohort summary"), ("small", "Small variants"),
    ("sv", "Structural variants"), ("tr", "Tandem repeats"),
    ("phasing", "Phasing"), ("assembly", "Assembly"),
    ("review", "Manual review"), ("annotation", "Annotation"),
    ("cohort", "Cohort"), ("benchmark", "Benchmark"),
)


def _write_text(path: Path, content: str, *, overwrite: bool, name: str) -> Path:
    path = Path(path)
    if path.exists() and not overwrite:
        raise OutputValidationError(f"{name} output already exists: '{path}'.")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8", newline="\n")
    return path


__all__ = ["FINAL_REPORT_SCHEMA_VERSION", "FinalRunReport", "FinalStatus", "ReportArtifact", "ToolRecord", "TrackReport"]
