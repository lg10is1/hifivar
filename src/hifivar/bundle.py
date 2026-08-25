"""Safe Phase 14 artifact and reproducibility bundle creation."""

from __future__ import annotations

import hashlib
import shutil
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Mapping

from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.report import FinalRunReport
from hifivar.serialization import redact_sensitive_data, utc_now_iso8601, write_json_atomic, write_yaml_atomic

_LARGE_SUFFIXES = (".bam", ".cram", ".sam", ".fastq", ".fastq.gz", ".fq", ".fq.gz", ".gfa", ".fasta", ".fasta.gz", ".fa", ".fa.gz")


@dataclass(frozen=True, slots=True)
class BundleItem:
    source: Path
    destination: PurePosixPath
    role: str
    checksum: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "source", Path(self.source))
        destination = PurePosixPath(self.destination)
        if destination.is_absolute() or ".." in destination.parts or not destination.parts:
            raise InputValidationError("Bundle destination must be a safe relative path.")
        object.__setattr__(self, "destination", destination)
        if not self.role.strip():
            raise InputValidationError("Bundle item role must be non-empty.")


@dataclass(frozen=True, slots=True)
class ReproducibilityRecord:
    software_versions: Mapping[str, object] = field(default_factory=dict)
    commands: tuple[tuple[str, ...], ...] = ()
    environment: Mapping[str, object] = field(default_factory=dict)
    reference: Mapping[str, object] = field(default_factory=dict)
    redact_values: tuple[str, ...] = ()
    sample_sheet: Path | None = None


@dataclass(frozen=True, slots=True)
class BundleResult:
    root: Path
    manifest: Path
    copied: tuple[Path, ...]
    pointers: tuple[Path, ...]


def create_release_bundle(
    report: FinalRunReport,
    destination: Path,
    *,
    items: tuple[BundleItem, ...] = (),
    reproducibility: ReproducibilityRecord | None = None,
    include_large: bool = False,
    overwrite: bool = False,
) -> BundleResult:
    """Create a directory bundle; large primary data remain pointers by default."""
    root = Path(destination)
    if root.exists() and any(root.iterdir()):
        if not overwrite:
            raise OutputValidationError(f"Release bundle directory is not empty: '{root}'.")
        if not root.joinpath("manifests", "bundle_manifest.json").is_file():
            raise OutputValidationError(
                f"Refusing overwrite because directory is not an owned release bundle: '{root}'."
            )
    root.mkdir(parents=True, exist_ok=True)
    report.write_json(root / "manifests" / "final_report.json", overwrite=overwrite)
    report.write_yaml(root / "manifests" / "final_report.yaml", overwrite=overwrite)
    report.write_markdown(root / "reports" / "final_report.md", overwrite=overwrite)
    report.write_html(root / "reports" / "final_report.html", overwrite=overwrite)
    copied: list[Path] = []
    pointers: list[Path] = []
    entries: list[dict[str, object]] = []
    seen: set[PurePosixPath] = set()
    for item in items:
        if item.destination in seen:
            raise InputValidationError(f"Duplicate bundle destination: '{item.destination}'.")
        seen.add(item.destination)
        if not item.source.is_file():
            raise InputValidationError(f"Bundle source does not exist: '{item.source}'.")
        is_large = _is_large_primary(item.source)
        target = root.joinpath(*item.destination.parts)
        if is_large and not include_large:
            pointer = target.with_suffix(target.suffix + ".pointer.json")
            write_json_atomic({"role": item.role, "source_path": str(item.source.absolute()), "size_bytes": item.source.stat().st_size, "copied": False}, pointer, overwrite=overwrite, artifact_name="Artifact pointer")
            pointers.append(pointer)
            entries.append({"role": item.role, "bundle_path": str(pointer.relative_to(root)), "copied": False, "sha256": None})
            continue
        if target.exists() and not overwrite:
            raise OutputValidationError(f"Bundle artifact already exists: '{target}'.")
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(item.source, target)
        checksum = _sha256(target) if item.checksum else None
        copied.append(target)
        entries.append({"role": item.role, "bundle_path": str(target.relative_to(root)), "copied": True, "size_bytes": target.stat().st_size, "sha256": checksum})
    repro = reproducibility or ReproducibilityRecord()
    _write_reproducibility(root, report, repro, overwrite=overwrite)
    manifest = root / "manifests" / "bundle_manifest.json"
    write_json_atomic({"schema_version": "1.0", "created_at": utc_now_iso8601(), "run_id": report.run_id, "report_status": report.status.value, "items": entries, "large_primary_data_included": include_large}, manifest, overwrite=overwrite, artifact_name="Bundle manifest")
    return BundleResult(root, manifest, tuple(copied), tuple(pointers))


def selected_report_items(report: FinalRunReport) -> tuple[BundleItem, ...]:
    """Map only explicitly selected report artifacts into stable result paths."""
    items=[]
    for track in report.tracks:
        for artifact in track.artifacts:
            if not artifact.selected_for_bundle:
                continue
            sample = artifact.sample_id or "run"
            destination = PurePosixPath("results", track.category, sample, artifact.path.name)
            items.append(BundleItem(artifact.path, destination, artifact.role, checksum=artifact.sha256 is not None))
    return tuple(items)


def _write_reproducibility(root: Path, report: FinalRunReport, record: ReproducibilityRecord, *, overwrite: bool) -> None:
    config = redact_sensitive_data(report.config)
    write_yaml_atomic(config, root / "configs" / "effective_config.yaml", overwrite=overwrite, artifact_name="Effective config")
    write_json_atomic(redact_sensitive_data(record.software_versions), root / "provenance" / "software_versions.json", overwrite=overwrite, artifact_name="Software versions")
    commands=[["***" if str(arg) in record.redact_values else str(arg) for arg in command] for command in record.commands]
    write_json_atomic(commands, root / "provenance" / "commands.json", overwrite=overwrite, artifact_name="Commands")
    write_json_atomic(redact_sensitive_data(record.environment), root / "provenance" / "environment.json", overwrite=overwrite, artifact_name="Environment summary")
    write_json_atomic(record.reference or report.reference, root / "provenance" / "reference.json", overwrite=overwrite, artifact_name="Reference metadata")
    git_path=root / "provenance" / "git_sha.txt"; git_path.parent.mkdir(parents=True,exist_ok=True)
    if git_path.exists() and not overwrite: raise OutputValidationError(f"Git provenance output already exists: '{git_path}'.")
    git_path.write_text((report.git_sha or "NOT_RECORDED") + "\n",encoding="utf-8",newline="\n")
    if record.sample_sheet is not None:
        source=Path(record.sample_sheet)
        if not source.is_file(): raise InputValidationError(f"Sample sheet does not exist: '{source}'.")
        target=root / "provenance" / "sample_sheet.tsv"
        if target.exists() and not overwrite: raise OutputValidationError(f"Sample-sheet bundle output exists: '{target}'.")
        target.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(source,target)


def _is_large_primary(path: Path) -> bool:
    return any(path.name.lower().endswith(suffix) for suffix in _LARGE_SUFFIXES)


def _sha256(path: Path) -> str:
    digest=hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""): digest.update(block)
    return digest.hexdigest()


__all__ = ["BundleItem", "BundleResult", "ReproducibilityRecord", "create_release_bundle", "selected_report_items"]
