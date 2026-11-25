"""Safe Phase 14 artifact and reproducibility bundle creation."""

from __future__ import annotations

import hashlib
import shutil
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath, PureWindowsPath
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
        if (
            destination.is_absolute()
            or ".." in destination.parts
            or not destination.parts
            or "\\" in str(destination)
            or any(":" in part for part in destination.parts)
            or any(
                part.endswith((".", " ")) or PureWindowsPath(part).is_reserved()
                or any(ord(character) < 32 or character in '<>"|?*' for character in part)
                for part in destination.parts
            )
        ):
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
    repro = reproducibility or ReproducibilityRecord()
    # Validate the complete write set before publishing even the first report.
    # In particular, pointer names are outputs, not the requested BAM names.
    _preflight_bundle(root, report, items, repro, include_large, overwrite)
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
    _write_reproducibility(root, report, repro, overwrite=overwrite)
    manifest = root / "manifests" / "bundle_manifest.json"
    write_json_atomic({"schema_version": "1.0", "created_at": utc_now_iso8601(), "run_id": report.run_id, "report_status": report.status.value, "items": entries, "large_primary_data_included": include_large}, manifest, overwrite=overwrite, artifact_name="Bundle manifest")
    return BundleResult(root, manifest, tuple(copied), tuple(pointers))


_GENERATED_PATHS = (
    "manifests/final_report.json", "manifests/final_report.yaml",
    "reports/final_report.md", "reports/final_report.html",
    "configs/effective_config.yaml", "provenance/software_versions.json",
    "provenance/commands.json", "provenance/environment.json",
    "provenance/reference.json", "provenance/git_sha.txt",
    "manifests/bundle_manifest.json",
)
_RESERVED_NAMESPACES = frozenset({"manifests", "reports", "configs", "provenance"})


def _preflight_bundle(
    root: Path,
    report: FinalRunReport,
    items: tuple[BundleItem, ...],
    record: ReproducibilityRecord,
    include_large: bool,
    overwrite: bool,
) -> None:
    """Reject static collisions/aliases without mutating a bundle or its sources.

    The bundle is a single-writer directory, not a defence against another
    process concurrently replacing its paths. I/O failure during publication
    is still reported as a failure, not as a successful transactional bundle.
    """
    if root.is_symlink() or getattr(root, "is_junction", lambda: False)():
        raise OutputValidationError(f"Release bundle root must not be a link: '{root}'.")
    if root.exists() and not root.is_dir():
        raise OutputValidationError(f"Release bundle root is not a directory: '{root}'.")
    outputs = [root / path for path in _GENERATED_PATHS]
    required_sources = [item.source for item in items]
    if record.sample_sheet is not None:
        required_sources.append(Path(record.sample_sheet))
        outputs.append(root / "provenance" / "sample_sheet.tsv")
    for source in required_sources:
        if not source.is_file():
            raise InputValidationError(f"Bundle source does not exist: '{source}'.")
    # Unselected/failed tracks may legitimately reference absent files; protect
    # their paths without requiring them to be deliverable artifacts.
    sources = required_sources + [
        artifact.path for track in report.tracks for artifact in track.artifacts
    ]
    for item in items:
        if item.destination.parts[0].casefold() in _RESERVED_NAMESPACES:
            raise InputValidationError(
                f"Bundle destination uses a reserved metadata namespace: '{item.destination}'."
            )
        target = root.joinpath(*item.destination.parts)
        if _is_large_primary(item.source) and not include_large:
            target = target.with_suffix(target.suffix + ".pointer.json")
        outputs.append(target)
    root_resolved = root.resolve()
    resolved_outputs: list[Path] = []
    for output in outputs:
        resolved = output.resolve()
        if not resolved.is_relative_to(root_resolved) or resolved == root_resolved:
            raise InputValidationError(f"Bundle output escapes its root: '{output}'.")
        for component in (output, *output.parents):
            if component == root.parent:
                break
            if component.is_symlink() or getattr(component, "is_junction", lambda: False)():
                raise OutputValidationError(f"Bundle output contains a link: '{component}'.")
            if component.exists() and (component == output) != component.is_file():
                raise OutputValidationError(f"Bundle file/directory conflict: '{component}'.")
        if output.exists() and not overwrite:
            raise OutputValidationError(f"Bundle output already exists: '{output}'.")
        for previous, previous_resolved in zip(outputs[:len(resolved_outputs)], resolved_outputs):
            if (
                resolved == previous_resolved
                or resolved.is_relative_to(previous_resolved)
                or previous_resolved.is_relative_to(resolved)
                or _same_existing_file(output, previous)
            ):
                raise InputValidationError(f"Conflicting bundle outputs: '{previous}' and '{output}'.")
        for source in sources:
            source_resolved = source.resolve()
            if (
                resolved == source_resolved
                or source_resolved.is_relative_to(resolved)
                or resolved.is_relative_to(source_resolved)
                or _same_existing_file(output, source)
            ):
                raise InputValidationError(f"Bundle output aliases a protected input: '{output}'.")
        resolved_outputs.append(resolved)


def _same_existing_file(first: Path, second: Path) -> bool:
    return first.exists() and second.exists() and first.samefile(second)


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
