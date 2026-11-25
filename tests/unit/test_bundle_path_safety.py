"""Bundle writes must not escape, corrupt inputs, or replace own metadata."""

from __future__ import annotations

import json
import os
from pathlib import Path, PurePosixPath

import pytest

from hifivar.bundle import BundleItem, ReproducibilityRecord, create_release_bundle
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.report import FinalRunReport, FinalStatus, ReportArtifact, TrackReport


def _report() -> FinalRunReport:
    return FinalRunReport("safe", "abc123", {}, (), (), {})


def _source(tmp_path: Path) -> Path:
    source = tmp_path / "source.tsv"
    source.write_bytes(b"protected input\n")
    return source


def _item(source: Path, destination: str) -> BundleItem:
    return BundleItem(source, PurePosixPath(destination), "table")


@pytest.mark.parametrize("destination", [r"..\outside.txt", r"a\..\outside.txt", "C:/outside.txt", "C:outside.txt", "//host/share/x", "results/NUL", "results/x.", "results/x ", "results/a\x00b", "results/a\x01b", "results/a\nb", "results/a\tb"])
def test_reject_nonportable_destinations(tmp_path: Path, destination: str) -> None:
    with pytest.raises(InputValidationError):
        _item(_source(tmp_path), destination)


@pytest.mark.parametrize("destination", ["manifests/final_report.json", "manifests/bundle_manifest.json", "reports/final_report.html", "configs/effective_config.yaml", "provenance/git_sha.txt"])
def test_reserve_generated_namespaces_before_writes(tmp_path: Path, destination: str) -> None:
    root = tmp_path / "bundle"
    with pytest.raises(InputValidationError):
        create_release_bundle(_report(), root, items=(_item(_source(tmp_path), destination),), overwrite=True)
    assert not root.exists()


def test_pointer_and_regular_item_collision_preflight(tmp_path: Path) -> None:
    bam = tmp_path / "input.bam"
    bam.write_bytes(b"input bam")
    root = tmp_path / "bundle"
    with pytest.raises(InputValidationError):
        create_release_bundle(_report(), root, items=(
            _item(bam, "results/a.bam"),
            _item(_source(tmp_path), "results/a.bam.pointer.json"),
        ), overwrite=True)
    assert not root.exists()
    assert bam.read_bytes() == b"input bam"


def test_missing_source_and_sample_sheet_do_not_publish_report(tmp_path: Path) -> None:
    root = tmp_path / "bundle"
    with pytest.raises(InputValidationError):
        create_release_bundle(_report(), root, items=(_item(tmp_path / "missing", "results/x"),))
    assert not root.exists()
    with pytest.raises(InputValidationError):
        create_release_bundle(_report(), root, reproducibility=ReproducibilityRecord(sample_sheet=tmp_path / "missing"))
    assert not root.exists()


def test_output_parent_collision_preflight(tmp_path: Path) -> None:
    source = _source(tmp_path)
    root = tmp_path / "bundle"
    with pytest.raises(InputValidationError):
        create_release_bundle(_report(), root, items=(
            _item(source, "results/x"), _item(source, "results/x/y"),
        ))
    assert not root.exists()


def test_existing_generated_file_hardlink_cannot_modify_source(tmp_path: Path) -> None:
    source = _source(tmp_path)
    root = tmp_path / "bundle"
    create_release_bundle(_report(), root)
    linked = root / "provenance" / "git_sha.txt"
    linked.unlink()  # Replace only this test's own previously generated file.
    os.link(source, linked)
    before = {str(p.relative_to(root)): p.read_bytes() for p in root.rglob("*") if p.is_file()}
    with pytest.raises((InputValidationError, OutputValidationError)):
        create_release_bundle(_report(), root, items=(_item(source, "results/x"),), overwrite=True)
    assert source.read_bytes() == b"protected input\n"
    assert before == {str(p.relative_to(root)): p.read_bytes() for p in root.rglob("*") if p.is_file()}


def test_target_hardlink_to_other_source_rejected_before_writes(tmp_path: Path) -> None:
    source = _source(tmp_path)
    other = tmp_path / "other.tsv"
    other.write_bytes(b"other protected input\n")
    root = tmp_path / "bundle"
    create_release_bundle(_report(), root)
    (root / "results").mkdir()
    os.link(other, root / "results" / "x.tsv")
    with pytest.raises((InputValidationError, OutputValidationError)):
        create_release_bundle(_report(), root, items=(
            _item(source, "results/x.tsv"), _item(other, "results/y.tsv"),
        ), overwrite=True)
    assert other.read_bytes() == b"other protected input\n"


def _symlink(link: Path, target: Path, *, directory: bool = False) -> None:
    try:
        link.symlink_to(target, target_is_directory=directory)
    except OSError as error:
        pytest.skip(f"Symlink fixture unavailable: {error}")


def test_symlink_parent_escape_preserves_outside(tmp_path: Path) -> None:
    source = _source(tmp_path)
    root = tmp_path / "bundle"
    create_release_bundle(_report(), root)
    outside = tmp_path / "outside"
    outside.mkdir()
    _symlink(root / "results", outside, directory=True)
    with pytest.raises((InputValidationError, OutputValidationError)):
        create_release_bundle(_report(), root, items=(_item(source, "results/x.tsv"),), overwrite=True)
    assert not (outside / "x.tsv").exists()


def test_same_source_target_rejected_with_inputs_intact(tmp_path: Path) -> None:
    root = tmp_path / "bundle"
    create_release_bundle(_report(), root)
    (root / "results").mkdir()
    source = root / "results" / "x.tsv"
    source.write_bytes(b"protected input\n")
    with pytest.raises((InputValidationError, OutputValidationError)):
        create_release_bundle(_report(), root, items=(_item(source, "results/x.tsv"),), overwrite=True)
    assert source.read_bytes() == b"protected input\n"


def test_valid_owned_bundle_overwrite_still_works(tmp_path: Path) -> None:
    source = _source(tmp_path)
    root = tmp_path / "bundle"
    item = _item(source, "results/S1/x.tsv")
    create_release_bundle(_report(), root, items=(item,))
    result = create_release_bundle(_report(), root, items=(item,), overwrite=True)
    assert result.copied[0].read_bytes() == source.read_bytes()
    assert json.loads(result.manifest.read_text())["run_id"] == "safe"


def test_unselected_absent_report_artifact_is_not_required(tmp_path: Path) -> None:
    track = TrackReport("small", "3", "small", FinalStatus.NOT_RUN, True,
                        (ReportArtifact(tmp_path / "absent.vcf.gz", "vcf", "vcf"),))
    report = FinalRunReport("not-run", None, {}, (), (track,), {})
    result = create_release_bundle(report, tmp_path / "bundle")
    assert json.loads(result.manifest.read_text())["report_status"] == "NOT_RUN"


def test_missing_report_artifact_cannot_become_an_output_directory(tmp_path: Path) -> None:
    root = tmp_path / "bundle"
    absent = root / "results" / "future.vcf.gz"
    track = TrackReport("small", "3", "small", FinalStatus.NOT_RUN, True,
                        (ReportArtifact(absent, "vcf", "vcf"),))
    report = FinalRunReport("not-run", None, {}, (), (track,), {})
    with pytest.raises(InputValidationError, match="protected input"):
        create_release_bundle(report, root, items=(
            _item(_source(tmp_path), "results/future.vcf.gz/nested.tsv"),
        ))
    assert not root.exists()
