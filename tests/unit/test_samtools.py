"""Tests for the minimal Phase 2.4 samtools indexing wrapper."""

from __future__ import annotations

from pathlib import Path

import pytest

from hifivar.alignment import AlignmentOutputFormat, AlignmentTool
from hifivar.alignment_postprocess import (
    AlignmentArtifact,
    AlignmentIndexFormat,
    AlignmentIndexRequest,
    AlignmentSortOrder,
    AlignmentSource,
)
from hifivar.command import CommandResult
from hifivar.exceptions import (
    CommandExecutionError,
    OutputValidationError,
    ToolNotFoundError,
    ToolVersionError,
)
from hifivar.reference import ReferenceGenome
from hifivar.samtools import IndexResultStatus, SamtoolsWrapper


class FakeCommandRunner:
    """Deterministic samtools CommandRunner test double."""

    def __init__(self) -> None:
        self.commands: list[tuple[tuple[str, ...], dict[str, object]]] = []
        self.require_calls: list[str] = []
        self.version_output = "samtools 1.22.1\nUsing htslib 1.22.1\n"
        self.create_output = True
        self.fail_index = False
        self.available = True

    def require_executable(self, executable: str) -> Path:
        self.require_calls.append(executable)
        if not self.available:
            raise ToolNotFoundError(f"Required executable was not found: {executable}")
        return Path("/opt/samtools")

    def run(self, command: object, **kwargs: object) -> CommandResult:
        args = tuple(str(arg) for arg in command)  # type: ignore[union-attr]
        self.commands.append((args, dict(kwargs)))
        if args[1:] == ("--version",):
            return CommandResult(
                args=args,
                returncode=0,
                stdout=self.version_output,
                stderr="",
                duration_seconds=0.01,
                cwd=None,
                executed=True,
            )
        if kwargs.get("dry_run") is True:
            return CommandResult(
                args=args,
                returncode=None,
                stdout=None,
                stderr=None,
                duration_seconds=0.0,
                cwd=None,
                executed=False,
            )
        if self.fail_index:
            raise CommandExecutionError("Command failed with return code 7: samtools")
        if self.create_output:
            output = Path(args[-1])
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_bytes(b"fake-index")
        return CommandResult(
            args=args,
            returncode=0,
            stdout="",
            stderr="",
            duration_seconds=1.25,
            cwd=None,
            executed=True,
        )


def write_reference(root: Path) -> ReferenceGenome:
    root.mkdir(parents=True, exist_ok=True)
    fasta = root / "reference.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    return ReferenceGenome.from_fasta(fasta, build="GRCh38")


def make_artifact(tmp_path: Path, *, cram: bool = False) -> AlignmentArtifact:
    path = tmp_path / ("S1.cram" if cram else "S1.bam")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"coordinate-sorted-alignment")
    return AlignmentArtifact(
        sample_id="S1",
        path=path,
        output_format=(
            AlignmentOutputFormat.CRAM if cram else AlignmentOutputFormat.BAM
        ),
        reference=write_reference(tmp_path / "reference"),
        source=AlignmentSource.GENERATED,
        sort_order=AlignmentSortOrder.COORDINATE,
        tool=AlignmentTool.PBMM2,
        tool_version="1.17.0",
    )


def make_wrapper(runner: FakeCommandRunner | None = None) -> SamtoolsWrapper:
    return SamtoolsWrapper(runner=runner or FakeCommandRunner())  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "index_format,flag,suffix",
    (
        (AlignmentIndexFormat.BAI, "-b", ".bai"),
        (AlignmentIndexFormat.CSI, "-c", ".csi"),
    ),
)
def test_build_bam_index_commands_are_explicit(
    tmp_path: Path,
    index_format: AlignmentIndexFormat,
    flag: str,
    suffix: str,
) -> None:
    request = AlignmentIndexRequest.create(
        make_artifact(tmp_path),
        index_format=index_format,
        threads=6,
    )
    command = make_wrapper().build_index_command(request)
    assert command[:4] == ["samtools", "index", "-@", "6"]
    assert flag in command
    assert command[-2] == str(request.artifact.path)
    assert command[-1].endswith(suffix)


def test_cram_index_command_uses_crai_default_without_bam_flag(
    tmp_path: Path,
) -> None:
    request = AlignmentIndexRequest.create(make_artifact(tmp_path, cram=True))
    command = make_wrapper().build_index_command(request)
    assert request.index_format is AlignmentIndexFormat.CRAI
    assert "-b" not in command and "-c" not in command
    assert command[-1].endswith(".crai")


def test_detect_version_and_missing_or_invalid_tool_output() -> None:
    runner = FakeCommandRunner()
    assert make_wrapper(runner).detect_version() == "1.22.1"
    assert runner.require_calls == ["samtools"]

    runner.version_output = "not samtools"
    with pytest.raises(ToolVersionError, match="parse samtools version"):
        make_wrapper(runner).detect_version()

    runner.available = False
    with pytest.raises(ToolNotFoundError, match="not found"):
        make_wrapper(runner).detect_version()


def test_dry_run_needs_no_samtools_and_creates_no_index(tmp_path: Path) -> None:
    runner = FakeCommandRunner()
    runner.available = False
    request = AlignmentIndexRequest.create(make_artifact(tmp_path), threads=4)
    result = make_wrapper(runner).run_index(request, dry_run=True)
    assert result.status is IndexResultStatus.PLANNED
    assert result.executed is False
    assert result.artifact.index_path is None
    assert runner.require_calls == []
    assert not request.output_path.exists()


def test_real_index_run_validates_and_attaches_index(tmp_path: Path) -> None:
    runner = FakeCommandRunner()
    request = AlignmentIndexRequest.create(make_artifact(tmp_path), threads=4)
    result = make_wrapper(runner).run_index(request)
    assert result.status is IndexResultStatus.COMPLETED
    assert result.executed is True
    assert result.tool_version == "1.22.1"
    assert result.duration_seconds == 1.25
    assert result.artifact.index_path == request.output_path
    assert request.output_path.read_bytes() == b"fake-index"


def test_missing_expected_index_is_output_validation_error(tmp_path: Path) -> None:
    runner = FakeCommandRunner()
    runner.create_output = False
    request = AlignmentIndexRequest.create(make_artifact(tmp_path))
    with pytest.raises(OutputValidationError, match="missing"):
        make_wrapper(runner).run_index(request)


def test_external_index_failure_propagates(tmp_path: Path) -> None:
    runner = FakeCommandRunner()
    runner.fail_index = True
    request = AlignmentIndexRequest.create(make_artifact(tmp_path))
    with pytest.raises(CommandExecutionError, match="return code 7"):
        make_wrapper(runner).run_index(request)


def test_index_output_race_is_not_silently_overwritten(tmp_path: Path) -> None:
    runner = FakeCommandRunner()
    request = AlignmentIndexRequest.create(make_artifact(tmp_path))
    request.output_path.write_bytes(b"race")
    with pytest.raises(OutputValidationError, match="already exists"):
        make_wrapper(runner).run_index(request)
    assert request.output_path.read_bytes() == b"race"


def test_explicit_overwrite_replaces_only_requested_index(tmp_path: Path) -> None:
    runner = FakeCommandRunner()
    request = AlignmentIndexRequest.create(make_artifact(tmp_path), overwrite=True)
    request.output_path.write_bytes(b"old-index")

    make_wrapper(runner).run_index(request)

    assert request.output_path.read_bytes() == b"fake-index"


def test_index_run_revalidates_alignment_path(tmp_path: Path) -> None:
    runner = FakeCommandRunner()
    request = AlignmentIndexRequest.create(make_artifact(tmp_path))
    request.artifact.path.unlink()
    with pytest.raises(OutputValidationError, match="missing"):
        make_wrapper(runner).run_index(request)
    assert runner.commands == []
