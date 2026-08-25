"""Tests for the Phase 2.3 pbmm2 CommandRunner wrapper."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from hifivar.alignment import (
    AlignmentOutputFormat,
    AlignmentRequest,
    AlignmentResources,
    AlignmentResultStatus,
    AlignmentTool,
)
from hifivar.command import CommandResult
from hifivar.exceptions import (
    CommandExecutionError,
    ConfigurationError,
    InputValidationError,
    OutputValidationError,
    ToolNotFoundError,
    ToolVersionError,
)
from hifivar.pbmm2 import Pbmm2Options, Pbmm2Wrapper
from hifivar.reference import ReferenceGenome
from hifivar.sample import InputDataset, Sample


class FakeCommandRunner:
    """Small deterministic CommandRunner test double."""

    def __init__(self) -> None:
        self.commands: list[tuple[tuple[str, ...], dict[str, object]]] = []
        self.require_calls: list[str] = []
        self.version_output = "pbmm2 1.17.0 (commit v1.17.0)\n"
        self.create_output = True
        self.fail_alignment = False
        self.available = True

    def require_executable(self, executable: str) -> Path:
        self.require_calls.append(executable)
        if not self.available:
            raise ToolNotFoundError(f"Required executable was not found: {executable}")
        return Path("/opt/pbmm2")

    def run(self, command: object, **kwargs: object) -> CommandResult:
        args = tuple(str(arg) for arg in command)  # type: ignore[union-attr]
        self.commands.append((args, dict(kwargs)))
        dry_run = kwargs.get("dry_run", False) is True
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
        if dry_run:
            return CommandResult(
                args=args,
                returncode=None,
                stdout=None,
                stderr=None,
                duration_seconds=0.0,
                cwd=None,
                executed=False,
                stderr_path=(
                    Path(kwargs["stderr_path"])
                    if kwargs.get("stderr_path") is not None
                    else None
                ),
            )
        if self.fail_alignment:
            raise CommandExecutionError("Command failed with return code 9: pbmm2")
        if self.create_output:
            output = Path(args[4])
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_bytes(b"fake-coordinate-sorted-bam")
        return CommandResult(
            args=args,
            returncode=0,
            stdout="",
            stderr="pbmm2 metrics",
            duration_seconds=12.5,
            cwd=None,
            executed=True,
        )


def write_reference(root: Path) -> ReferenceGenome:
    root.mkdir(parents=True, exist_ok=True)
    fasta = root / "reference.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    return ReferenceGenome.from_fasta(fasta, build="GRCh38")


def write_fastq(root: Path, name: str) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    path = root / name
    path.write_text("@read1\nACGT\n+\nIIII\n", encoding="utf-8")
    return path


def make_request(
    tmp_path: Path,
    *,
    files: tuple[Path, ...] | None = None,
    tool: AlignmentTool = AlignmentTool.PBMM2,
    output_format: AlignmentOutputFormat = AlignmentOutputFormat.BAM,
    overwrite: bool = False,
) -> AlignmentRequest:
    if files is None:
        files = (write_fastq(tmp_path / "reads", "S1.fastq"),)
    return AlignmentRequest(
        sample=Sample("S1", InputDataset.from_files(files)),
        reference=write_reference(tmp_path / "reference"),
        output_path=tmp_path / "results" / f"S1.aligned.{output_format.value}",
        tool=tool,
        output_format=output_format,
        resources=AlignmentResources(
            threads=24,
            memory_mb=64_000,
            runtime_minutes=720,
        ),
        overwrite=overwrite,
    )


def make_wrapper(runner: FakeCommandRunner | None = None) -> Pbmm2Wrapper:
    return Pbmm2Wrapper(
        runner=runner or FakeCommandRunner(),  # type: ignore[arg-type]
        options=Pbmm2Options(preset="CCS", log_level="INFO"),
    )


def test_pbmm2_options_normalize_and_load_from_config() -> None:
    options = Pbmm2Options.from_config(
        {
            "alignment": {
                "pbmm2_preset": "hifi",
                "pbmm2_log_level": "debug",
            }
        }
    )
    assert options.to_dict() == {"preset": "HIFI", "log_level": "DEBUG"}


@pytest.mark.parametrize(
    "kwargs",
    (
        {"preset": "SUBREAD"},
        {"preset": 1},
        {"log_level": "VERBOSE"},
        {"log_level": None},
    ),
)
def test_pbmm2_options_reject_non_hifi_settings(
    kwargs: dict[str, object],
) -> None:
    with pytest.raises(ConfigurationError):
        Pbmm2Options(**kwargs)  # type: ignore[arg-type]


def test_build_command_matches_official_sorted_hifi_shape(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    command = make_wrapper().build_command(request)

    assert command[:5] == [
        "pbmm2",
        "align",
        str(request.reference.fasta),
        str(request.input_paths[0]),
        str(request.output_path),
    ]
    assert command[command.index("--preset") + 1] == "CCS"
    assert "--sort" in command
    assert command[command.index("--bam-index") + 1] == "NONE"
    assert command[command.index("-j") + 1] == "24"
    read_group = command[command.index("--rg") + 1]
    assert read_group == "@RG\tID:S1\tSM:S1\tPL:PACBIO"


def test_multiple_fastqs_plan_deterministic_fofn_without_writing(
    tmp_path: Path,
) -> None:
    files = (
        write_fastq(tmp_path / "reads", "movie1.fastq"),
        write_fastq(tmp_path / "reads", "movie2.fastq"),
    )
    request = make_request(tmp_path, files=files)
    wrapper = make_wrapper()

    input_argument = wrapper.input_argument(request)
    command = wrapper.plan_command(request)

    assert input_argument.name == "S1.fastq.fofn"
    assert str(input_argument) in command.args
    assert not input_argument.exists()


def test_detect_version_uses_executable_check_and_parses_output(
    tmp_path: Path,
) -> None:
    runner = FakeCommandRunner()
    wrapper = make_wrapper(runner)
    assert wrapper.detect_version() == "1.17.0"
    assert runner.require_calls == ["pbmm2"]
    assert runner.commands[0][0] == ("pbmm2", "--version")


def test_detect_version_rejects_unparseable_output() -> None:
    runner = FakeCommandRunner()
    runner.version_output = "unknown tool output"
    with pytest.raises(ToolVersionError, match="parse pbmm2 version"):
        make_wrapper(runner).detect_version()


def test_missing_pbmm2_is_reported_by_command_runner_boundary() -> None:
    runner = FakeCommandRunner()
    runner.available = False
    with pytest.raises(ToolNotFoundError, match="not found"):
        make_wrapper(runner).detect_version()


def test_dry_run_does_not_require_executable_or_create_outputs(
    tmp_path: Path,
) -> None:
    runner = FakeCommandRunner()
    runner.available = False
    request = make_request(tmp_path)
    result = make_wrapper(runner).run(request, dry_run=True)

    assert result.status is AlignmentResultStatus.PLANNED
    assert result.executed is False
    assert result.tool_version is None
    assert runner.require_calls == []
    assert runner.commands[-1][1]["dry_run"] is True
    assert not request.output_path.exists()
    assert not request.output_path.parent.exists()


def test_real_run_validates_output_and_records_version_runtime(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    runner = FakeCommandRunner()
    request = make_request(tmp_path)
    with caplog.at_level(logging.INFO, logger="hifivar"):
        result = make_wrapper(runner).run(request)

    assert result.status is AlignmentResultStatus.COMPLETED
    assert result.executed is True
    assert result.tool_version == "1.17.0"
    assert result.duration_seconds == 12.5
    assert request.output_path.read_bytes() == b"fake-coordinate-sorted-bam"
    assert "sample=S1" in caplog.text
    assert "version=1.17.0" in caplog.text


def test_real_multi_fastq_run_writes_ordered_absolute_fofn(
    tmp_path: Path,
) -> None:
    files = (
        write_fastq(tmp_path / "reads", "movie2.fastq"),
        write_fastq(tmp_path / "reads", "movie1.fastq"),
    )
    request = make_request(tmp_path, files=files)
    wrapper = make_wrapper()
    wrapper.run(request)

    fofn = wrapper.input_argument(request)
    assert fofn.read_text(encoding="utf-8").splitlines() == [
        str(path.absolute()) for path in files
    ]


def test_missing_expected_bam_is_output_validation_error(tmp_path: Path) -> None:
    runner = FakeCommandRunner()
    runner.create_output = False
    with pytest.raises(OutputValidationError, match="missing"):
        make_wrapper(runner).run(make_request(tmp_path))


def test_external_failure_propagates_command_execution_error(tmp_path: Path) -> None:
    runner = FakeCommandRunner()
    runner.fail_alignment = True
    with pytest.raises(CommandExecutionError, match="return code 9"):
        make_wrapper(runner).run(make_request(tmp_path))


def test_disappeared_fastq_fails_before_external_execution(tmp_path: Path) -> None:
    runner = FakeCommandRunner()
    request = make_request(tmp_path)
    request.input_paths[0].unlink()

    with pytest.raises(InputValidationError, match="missing"):
        make_wrapper(runner).run(request)
    assert runner.commands == []


def test_wrapper_rejects_non_pbmm2_and_cram_requests(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="cannot handle"):
        make_wrapper().build_command(
            make_request(tmp_path / "wrong-tool", tool=AlignmentTool.MINIMAP2)
        )
    with pytest.raises(InputValidationError, match="must be BAM"):
        make_wrapper().build_command(
            make_request(
                tmp_path / "cram",
                output_format=AlignmentOutputFormat.CRAM,
            )
        )


def test_output_created_after_planning_is_not_silently_overwritten(
    tmp_path: Path,
) -> None:
    runner = FakeCommandRunner()
    request = make_request(tmp_path)
    request.output_path.parent.mkdir(parents=True)
    request.output_path.write_bytes(b"race")

    with pytest.raises(OutputValidationError, match="already exists"):
        make_wrapper(runner).run(request)
    assert request.output_path.read_bytes() == b"race"


def test_explicit_overwrite_replaces_only_requested_output(tmp_path: Path) -> None:
    runner = FakeCommandRunner()
    request = make_request(tmp_path, overwrite=True)
    request.output_path.parent.mkdir(parents=True)
    request.output_path.write_bytes(b"old-alignment")

    make_wrapper(runner).run(request)

    assert request.output_path.read_bytes() == b"fake-coordinate-sorted-bam"


def test_run_forwards_redactions_timeout_and_tool_log_path(tmp_path: Path) -> None:
    runner = FakeCommandRunner()
    request = make_request(tmp_path)
    log_path = tmp_path / "logs" / "pbmm2.log"
    make_wrapper(runner).run(
        request,
        dry_run=True,
        timeout=30,
        redact_values={"secret"},
        stderr_path=log_path,
    )

    kwargs = runner.commands[-1][1]
    assert kwargs["timeout"] == 30
    assert kwargs["redact_values"] == {"secret"}
    assert kwargs["stderr_path"] == log_path
