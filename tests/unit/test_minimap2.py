"""Tests for the Minimap2Wrapper.

R2-MAJOR-003: Tests atomic output replacement with .partial paths.
R2-MAJOR-004: Tests comprehensive validation after completion.
R2-MAJOR-006: Tests structured provenance for each step.
R2-MAJOR-007: Tests separate logs, unique intermediates.
R2-MAJOR-008: Includes non-dry-run fake-executable E2E test.
"""

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
from hifivar.minimap2 import (
    AlignmentChainProvenance,
    Minimap2Options,
    Minimap2Wrapper,
    StepProvenance,
)
from hifivar.reference import ReferenceGenome
from hifivar.sample import InputDataset, Sample


class FakeCommandRunner:
    """Small deterministic CommandRunner test double."""

    def __init__(self) -> None:
        self.commands: list[tuple[tuple[str, ...], dict[str, object]]] = []
        self.require_calls: list[str] = []
        self.version_output = "minimap2 2.28-r1209\n"
        self.samtools_version_output = "samtools 1.20\n"
        self.create_output = True
        self.fail_alignment = False
        self.fail_sort = False
        self.fail_index = False
        self.fail_quickcheck = False
        self.available = True

    def require_executable(self, executable: str) -> Path:
        self.require_calls.append(executable)
        if not self.available:
            raise ToolNotFoundError(f"Required executable was not found: {executable}")
        return Path("/opt/minimap2")

    def run(self, command: object, **kwargs: object) -> CommandResult:
        args = tuple(str(arg) for arg in command)
        self.commands.append((args, dict(kwargs)))
        dry_run = kwargs.get("dry_run", False) is True

        # Handle version detection
        if args[1:] == ("--version",) and args[0] == "minimap2":
            return CommandResult(
                args=args, returncode=0, stdout=self.version_output,
                stderr="", duration_seconds=0.01, cwd=None, executed=True,
            )
        if args[1:] == ("--version",) and args[0] == "samtools":
            return CommandResult(
                args=args, returncode=0, stdout=self.samtools_version_output,
                stderr="", duration_seconds=0.01, cwd=None, executed=True,
            )

        # Handle samtools view -H for header parsing
        if args[0] == "samtools" and len(args) >= 3 and args[1] == "view" and args[2] == "-H":
            header = (
                "@HD\tVN:1.6\tSO:coordinate\n"
                "@SQ\tSN:chr1\tLN:1000\n"
                "@RG\tID:RG1\tSM:S1\tPL:PACBIO\n"
            )
            return CommandResult(
                args=args, returncode=0, stdout=header,
                stderr="", duration_seconds=0.01, cwd=None, executed=True,
            )

        # Handle samtools quickcheck
        if args[0] == "samtools" and len(args) >= 2 and args[1] == "quickcheck":
            if self.fail_quickcheck:
                return CommandResult(
                    args=args, returncode=1, stdout="",
                    stderr="fail", duration_seconds=0.01, cwd=None, executed=True,
                )
            return CommandResult(
                args=args, returncode=0, stdout="",
                stderr="", duration_seconds=0.01, cwd=None, executed=True,
            )

        # Handle samtools index
        if args[0] == "samtools" and len(args) >= 2 and args[1] == "index":
            if self.fail_index:
                raise CommandExecutionError("Command failed with return code 1: samtools index")
            if self.create_output and len(args) >= 4:
                bai_path = Path(args[3])
                bai_path.parent.mkdir(parents=True, exist_ok=True)
                bai_path.write_bytes(b"fake-bai")
            return CommandResult(
                args=args, returncode=0, stdout="",
                stderr="", duration_seconds=0.01, cwd=None, executed=True,
            )

        if dry_run:
            return CommandResult(
                args=args, returncode=None, stdout=None, stderr=None,
                duration_seconds=0.0, cwd=None, executed=False,
                stderr_path=(
                    Path(kwargs["stderr_path"])
                    if kwargs.get("stderr_path") is not None
                    else None
                ),
            )

        if self.fail_alignment:
            raise CommandExecutionError("Command failed with return code 1: minimap2")

        # Handle minimap2 alignment command
        if args[0] == "minimap2" and "-a" in args:
            rg_idx = args.index("-R")
            ref_idx = rg_idx + 2
            fastq_args = list(args[ref_idx + 1:])
            if not fastq_args or all(f == "-a" or f == "-x" or f.startswith("-") for f in fastq_args):
                raise ValueError("No FASTQ paths found in minimap2 command")

            if self.create_output:
                stdout_path = kwargs.get("stdout_path")
                if stdout_path:
                    out = Path(str(stdout_path))
                    out.parent.mkdir(parents=True, exist_ok=True)
                    out.write_bytes(b"fake-sam-output")

            return CommandResult(
                args=args, returncode=0, stdout="",
                stderr="minimap2 metrics", duration_seconds=12.5,
                cwd=None, executed=True,
            )

        # Handle samtools sort
        if args[0] == "samtools" and len(args) >= 2 and args[1] == "sort":
            if self.fail_sort:
                raise CommandExecutionError("Command failed with return code 1: samtools sort")
            if self.create_output:
                for i, arg in enumerate(args):
                    if arg == "-o" and i + 1 < len(args):
                        output = Path(args[i + 1])
                        output.parent.mkdir(parents=True, exist_ok=True)
                        output.write_bytes(b"fake-coordinate-sorted-bam")
                        break

            return CommandResult(
                args=args, returncode=0, stdout="",
                stderr="", duration_seconds=10.0,
                cwd=None, executed=True,
            )

        # Default success
        return CommandResult(
            args=args, returncode=0, stdout="",
            stderr="", duration_seconds=0.01,
            cwd=None, executed=True,
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
    tool: AlignmentTool = AlignmentTool.MINIMAP2,
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


def make_wrapper(runner: FakeCommandRunner | None = None) -> Minimap2Wrapper:
    return Minimap2Wrapper(
        runner=runner or FakeCommandRunner(),
        options=Minimap2Options(preset="map-hifi", output_format="bam", log_level="INFO"),
    )


class TestMinimap2Options:
    def test_options_normalize_and_load_from_config(self) -> None:
        options = Minimap2Options.from_config(
            {
                "alignment": {
                    "minimap2_preset": "ASM20",
                    "minimap2_output_format": "CRAM",
                    "minimap2_log_level": "debug",
                }
            }
        )
        assert options.to_dict() == {
            "preset": "asm20",
            "output_format": "cram",
            "log_level": "DEBUG",
            "extra_flags": [],
        }

    def test_options_with_extra_flags(self) -> None:
        options = Minimap2Options(
            preset="map-hifi",
            output_format="bam",
            extra_flags=("--eqx", "--MD"),
        )
        assert options.extra_flags == ("--eqx", "--MD")

    def test_options_from_config_with_string_flags(self) -> None:
        options = Minimap2Options.from_config(
            {
                "alignment": {
                    "minimap2_extra_flags": "--eqx --MD",
                }
            }
        )
        assert options.extra_flags == ("--eqx", "--MD")

    @pytest.mark.parametrize(
        "kwargs",
        (
            {"preset": "SUBREAD"},
            {"preset": 1},
            {"output_format": "sam"},
            {"output_format": None},
            {"log_level": "VERBOSE"},
            {"log_level": None},
            {"extra_flags": "not-a-tuple"},
        ),
    )
    def test_options_reject_invalid_settings(
        self, kwargs: dict[str, object]
    ) -> None:
        with pytest.raises(ConfigurationError):
            Minimap2Options(**kwargs)  # type: ignore[arg-type]


class TestMinimap2Wrapper:
    def test_build_command_matches_official_hifi_shape(self, tmp_path: Path) -> None:
        request = make_request(tmp_path)
        command = make_wrapper().build_command(request)

        assert command[0] == "minimap2"
        assert "-a" in command
        assert "-x" in command
        assert command[command.index("-x") + 1] == "map-hifi"
        assert "-t" in command
        assert command[command.index("-t") + 1] == "24"
        assert "--eqx" in command
        assert "-Y" in command
        assert "--MD" in command
        assert "-R" in command
        rg_idx = command.index("-R")
        read_group = command[rg_idx + 1]
        assert read_group == "@RG\tID:S1\tSM:S1\tPL:PACBIO"
        assert str(request.reference.fasta) in command

    def test_build_command_includes_fastq_paths(self, tmp_path: Path) -> None:
        fastq1 = write_fastq(tmp_path / "reads", "S1_R1.fastq")
        fastq2 = write_fastq(tmp_path / "reads", "S1_R2.fastq")
        request = make_request(tmp_path, files=(fastq1, fastq2))
        command = make_wrapper().build_command(request)

        fastq_in_command = [
            arg for arg in command
            if arg.endswith(".fastq") or arg.endswith(".fq")
        ]
        assert len(fastq_in_command) == 2
        assert str(fastq1) in command
        assert str(fastq2) in command

    def test_build_command_single_fastq(self, tmp_path: Path) -> None:
        request = make_request(tmp_path)
        command = make_wrapper().build_command(request)

        fastq_in_command = [
            arg for arg in command
            if arg.endswith(".fastq") or arg.endswith(".fq")
        ]
        assert len(fastq_in_command) == 1

    def test_build_samtools_sort_command(self, tmp_path: Path) -> None:
        request = make_request(tmp_path)
        wrapper = make_wrapper()
        samtools_cmd = wrapper.build_samtools_sort_command(
            request, input_mode="file", input_path=tmp_path / "intermediate.sam"
        )

        assert samtools_cmd[0] == "samtools"
        assert "sort" in samtools_cmd
        assert "-@" in samtools_cmd
        assert "-o" in samtools_cmd
        assert str(request.output_path) in samtools_cmd
        assert str(tmp_path / "intermediate.sam") in samtools_cmd

    def test_build_samtools_sort_command_requires_input_path_for_file_mode(
        self, tmp_path: Path
    ) -> None:
        request = make_request(tmp_path)
        wrapper = make_wrapper()
        with pytest.raises(InputValidationError, match="requires input_path"):
            wrapper.build_samtools_sort_command(request, input_mode="file")

    def test_build_samtools_sort_command_pipe_mode(self, tmp_path: Path) -> None:
        request = make_request(tmp_path)
        wrapper = make_wrapper()
        samtools_cmd = wrapper.build_samtools_sort_command(request, input_mode="pipe")

        assert samtools_cmd[-1] == "-"

    def test_detect_version_uses_executable_check_and_parses_output(
        self, tmp_path: Path
    ) -> None:
        runner = FakeCommandRunner()
        wrapper = make_wrapper(runner)
        assert wrapper.detect_version() == "2.28-r1209"
        assert runner.require_calls == ["minimap2"]
        assert runner.commands[0][0] == ("minimap2", "--version")

    def test_detect_version_rejects_unparseable_output(self) -> None:
        runner = FakeCommandRunner()
        runner.version_output = "unknown tool output"
        with pytest.raises(ToolVersionError, match="parse minimap2 version"):
            make_wrapper(runner).detect_version()

    def test_missing_minimap2_is_reported_by_command_runner_boundary(self) -> None:
        runner = FakeCommandRunner()
        runner.available = False
        with pytest.raises(ToolNotFoundError, match="not found"):
            make_wrapper(runner).detect_version()

    def test_dry_run_does_not_require_executable_or_create_outputs(
        self, tmp_path: Path
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

    def test_output_created_after_planning_is_not_silently_overwritten(
        self, tmp_path: Path
    ) -> None:
        runner = FakeCommandRunner()
        request = make_request(tmp_path)
        request.output_path.parent.mkdir(parents=True)
        request.output_path.write_bytes(b"race")

        with pytest.raises(OutputValidationError, match="already exists"):
            make_wrapper(runner).run(request)
        assert request.output_path.read_bytes() == b"race"

    def test_explicit_overwrite_replaces_only_requested_output(
        self, tmp_path: Path
    ) -> None:
        runner = FakeCommandRunner()
        request = make_request(tmp_path, overwrite=True)
        request.output_path.parent.mkdir(parents=True)
        request.output_path.write_bytes(b"old-alignment")

        result = make_wrapper(runner).run(request, dry_run=True)

        assert result.status is AlignmentResultStatus.PLANNED
        assert request.output_path.read_bytes() == b"old-alignment"

    def test_wrapper_rejects_non_minimap2_requests(self, tmp_path: Path) -> None:
        with pytest.raises(InputValidationError, match="cannot handle"):
            make_wrapper().build_command(
                make_request(tmp_path / "wrong-tool", tool=AlignmentTool.PBMM2)
            )

    def test_wrapper_rejects_invalid_output_format(self, tmp_path: Path) -> None:
        from hifivar.alignment import AlignmentOutputFormat
        with pytest.raises(ConfigurationError, match="must be one of"):
            Minimap2Options(output_format="sam")

    def test_run_forwards_redactions_timeout_and_tool_log_path(
        self, tmp_path: Path
    ) -> None:
        runner = FakeCommandRunner()
        request = make_request(tmp_path)
        log_path = tmp_path / "logs" / "minimap2.log"
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

    def test_wrapper_rejects_empty_executable(self) -> None:
        with pytest.raises(InputValidationError, match="non-empty string"):
            Minimap2Wrapper(executable="")

    def test_wrapper_rejects_whitespace_executable(self) -> None:
        with pytest.raises(InputValidationError, match="non-empty string"):
            Minimap2Wrapper(executable="   ")

    def test_wrapper_accepts_path_executable(self) -> None:
        wrapper = Minimap2Wrapper(executable=Path("/usr/bin/minimap2"))
        assert "minimap2" in wrapper.executable

    def test_plan_command_returns_alignment_command_plan(self, tmp_path: Path) -> None:
        request = make_request(tmp_path)
        plan = make_wrapper().plan_command(request)
        assert plan.tool is AlignmentTool.MINIMAP2
        assert plan.args[0] == "minimap2"

    def test_request_is_frozen(self, tmp_path: Path) -> None:
        from dataclasses import FrozenInstanceError
        request = make_request(tmp_path)
        with pytest.raises(FrozenInstanceError):
            request.threads = 2  # type: ignore[misc]

    def test_cram_output_format_accepted(self, tmp_path: Path) -> None:
        request = make_request(tmp_path, output_format=AlignmentOutputFormat.CRAM)
        command = make_wrapper().build_command(request)
        assert command[0] == "minimap2"

    def test_different_presets_accepted(self, tmp_path: Path) -> None:
        for preset in ("map-hifi", "asm20", "asm10", "asm5"):
            wrapper = Minimap2Wrapper(
                options=Minimap2Options(preset=preset, output_format="bam")
            )
            request = make_request(tmp_path)
            command = wrapper.build_command(request)
            assert command[command.index("-x") + 1] == preset

    def test_wrapper_rejects_empty_samtools_executable(self) -> None:
        with pytest.raises(InputValidationError, match="non-empty string"):
            Minimap2Wrapper(samtools_executable="")

    def test_wrapper_rejects_whitespace_samtools_executable(self) -> None:
        with pytest.raises(InputValidationError, match="non-empty string"):
            Minimap2Wrapper(samtools_executable="   ")

    def test_detect_samtools_version(self, tmp_path: Path) -> None:
        runner = FakeCommandRunner()
        wrapper = make_wrapper(runner)
        assert wrapper.detect_samtools_version() == "1.20"


class TestMinimap2WrapperE2E:
    """R2-MAJOR-008: Non-dry-run fake-executable E2E test.

    Tests the full minimap2 -> samtools sort -> samtools index chain.
    """

    def test_e2e_alignment_chain_succeeds(self, tmp_path: Path) -> None:
        runner = FakeCommandRunner()
        request = make_request(tmp_path)
        wrapper = make_wrapper(runner)

        result = wrapper.run(request)

        assert result.status is AlignmentResultStatus.COMPLETED
        assert result.executed is True
        assert result.tool_version is not None
        assert "minimap2=" in result.tool_version
        assert "samtools=" in result.tool_version

        # Verify the output file exists
        assert request.output_path.exists()
        assert request.output_path.stat().st_size > 0

        # Verify the BAI index exists
        bai_path = Path(f"{request.output_path}.bai")
        assert bai_path.exists()

        # Verify provenance commands were called
        command_names = [cmd[0][0] for cmd in runner.commands]
        assert "minimap2" in command_names
        assert "samtools" in command_names

    def test_e2e_alignment_chain_minimap2_failure(self, tmp_path: Path) -> None:
        runner = FakeCommandRunner()
        runner.fail_alignment = True
        request = make_request(tmp_path)
        wrapper = make_wrapper(runner)

        with pytest.raises(CommandExecutionError):
            wrapper.run(request)

        # Output should not exist
        assert not request.output_path.exists()

    def test_e2e_alignment_chain_sort_failure(self, tmp_path: Path) -> None:
        runner = FakeCommandRunner()
        runner.fail_sort = True
        request = make_request(tmp_path)
        wrapper = make_wrapper(runner)

        with pytest.raises(CommandExecutionError):
            wrapper.run(request)

        # Output should not exist
        assert not request.output_path.exists()

    def test_e2e_alignment_chain_index_failure(self, tmp_path: Path) -> None:
        runner = FakeCommandRunner()
        runner.fail_index = True
        request = make_request(tmp_path)
        wrapper = make_wrapper(runner)

        with pytest.raises(CommandExecutionError):
            wrapper.run(request)

        # Output should not exist (atomic replacement failed)
        assert not request.output_path.exists()

    def test_e2e_alignment_chain_quickcheck_failure(self, tmp_path: Path) -> None:
        runner = FakeCommandRunner()
        runner.fail_quickcheck = True
        request = make_request(tmp_path)
        wrapper = make_wrapper(runner)

        with pytest.raises(OutputValidationError, match="quickcheck"):
            wrapper.run(request)

        # Output should not exist (validation failed)
        assert not request.output_path.exists()

    def test_e2e_alignment_creates_unique_intermediates(self, tmp_path: Path) -> None:
        """R2-MAJOR-007: Intermediate paths are unique and task-owned."""
        runner = FakeCommandRunner()
        request = make_request(tmp_path)
        wrapper = make_wrapper(runner)

        result = wrapper.run(request)

        assert result.status is AlignmentResultStatus.COMPLETED
        # No intermediate files should remain
        intermediate_files = list(tmp_path.glob(".S1_*.intermediate.*"))
        assert len(intermediate_files) == 0

    def test_e2e_alignment_creates_separate_logs(self, tmp_path: Path) -> None:
        """R2-MAJOR-007: Separate stderr logs per step."""
        runner = FakeCommandRunner()
        request = make_request(tmp_path)
        wrapper = make_wrapper(runner)

        result = wrapper.run(request)

        assert result.status is AlignmentResultStatus.COMPLETED
        # Check that separate log files were created
        logs_dir = request.output_path.parent / ".logs"
        assert logs_dir.exists()

    def test_e2e_alignment_atomic_replacement(self, tmp_path: Path) -> None:
        """R2-MAJOR-003: Atomic replacement preserves existing final on failure."""
        runner = FakeCommandRunner()
        request = make_request(tmp_path, overwrite=True)
        request.output_path.parent.mkdir(parents=True)
        request.output_path.write_bytes(b"existing-good-bam")
        bai_path = Path(f"{request.output_path}.bai")
        bai_path.write_bytes(b"existing-good-bai")

        # Now make sort fail
        runner.fail_sort = True
        wrapper = make_wrapper(runner)

        with pytest.raises(CommandExecutionError):
            wrapper.run(request)

        # Existing final artifacts should be preserved
        assert request.output_path.read_bytes() == b"existing-good-bam"
        assert bai_path.read_bytes() == b"existing-good-bai"

    def test_e2e_alignment_overwrite_replaces_old_output(self, tmp_path: Path) -> None:
        """R2-MAJOR-003: Overwrite mode replaces old output."""
        runner = FakeCommandRunner()
        request = make_request(tmp_path, overwrite=True)
        request.output_path.parent.mkdir(parents=True)
        request.output_path.write_bytes(b"old-alignment")
        bai_path = Path(f"{request.output_path}.bai")
        bai_path.write_bytes(b"old-bai")

        wrapper = make_wrapper(runner)
        result = wrapper.run(request)

        assert result.status is AlignmentResultStatus.COMPLETED
        assert request.output_path.read_bytes() != b"old-alignment"
        assert bai_path.read_bytes() != b"old-bai"


class TestStepProvenance:
    def test_step_provenance_serialization(self) -> None:
        prov = StepProvenance(
            step_name="minimap2",
            command=("minimap2", "-a", "-x", "map-hifi"),
            returncode=0,
            tool="minimap2",
            tool_version="2.28-r1209",
            duration_seconds=12.5,
            stdout_path=Path("/tmp/output.sam"),
            stderr_path=Path("/tmp/minimap2.stderr"),
        )
        payload = prov.to_dict()
        assert payload["step_name"] == "minimap2"
        assert payload["returncode"] == 0
        assert payload["tool_version"] == "2.28-r1209"
        assert payload["duration_seconds"] == 12.5

    def test_alignment_chain_provenance_serialization(self) -> None:
        minimap2 = StepProvenance(
            step_name="minimap2", command=("minimap2",), returncode=0,
            tool="minimap2", tool_version="2.28", duration_seconds=1.0,
            stdout_path=None, stderr_path=None,
        )
        sort = StepProvenance(
            step_name="samtools_sort", command=("samtools", "sort"), returncode=0,
            tool="samtools", tool_version="1.20", duration_seconds=2.0,
            stdout_path=None, stderr_path=None,
        )
        index = StepProvenance(
            step_name="samtools_index", command=("samtools", "index"), returncode=0,
            tool="samtools", tool_version="1.20", duration_seconds=0.5,
            stdout_path=None, stderr_path=None,
        )
        chain = AlignmentChainProvenance(
            sample_id="S1",
            minimap2=minimap2,
            samtools_sort=sort,
            samtools_index=index,
        )
        payload = chain.to_dict()
        assert payload["sample_id"] == "S1"
        assert "minimap2" in payload["chain"]
        assert payload["minimap2"]["tool"] == "minimap2"
        assert payload["samtools_sort"]["tool"] == "samtools"
        assert payload["samtools_index"]["tool"] == "samtools"
