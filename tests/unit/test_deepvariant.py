"""Tests for the Phase 3.1 DeepVariant CommandRunner wrapper."""

from __future__ import annotations

from pathlib import Path
import struct
import zlib

import pytest

import hifivar.deepvariant as deepvariant_module

from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import (
    AlignmentArtifact,
    AlignmentSortOrder,
    AlignmentSource,
)
from hifivar.command import CommandResult
from hifivar.deepvariant import (
    DeepVariantExecutionMode,
    DeepVariantRuntime,
    DeepVariantWrapper,
)
from hifivar.exceptions import (
    CommandExecutionError,
    ConfigurationError,
    OutputValidationError,
    ToolNotFoundError,
    ToolVersionError,
)
from hifivar.reference import ReferenceGenome
from hifivar.small import DeepVariantRequest, SmallVariantResultStatus


def write_bgzf(path: Path, payload: bytes) -> None:
    compressor = zlib.compressobj(level=6, wbits=-15)
    compressed = compressor.compress(payload) + compressor.flush()
    total_size = 18 + len(compressed) + 8
    header = (
        b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
        + struct.pack("<H", total_size - 1)
    )
    footer = struct.pack("<II", zlib.crc32(payload), len(payload))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(header + compressed + footer)


def write_deepvariant_outputs(vcf: Path, gvcf: Path, sample_id: str) -> None:
    base = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=4>\n"
    columns = f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}\n"
    markers = (
        "##DeepVariant_version=1.10.0\n"
        "##FILTER=<ID=RefCall,Description=\"Reference call\">\n"
        "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP\">\n"
    )
    write_bgzf(vcf, (base + columns).encode())
    write_bgzf(gvcf, (base + markers + columns).encode())
    write_bgzf(Path(f"{vcf}.tbi"), b"TBI\x01")
    write_bgzf(Path(f"{gvcf}.tbi"), b"TBI\x01")


class FakeCommandRunner:
    def __init__(self) -> None:
        self.commands: list[tuple[tuple[str, ...], dict[str, object]]] = []
        self.require_calls: list[str] = []
        self.version_output = "DeepVariant version 1.10.0\n"
        self.available = True
        self.fail = False
        self.create_outputs = True

    def require_executable(self, executable: str) -> Path:
        self.require_calls.append(executable)
        if not self.available:
            raise ToolNotFoundError(f"Required executable was not found: {executable}")
        return Path(f"/opt/{executable}")

    def run(self, command: object, **kwargs: object) -> CommandResult:
        args = tuple(str(value) for value in command)  # type: ignore[union-attr]
        self.commands.append((args, dict(kwargs)))
        if args[-1] == "--version":
            return CommandResult(args, 0, self.version_output, "", 0.1, None, True)
        if kwargs.get("dry_run") is True:
            return CommandResult(args, None, None, None, 0.0, None, False)
        if self.fail:
            raise CommandExecutionError("Command failed with return code 17: DeepVariant")
        if self.create_outputs:
            outputs: dict[str, Path] = {}
            for argument in args:
                if argument.startswith("--output_vcf="):
                    outputs["vcf"] = Path(argument.split("=", 1)[1])
                elif argument.startswith("--output_gvcf="):
                    outputs["gvcf"] = Path(argument.split("=", 1)[1])
            sample = next(
                argument.split("=", 1)[1]
                for argument in args
                if argument.startswith("--sample_name=")
            )
            write_deepvariant_outputs(outputs["vcf"], outputs["gvcf"], sample)
        return CommandResult(args, 0, "", "metrics", 12.0, None, True)


def make_request(tmp_path: Path, *, overwrite: bool = False) -> DeepVariantRequest:
    tmp_path.mkdir(parents=True, exist_ok=True)
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    bam = tmp_path / "S1.bam"
    bam.write_bytes(b"alignment")
    bai = Path(f"{bam}.bai")
    bai.write_bytes(b"index")
    artifact = AlignmentArtifact(
        sample_id="S1",
        path=bam,
        output_format=AlignmentOutputFormat.BAM,
        reference=ReferenceGenome.from_fasta(fasta, build="GRCh38"),
        source=AlignmentSource.EXISTING,
        sort_order=AlignmentSortOrder.UNKNOWN,
        index_path=bai,
    )
    return DeepVariantRequest.create(
        artifact,
        tmp_path / "results" / "small",
        overwrite=overwrite,
    )


def test_native_command_is_deterministic_hifi_command(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    command = DeepVariantWrapper().build_command(request)
    assert command[0] == "run_deepvariant"
    assert "--model_type=PACBIO" in command
    assert f"--reads={request.alignment_path.absolute()}" in command
    assert f"--ref={request.reference_fasta.absolute()}" in command
    assert "--sample_name=S1" in command
    assert f"--output_vcf={request.output_vcf.absolute()}" in command
    assert f"--output_gvcf={request.output_gvcf.absolute()}" in command
    assert "--num_shards=8" in command
    assert request.temporary_directory == request.intermediate_directory / "tmp"
    assert request.to_dict()["temporary_directory"] == str(
        request.temporary_directory
    )


@pytest.mark.parametrize(
    "mode,launcher",
    ((DeepVariantExecutionMode.DOCKER, "docker"), (DeepVariantExecutionMode.APPTAINER, "apptainer")),
)
def test_container_runtime_isolated_from_caller_arguments(
    tmp_path: Path,
    mode: DeepVariantExecutionMode,
    launcher: str,
) -> None:
    request = make_request(tmp_path)
    runtime = DeepVariantRuntime(mode=mode, image="deepvariant:test")
    command = DeepVariantWrapper(runtime=runtime).build_command(request)
    assert command[0] == launcher
    assert "/opt/deepvariant/bin/run_deepvariant" in command
    assert "--model_type=PACBIO" in command
    assert "--env" in command
    assert f"TMPDIR={request.temporary_directory.absolute()}" in command
    if mode is DeepVariantExecutionMode.DOCKER:
        assert any(
            str(request.temporary_directory.absolute()) in argument
            for argument in command
            if argument.startswith("type=bind")
        )
    else:
        assert (
            f"{request.temporary_directory.absolute()}:"
            f"{request.temporary_directory.absolute()}"
        ) in command
    assert runtime.to_dict()["image"] == "deepvariant:test"


def test_container_command_precreates_all_writable_bind_sources(
    tmp_path: Path,
) -> None:
    request = make_request(tmp_path)
    writable = {
        request.output_vcf.parent,
        request.output_gvcf.parent,
        request.intermediate_directory,
        request.logging_directory,
        request.temporary_directory,
    }
    assert all(not path.exists() for path in writable)

    DeepVariantWrapper(
        runtime=DeepVariantRuntime(
            mode=DeepVariantExecutionMode.APPTAINER,
            image="deepvariant_1.10.0.sif",
        )
    ).build_command(request, create_writable_mounts=True)

    assert all(path.is_dir() for path in writable)


def test_container_dry_run_does_not_create_mount_sources(tmp_path: Path) -> None:
    runner = FakeCommandRunner()
    runner.available = False
    request = make_request(tmp_path)
    wrapper = DeepVariantWrapper(
        runner=runner,  # type: ignore[arg-type]
        runtime=DeepVariantRuntime(
            mode=DeepVariantExecutionMode.APPTAINER,
            image="deepvariant_1.10.0.sif",
        ),
    )

    result = wrapper.run(request, dry_run=True)

    assert result.status is SmallVariantResultStatus.PLANNED
    assert not request.output_vcf.parent.exists()
    assert not request.intermediate_directory.exists()
    assert not request.logging_directory.exists()
    assert not request.temporary_directory.exists()


def test_runtime_config_and_invalid_combinations() -> None:
    runtime = DeepVariantRuntime.from_config(
        {"small": {"execution_mode": "docker", "deepvariant_image": "google/deepvariant:1.10.0"}}
    )
    assert runtime.mode is DeepVariantExecutionMode.DOCKER
    assert runtime.launcher == "docker"
    with pytest.raises(ConfigurationError, match="requires an image"):
        DeepVariantRuntime(mode=DeepVariantExecutionMode.APPTAINER)
    with pytest.raises(ConfigurationError, match="cannot set an image"):
        DeepVariantRuntime(image="unexpected")


def test_version_detection_and_errors() -> None:
    runner = FakeCommandRunner()
    wrapper = DeepVariantWrapper(runner=runner)  # type: ignore[arg-type]
    assert wrapper.detect_version() == "1.10.0"
    assert runner.require_calls == ["run_deepvariant"]
    runner.version_output = "unknown"
    with pytest.raises(ToolVersionError, match="parse DeepVariant"):
        wrapper.detect_version()
    runner.available = False
    with pytest.raises(ToolNotFoundError):
        wrapper.detect_version()


def test_dry_run_requires_inputs_but_not_deepvariant(tmp_path: Path) -> None:
    runner = FakeCommandRunner()
    runner.available = False
    request = make_request(tmp_path)
    result = DeepVariantWrapper(runner=runner).run(request, dry_run=True)  # type: ignore[arg-type]
    assert result.status is SmallVariantResultStatus.PLANNED
    assert result.executed is False
    assert runner.require_calls == []
    assert not request.output_vcf.exists()
    assert not request.output_vcf.parent.exists()


def test_real_run_records_version_runtime_and_outputs(tmp_path: Path) -> None:
    runner = FakeCommandRunner()
    request = make_request(tmp_path)
    result = DeepVariantWrapper(runner=runner).run(request)  # type: ignore[arg-type]
    assert result.status is SmallVariantResultStatus.COMPLETED
    assert result.tool_version == "1.10.0"
    assert result.duration_seconds == 12.0
    assert request.output_vcf.is_file()
    assert request.output_gvcf.is_file()
    assert request.output_vcf_index.is_file()
    assert request.output_gvcf_index.is_file()
    assert request.temporary_directory.is_dir()
    assert runner.commands[-1][1]["env"] == {
        "TMPDIR": str(request.temporary_directory.absolute())
    }


def test_native_dry_run_records_sample_specific_tmpdir_without_creating_it(
    tmp_path: Path,
) -> None:
    runner = FakeCommandRunner()
    request = make_request(tmp_path)

    DeepVariantWrapper(runner=runner).run(request, dry_run=True)  # type: ignore[arg-type]

    assert runner.commands[-1][1]["env"] == {
        "TMPDIR": str(request.temporary_directory.absolute())
    }
    assert not request.temporary_directory.exists()


def test_command_failure_and_missing_outputs_are_clear(tmp_path: Path) -> None:
    runner = FakeCommandRunner()
    runner.fail = True
    with pytest.raises(CommandExecutionError, match="return code 17"):
        DeepVariantWrapper(runner=runner).run(make_request(tmp_path / "fail"))  # type: ignore[arg-type]
    runner.fail = False
    runner.create_outputs = False
    with pytest.raises(OutputValidationError, match="missing"):
        DeepVariantWrapper(runner=runner).run(make_request(tmp_path / "missing"))  # type: ignore[arg-type]


def test_output_race_refused_and_explicit_overwrite_replaces(tmp_path: Path) -> None:
    runner = FakeCommandRunner()
    request = make_request(tmp_path / "race")
    request.output_vcf.parent.mkdir(parents=True)
    request.output_vcf.write_bytes(b"race")
    with pytest.raises(OutputValidationError, match="already exists"):
        DeepVariantWrapper(runner=runner).run(request)  # type: ignore[arg-type]
    assert request.output_vcf.read_bytes() == b"race"

    overwrite_request = make_request(tmp_path / "overwrite", overwrite=True)
    overwrite_request.output_vcf.parent.mkdir(parents=True)
    overwrite_request.output_vcf.write_bytes(b"old")
    DeepVariantWrapper(runner=runner).run(overwrite_request)  # type: ignore[arg-type]
    assert overwrite_request.output_vcf.read_bytes() != b"old"


def test_missing_alignment_index_fails_before_execution(tmp_path: Path) -> None:
    runner = FakeCommandRunner()
    request = make_request(tmp_path)
    request.artifact.index_path.unlink()  # type: ignore[union-attr]
    with pytest.raises(OutputValidationError, match="missing"):
        DeepVariantWrapper(runner=runner).run(request)  # type: ignore[arg-type]
    assert runner.commands == []


def test_low_file_descriptor_limit_fails_before_tool_launch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        deepvariant_module,
        "_get_file_descriptor_limits",
        lambda: (1024, 4096),
    )

    with pytest.raises(ConfigurationError, match="open-file soft limit is 1024"):
        DeepVariantWrapper().run(make_request(tmp_path))


def test_validation_failure_quarantines_completed_outputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runner = FakeCommandRunner()
    request = make_request(tmp_path)

    def reject_outputs(*args: object, **kwargs: object) -> object:
        raise OutputValidationError("synthetic validation failure")

    monkeypatch.setattr(
        deepvariant_module,
        "validate_small_variant_outputs",
        reject_outputs,
    )

    with pytest.raises(OutputValidationError, match="moved to quarantine"):
        DeepVariantWrapper(runner=runner).run(request)  # type: ignore[arg-type]

    assert not request.output_vcf.exists()
    assert not request.output_gvcf.exists()
    runs = list((request.output_vcf.parent / "quarantine").iterdir())
    assert len(runs) == 1
    quarantined = runs[0]
    assert (quarantined / request.output_vcf.name).is_file()
    assert (quarantined / request.output_gvcf.name).is_file()
    assert (quarantined / request.output_vcf_index.name).is_file()
    assert (quarantined / request.output_gvcf_index.name).is_file()
    assert "synthetic validation failure" in (
        quarantined / "VALIDATION_ERROR.txt"
    ).read_text(encoding="utf-8")
