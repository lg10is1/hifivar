from pathlib import Path

import pytest

from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import AlignmentArtifact, AlignmentSortOrder, AlignmentSource
from hifivar.command import CommandResult
from hifivar.exceptions import CommandExecutionError, OutputValidationError, ToolNotFoundError, ToolVersionError
from hifivar.reference import ReferenceGenome
from hifivar.sawfish import SawfishRequest, SawfishResources, SawfishResultStatus, SawfishWrapper


class FakeRunner:
    def __init__(self, *, version="sawfish 2.2.1", missing=False, fail_step=None, materialize=False):
        self.version = version
        self.missing = missing
        self.fail_step = fail_step
        self.materialize = materialize
        self.calls = []
        self.request = None

    def require_executable(self, executable):
        if self.missing:
            raise ToolNotFoundError("missing sawfish")
        return Path(executable)

    def run(self, command, **kwargs):
        args = tuple(str(item) for item in command)
        self.calls.append((args, kwargs))
        if "--version" in args:
            return CommandResult(args, 0, self.version, "", 0.1, None, True)
        step = args[1]
        if self.fail_step == step:
            raise CommandExecutionError("sawfish failed")
        if self.materialize and step == "joint-call":
            self.request.native_vcf.parent.mkdir(parents=True, exist_ok=True)
            self.request.native_vcf.write_bytes(b"vcf")
            self.request.native_index.write_bytes(b"tbi")
        return CommandResult(args, 0, "", "", 0.2, None, not kwargs.get("dry_run", False))


@pytest.fixture
def sawfish_request(tmp_path):
    fasta = tmp_path / "参考.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    bam = tmp_path / "样本.bam"
    bam.write_bytes(b"BAM")
    bai = Path(f"{bam}.bai")
    bai.write_bytes(b"BAI")
    artifact = AlignmentArtifact(
        sample_id="S1",
        path=bam,
        output_format=AlignmentOutputFormat.BAM,
        reference=ReferenceGenome.from_fasta(fasta, build="GRCh38"),
        source=AlignmentSource.EXISTING,
        sort_order=AlignmentSortOrder.UNKNOWN,
        index_path=bai,
    )
    return SawfishRequest.create(artifact, tmp_path / "结果", tmp_path / "工作", resources=SawfishResources(threads=8))


def test_command_construction_is_deterministic(sawfish_request):
    wrapper = SawfishWrapper(executable="/opt/sawfish")
    first = wrapper.build_commands(sawfish_request)
    assert first == wrapper.build_commands(sawfish_request)
    assert first[0][1:4] == ["discover", "--threads", "8"]
    assert "--bam" in first[0] and "--ref" in first[0]
    assert first[1][1] == "joint-call" and "--sample" in first[1]
    assert all(isinstance(item, str) for command in first for item in command)


def test_version_detection(sawfish_request):
    assert SawfishWrapper(runner=FakeRunner()).detect_version() == "2.2.1"
    with pytest.raises(ToolVersionError):
        SawfishWrapper(runner=FakeRunner(version="unknown")).detect_version()
    with pytest.raises(ToolNotFoundError):
        SawfishWrapper(runner=FakeRunner(missing=True)).detect_version()


def test_dry_run_needs_no_executable_and_writes_nothing(sawfish_request):
    runner = FakeRunner(missing=True)
    result = SawfishWrapper(runner=runner).run(sawfish_request, dry_run=True)
    assert result.status is SawfishResultStatus.PLANNED
    assert len(runner.calls) == 2
    assert not sawfish_request.output_vcf.exists()
    assert not sawfish_request.work_directory.exists()


def test_execution_materializes_named_output(sawfish_request):
    runner = FakeRunner(materialize=True)
    runner.request = sawfish_request
    result = SawfishWrapper(runner=runner).run(sawfish_request)
    assert result.status is SawfishResultStatus.COMPLETED
    assert sawfish_request.output_vcf.read_bytes() == b"vcf"
    assert sawfish_request.output_index.read_bytes() == b"tbi"
    assert result.tool_version == "2.2.1"


@pytest.mark.parametrize("step", ["discover", "joint-call"])
def test_each_step_failure_is_propagated(sawfish_request, step):
    runner = FakeRunner(fail_step=step)
    with pytest.raises(CommandExecutionError):
        SawfishWrapper(runner=runner).run(sawfish_request)


def test_no_silent_overwrite(sawfish_request):
    sawfish_request.output_vcf.parent.mkdir(parents=True)
    sawfish_request.output_vcf.write_bytes(b"old")
    with pytest.raises(OutputValidationError):
        SawfishWrapper(runner=FakeRunner()).run(sawfish_request)


def test_missing_native_output_is_rejected(sawfish_request):
    with pytest.raises(OutputValidationError):
        SawfishWrapper(runner=FakeRunner()).run(sawfish_request)
