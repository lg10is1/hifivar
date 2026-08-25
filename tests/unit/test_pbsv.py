from pathlib import Path

import pytest

from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import AlignmentArtifact, AlignmentSortOrder, AlignmentSource
from hifivar.command import CommandResult
from hifivar.exceptions import CommandExecutionError, InputValidationError, OutputValidationError, ToolNotFoundError, ToolVersionError
from hifivar.pbsv import PbsvRequest, PbsvResources, PbsvResultStatus, PbsvWrapper
from hifivar.reference import ReferenceGenome


class FakeRunner:
    def __init__(self, *, version="pbsv 2.11.0", missing=False, fail_step=None, materialize=False):
        self.version, self.missing, self.fail_step, self.materialize = version, missing, fail_step, materialize
        self.calls = []
        self.request = None

    def require_executable(self, executable):
        if self.missing:
            raise ToolNotFoundError("missing pbsv")
        return Path(executable)

    def run(self, command, **kwargs):
        args = tuple(str(item) for item in command)
        self.calls.append((args, kwargs))
        if "--version" in args:
            return CommandResult(args, 0, self.version, "", 0.1, None, True)
        step = args[1]
        if self.fail_step == step:
            raise CommandExecutionError("pbsv failed")
        if self.materialize and step == "discover":
            self.request.signatures_path.parent.mkdir(parents=True, exist_ok=True)
            self.request.signatures_path.write_bytes(b"sig")
        if self.materialize and step == "call":
            self.request.raw_vcf.parent.mkdir(parents=True, exist_ok=True)
            self.request.raw_vcf.write_text("##fileformat=VCFv4.3\n", encoding="utf-8")
        return CommandResult(args, 0, "", "", 0.2, None, not kwargs.get("dry_run", False))


@pytest.fixture
def pbsv_request(tmp_path):
    fasta = tmp_path / "参考.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    bam = tmp_path / "样本.bam"
    bam.write_bytes(b"BAM")
    bai = Path(f"{bam}.bai")
    bai.write_bytes(b"BAI")
    artifact = AlignmentArtifact("S1", bam, AlignmentOutputFormat.BAM, ReferenceGenome.from_fasta(fasta, build="GRCh38"), AlignmentSource.EXISTING, AlignmentSortOrder.UNKNOWN, bai)
    return PbsvRequest.create(artifact, tmp_path / "结果", tmp_path / "工作", resources=PbsvResources(threads=12))


def test_commands_are_separate_and_deterministic(pbsv_request):
    wrapper = PbsvWrapper(executable="/opt/pbsv")
    commands = wrapper.build_commands(pbsv_request)
    assert commands == wrapper.build_commands(pbsv_request)
    assert commands[0][1] == "discover"
    assert commands[1][1:5] == ["call", "--ccs", "-j", "12"]
    assert not any("|" in item for command in commands for item in command)


def test_version_detection_and_missing_executable():
    assert PbsvWrapper(runner=FakeRunner()).detect_version() == "2.11.0"
    with pytest.raises(ToolVersionError):
        PbsvWrapper(runner=FakeRunner(version="unknown")).detect_version()
    with pytest.raises(ToolNotFoundError):
        PbsvWrapper(runner=FakeRunner(missing=True)).detect_version()


def test_dry_run_plans_both_steps_without_writes(pbsv_request):
    runner = FakeRunner(missing=True)
    result = PbsvWrapper(runner=runner).run(pbsv_request, dry_run=True)
    assert result.status is PbsvResultStatus.PLANNED
    assert len(runner.calls) == 2
    assert not pbsv_request.signatures_path.exists() and not pbsv_request.raw_vcf.exists()


def test_execution_validates_signature_and_raw_vcf(pbsv_request):
    runner = FakeRunner(materialize=True)
    runner.request = pbsv_request
    result = PbsvWrapper(runner=runner).run(pbsv_request)
    assert result.status is PbsvResultStatus.COMPLETED
    assert pbsv_request.signatures_path.exists() and pbsv_request.raw_vcf.exists()


@pytest.mark.parametrize("step", ["discover", "call"])
def test_each_step_failure_is_propagated(pbsv_request, step):
    runner = FakeRunner(fail_step=step, materialize=True)
    runner.request = pbsv_request
    with pytest.raises(CommandExecutionError):
        PbsvWrapper(runner=runner).run(pbsv_request)


def test_missing_discover_output_stops_before_call(pbsv_request):
    runner = FakeRunner()
    with pytest.raises(OutputValidationError):
        PbsvWrapper(runner=runner).run(pbsv_request)
    assert [call[0][1] for call in runner.calls] == ["--version", "discover"]


def test_no_silent_overwrite(pbsv_request):
    pbsv_request.signatures_path.parent.mkdir(parents=True)
    pbsv_request.signatures_path.write_bytes(b"old")
    with pytest.raises(OutputValidationError):
        PbsvWrapper(runner=FakeRunner()).run(pbsv_request)


def test_cram_is_rejected_explicitly(pbsv_request):
    artifact = pbsv_request.artifact
    cram_index = artifact.path.with_suffix(".cram.crai")
    with pytest.raises(InputValidationError, match="requires BAM"):
        PbsvRequest(AlignmentArtifact(artifact.sample_id, artifact.path.with_suffix(".cram"), AlignmentOutputFormat.CRAM, artifact.reference, artifact.source, artifact.sort_order, cram_index), pbsv_request.signatures_path, pbsv_request.raw_vcf)
