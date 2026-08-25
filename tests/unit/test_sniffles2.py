from pathlib import Path

import pytest

from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import AlignmentArtifact, AlignmentSortOrder, AlignmentSource
from hifivar.command import CommandResult
from hifivar.exceptions import CommandExecutionError, OutputValidationError, ToolNotFoundError, ToolVersionError
from hifivar.reference import ReferenceGenome
from hifivar.sniffles2 import Sniffles2Request, Sniffles2Resources, Sniffles2ResultStatus, Sniffles2Wrapper


class FakeRunner:
    def __init__(self, *, version="Sniffles2 Version 2.8.0", missing=False, fail=False, materialize=False):
        self.version, self.missing, self.fail, self.materialize = version, missing, fail, materialize
        self.calls = []
        self.request = None

    def require_executable(self, executable):
        if self.missing:
            raise ToolNotFoundError("missing sniffles")
        return Path(executable)

    def run(self, command, **kwargs):
        args = tuple(str(item) for item in command)
        self.calls.append((args, kwargs))
        if "--version" in args:
            return CommandResult(args, 0, self.version, "", 0.1, None, True)
        if self.fail:
            raise CommandExecutionError("sniffles failed")
        if self.materialize:
            self.request.output_vcf.parent.mkdir(parents=True, exist_ok=True)
            self.request.output_vcf.write_bytes(b"vcf")
            self.request.output_index.write_bytes(b"tbi")
        return CommandResult(args, 0, "", "", 0.2, None, not kwargs.get("dry_run", False))


@pytest.fixture
def sniffles_request(tmp_path):
    fasta = tmp_path / "参考.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    cram = tmp_path / "样本.cram"
    cram.write_bytes(b"CRAM")
    crai = Path(f"{cram}.crai")
    crai.write_bytes(b"CRAI")
    artifact = AlignmentArtifact("S1", cram, AlignmentOutputFormat.CRAM, ReferenceGenome.from_fasta(fasta, build="GRCh38"), AlignmentSource.EXISTING, AlignmentSortOrder.UNKNOWN, crai)
    return Sniffles2Request.create(artifact, tmp_path / "结果", resources=Sniffles2Resources(threads=6), minimum_support=4, minimum_sv_length=35)


def test_command_is_deterministic_and_caller_parameters_are_explicit(sniffles_request):
    wrapper = Sniffles2Wrapper(executable="/opt/sniffles")
    command = wrapper.build_command(sniffles_request)
    assert command == wrapper.build_command(sniffles_request)
    assert command[:3] == ["/opt/sniffles", "--input", str(sniffles_request.artifact.path.absolute())]
    assert command[command.index("--threads") + 1] == "6"
    assert command[command.index("--sample-id") + 1] == "S1"
    assert command[command.index("--minsupport") + 1] == "4"
    assert command[command.index("--minsvlen") + 1] == "35"


def test_version_detection_and_missing_executable():
    assert Sniffles2Wrapper(runner=FakeRunner()).detect_version() == "2.8.0"
    with pytest.raises(ToolVersionError):
        Sniffles2Wrapper(runner=FakeRunner(version="unknown")).detect_version()
    with pytest.raises(ToolNotFoundError):
        Sniffles2Wrapper(runner=FakeRunner(missing=True)).detect_version()


def test_dry_run_does_not_require_tool_or_write(sniffles_request):
    runner = FakeRunner(missing=True)
    result = Sniffles2Wrapper(runner=runner).run(sniffles_request, dry_run=True)
    assert result.status is Sniffles2ResultStatus.PLANNED
    assert len(runner.calls) == 1 and runner.calls[0][1]["dry_run"] is True
    assert not sniffles_request.output_vcf.exists()


def test_execution_validates_vcf_and_index(sniffles_request):
    runner = FakeRunner(materialize=True)
    runner.request = sniffles_request
    result = Sniffles2Wrapper(runner=runner).run(sniffles_request)
    assert result.status is Sniffles2ResultStatus.COMPLETED
    assert sniffles_request.output_vcf.exists() and sniffles_request.output_index.exists()


def test_failure_and_missing_outputs_propagate(sniffles_request):
    with pytest.raises(CommandExecutionError):
        Sniffles2Wrapper(runner=FakeRunner(fail=True)).run(sniffles_request)
    with pytest.raises(OutputValidationError):
        Sniffles2Wrapper(runner=FakeRunner()).run(sniffles_request)


def test_no_silent_overwrite(sniffles_request):
    sniffles_request.output_vcf.parent.mkdir(parents=True)
    sniffles_request.output_vcf.write_bytes(b"old")
    with pytest.raises(OutputValidationError):
        Sniffles2Wrapper(runner=FakeRunner()).run(sniffles_request)
