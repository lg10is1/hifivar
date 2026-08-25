from pathlib import Path
from dataclasses import replace

import pytest

from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import AlignmentArtifact, AlignmentSortOrder, AlignmentSource
from hifivar.command import CommandResult
from hifivar.cutesv import CuteSvRequest, CuteSvResources, CuteSvResultStatus, CuteSvWrapper
from hifivar.exceptions import CommandExecutionError, OutputValidationError, ToolNotFoundError, ToolVersionError
from hifivar.reference import ReferenceGenome


class FakeRunner:
    def __init__(self, *, version="cuteSV 2.1.4", missing=False, fail=False, materialize=False):
        self.version, self.missing, self.fail, self.materialize = version, missing, fail, materialize
        self.calls = []
        self.request = None

    def require_executable(self, executable):
        if self.missing:
            raise ToolNotFoundError("missing cuteSV")
        return Path(executable)

    def run(self, command, **kwargs):
        args = tuple(str(item) for item in command)
        self.calls.append((args, kwargs))
        if "--version" in args:
            return CommandResult(args, 0, self.version, "", 0.1, None, True)
        if self.fail:
            raise CommandExecutionError("cuteSV failed")
        if self.materialize:
            self.request.raw_vcf.write_text("##fileformat=VCFv4.3\n", encoding="utf-8")
        return CommandResult(args, 0, "", "", 0.2, None, not kwargs.get("dry_run", False))


@pytest.fixture
def cutesv_request(tmp_path):
    fasta = tmp_path / "参考.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    bam = tmp_path / "样本.bam"
    bam.write_bytes(b"BAM")
    bai = Path(f"{bam}.bai")
    bai.write_bytes(b"BAI")
    artifact = AlignmentArtifact("S1", bam, AlignmentOutputFormat.BAM, ReferenceGenome.from_fasta(fasta, build="GRCh38"), AlignmentSource.EXISTING, AlignmentSortOrder.UNKNOWN, bai)
    return CuteSvRequest.create(artifact, tmp_path / "结果", tmp_path / "工作", resources=CuteSvResources(threads=10), minimum_support=5)


def test_command_is_deterministic_and_hifi_parameters_are_separate(cutesv_request):
    wrapper = CuteSvWrapper(executable="/opt/cuteSV")
    command = wrapper.build_command(cutesv_request)
    assert command == wrapper.build_command(cutesv_request)
    assert command[:5] == ["/opt/cuteSV", str(cutesv_request.artifact.path.absolute()), str(cutesv_request.artifact.reference.fasta.absolute()), str(cutesv_request.raw_vcf.absolute()), str(cutesv_request.work_directory.absolute())]
    assert command[command.index("--threads") + 1] == "10"
    assert command[command.index("--min_support") + 1] == "5"
    assert command[command.index("--max_cluster_bias_INS") + 1] == "1000"
    assert "--genotype" in command


def test_genotype_can_be_disabled_explicitly(cutesv_request):
    request = replace(cutesv_request, genotype=False)
    assert "--genotype" not in CuteSvWrapper().build_command(request)
    assert request.to_dict()["genotype"] is False


def test_version_detection_and_missing_executable():
    assert CuteSvWrapper(runner=FakeRunner()).detect_version() == "2.1.4"
    with pytest.raises(ToolVersionError):
        CuteSvWrapper(runner=FakeRunner(version="unknown")).detect_version()
    with pytest.raises(ToolNotFoundError):
        CuteSvWrapper(runner=FakeRunner(missing=True)).detect_version()


def test_dry_run_writes_nothing(cutesv_request):
    runner = FakeRunner(missing=True)
    result = CuteSvWrapper(runner=runner).run(cutesv_request, dry_run=True)
    assert result.status is CuteSvResultStatus.PLANNED
    assert not cutesv_request.work_directory.exists() and not cutesv_request.raw_vcf.exists()


def test_execution_validates_native_vcf(cutesv_request):
    runner = FakeRunner(materialize=True)
    runner.request = cutesv_request
    result = CuteSvWrapper(runner=runner).run(cutesv_request)
    assert result.status is CuteSvResultStatus.COMPLETED
    assert cutesv_request.raw_vcf.exists()


def test_failure_and_missing_output_propagate(cutesv_request):
    with pytest.raises(CommandExecutionError):
        CuteSvWrapper(runner=FakeRunner(fail=True)).run(cutesv_request)
    with pytest.raises(OutputValidationError):
        CuteSvWrapper(runner=FakeRunner()).run(cutesv_request)


def test_no_silent_overwrite(cutesv_request):
    cutesv_request.raw_vcf.parent.mkdir(parents=True)
    cutesv_request.raw_vcf.write_bytes(b"old")
    with pytest.raises(OutputValidationError):
        CuteSvWrapper(runner=FakeRunner()).run(cutesv_request)


def test_overwrite_cleans_only_owned_work_directory(cutesv_request):
    runner = FakeRunner(materialize=True)
    runner.request = cutesv_request
    CuteSvWrapper(runner=runner).run(cutesv_request)
    stale = cutesv_request.work_directory / "signatures" / "stale.pickle"
    stale.parent.mkdir()
    stale.write_bytes(b"stale")

    replacement = replace(cutesv_request, overwrite=True)
    runner.request = replacement
    CuteSvWrapper(runner=runner).run(replacement)
    assert not stale.exists()


def test_overwrite_refuses_unowned_work_directory(cutesv_request):
    cutesv_request.work_directory.mkdir(parents=True)
    foreign = cutesv_request.work_directory / "user-data.txt"
    foreign.write_text("keep", encoding="utf-8")
    request = replace(cutesv_request, overwrite=True)
    with pytest.raises(OutputValidationError, match="not marked as HiFiVar-owned"):
        CuteSvWrapper(runner=FakeRunner()).run(request)
    assert foreign.read_text(encoding="utf-8") == "keep"


def test_execution_refuses_foreign_ownership_marker(cutesv_request):
    marker = cutesv_request.work_directory.parent / ".S1.hifivar-cutesv-owned"
    marker.parent.mkdir(parents=True)
    marker.write_text("foreign", encoding="utf-8")
    with pytest.raises(OutputValidationError, match="unrecognized"):
        CuteSvWrapper(runner=FakeRunner()).run(cutesv_request)
    assert marker.read_text(encoding="utf-8") == "foreign"
