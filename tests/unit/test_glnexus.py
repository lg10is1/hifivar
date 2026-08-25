from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from hifivar.cohort import CohortDefinition, CohortSampleInput, SampleCallState
from hifivar.command import CommandResult
from hifivar.exceptions import CommandExecutionError, InputValidationError, OutputValidationError, ToolNotFoundError
from hifivar.glnexus import GLnexusRequest, GLnexusResources, GLnexusRunStatus, GLnexusWrapper
from hifivar.reference import Contig, ReferenceGenome


def write_gvcf(path: Path, sample: str, contig: str = "chr1") -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(f"##fileformat=VCFv4.2\n##contig=<ID={contig},length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
    Path(f"{path}.tbi").write_bytes(b"index")


def make_request(tmp_path: Path) -> GLnexusRequest:
    tmp_path.mkdir(parents=True, exist_ok=True)
    ref = ReferenceGenome(tmp_path / "ref.fa", tmp_path / "ref.fa.fai", "GRCh38", (Contig("chr1", 1000),))
    inputs = []
    for sample in ("S2", "S1"):
        path = tmp_path / f"{sample}.g.vcf.gz"
        write_gvcf(path, sample)
        inputs.append(CohortSampleInput(sample, SampleCallState.CALLED, path, Path(f"{path}.tbi"), "deepvariant", "1.10.0", "GRCh38"))
    return GLnexusRequest(CohortDefinition("C1", ("S2", "S1"), ref), tuple(inputs), tmp_path / "db", tmp_path / "out" / "C1.small.bcf", tmp_path / "out" / "C1.small.vcf.gz", resources=GLnexusResources(4, 8))


class FakeRunner:
    def __init__(self, request: GLnexusRequest, *, fail: bool = False) -> None:
        self.request, self.fail, self.calls = request, fail, []

    def require_executable(self, executable, **kwargs):
        return Path(str(executable))

    def run(self, command, **kwargs):
        args = tuple(str(item) for item in command)
        self.calls.append((args, kwargs))
        if self.fail and not kwargs.get("dry_run"):
            raise CommandExecutionError("fake failure")
        stdout = ""
        if args[-1:] == ("--help",):
            stdout = "glnexus_cli release v1.4.1\n"
        elif args[-1:] == ("--version",):
            stdout = "bcftools 1.21\n"
        elif not kwargs.get("dry_run") and args[0] == "glnexus_cli":
            Path(kwargs["stdout_path"]).parent.mkdir(parents=True, exist_ok=True)
            Path(kwargs["stdout_path"]).write_bytes(b"BCF")
        elif not kwargs.get("dry_run") and args[1] == "view":
            output = Path(args[args.index("-o") + 1])
            with gzip.open(output, "wt", encoding="utf-8") as handle:
                handle.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\tS1\nchr1\t10\tv1\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\n")
        elif not kwargs.get("dry_run") and args[1] == "index":
            target = Path(args[-1])
            Path(f"{target}.tbi" if "--tbi" in args else f"{target}.csi").write_bytes(b"index")
        return CommandResult(args, 0, stdout, "", 0.1, None, not kwargs.get("dry_run"), kwargs.get("stdout_path"), kwargs.get("stderr_path"))


class MissingRunner(FakeRunner):
    def require_executable(self, executable, **kwargs):
        raise ToolNotFoundError(f"missing {executable}")


def test_command_is_deterministic_shell_free_and_preserves_input_order(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    command = GLnexusWrapper().plan_commands(request)[0]
    assert command.args[:9] == ("glnexus_cli", "--dir", str(request.work_directory.absolute()), "--config", "DeepVariantWGS", "--threads", "4", "--mem-gbytes", "8")
    assert command.args[-2:] == tuple(str(item.source_path.absolute()) for item in request.inputs)
    assert command.stdout_path == request.output_bcf
    assert "|" not in command.args


def test_dry_run_creates_no_outputs(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    runner = FakeRunner(request)
    result = GLnexusWrapper(runner=runner).run(request, dry_run=True)
    assert result.status is GLnexusRunStatus.PLANNED
    assert not request.output_bcf.exists()


def test_fake_end_to_end_validates_sample_set_and_qc(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    runner = FakeRunner(request)
    result = GLnexusWrapper(runner=runner).run(request)
    assert result.status is GLnexusRunStatus.COMPLETED
    assert result.versions == {"glnexus": "1.4.1", "bcftools": "1.21"}
    assert result.qc["sample_count"] == 2
    assert all(path.exists() for path in request.expected_outputs)
    assert result.as_track_result().tool == "glnexus"


def test_missing_index_and_reference_contig_mismatch(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    request.inputs[0].index_path.unlink()
    with pytest.raises(InputValidationError, match="index"):
        GLnexusWrapper(runner=FakeRunner(request)).run(request)
    request = make_request(tmp_path / "other")
    write_gvcf(request.inputs[0].source_path, "S2", "1")
    with pytest.raises(InputValidationError, match="REFERENCE_CONTIG_MISMATCH"):
        GLnexusWrapper(runner=FakeRunner(request)).run(request)


def test_no_silent_partial_cohort_or_overwrite(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    bad = (request.inputs[0], CohortSampleInput("S1", SampleCallState.FAILED))
    with pytest.raises(InputValidationError, match="every cohort sample"):
        GLnexusRequest(request.cohort, bad, request.work_directory, request.output_bcf, request.output_vcf)
    request.output_bcf.parent.mkdir(parents=True, exist_ok=True)
    request.output_bcf.write_bytes(b"old")
    with pytest.raises(OutputValidationError, match="already exists"):
        GLnexusRequest(request.cohort, request.inputs, request.work_directory, request.output_bcf, request.output_vcf)


def test_external_failure_is_not_hidden(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    with pytest.raises(CommandExecutionError, match="fake failure"):
        GLnexusWrapper(runner=FakeRunner(request, fail=True)).run(request)


def test_missing_executable_is_clear(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    with pytest.raises(ToolNotFoundError, match="glnexus_cli"):
        GLnexusWrapper(runner=MissingRunner(request)).run(request)


def test_wrapper_rejects_output_sample_mismatch(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    runner = FakeRunner(request)
    runner.request = request
    original_run = runner.run
    def reversed_output(command, **kwargs):
        result = original_run(command, **kwargs)
        args = tuple(str(item) for item in command)
        if len(args) > 1 and args[1] == "view":
            with gzip.open(request.output_vcf, "wt", encoding="utf-8") as handle:
                handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n")
        return result
    runner.run = reversed_output
    with pytest.raises(OutputValidationError, match="order/set mismatch"):
        GLnexusWrapper(runner=runner).run(request)
