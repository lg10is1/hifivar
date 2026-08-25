from __future__ import annotations

from pathlib import Path

import pytest

from hifivar.assembly import (
    AssemblyRequest,
    AssemblyResources,
    AssemblyRole,
    convert_gfa_to_fasta,
)
from hifivar.command import CommandResult
from hifivar.exceptions import (
    CommandExecutionError,
    InputValidationError,
    OutputValidationError,
    ToolNotFoundError,
    ToolVersionError,
)
from hifivar.hifiasm import AssemblyResultStatus, HifiasmWrapper
from hifivar.sample import InputDataset, InputType, Sample


def write_fastq(path: Path, name: str = "read") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"@{name}\nACGT\n+\n!!!!\n", encoding="utf-8")


def make_request(tmp_path: Path, *, multiple: bool = False) -> AssemblyRequest:
    first = tmp_path / "reads1.fastq"
    write_fastq(first, "one")
    files = [first]
    if multiple:
        second = tmp_path / "reads2.fastq.gz"
        import gzip

        with gzip.open(second, "wt", encoding="utf-8") as handle:
            handle.write("@two\nTGCA\n+\n!!!!\n")
        files.append(second)
    sample = Sample("S1", InputDataset.from_files(files))
    return AssemblyRequest(
        sample,
        tmp_path / "work" / "hifiasm" / "S1" / "S1.asm",
        tmp_path / "results" / "assembly" / "S1",
        AssemblyResources(threads=12),
    )


class FakeRunner:
    def __init__(
        self,
        *,
        missing: bool = False,
        fail: bool = False,
        materialize: bool = False,
        bad_version: bool = False,
    ) -> None:
        self.missing = missing
        self.fail = fail
        self.materialize = materialize
        self.bad_version = bad_version
        self.calls: list[tuple[tuple[str, ...], dict[str, object]]] = []

    def require_executable(self, executable: str) -> Path:
        if self.missing:
            raise ToolNotFoundError(f"missing {executable}")
        return Path(executable)

    def run(self, command, **kwargs):
        args = tuple(str(value) for value in command)
        self.calls.append((args, kwargs))
        if "--version" in args:
            value = "unknown" if self.bad_version else "hifiasm version 0.25.0-r726"
            return CommandResult(args, 0, value, "", 0.01, None, True)
        if kwargs.get("dry_run"):
            return CommandResult(args, None, None, None, 0.0, None, False)
        if self.fail:
            raise CommandExecutionError("synthetic hifiasm failure")
        if self.materialize:
            prefix = Path(args[args.index("-o") + 1])
            prefix.parent.mkdir(parents=True, exist_ok=True)
            for suffix, sequence in (
                (".bp.p_ctg.gfa", "ACGT"),
                (".bp.hap1.p_ctg.gfa", "AAAA"),
                (".bp.hap2.p_ctg.gfa", "TTTT"),
            ):
                Path(f"{prefix}{suffix}").write_text(
                    f"H\tVN:Z:1.0\nS\tctg1\t{sequence}\n",
                    encoding="utf-8",
                )
        return CommandResult(args, 0, "", "", 1.5, None, True)


def test_single_and_multiple_fastq_command_order(tmp_path: Path) -> None:
    single = make_request(tmp_path / "single")
    multiple = make_request(tmp_path / "multiple", multiple=True)
    wrapper = HifiasmWrapper(executable="/opt/hifiasm")
    one = wrapper.plan_command(single).args
    many = wrapper.plan_command(multiple).args
    assert one[:5] == (
        "/opt/hifiasm",
        "-o",
        str(single.output_prefix.absolute()),
        "-t",
        "12",
    )
    assert one[5:] == tuple(str(path.absolute()) for path in single.fastq_files)
    assert many[5:] == tuple(str(path.absolute()) for path in multiple.fastq_files)


def test_bam_is_explicitly_unsupported(tmp_path: Path) -> None:
    bam = tmp_path / "S1.bam"
    bam.write_bytes(b"BAM")
    sample = Sample("S1", InputDataset(InputType.BAM, (bam,)))
    with pytest.raises(InputValidationError, match="not converted"):
        AssemblyRequest(
            sample,
            tmp_path / "work" / "S1.asm",
            tmp_path / "assembly",
        )


def test_version_missing_and_parse_failure() -> None:
    with pytest.raises(ToolNotFoundError):
        HifiasmWrapper(runner=FakeRunner(missing=True)).detect_version()  # type: ignore[arg-type]
    with pytest.raises(ToolVersionError):
        HifiasmWrapper(runner=FakeRunner(bad_version=True)).detect_version()  # type: ignore[arg-type]


def test_dry_run_writes_nothing_and_needs_no_tool(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    runner = FakeRunner(missing=True)
    result = HifiasmWrapper(runner=runner).run(request, dry_run=True)  # type: ignore[arg-type]
    assert result.status is AssemblyResultStatus.PLANNED
    assert len(runner.calls) == 1
    assert not request.output_prefix.parent.exists()
    assert not request.assembly_directory.exists()


def test_fake_execution_preserves_gfa_and_builds_haplotype_fastas(tmp_path: Path) -> None:
    request = make_request(tmp_path, multiple=True)
    result = HifiasmWrapper(runner=FakeRunner(materialize=True)).run(request)  # type: ignore[arg-type]
    assert result.status is AssemblyResultStatus.COMPLETED
    assert result.hifiasm_version == "0.25.0-r726"
    assert result.artifact is not None
    assert {item.role for item in result.artifact.assemblies} == {
        AssemblyRole.PRIMARY,
        AssemblyRole.HAPLOTYPE1,
        AssemblyRole.HAPLOTYPE2,
    }
    for item in result.artifact.assemblies:
        assert item.path.is_file()
        assert item.source_gfa.is_file()
        assert item.reference_independent is True
        assert item.file_size > 0
    assert result.artifact.assemblies[1].path.read_text(encoding="utf-8") == ">ctg1\nAAAA\n"


def test_explicit_conversion_rejects_sequence_free_gfa(tmp_path: Path) -> None:
    source = tmp_path / "empty.gfa"
    source.write_text("H\tVN:Z:1.0\nS\tctg1\t*\n", encoding="utf-8")
    with pytest.raises(OutputValidationError, match="lacks an embedded sequence"):
        convert_gfa_to_fasta(source, tmp_path / "output.fa")
    assert source.is_file()
    assert not (tmp_path / "output.fa").exists()


def test_failure_and_overwrite_policy(tmp_path: Path) -> None:
    request = make_request(tmp_path)
    with pytest.raises(CommandExecutionError):
        HifiasmWrapper(runner=FakeRunner(fail=True)).run(request)  # type: ignore[arg-type]
    raw = request.raw_gfa_paths[AssemblyRole.PRIMARY]
    raw.parent.mkdir(parents=True, exist_ok=True)
    raw.write_text("old", encoding="utf-8")
    with pytest.raises(OutputValidationError, match="already exists"):
        make_request(tmp_path)
