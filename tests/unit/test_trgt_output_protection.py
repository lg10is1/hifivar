"""TRGT must reject output/log aliases before probing tools or mutating inputs."""
from __future__ import annotations

from dataclasses import replace
import hashlib
from pathlib import Path

import pytest

from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import AlignmentArtifact, AlignmentSortOrder, AlignmentSource
from hifivar.command import _default_stdout_log_path
from hifivar.exceptions import OutputValidationError
from hifivar.reference import ReferenceGenome
from hifivar.tr import TandemRepeatCatalog
from hifivar.trgt import TrgtRequest, TrgtWrapper, _command_log_path


class NoExecution:
    def require_executable(self, executable):
        pytest.fail("Output protection must precede tool discovery/version execution")

    def run(self, *args, **kwargs):
        pytest.fail("Unsafe paths must fail before any external command")


@pytest.fixture
def case(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    ref = inputs / "reference.fa"
    ref.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
    fai = Path(str(ref) + ".fai")
    fai.write_text("chr1\t8\t6\t8\t9\n", encoding="utf-8")
    bam = inputs / "S1.tr.spanning.bam"
    bam.write_bytes(b"immutable alignment")
    bai = Path(str(bam) + ".bai")
    bai.write_bytes(b"immutable index")
    catalog = inputs / "catalog.bed"
    catalog.write_text("chr1\t1\t7\tID=L1;MOTIFS=CG;STRUC=(CG)n\n", encoding="utf-8")
    resources = {"bam": bam, "index": bai, "reference": ref, "fai": fai, "catalog": catalog}
    for name in ("trgt", "bcftools", "samtools"):
        path = inputs / name
        path.write_bytes(b"immutable tool " + name.encode())
        resources[name] = path
    artifact = AlignmentArtifact("S1", bam, AlignmentOutputFormat.BAM,
                                 ReferenceGenome.from_fasta(ref, build="GRCh38"),
                                 AlignmentSource.EXISTING, AlignmentSortOrder.UNKNOWN, bai)
    request = TrgtRequest(artifact, TandemRepeatCatalog(catalog, "GRCh38"),
                          tmp_path / "work" / "S1.trgt",
                          tmp_path / "outputs" / "S1.tr.vcf.gz",
                          tmp_path / "outputs" / "S1.tr.spanning.bam", "XY", overwrite=True)
    wrapper = TrgtWrapper(executable=str(resources["trgt"]),
                          bcftools_executable=str(resources["bcftools"]),
                          samtools_executable=str(resources["samtools"]), runner=NoExecution())
    return request, wrapper, resources


def fingerprints(resources):
    return {key: (path.read_bytes(), hashlib.sha256(path.read_bytes()).hexdigest())
            for key, path in resources.items()}


def reject_unchanged(case, request=None, *, stderr_path=None):
    original, wrapper, resources = case
    before = fingerprints(resources)
    with pytest.raises(OutputValidationError, match="TRGT_OUTPUT_PROTECTION"):
        wrapper.run(request or original, stderr_path=stderr_path)
    assert fingerprints(resources) == before


@pytest.mark.parametrize("overwrite", [False, True])
def test_direct_input_bam_and_bai_alias_never_deleted(case, overwrite):
    request, _, _ = case
    reject_unchanged(case, replace(request, final_spanning_bam=request.artifact.path,
                                   overwrite=overwrite))


@pytest.mark.parametrize("resource", ["bam", "index", "reference", "fai", "catalog", "trgt", "bcftools", "samtools"])
@pytest.mark.parametrize("role", ["raw_vcf", "final_vcf", "final_spanning_bam_index"])
def test_hardlinked_outputs_preserve_every_protected_resource(case, resource, role):
    request, _, resources = case
    output = getattr(request, role)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.hardlink_to(resources[resource])
    reject_unchanged(case)
    assert output.samefile(resources[resource])


@pytest.mark.parametrize("command_index", range(5))
@pytest.mark.parametrize("stream", ["stderr", "implicit_stdout"])
def test_each_actual_log_target_is_protected_before_any_command(case, command_index, stream):
    request, _, resources = case
    base = request.final_vcf.parent / "tool.log"
    log = _command_log_path(base, command_index)
    assert log is not None
    if stream == "implicit_stdout":
        log = _default_stdout_log_path(log)
    log.parent.mkdir(parents=True, exist_ok=True)
    log.hardlink_to(resources["catalog"])
    reject_unchanged(case, stderr_path=base)


def test_raw_and_final_outputs_must_not_alias(case):
    request, _, _ = case
    reject_unchanged(case, replace(request, raw_output_prefix=request.final_vcf.parent / "S1.tr"))


def test_hardlinked_outputs_must_not_alias_each_other(case):
    request, _, _ = case
    request.raw_vcf.parent.mkdir(parents=True, exist_ok=True)
    request.final_vcf.parent.mkdir(parents=True, exist_ok=True)
    request.raw_vcf.write_bytes(b"old output")
    request.final_vcf.hardlink_to(request.raw_vcf)
    reject_unchanged(case)
    assert request.raw_vcf.read_bytes() == b"old output"
    assert request.final_vcf.read_bytes() == b"old output"


def test_output_and_log_must_not_alias(case):
    request, _, _ = case
    base = request.final_vcf.parent / "tool.log"
    log = _command_log_path(base, 0)
    assert log is not None
    log.parent.mkdir(parents=True, exist_ok=True)
    request.final_vcf.write_bytes(b"old output")
    log.hardlink_to(request.final_vcf)
    reject_unchanged(case, stderr_path=base)
    assert request.final_vcf.read_bytes() == b"old output"


def test_output_parent_cannot_be_an_immutable_file(case):
    request, _, _ = case
    reject_unchanged(case, replace(request, raw_output_prefix=request.artifact.path / "child"))


def test_output_parent_cannot_be_another_output_file(case):
    request, _, _ = case
    reject_unchanged(case, replace(request, raw_output_prefix=request.final_vcf / "child"))


def test_output_directory_is_rejected_without_deleting_contents(case):
    request, _, _ = case
    request.final_vcf.mkdir(parents=True)
    marker = request.final_vcf / "keep"
    marker.write_bytes(b"user directory")
    reject_unchanged(case)
    assert marker.read_bytes() == b"user directory"


@pytest.mark.parametrize("directory_link", [False, True])
def test_resolved_symlink_alias_preserves_input(case, directory_link):
    request, _, resources = case
    link = request.final_vcf.parent
    try:
        if directory_link:
            link.symlink_to(resources["bam"].parent, target_is_directory=True)
            request = replace(request, final_spanning_bam=link / resources["bam"].name)
        else:
            link.mkdir(parents=True)
            request.final_vcf.symlink_to(resources["catalog"])
    except OSError as error:
        pytest.skip(f"Host does not permit creating the symlink fixture: {error}")
    reject_unchanged(case, request)


def test_dangling_output_symlink_fails_closed(case):
    request, _, _ = case
    request.final_vcf.parent.mkdir(parents=True)
    try:
        request.final_vcf.symlink_to(request.final_vcf.parent / "missing")
    except OSError as error:
        pytest.skip(f"Host does not permit creating the symlink fixture: {error}")
    reject_unchanged(case)
    assert request.final_vcf.is_symlink()


def test_noncolliding_siblings_pass_read_only_guard_without_execution(case):
    request, wrapper, resources = case
    before = fingerprints(resources)
    wrapper._protect_write_paths(request, stderr_path=request.raw_output_prefix.parent / "tool.log")
    assert fingerprints(resources) == before
    assert not request.raw_output_prefix.parent.exists()
    assert not request.final_vcf.parent.exists()


def test_bare_executable_discovery_is_also_protected(case, monkeypatch):
    request, wrapper, resources = case
    monkeypatch.setattr("hifivar.trgt.shutil.which", lambda value: str(resources["trgt"]) if value == "trgt" else None)
    wrapper.executable = "trgt"
    request.final_vcf.parent.mkdir(parents=True)
    request.final_vcf.hardlink_to(resources["trgt"])
    reject_unchanged(case)
