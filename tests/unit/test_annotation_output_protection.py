"""Immutable resource aliases must fail before any annotation side effect."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import os

import pytest

from hifivar.annotation import AnnotationInput, VariantCategory
from hifivar.annovar import AnnovarRequest, AnnovarWrapper
from hifivar.exceptions import OutputValidationError
from hifivar.reference import ReferenceGenome
from hifivar.vep import VepRequest, VepWrapper


class ProbeRunner:
    def __init__(self):
        self.calls = []

    def require_executable(self, executable):
        self.calls.append(("require_executable", executable))
        pytest.fail("Unsafe paths reached executable detection")

    def run(self, command, **kwargs):
        self.calls.append(("run", tuple(command), kwargs))
        if not kwargs.get("dry_run"):
            pytest.fail("Unsafe paths reached execution")
        return None


@pytest.fixture
def request_files(tmp_path):
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\n" + "A" * 1000 + "\n", encoding="utf-8")
    fai = Path(f"{fasta}.fai")
    fai.write_text("chr1\t1000\t6\t1000\t1001\n", encoding="utf-8")
    source = tmp_path / "caller.hg38_multianno.vcf"
    source.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t100\tv1\tA\tG\t30\tPASS\t.\tGT\t0/1\n", encoding="utf-8")
    tbi, csi = Path(f"{source}.tbi"), Path(f"{source}.csi")
    tbi.write_bytes(b"SOURCE INDEX TBI")
    csi.write_bytes(b"SOURCE INDEX CSI")
    database = tmp_path / "database"
    cache = tmp_path / "cache"
    database.mkdir()
    cache.mkdir()
    (database / "nested").mkdir()
    db_file = database / "nested/hg38_refGene.txt"
    db_file.write_text("immutable gene data", encoding="utf-8")
    cache_file = cache / "species.cache.gz"
    cache_file.write_bytes(b"immutable cache data")
    genome = ReferenceGenome.from_fasta(fasta, build="GRCh38")
    item = AnnotationInput("HG002", VariantCategory.SMALL, source, "DeepVariant", genome)
    annovar = AnnovarRequest(item, database, "snapshot1", ("refGene",), ("g",),
                             tmp_path / "outputs/result", "2020-06-08")
    vep = VepRequest(item, cache, "116", "homo_sapiens", "GRCh38", tmp_path / "outputs/result.vep.tsv")
    return {"annovar": annovar, "vep": vep, "source": source, "tbi": tbi, "csi": csi,
            "fasta": fasta, "fai": fai, "database": database, "cache": cache,
            "db_file": db_file, "cache_file": cache_file}


def snapshot(request_files):
    return {str(request_files[name]): request_files[name].read_bytes()
            for name in ("source", "tbi", "csi", "fasta", "fai", "db_file", "cache_file")}


def check_unchanged(request_files, before):
    assert snapshot(request_files) == before


@pytest.mark.parametrize("overwrite", [False, True])
def test_annovar_native_vcf_must_not_be_source_even_with_overwrite(tmp_path, request_files, overwrite):
    request = replace(request_files["annovar"], output_prefix=tmp_path / "caller", overwrite=overwrite)
    assert request.output_vcf == request_files["source"]
    runner = ProbeRunner()
    before = snapshot(request_files)
    with pytest.raises(OutputValidationError, match="ANNOTATION_OUTPUT_PROTECTION"):
        AnnovarWrapper(runner=runner).run(request)
    assert runner.calls == []
    check_unchanged(request_files, before)


@pytest.mark.parametrize("overwrite", [False, True])
@pytest.mark.parametrize("protected", ["source", "tbi", "csi", "fasta", "fai", "cache_file"])
def test_vep_output_direct_alias_fails_before_unlink_or_probe(request_files, protected, overwrite):
    request = replace(request_files["vep"], output_tsv=request_files[protected], overwrite=overwrite)
    runner = ProbeRunner()
    before = snapshot(request_files)
    with pytest.raises(OutputValidationError, match="ANNOTATION_OUTPUT_PROTECTION"):
        VepWrapper(runner=runner).run(request)
    assert runner.calls == []
    check_unchanged(request_files, before)


@pytest.mark.parametrize("tool", ["annovar", "vep"])
@pytest.mark.parametrize("log_kind", ["stdout_path", "stderr_path"])
@pytest.mark.parametrize("protected", ["source", "tbi", "fasta", "fai"])
def test_logs_cannot_truncate_immutable_inputs(request_files, tool, log_kind, protected):
    request = replace(request_files[tool], overwrite=True)
    wrapper = AnnovarWrapper if tool == "annovar" else VepWrapper
    runner = ProbeRunner()
    before = snapshot(request_files)
    with pytest.raises(OutputValidationError, match="ANNOTATION_OUTPUT_PROTECTION"):
        wrapper(runner=runner).run(request, **{log_kind: request_files[protected]})
    assert runner.calls == []
    assert not request.output_tsv.parent.exists()
    check_unchanged(request_files, before)


@pytest.mark.parametrize("tool,resource", [("annovar", "database"), ("vep", "cache")])
@pytest.mark.parametrize("overwrite", [False, True])
def test_new_output_cannot_be_written_inside_database_or_cache(request_files, tool, resource, overwrite):
    directory = request_files[resource]
    before_names = sorted(str(x.relative_to(directory)) for x in directory.rglob("*"))
    if tool == "annovar":
        request = replace(request_files[tool], output_prefix=directory / "new_result", overwrite=overwrite)
        wrapper = AnnovarWrapper
    else:
        request = replace(request_files[tool], output_tsv=directory / "new/result.tsv", overwrite=overwrite)
        wrapper = VepWrapper
    runner = ProbeRunner()
    with pytest.raises(OutputValidationError, match="ANNOTATION_OUTPUT_PROTECTION"):
        wrapper(runner=runner).run(request)
    assert runner.calls == []
    assert sorted(str(x.relative_to(directory)) for x in directory.rglob("*")) == before_names


@pytest.mark.parametrize("tool,resource", [("annovar", "database"), ("vep", "cache")])
def test_log_target_cannot_be_resource_root_or_ancestor(tmp_path, request_files, tool, resource):
    wrapper = AnnovarWrapper if tool == "annovar" else VepWrapper
    for path in (request_files[resource], tmp_path):
        runner = ProbeRunner()
        with pytest.raises(OutputValidationError, match="ANNOTATION_OUTPUT_PROTECTION"):
            wrapper(runner=runner).run(request_files[tool], stderr_path=path)
        assert runner.calls == []


def link_or_skip(source: Path, destination: Path, kind: str):
    try:
        if kind == "hardlink":
            os.link(source, destination)
        else:
            destination.symlink_to(source, target_is_directory=source.is_dir())
    except (OSError, NotImplementedError) as error:
        pytest.skip(f"Filesystem {kind} privilege unavailable: {error}")


@pytest.mark.parametrize("tool,protected", [("annovar", "db_file"), ("vep", "cache_file")])
@pytest.mark.parametrize("kind", ["hardlink", "symlink"])
@pytest.mark.parametrize("overwrite", [False, True])
def test_external_output_alias_of_internal_resource_rejected(tmp_path, request_files, tool, protected, kind, overwrite):
    request = replace(request_files[tool], overwrite=overwrite)
    request.output_tsv.parent.mkdir()
    link_or_skip(request_files[protected], request.output_tsv, kind)
    runner = ProbeRunner()
    before = snapshot(request_files)
    wrapper = AnnovarWrapper if tool == "annovar" else VepWrapper
    with pytest.raises(OutputValidationError, match="ANNOTATION_OUTPUT_PROTECTION"):
        wrapper(runner=runner).run(request)
    assert runner.calls == []
    check_unchanged(request_files, before)


@pytest.mark.parametrize("protected", ["source", "tbi", "csi", "fasta", "fai"])
@pytest.mark.parametrize("kind", ["hardlink", "symlink"])
def test_vep_external_log_alias_of_input_rejected(tmp_path, request_files, protected, kind):
    log = tmp_path / "aliased.log"
    link_or_skip(request_files[protected], log, kind)
    runner = ProbeRunner()
    before = snapshot(request_files)
    with pytest.raises(OutputValidationError, match="ANNOTATION_OUTPUT_PROTECTION"):
        VepWrapper(runner=runner).run(request_files["vep"], stdout_path=log)
    assert runner.calls == []
    check_unchanged(request_files, before)


@pytest.mark.parametrize("tool,resource", [("annovar", "database"), ("vep", "cache")])
def test_nested_resource_directory_symlink_protects_external_target(tmp_path, request_files, tool, resource):
    outside = tmp_path / "linked_resource"
    outside.mkdir()
    link_or_skip(outside, request_files[resource] / "link", "symlink")
    wrapper = AnnovarWrapper if tool == "annovar" else VepWrapper
    runner = ProbeRunner()
    with pytest.raises(OutputValidationError, match="ANNOTATION_OUTPUT_PROTECTION"):
        wrapper(runner=runner).run(request_files[tool], stderr_path=outside / "new.log")
    assert runner.calls == []
    assert not (outside / "new.log").exists()


def test_annovar_existing_helper_alias_is_protected(tmp_path, request_files):
    request = request_files["annovar"]
    request.output_prefix.parent.mkdir()
    helper = Path(f"{request.output_prefix}.refGene.variant_function")
    link_or_skip(request_files["db_file"], helper, "hardlink")
    runner = ProbeRunner()
    before = snapshot(request_files)
    with pytest.raises(OutputValidationError, match="ANNOTATION_OUTPUT_PROTECTION"):
        AnnovarWrapper(runner=runner).run(request)
    assert runner.calls == []
    check_unchanged(request_files, before)


def test_owned_annovar_helpers_and_sibling_outputs_are_allowed(tmp_path, request_files):
    request = request_files["annovar"]
    request.output_prefix.parent.mkdir()
    helper = Path(f"{request.output_prefix}.refGene.variant_function")
    helper.write_text("own old intermediate", encoding="utf-8")
    runner = ProbeRunner()
    before = snapshot(request_files)
    result = AnnovarWrapper(runner=runner).run(request, dry_run=True,
                                             stdout_path=tmp_path / "newlogs/stdout.log",
                                             stderr_path=tmp_path / "newlogs/stderr.log")
    assert result.status.value == "planned"
    assert len(runner.calls) == 1
    assert helper.read_text() == "own old intermediate"
    assert not (tmp_path / "newlogs").exists()
    check_unchanged(request_files, before)


def test_output_and_log_must_not_collide(request_files):
    request = request_files["vep"]
    runner = ProbeRunner()
    with pytest.raises(OutputValidationError, match="aliases collide"):
        VepWrapper(runner=runner).run(request, stderr_path=request.output_tsv)
    assert runner.calls == []


@pytest.mark.parametrize("tool,protected", [("annovar", "db_file"), ("vep", "cache_file")])
@pytest.mark.parametrize("log_kind", ["stdout_path", "stderr_path"])
def test_logs_inside_or_hardlinked_to_database_are_protected(tmp_path, request_files, tool, protected, log_kind):
    alias = tmp_path / "database_alias.log"
    link_or_skip(request_files[protected], alias, "hardlink")
    wrapper = AnnovarWrapper if tool == "annovar" else VepWrapper
    for target in (request_files[protected], alias):
        runner = ProbeRunner()
        before = snapshot(request_files)
        with pytest.raises(OutputValidationError, match="ANNOTATION_OUTPUT_PROTECTION"):
            wrapper(runner=runner).run(request_files[tool], **{log_kind: target})
        assert runner.calls == []
        check_unchanged(request_files, before)


def test_unreadable_database_tree_fails_closed_before_execution(monkeypatch, request_files):
    from hifivar import annotation_paths

    def unreadable(_):
        raise PermissionError("simulated resource metadata denial")

    monkeypatch.setattr(annotation_paths.os, "scandir", unreadable)
    runner = ProbeRunner()
    with pytest.raises(OutputValidationError, match="cannot inspect resource tree"):
        VepWrapper(runner=runner).run(request_files["vep"])
    assert runner.calls == []
    assert not request_files["vep"].output_tsv.parent.exists()


@pytest.mark.parametrize("tool", ["annovar", "vep"])
@pytest.mark.parametrize("log_kind", ["stdout_path", "stderr_path"])
def test_program_file_cannot_be_a_log(tmp_path, request_files, tool, log_kind):
    program = tmp_path / "private_program.pl"
    program.write_text("# immutable executable\n", encoding="utf-8")
    runner = ProbeRunner()
    wrapper = AnnovarWrapper if tool == "annovar" else VepWrapper
    with pytest.raises(OutputValidationError, match="ANNOTATION_OUTPUT_PROTECTION"):
        wrapper(executable=str(program), runner=runner).run(request_files[tool], **{log_kind: program})
    assert runner.calls == [] and program.read_text() == "# immutable executable\n"


@pytest.mark.parametrize("tool", ["annovar", "vep"])
def test_program_file_hardlinked_to_overwritten_output_is_protected(tmp_path, request_files, tool):
    program = tmp_path / "private_program.pl"
    program.write_text("# immutable executable\n", encoding="utf-8")
    request = replace(request_files[tool], overwrite=True)
    request.output_tsv.parent.mkdir()
    link_or_skip(program, request.output_tsv, "hardlink")
    runner = ProbeRunner()
    wrapper = AnnovarWrapper if tool == "annovar" else VepWrapper
    with pytest.raises(OutputValidationError, match="ANNOTATION_OUTPUT_PROTECTION"):
        wrapper(executable=str(program), runner=runner).run(request)
    assert runner.calls == [] and program.read_text() == "# immutable executable\n"


def test_annovar_sibling_perl_helpers_protected_without_system_tree_walk(tmp_path, request_files):
    program = tmp_path / "software/table_annovar.pl"
    program.parent.mkdir()
    program.write_text("# table script\n", encoding="utf-8")
    helper = program.parent / "annotate_variation.pl"
    helper.write_text("# immutable helper\n", encoding="utf-8")
    runner = ProbeRunner()
    with pytest.raises(OutputValidationError, match="ANNOTATION_OUTPUT_PROTECTION"):
        AnnovarWrapper(executable=str(program), runner=runner).run(request_files["annovar"], stdout_path=helper)
    assert runner.calls == [] and helper.read_text() == "# immutable helper\n"


def test_bare_program_path_is_discovered_read_only(tmp_path, monkeypatch, request_files):
    from hifivar import annotation_paths

    program = tmp_path / "real_vep"
    program.write_text("# existing program\n", encoding="utf-8")
    monkeypatch.setattr(annotation_paths.shutil, "which", lambda name: str(program) if name == "vep" else None)
    runner = ProbeRunner()
    with pytest.raises(OutputValidationError, match="ANNOTATION_OUTPUT_PROTECTION"):
        VepWrapper(runner=runner).run(request_files["vep"], stderr_path=program)
    assert runner.calls == []
