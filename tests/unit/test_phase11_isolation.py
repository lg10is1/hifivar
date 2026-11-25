"""Expected backend failures must not erase independent successful evidence."""

from pathlib import Path
from unittest.mock import Mock

import pytest

from hifivar.annotation import AnnotationInput, AnnotationResult, AnnotationRunStatus, AnnotationSource, VariantCategory
from hifivar.annovar import AnnovarRequest
from hifivar.exceptions import CommandExecutionError, InputValidationError
from hifivar.phase11 import AnnotationJob, run_phase11
from hifivar.reference import ReferenceGenome
from hifivar.vep import VepRequest


@pytest.fixture
def job(tmp_path):
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\n" + "A" * 100 + "\n", encoding="utf-8")
    Path(str(fasta) + ".fai").write_text("chr1\t100\t6\t100\t101\n", encoding="utf-8")
    source = tmp_path / "source.vcf"
    source.write_text(
        "##fileformat=VCFv4.3\n##reference=GRCh38\n##contig=<ID=chr1,length=100>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t10\tv1\tA\tT\t.\tPASS\t.\tGT\t0/1\n", encoding="utf-8",
    )
    database = tmp_path / "database"
    database.mkdir()
    cache = tmp_path / "cache"
    cache.mkdir()
    item = AnnotationInput("S1", VariantCategory.SMALL, source, "deepvariant",
                           ReferenceGenome.from_fasta(fasta, build="GRCh38"))
    return AnnotationJob(
        item,
        AnnovarRequest(item, database, "snapshot", ("refGene",), ("g",),
                       tmp_path / "annovar", "2020-06-08"),
        VepRequest(item, cache, "116", "homo_sapiens", "GRCh38", tmp_path / "vep.tsv"),
    )


def success(job, backend):
    return AnnotationResult(job.input, AnnotationSource(backend), AnnotationRunStatus.COMPLETED,
                            (backend, "source.vcf"), "1.0")


def test_default_remains_fail_fast(job):
    annovar = Mock()
    annovar.run.side_effect = CommandExecutionError("failed with exit 2")
    vep = Mock()
    with pytest.raises(CommandExecutionError):
        run_phase11([job], annovar_wrapper=annovar, vep_wrapper=vep)
    vep.run.assert_not_called()


def test_isolated_failure_keeps_later_backend_and_source(job, tmp_path):
    before = job.input.source_vcf.read_bytes()
    annovar = Mock()
    annovar.run.side_effect = CommandExecutionError("failed with exit 2")
    vep = Mock()
    marker = tmp_path / "vep.completed.tsv"

    def run(*args, **kwargs):
        marker.write_text("independent result\n", encoding="utf-8")
        return success(job, "vep")

    vep.run.side_effect = run
    report = run_phase11([job], annovar_wrapper=annovar, vep_wrapper=vep,
                         isolate_failures=True, log_directory=tmp_path / "logs")
    payload = report.to_dict()
    assert payload["status"] == "partial"
    assert payload["backend_failures"][0]["backend"] == "annovar"
    assert payload["backend_failures"][0]["error_type"] == "CommandExecutionError"
    assert payload["failure_isolation_enabled"] is True
    assert payload["tool_annotation_scope"] == "whole_source_vcf"
    assert len(report.annotation_results) == 1
    assert marker.read_text(encoding="utf-8") == "independent result\n"
    assert job.input.source_vcf.read_bytes() == before
    report.write_json(tmp_path / "report.json")
    report.write_yaml(tmp_path / "report.yaml")


def test_isolation_preserves_earlier_success(job):
    annovar, vep = Mock(), Mock()
    annovar.run.return_value = success(job, "annovar")
    vep.run.side_effect = CommandExecutionError("vep failed")
    report = run_phase11([job], annovar_wrapper=annovar, vep_wrapper=vep, isolate_failures=True)
    assert report.to_dict()["status"] == "partial"
    assert report.annotation_results[0].source is AnnotationSource.ANNOVAR
    assert report.backend_failures[0].backend == "vep"


def test_all_failed_is_not_complete(job):
    annovar, vep = Mock(), Mock()
    annovar.run.side_effect = CommandExecutionError("A failed")
    vep.run.side_effect = CommandExecutionError("V failed")
    report = run_phase11([job], annovar_wrapper=annovar, vep_wrapper=vep, isolate_failures=True)
    assert report.to_dict()["status"] == "failed"
    assert len(report.backend_failures) == 2
    assert not report.annotation_results


def test_incomplete_annotation_does_not_trigger_functional_prediction(job):
    annovar, vep, functional = Mock(), Mock(), Mock()
    annovar.run.side_effect = CommandExecutionError("A failed")
    vep.run.return_value = success(job, "vep")
    report = run_phase11([job], annovar_wrapper=annovar, vep_wrapper=vep,
                         functional_request=Mock(), functional_backend=functional,
                         isolate_failures=True)
    functional.predict.assert_not_called()
    assert report.functional_result is None
    assert report.to_dict()["functional_status"] == "not_run_annotation_failure"
    assert report.to_dict()["status"] == "partial"


def test_unexpected_programming_error_propagates(job):
    annovar, vep = Mock(), Mock()
    annovar.run.side_effect = TypeError("programming bug")
    with pytest.raises(TypeError, match="programming bug"):
        run_phase11([job], annovar_wrapper=annovar, vep_wrapper=vep, isolate_failures=True)
    vep.run.assert_not_called()


@pytest.mark.parametrize("isolation", [False, True])
def test_dry_run_remains_planned(job, isolation):
    annovar, vep = Mock(), Mock()
    annovar.run.return_value = AnnotationResult(job.input, AnnotationSource.ANNOVAR,
                                              AnnotationRunStatus.PLANNED, ("annovar",))
    vep.run.return_value = AnnotationResult(job.input, AnnotationSource.VEP,
                                          AnnotationRunStatus.PLANNED, ("vep",))
    report = run_phase11([job], annovar_wrapper=annovar, vep_wrapper=vep,
                         dry_run=True, isolate_failures=isolation)
    assert report.to_dict()["status"] == "planned"
    assert report.to_dict()["functional_status"] == "disabled"
    assert not job.vep.output_tsv.exists()


def test_isolation_flag_is_strict_boolean(job):
    with pytest.raises(InputValidationError, match="boolean"):
        run_phase11([job], isolate_failures="yes")
