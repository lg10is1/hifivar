from pathlib import Path

import pytest

from hifivar.command import CommandResult
from hifivar.config import validate_config
from hifivar.exceptions import ConfigurationError
from hifivar.exceptions import CommandExecutionError, InputValidationError, OutputValidationError
from hifivar.igv import IgvRunStatus, IgvWrapper
from hifivar.phase10 import Phase10Status, run_phase10
from hifivar.review import (
    EvidenceClass,
    ReviewManifest,
    ReviewResult,
    ReviewStatus,
    ReviewTarget,
    VariantClass,
)


def inputs(tmp_path: Path):
    reference = tmp_path / "参考 genome.fa"
    reference.write_text(">chr1\n" + "A" * 1000 + "\n>chr2\n" + "C" * 1000 + "\n", encoding="utf-8")
    Path(f"{reference}.fai").write_text(
        "chr1\t1000\t6\t1000\t1001\nchr2\t1000\t1013\t1000\t1001\n",
        encoding="utf-8",
    )
    bam = tmp_path / "S1.bam"
    bam.write_bytes(b"BAM")
    Path(f"{bam}.bai").write_bytes(b"BAI")
    vcf = tmp_path / "source.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\tv1\tA\tT\t.\tPASS\t.\tGT\t0/1\n",
        encoding="utf-8",
    )
    return reference, bam, vcf


def target(tmp_path: Path, variant_class: VariantClass, **overrides) -> ReviewTarget:
    reference, bam, vcf = inputs(tmp_path)
    values = {
        "review_id": f"R-{variant_class.value}",
        "sample_id": "S1",
        "variant_id": "source-v1",
        "variant_class": variant_class,
        "contig": "chr1",
        "start": 100,
        "end": 105,
        "source_vcf": vcf,
        "source_caller": "explicit",
        "evidence_class": EvidenceClass.EXPLICIT,
        "alignment_path": bam,
        "reference_fasta": reference,
        "output_directory": tmp_path / "审阅 output",
        "flank_bp": 10,
    }
    values.update(overrides)
    return ReviewTarget(**values)


@pytest.mark.parametrize(
    ("variant_class", "expected_end"),
    ((VariantClass.SNV, 110), (VariantClass.INS, 110), (VariantClass.INDEL, 115),
     (VariantClass.DEL, 115), (VariantClass.DUP, 115), (VariantClass.INV, 115),
     (VariantClass.TR, 115)),
)
def test_variant_centered_window_rules(tmp_path, variant_class, expected_end):
    item = target(tmp_path, variant_class)
    locus = item.loci[0]
    assert locus.window_start == 90 and locus.window_end == expected_end
    if variant_class is VariantClass.INS:
        assert locus.variant_end == 100


def test_bnd_has_two_independent_loci(tmp_path):
    item = target(
        tmp_path,
        VariantClass.BND,
        mate_contig="chr2",
        mate_position=500,
    )
    assert [locus.igv_locus for locus in item.loci] == ["chr1:90-110", "chr2:490-510"]
    assert [locus.variant_end for locus in item.loci] == [100, 500]
    assert [path.name for path in item.screenshot_paths] == ["R-BND.locus01.png", "R-BND.locus02.png"]


def test_trgt_visualization_is_metadata_only(tmp_path):
    plot = tmp_path / "trgt.plot.svg"
    plot.write_text("<svg/>\n", encoding="utf-8")
    item = target(
        tmp_path,
        VariantClass.TR,
        evidence_class=EvidenceClass.TANDEM_REPEAT,
        source_caller="trgt",
        trgt_visualization_path=plot,
    )
    assert item.to_dict()["trgt_visualization"] == str(plot)
    with pytest.raises(InputValidationError, match="only for TR"):
        target(tmp_path, VariantClass.SNV, trgt_visualization_path=plot)


class FakeIgvRunner:
    def __init__(self):
        self.calls = []

    def require_executable(self, executable):
        return Path(executable)

    def run(self, command, **kwargs):
        args = tuple(str(value) for value in command)
        self.calls.append((args, kwargs))
        if "--version" in args:
            return CommandResult(args, 0, "IGV Version 2.19.8", "", 0.1, None, True)
        if kwargs.get("dry_run"):
            return CommandResult(args, 0, None, None, 0.0, None, False)
        batch = Path(args[-1])
        snapshot_root = None
        for line in batch.read_text(encoding="utf-8").splitlines():
            if line.startswith("snapshotDirectory "):
                snapshot_root = Path(line.split(" ", 1)[1].strip('"'))
                snapshot_root.mkdir(parents=True, exist_ok=True)
            elif line.startswith("snapshot "):
                assert snapshot_root is not None
                (snapshot_root / line.split(" ", 1)[1]).write_bytes(b"PNG")
        return CommandResult(args, 0, "", "", 0.2, None, True)


class FailedIgvRunner(FakeIgvRunner):
    def run(self, command, **kwargs):
        args = tuple(str(value) for value in command)
        if "--batch" in args:
            raise CommandExecutionError("IGV failed")
        return super().run(command, **kwargs)


def test_batch_generation_unicode_paths_and_evidence(tmp_path):
    item = target(tmp_path, VariantClass.SNV)
    runner = FakeIgvRunner()
    wrapper = IgvWrapper(runner=runner)
    report = run_phase10((item,), output_directory=item.output_directory, igv_wrapper=wrapper)
    assert report.status is Phase10Status.COMPLETED
    assert report.igv.status is IgvRunStatus.COMPLETED
    text = (item.output_directory / "review.igv.batch").read_text(encoding="utf-8")
    assert 'genome "' in text and 'load "' in text
    assert "goto chr1:90-110" in text
    assert "snapshot R-SNV.locus01.png" in text
    assert item.screenshot_paths[0].read_bytes() == b"PNG"
    assert report.manifest.results[0].evidence.target.source_vcf == item.source_vcf


def test_manual_status_is_not_truth_and_uncertain_is_retained(tmp_path):
    item = target(tmp_path, VariantClass.DEL)
    wrapper = IgvWrapper(runner=FakeIgvRunner())
    report = run_phase10((item,), output_directory=item.output_directory, igv_wrapper=wrapper)
    support = report.manifest.results[0].with_manual_review(ReviewStatus.SUPPORT, "read support")
    uncertain = support.with_manual_review(ReviewStatus.UNCERTAIN, "ambiguous")
    manifest = ReviewManifest((uncertain,), created_at="2026-08-23T00:00:00Z")
    assert support.is_truth is False
    assert uncertain.status is ReviewStatus.UNCERTAIN
    assert [row["status"] for row in manifest.to_dict()["results"]] == ["UNCERTAIN"]
    assert manifest.to_dict()["scientific_policy"]["manual_status_is_truth"] is False


def test_manifest_formats_and_overwrite_protection(tmp_path):
    item = target(tmp_path, VariantClass.INDEL)
    report = run_phase10((item,), output_directory=item.output_directory, igv_wrapper=IgvWrapper(runner=FakeIgvRunner()))
    assert (item.output_directory / "review_manifest.json").exists()
    assert "NOT_REVIEWED" in (item.output_directory / "review_manifest.tsv").read_text(encoding="utf-8")
    assert "clinical_classification" not in (item.output_directory / "review_manifest.yaml").read_text(encoding="utf-8")
    with pytest.raises(OutputValidationError, match="already exists"):
        run_phase10((item,), output_directory=item.output_directory, igv_wrapper=IgvWrapper(runner=FakeIgvRunner()))


def test_dry_run_and_empty_target_list(tmp_path):
    item = target(tmp_path, VariantClass.SNV)
    wrapper = IgvWrapper(runner=FakeIgvRunner())
    planned = run_phase10((item,), output_directory=item.output_directory, igv_wrapper=wrapper, dry_run=True)
    assert planned.status is Phase10Status.PLANNED
    assert not item.output_directory.exists()
    empty_output = tmp_path / "empty"
    empty = run_phase10((), output_directory=empty_output, igv_wrapper=wrapper)
    assert empty.status is Phase10Status.COMPLETED and not empty.manifest.results
    assert (empty_output / "review_manifest.tsv").read_text().count("\n") == 1
    assert (empty_output / "screenshots").is_dir()


def test_missing_bam_or_reference_is_rejected(tmp_path):
    reference, bam, vcf = inputs(tmp_path)
    bam.unlink()
    with pytest.raises(InputValidationError, match="missing"):
        ReviewTarget("R1", "S1", "v1", VariantClass.SNV, "chr1", 10, 10, vcf,
                     "deepvariant", EvidenceClass.SMALL_VARIANT, bam, reference, tmp_path / "out")
    bam.write_bytes(b"BAM")
    reference.unlink()
    with pytest.raises(InputValidationError, match="missing"):
        ReviewTarget("R2", "S1", "v2", VariantClass.SNV, "chr1", 10, 10, vcf,
                     "deepvariant", EvidenceClass.SMALL_VARIANT, tmp_path / "S1.bam", reference, tmp_path / "out")


def test_external_failure_does_not_modify_source_artifacts(tmp_path):
    item = target(tmp_path, VariantClass.SNV)
    before = item.source_vcf.read_bytes()
    with pytest.raises(CommandExecutionError, match="IGV failed"):
        run_phase10(
            (item,), output_directory=item.output_directory,
            igv_wrapper=IgvWrapper(runner=FailedIgvRunner()),
        )
    assert item.source_vcf.read_bytes() == before
    assert not item.screenshot_paths[0].exists()


def test_review_configuration_is_explicit_and_typed():
    validate_config({
        "review": {
            "enabled": True, "selection_file": "selected.tsv",
            "igv_executable": "igv.sh", "flank_bp": 0,
            "threads": 1, "memory_mb": 8000, "runtime_minutes": 240,
            "overwrite": False,
        }
    })
    with pytest.raises(ConfigurationError, match="selection_file"):
        validate_config({"review": {"enabled": True, "selection_file": None}})
    with pytest.raises(ConfigurationError, match="flank_bp"):
        validate_config({"review": {"flank_bp": -1}})
