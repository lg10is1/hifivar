from pathlib import Path

from hifivar.annotation import AnnotationInput, VariantCategory
from hifivar.annovar import AnnovarRequest, AnnovarWrapper
from hifivar.command import CommandResult
from hifivar.functional import (
    FunctionalPrediction,
    FunctionalPrioritizationRequest,
    FunctionalPrioritizationResult,
    PrioritizedVariant,
)
from hifivar.phase11 import AnnotationJob, run_phase11
from hifivar.reference import ReferenceGenome
from hifivar.vep import VepRequest, VepWrapper


class FakeAnnotationRunner:
    def require_executable(self, executable):
        return Path(executable)

    def run(self, command, **kwargs):
        args = tuple(str(item) for item in command)
        if "--help" in args:
            return CommandResult(args, 0, "Versions:\nensembl-vep : 115", "", 0.1, None, True)
        if "-out" in args:
            prefix = Path(args[args.index("-out") + 1])
            build = args[args.index("-buildver") + 1]
            Path(f"{prefix}.{build}_multianno.txt").write_text("Chr\tStart\nchr1\t10\n", encoding="utf-8")
            Path(f"{prefix}.{build}_multianno.vcf").write_text(
                "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t10\tv1\tA\tT\t.\tPASS\tANNOVAR_DATE=20260823\n",
                encoding="utf-8",
            )
        if "--output_file" in args:
            Path(args[args.index("--output_file") + 1]).write_text(
                "#Uploaded_variation\tConsequence\nv1\tfeature\n", encoding="utf-8"
            )
        return CommandResult(args, 0, "", "", 0.2, None, True)


class FakeAlphaGenome:
    def predict(self, request):
        return FunctionalPrioritizationResult(
            request,
            tuple(FunctionalPrediction(item, {"RNA_SEQ": 0.5}) for item in request.selected_variants),
            "mock-cloud-1",
        )


def test_phase11_small_sv_tr_annotation_and_prioritization(tmp_path):
    reference_path = tmp_path / "reference.fa"
    reference_path.write_text(">chr1\n" + "A" * 1000 + "\n", encoding="utf-8")
    Path(f"{reference_path}.fai").write_text("chr1\t1000\t6\t1000\t1001\n", encoding="utf-8")
    reference = ReferenceGenome.from_fasta(reference_path, build="GRCh38")
    database = tmp_path / "humandb"
    cache = tmp_path / "vep-cache"
    database.mkdir()
    cache.mkdir()
    output = tmp_path / "annotation"
    originals = {}
    jobs = []
    for category, caller in (
        (VariantCategory.SMALL, "deepvariant"),
        (VariantCategory.SV, "jasmine"),
        (VariantCategory.TR, "trgt"),
    ):
        vcf = tmp_path / f"S1.{category.value}.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t10\tv1\tA\tT\t.\tPASS\t.\n",
            encoding="utf-8",
        )
        originals[vcf] = vcf.read_bytes()
        item = AnnotationInput("S1", category, vcf, caller, reference, ("v1",))
        root = output / category.value
        annovar = AnnovarRequest(
            item, database, "db-2026", ("refGeneWithVer",), ("g",),
            root / f"S1.{category.value}.annovar", "2025-03-02",
        )
        vep = VepRequest(
            item, cache, "115", "homo_sapiens", "GRCh38",
            root / f"S1.{category.value}.vep.tsv", 2,
        )
        jobs.append(AnnotationJob(item, annovar, vep))

    selection_source = tmp_path / "explicit.annotation.tsv"
    selection_source.write_text("variant\timpact\nv1\tmodifier\n", encoding="utf-8")
    prioritized = PrioritizedVariant(
        "S1", "v1", VariantCategory.SMALL, "chr1", 10, "A", "T",
        selection_source, "explicit candidate",
    )
    functional = FunctionalPrioritizationRequest(
        (prioritized,), "AlphaGenome", "API_VERSION_PENDING_LINUX_VERIFICATION", ("RNA_SEQ",),
    )
    runner = FakeAnnotationRunner()
    report = run_phase11(
        jobs,
        annovar_wrapper=AnnovarWrapper(runner=runner),
        vep_wrapper=VepWrapper(runner=runner),
        functional_request=functional,
        functional_backend=FakeAlphaGenome(),
        log_directory=output / "logs",
    )
    assert len(report.annotation_results) == 6
    assert {result.input.variant_category for result in report.annotation_results} == set(VariantCategory)
    assert report.functional_result is not None
    assert report.to_dict()["scientific_policy"]["variant_classes_remain_separate"] is True
    assert report.to_dict()["scientific_policy"]["functional_impact_is_variant_call_confidence"] is False
    assert all(path.read_bytes() == contents for path, contents in originals.items())
    report.write_json(output / "phase11.provenance.json")
    report.write_yaml(output / "phase11.provenance.yaml")
