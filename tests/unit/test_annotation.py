from pathlib import Path

import pytest

from hifivar.annotation import (
    AnnotationInput,
    AnnotationRunStatus,
    AnnotationSource,
    RegionCategory,
    RegionDatabase,
    VariantCategory,
    VariantLocus,
    annotate_region_overlaps,
    read_annotation_inputs,
    read_selected_variant_loci,
)
from hifivar.annovar import AnnovarRequest, AnnovarWrapper
from hifivar.command import CommandResult
from hifivar.config import validate_config
from hifivar.exceptions import (
    CommandExecutionError,
    ConfigurationError,
    InputValidationError,
    OutputValidationError,
)
from hifivar.functional import (
    FunctionalPrediction,
    FunctionalPrioritizationRequest,
    FunctionalPrioritizationResult,
    PrioritizedVariant,
    read_functional_selection,
    run_functional_prioritization,
)
from hifivar.phase11 import AnnotationJob, run_phase11
from hifivar.reference import ReferenceGenome
from hifivar.vep import VepRequest, VepWrapper


def files(tmp_path: Path):
    tmp_path.mkdir(parents=True, exist_ok=True)
    reference = tmp_path / "参考 genome.fa"
    reference.write_text(">chr1\n" + "A" * 1000 + "\n", encoding="utf-8")
    Path(f"{reference}.fai").write_text("chr1\t1000\t6\t1000\t1001\n", encoding="utf-8")
    genome = ReferenceGenome.from_fasta(reference, build="GRCh38")
    vcf = tmp_path / "variants.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\tsv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=150\n",
        encoding="utf-8",
    )
    database = tmp_path / "数据库"
    database.mkdir()
    cache = tmp_path / "vep cache"
    cache.mkdir()
    return genome, vcf, database, cache


def annotation_input(tmp_path: Path, category=VariantCategory.SMALL, ids=()):
    genome, vcf, database, cache = files(tmp_path)
    return AnnotationInput("S1", category, vcf, "source-caller", genome, tuple(ids)), database, cache


class FakeToolRunner:
    def __init__(self, *, fail=False):
        self.calls = []
        self.fail = fail

    def require_executable(self, executable):
        return Path(executable)

    def run(self, command, **kwargs):
        args = tuple(str(item) for item in command)
        self.calls.append((args, kwargs))
        if kwargs.get("dry_run"):
            return CommandResult(args, 0, None, None, 0.0, None, False)
        if "--help" in args:
            return CommandResult(args, 0, "Versions:\nensembl-vep : 115", "", 0.1, None, True)
        if self.fail:
            raise CommandExecutionError("annotation tool failed")
        if "-out" in args:
            prefix = Path(args[args.index("-out") + 1])
            build = args[args.index("-buildver") + 1]
            Path(f"{prefix}.{build}_multianno.txt").write_text("Chr\tStart\tFunc\nchr1\t100\texonic\n", encoding="utf-8")
            Path(f"{prefix}.{build}_multianno.vcf").write_text(
                "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t100\tsv1\tN\t<DEL>\t.\tPASS\tANNOVAR_DATE=20260823\n",
                encoding="utf-8",
            )
        if "--output_file" in args:
            Path(args[args.index("--output_file") + 1]).write_text(
                "#Uploaded_variation\tLocation\tConsequence\nsv1\tchr1:100\tfeature\n",
                encoding="utf-8",
            )
        return CommandResult(args, 0, "", "", 0.2, None, True)


def requests(tmp_path: Path, category=VariantCategory.SMALL):
    item, database, cache = annotation_input(tmp_path, category)
    annovar = AnnovarRequest(
        item, database, "refGeneWithVer-2026", ("refGeneWithVer",), ("g",),
        tmp_path / "结果" / "S1.small.annovar", "2025-03-02", False,
    )
    vep = VepRequest(
        item, cache, "115", "homo_sapiens", "GRCh38",
        tmp_path / "结果" / "S1.small.vep.tsv", 4, False,
    )
    return item, annovar, vep


def test_annovar_command_and_provenance(tmp_path):
    item, request, _ = requests(tmp_path)
    runner = FakeToolRunner()
    wrapper = AnnovarWrapper(runner=runner)
    command = wrapper.build_command(request)
    assert command[:3] == ("table_annovar.pl", str(item.source_vcf.absolute()), str(request.database_root.absolute()))
    assert command[command.index("-protocol") + 1] == "refGeneWithVer"
    assert command[command.index("-operation") + 1] == "g"
    assert "-vcfinput" in command and "-remove" not in command
    before = item.source_vcf.read_bytes()
    result = wrapper.run(request)
    assert result.status is AnnotationRunStatus.COMPLETED
    assert result.artifact is not None
    assert result.artifact.source is AnnotationSource.ANNOVAR
    assert result.artifact.databases[0].version == "refGeneWithVer-2026"
    assert item.source_vcf.read_bytes() == before


def test_vep_offline_command_version_and_provenance(tmp_path):
    item, _, request = requests(tmp_path)
    runner = FakeToolRunner()
    wrapper = VepWrapper(runner=runner)
    command = wrapper.build_command(request)
    assert "--offline" in command and "--cache" in command and "--tab" in command
    assert command[command.index("--fasta") + 1] == str(item.reference.fasta.absolute())
    result = wrapper.run(request)
    assert result.tool_version == "115"
    assert result.artifact.source is AnnotationSource.VEP
    assert result.artifact.databases[0].version == "115"


def test_dry_run_does_not_require_executable_or_write(tmp_path):
    _, annovar, vep = requests(tmp_path)
    runner = FakeToolRunner()
    assert AnnovarWrapper(runner=runner).run(annovar, dry_run=True).status is AnnotationRunStatus.PLANNED
    assert VepWrapper(runner=runner).run(vep, dry_run=True).status is AnnotationRunStatus.PLANNED
    assert not annovar.output_tsv.exists() and not vep.output_tsv.exists()


def test_missing_database_and_cache_are_rejected(tmp_path):
    item, database, cache = annotation_input(tmp_path)
    database.rmdir()
    with pytest.raises(InputValidationError, match="database root"):
        AnnovarRequest(item, database, "db1", ("refGene",), ("g",), tmp_path / "out", "v1")
    cache.rmdir()
    with pytest.raises(InputValidationError, match="cache directory"):
        VepRequest(item, cache, "115", "homo_sapiens", "GRCh38", tmp_path / "out.tsv")


def test_overwrite_and_external_failure_preserve_source(tmp_path):
    item, annovar, vep = requests(tmp_path)
    annovar.output_tsv.parent.mkdir(parents=True)
    annovar.output_tsv.write_text("existing\n", encoding="utf-8")
    with pytest.raises(OutputValidationError, match="already exists"):
        AnnovarWrapper(runner=FakeToolRunner()).run(annovar)
    before = item.source_vcf.read_bytes()
    with pytest.raises(CommandExecutionError, match="failed"):
        VepWrapper(runner=FakeToolRunner(fail=True)).run(vep)
    assert item.source_vcf.read_bytes() == before


def test_region_overlap_categories_and_breakpoints_are_immutable(tmp_path):
    item, _, _ = annotation_input(tmp_path, VariantCategory.SV, ("sv1",))
    variants = read_selected_variant_loci(item)
    databases = []
    for category in RegionCategory:
        bed = tmp_path / f"{category.value}.bed"
        bed.write_text(f"chr1\t90\t120\t{category.value}-feature\n", encoding="utf-8")
        databases.append(RegionDatabase(category, bed, "2026.1", "GRCh38"))
    output = tmp_path / "overlap.tsv"
    result = annotate_region_overlaps(variants, databases, reference=item.reference, output_tsv=output)
    assert {row.category for row in result.overlaps} == set(RegionCategory)
    assert all((row.variant.start, row.variant.end) == (100, 150) for row in result.overlaps)
    assert "breakpoint_modified\tfunctional_impact_is_call_confidence" in output.read_text(encoding="utf-8")


class MockAlphaGenomeBackend:
    def predict(self, request):
        return FunctionalPrioritizationResult(
            request,
            tuple(FunctionalPrediction(item, {"RNA_SEQ": 0.75}) for item in request.selected_variants),
            "mock-api-1",
        )


def test_alphagenome_requires_explicit_selection_and_impact_is_not_confidence(tmp_path):
    annotation = tmp_path / "annotation.tsv"
    annotation.write_text("variant\timpact\nsv1\tmodifier\n", encoding="utf-8")
    selected = PrioritizedVariant(
        "S1", "sv1", VariantCategory.SV, "chr1", 100, "N", "<DEL>",
        annotation, "explicit research candidate",
    )
    request = FunctionalPrioritizationRequest(
        (selected,), "AlphaGenome", "API_VERSION_PENDING_LINUX_VERIFICATION", ("RNA_SEQ",),
    )
    result = run_functional_prioritization(request, backend=MockAlphaGenomeBackend())
    payload = result.to_dict()
    assert payload["request"]["whole_genome_unselected_execution"] is False
    assert payload["predictions"][0]["functional_impact_is_call_confidence"] is False
    with pytest.raises(InputValidationError, match="explicit non-empty selection"):
        FunctionalPrioritizationRequest((), "AlphaGenome", "v1", ("RNA_SEQ",))


def test_functional_selection_file_is_explicit_and_traceable(tmp_path):
    annotation = tmp_path / "S1.small.vep.tsv"
    annotation.write_text("variant\timpact\nv1\tmodifier\n", encoding="utf-8")
    selection = tmp_path / "selected.tsv"
    selection.write_text(
        "sample\tsource_variant_id\tvariant_category\tcontig\tposition\t"
        "reference_bases\talternate_bases\tsource_annotation\tselection_reason\n"
        "S1\tv1\tsmall\tchr1\t10\tA\tT\tS1.small.vep.tsv\texplicit candidate\n",
        encoding="utf-8",
    )
    selected = read_functional_selection(selection)
    assert selected[0].source_annotation == annotation
    assert selected[0].to_dict()["selection_is_explicit"] is True


def test_annotation_manifest_and_phase11_report(tmp_path):
    genome, vcf, _, _ = files(tmp_path)
    manifest = tmp_path / "annotation_inputs.tsv"
    manifest.write_text(
        "sample\tvariant_category\tsource_vcf\tsource_tool\tsource_variant_ids\n"
        f"S1\tsmall\t{vcf.name}\tdeepvariant\tsv1\n",
        encoding="utf-8",
    )
    inputs = read_annotation_inputs(manifest, reference=genome)
    assert inputs[0].source_variant_ids == ("sv1",)
    _, annovar, vep = requests(tmp_path / "run")
    report = run_phase11(
        (AnnotationJob(annovar.input, annovar, vep),),
        annovar_wrapper=AnnovarWrapper(runner=FakeToolRunner()),
        vep_wrapper=VepWrapper(runner=FakeToolRunner()),
    )
    payload = report.to_dict()
    assert len(payload["annotation_results"]) == 2
    assert payload["scientific_policy"]["functional_impact_is_variant_call_confidence"] is False
    report.write_json(tmp_path / "phase11.json")
    report.write_yaml(tmp_path / "phase11.yaml")


def test_annotation_config_validation():
    validate_config({
        "annotation": {
            "enabled": True, "input_manifest": "inputs.tsv", "overwrite": False,
            "threads": 4, "memory_mb": 16000, "runtime_minutes": 480,
            "annovar_enabled": True, "annovar_executable": "table_annovar.pl",
            "annovar_version": "2025-03-02", "annovar_database_root": "humandb",
            "annovar_database_version": "2026.1", "annovar_protocols": ["refGeneWithVer"],
            "annovar_operations": ["g"], "vep_enabled": False,
            "overlap_enabled": False, "functional_enabled": False,
        }
    })
    with pytest.raises(ConfigurationError, match="input_manifest"):
        validate_config({"annotation": {"enabled": True, "annovar_enabled": True}})
    with pytest.raises(ConfigurationError, match="equal length"):
        validate_config({"annotation": {
            "annovar_enabled": True, "annovar_executable": "a", "annovar_version": "v",
            "annovar_database_root": "d", "annovar_database_version": "dv",
            "annovar_protocols": ["p"], "annovar_operations": ["g", "r"],
        }})
