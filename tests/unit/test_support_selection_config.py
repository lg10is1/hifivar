"""User-configured support views preserve native Phase 9 output authority.

All VCFs and tool results here are synthetic engineering fixtures, not accuracy
evidence. Native benchmark results are documented separately.
"""
from __future__ import annotations

from contextlib import AbstractContextManager
from copy import deepcopy
from dataclasses import replace
import gzip
import hashlib
from pathlib import Path
import runpy
import struct
from types import SimpleNamespace
import zlib

import pytest
import yaml

import hifivar
from hifivar import launcher
from hifivar.assembly_sv import SVEvidenceSource
from hifivar.command import CommandResult
from hifivar.config import HiFiVarConfig, load_config, validate_config
from hifivar.context import AnalysisContext
from hifivar.exceptions import ConfigurationError, InputValidationError, OutputValidationError
from hifivar.harmonization import EvidenceRunStatus, SVEvidenceSourceArtifact
from hifivar.jasmine import JasmineWrapper
from hifivar.phase9 import Phase9Settings, run_phase9
from hifivar.reference import ReferenceGenome
from hifivar.sample import InputDataset, Sample
from hifivar.sample_sheet import SampleRecord
from hifivar.sv_support_selection import SupportSelectionWrapper, support_selection_paths
from hifivar.truvari import TruvariWrapper


ROOT = Path(hifivar.__file__).resolve().parents[2]
DEFAULT = ROOT / "src/hifivar/resources/configs/default.yaml"
CALLERS = (
    ("sawfish", SVEvidenceSource.READ),
    ("sniffles2", SVEvidenceSource.READ),
    ("pbsv", SVEvidenceSource.READ),
    ("cutesv", SVEvidenceSource.READ),
    ("pav", SVEvidenceSource.ASSEMBLY),
    ("svim_asm", SVEvidenceSource.ASSEMBLY),
)


def _bgzf(path: Path, payload: bytes) -> None:
    compressor = zlib.compressobj(level=6, wbits=-15)
    data = compressor.compress(payload) + compressor.flush()
    size = 18 + len(data) + 8
    header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", size - 1)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(header + data + struct.pack("<II", zlib.crc32(payload), len(payload)))


def _vcf(sample: str, records: str) -> bytes:
    return (
        "##fileformat=VCFv4.3\n"
        "##contig=<ID=chr1,length=1000>\n"
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
        '##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Source vector">\n'
        '##INFO=<ID=SUPP,Number=1,Type=Integer,Description="Source count">\n'
        '##INFO=<ID=IDLIST,Number=.,Type=String,Description="Source IDs">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
        + records
    ).encode("utf-8")


def _settings(minimum=None, *, overwrite=False):
    section = {
        "enabled": True, "backend": "jasmine", "jasmine_executable": "jasmine",
        "truvari_executable": "truvari", "bgzip_executable": "bgzip", "tabix_executable": "tabix",
        "threads": 2, "memory_mb": 1000, "runtime_minutes": 30,
        "max_dist": 1000, "distance_type": "linear", "overwrite": overwrite,
    }
    if minimum is not None:
        section["minimum_caller_support"] = minimum
    return {"sv": {"harmonization": section}}


@pytest.mark.parametrize("minimum", [1, 2, 3, 4, 5, 6])
def test_support_config_accepts_partial_integer_override(minimum):
    # A user layer may supply only this value; merged settings own enablement.
    validate_config({"sv": {"harmonization": {"minimum_caller_support": minimum}}})


@pytest.mark.parametrize("minimum", [True, False, 2.0, "2", 0, -1, 7, None])
def test_support_config_rejects_invalid_raw_values(minimum):
    with pytest.raises(ConfigurationError, match="minimum_caller_support"):
        validate_config({"sv": {"harmonization": {"minimum_caller_support": minimum}}})


def test_explicitly_disabled_selection_is_rejected():
    with pytest.raises(ConfigurationError, match="requires.*enabled"):
        validate_config({"sv": {"harmonization": {"enabled": False, "minimum_caller_support": 2}}})


def test_phase9_settings_default_is_unfiltered_and_serialized():
    settings = Phase9Settings.from_config(_settings())
    assert settings.minimum_caller_support == 1
    assert settings.to_dict()["minimum_caller_support"] == 1


@pytest.mark.parametrize("minimum", [True, False, 2.0, "2", 0, -1, 7])
def test_phase9_settings_does_not_coerce_support_threshold(minimum):
    with pytest.raises(ConfigurationError, match="minimum_caller_support"):
        Phase9Settings.from_config(_settings(minimum))


@pytest.fixture
def launch_case(tmp_path):
    (tmp_path / "ref.fa").write_text(">chr1\n" + "A" * 1000 + "\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t1000\t6\t1000\t1001\n", encoding="utf-8")
    for sample in ("S2", "S1"):
        (tmp_path / f"{sample}.bam").write_bytes(b"SYNTHETIC_BAM_" + sample.encode())
        (tmp_path / f"{sample}.bam.bai").write_bytes(b"SYNTHETIC_INDEX")
    (tmp_path / "samples.tsv").write_text(
        "sample_id\tinput\tinput_type\nS2\tS2.bam\tbam\nS1\tS1.bam\tbam\n", encoding="utf-8")
    user = tmp_path / "config.yaml"
    user.write_text(yaml.safe_dump({
        "reference": {"fasta": "ref.fa", "build": "CUSTOM_TEST"},
        "samples": {"sheet": "samples.tsv"}, "runtime": {"threads": 4},
        "workflow": {"memory_mb": 2000},
        "sv": {"enabled": True,
            "sawfish": {"enabled": True, "threads": 2, "memory_mb": 1000},
            "sniffles2": {"enabled": True, "threads": 2, "memory_mb": 1000},
            "pbsv": {"enabled": False}, "cutesv": {"enabled": False},
            "harmonization": {"enabled": True, "threads": 2, "memory_mb": 1000,
                "minimum_caller_support": 2}}}), encoding="utf-8")
    return load_config(DEFAULT, user_config=user), tmp_path


def test_launcher_preserves_two_source_threshold_without_writes(launch_case):
    cfg, root = launch_case
    before = cfg.to_dict()
    planned = root / "planned"
    effective, context, sheet, _ = launcher.prepare_launch(cfg, planned)
    assert context.sample_ids == ("S2", "S1")
    assert effective["sv"]["harmonization"]["minimum_caller_support"] == 2
    assert "sv.harmonization" in dict(launcher._enabled_tracks(effective))
    assert sheet.index("S2\t") < sheet.index("S1\t")
    assert cfg.to_dict() == before
    assert not planned.exists()


def test_launcher_rejects_threshold_above_actual_callers_before_writes(launch_case):
    cfg, root = launch_case
    values = cfg.to_dict()
    values["sv"]["harmonization"]["minimum_caller_support"] = 3
    planned = root / "planned"
    with pytest.raises(ConfigurationError, match="exceeds.*actual input callers"):
        launcher.prepare_launch(HiFiVarConfig(values, cfg.sources), planned)
    assert not planned.exists()


class _FixtureRunner:
    """Write parseable synthetic products and retain every argv for assertions."""
    def __init__(self):
        self.commands = []

    def require_executable(self, executable):
        return Path(executable)

    def run(self, command, **kwargs):
        args = tuple(str(value) for value in command)
        self.commands.append(args)
        dry = kwargs.get("dry_run", False)
        if "--version" in args or args[-1:] == ("version",):
            versions = {"jasmine": "1.1.5", "bgzip": "1.24", "tabix": "1.24", "truvari": "5.4.0"}
            return CommandResult(args, 0, f"{args[0]} {versions[args[0]]}", "", 0.01, None, True)
        if not dry and args[0] == "jasmine":
            raw = Path(next(value.split("=", 1)[1] for value in args if value.startswith("out_file=")))
            raw.write_bytes(_vcf("S1",
                "chr1\t10\tJ1\tA\t<DEL>\t17\tPASS\tSVTYPE=DEL;END=70;SUPP=1;SUPP_VEC=100000;IDLIST=r1\tGT\t0|1\n"
                "chr1\t30\tJ2\tA\t<DEL>\t25\tPASS\tSVTYPE=DEL;END=90;SUPP=2;SUPP_VEC=110000;IDLIST=r1,r2\tGT\t1|0\n"
                "chr1\t60\tJ6\tA\t<DEL>\t38\tPASS\tSVTYPE=DEL;END=120;SUPP=6;SUPP_VEC=111111;IDLIST=r1,r2,r3,r4,r5,r6\tGT\t1|1\n"))
        elif not dry and args[0] == "bgzip" and "-c" in args:
            _bgzf(Path(kwargs["stdout_path"]), Path(args[-1]).read_bytes())
        elif not dry and args[0] == "bgzip" and "-t" in args:
            with gzip.open(args[-1], "rb") as handle:
                handle.read()
        elif not dry and args[0] == "tabix" and "-p" in args:
            _bgzf(Path(args[-1] + ".tbi"), b"TBI\x01")
        elif not dry and args[:2] == ("truvari", "bench"):
            directory = Path(args[args.index("-o") + 1])
            directory.mkdir(parents=True)
            (directory / "summary.json").write_text('{"comparison_only": true}\n', encoding="utf-8")
        return CommandResult(args, 0, "", "", 0.01, None, not dry)


def _case(root: Path, samples=("S1",)):
    fasta = root / "ref.fa"
    fasta.write_text(">chr1\n" + "A" * 1000 + "\n", encoding="utf-8")
    Path(str(fasta) + ".fai").write_text("chr1\t1000\t6\t1000\t1001\n", encoding="utf-8")
    records, all_sources = [], {}
    for sample in samples:
        bam = root / f"{sample}.bam"
        bam.write_bytes(b"SYNTHETIC_" + sample.encode())
        records.append(SampleRecord(Sample(sample, InputDataset.from_files([bam]))))
        sources = []
        for number, (caller, kind) in enumerate(CALLERS, 1):
            path = root / f"{sample}.{caller}.vcf.gz"
            _bgzf(path, _vcf(sample,
                f"chr1\t{number}\tr{number}\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=70\tGT\t0/1\n"))
            index = Path(str(path) + ".tbi")
            _bgzf(index, b"TBI\x01")
            sources.append(SVEvidenceSourceArtifact(sample, kind, caller, path, index,
                EvidenceRunStatus.COMPLETED, ("haplotype1", "haplotype2") if kind is SVEvidenceSource.ASSEMBLY else ()))
        all_sources[sample] = tuple(sources)
    context = AnalysisContext(ReferenceGenome.from_fasta(fasta, build="GRCh38"), tuple(records), {})
    return context, all_sources


def _data_ids(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.split("\t")[2] for line in handle if not line.startswith("#")]


@pytest.mark.parametrize("minimum,selected_ids", [(2, ["J2", "J6"]), (6, ["J6"])])
def test_phase9_preserves_native_and_separate_selected_outputs(tmp_path, minimum, selected_ids):
    context, sources = _case(tmp_path)
    protected = {path: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in tmp_path.iterdir() if path.is_file()}
    runner = _FixtureRunner()
    report = run_phase9(context, sources, output_directory=tmp_path / "out", work_directory=tmp_path / "work",
        config=_settings(minimum), jasmine_wrapper=JasmineWrapper(runner=runner),
        truvari_wrapper=TruvariWrapper(runner=runner), support_wrapper=SupportSelectionWrapper(runner=runner))
    sample = report.sample_results[0]
    native = sample.jasmine.artifact
    selected = sample.support_selection
    assert native is not None and selected is not None
    assert _data_ids(native.vcf_path) == ["J1", "J2", "J6"]
    assert _data_ids(selected.vcf_path) == selected_ids
    assert selected.selected_records == len(selected_ids) and selected.source_records == 3
    assert selected.vcf_path != native.vcf_path
    assert selected.index_path.is_file() and selected.evidence_table.is_file()
    assert len(selected.evidence_table.read_text(encoding="utf-8").splitlines()) == len(selected_ids) + 1
    payload = report.to_dict()
    assert payload["settings"]["minimum_caller_support"] == minimum
    assert payload["sample_results"][0]["support_selection"]["request"]["minimum_caller_support"] == minimum
    assert payload["sample_results"][0]["support_selection"]["canonical_output_role_unchanged"] is True
    comparisons = [command for command in runner.commands if command[:2] == ("truvari", "bench")]
    assert len(comparisons) == 6
    assert {command[command.index("-b") + 1] for command in comparisons} == {str(native.vcf_path)}
    native_commands = [command for command in runner.commands if command[0] == "jasmine" and any(x.startswith("out_file=") for x in command)]
    assert len(native_commands) == 1 and not any(x.startswith("min_support=") for x in native_commands[0])
    assert {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in protected} == protected


def test_phase9_default_one_creates_no_selection_bundle(tmp_path):
    context, sources = _case(tmp_path)
    runner = _FixtureRunner()
    report = run_phase9(context, sources, output_directory=tmp_path / "out", work_directory=tmp_path / "work",
        config=_settings(), jasmine_wrapper=JasmineWrapper(runner=runner), truvari_wrapper=TruvariWrapper(runner=runner))
    assert report.settings.minimum_caller_support == 1
    assert report.sample_results[0].support_selection is None
    assert report.to_dict()["sample_results"][0]["support_selection"] is None
    assert not list((tmp_path / "out").rglob("support-ge-*"))


class _NeverRun:
    def run(self, *args, **kwargs):
        pytest.fail("No wrapper may execute before every sample passes preflight")


@pytest.mark.parametrize("available,minimum", [(1, 2), (2, 3), (5, 6)])
def test_phase9_actual_source_ceiling_fails_before_tools(tmp_path, available, minimum):
    context, sources = _case(tmp_path)
    sources["S1"] = tuple(source if number < available else replace(source,
        status=EvidenceRunStatus.DISABLED, vcf_path=None, index_path=None)
        for number, source in enumerate(sources["S1"]))
    with pytest.raises(ConfigurationError, match="exceeds.*actual input callers"):
        run_phase9(context, sources, output_directory=tmp_path / "out", work_directory=tmp_path / "work",
            config=_settings(minimum), jasmine_wrapper=_NeverRun(), truvari_wrapper=_NeverRun())
    assert not (tmp_path / "out").exists() and not (tmp_path / "work").exists()


@pytest.mark.parametrize("second_sources", ["missing", "insufficient"])
def test_phase9_validates_all_samples_before_first_command(tmp_path, second_sources):
    context, sources = _case(tmp_path, samples=("S1", "S2"))
    if second_sources == "missing":
        del sources["S2"]
        error = InputValidationError
    else:
        sources["S2"] = sources["S2"][:1]
        error = ConfigurationError
    with pytest.raises(error):
        run_phase9(context, sources, output_directory=tmp_path / "out", work_directory=tmp_path / "work",
            config=_settings(2), jasmine_wrapper=_NeverRun(), truvari_wrapper=_NeverRun())
    assert not (tmp_path / "out").exists() and not (tmp_path / "work").exists()


def test_phase9_protects_existing_selection_for_later_sample_even_overwrite(tmp_path):
    context, sources = _case(tmp_path, samples=("S1", "S2"))
    bundle = support_selection_paths(tmp_path / "out", "S2", 2)["directory"]
    bundle.mkdir(parents=True)
    marker = bundle / "protected.txt"
    marker.write_bytes(b"USER_OUTPUT_MUST_SURVIVE")
    with pytest.raises(OutputValidationError, match="fresh per-sample bundle"):
        run_phase9(context, sources, output_directory=tmp_path / "out", work_directory=tmp_path / "work",
            config=_settings(2, overwrite=True), jasmine_wrapper=_NeverRun(), truvari_wrapper=_NeverRun())
    assert marker.read_bytes() == b"USER_OUTPUT_MUST_SURVIVE"
    assert not (tmp_path / "out/S1").exists() and not (tmp_path / "work").exists()


class _Receipt(AbstractContextManager):
    def __init__(self, observed):
        self.observed = observed
    def __exit__(self, *args):
        return False
    def complete(self, **values):
        self.observed["complete"] = values


@pytest.mark.parametrize("mode", ["correct", "wrong_bundle", "missing_result"])
def test_rule_bridge_binds_selected_bundle_and_receipt(launch_case, monkeypatch, mode):
    from hifivar import attempt_receipts, phase9
    cfg, root = launch_case
    config, _, _, _ = launcher.prepare_launch(cfg, root / "planned")
    config["samples"]["sheet"] = str(root / "samples.tsv")
    monkeypatch.delenv(attempt_receipts.ATTEMPT_CONTEXT_ENV, raising=False)
    output = root / "results/sv_harmonized/S1"
    paths = support_selection_paths(output.parent, "S1", 2)
    observed = {}

    def attempt_step(track, samples, **values):
        observed["step"] = {"track": track, "samples": samples, **values}
        return _Receipt(observed)

    def fake_phase9(context, evidence, **kwargs):
        assert context.sample_ids == ("S1",) and set(evidence) == {"S1"}
        assert kwargs["config"]["sv"]["harmonization"]["minimum_caller_support"] == 2
        observed["phase9"] = True
        selected = None if mode == "missing_result" else SimpleNamespace(
            to_dict=lambda: {"status": "fixture_only", "minimum_caller_support": 2}, tool_versions={})
        return SimpleNamespace(write_json=lambda *args, **kwargs: None,
            sample_results=(SimpleNamespace(jasmine=SimpleNamespace(tool_versions={"jasmine": "fixture"}),
                truvari=(SimpleNamespace(version="fixture"), SimpleNamespace(version="fixture")),
                support_selection=selected),))

    monkeypatch.setattr(attempt_receipts, "attempt_step", attempt_step)
    monkeypatch.setattr(phase9, "run_phase9", fake_phase9)
    inputs = []
    for caller in ("sawfish", "sniffles2"):
        path = root / "results/sv" / f"S1.{caller}.sv.vcf.gz"
        inputs.extend((str(path), str(path) + ".tbi"))
    declared = root / "wrong-bundle" if mode == "wrong_bundle" else paths["directory"]
    snakemake = SimpleNamespace(wildcards=SimpleNamespace(sample="S1"), config=config,
        input=SimpleNamespace(sources=inputs),
        output=SimpleNamespace(vcf=str(output / "S1.harmonized.sv.vcf.gz"),
            index=str(output / "S1.harmonized.sv.vcf.gz.tbi"), evidence=str(output / "S1.sv.evidence.tsv"),
            provenance=str(output / "S1.phase9.provenance.json"), support_bundle=str(declared)),
        params=SimpleNamespace(work_root=str(root / "work")))
    if mode == "correct":
        runpy.run_path(str(ROOT / "workflow/scripts/run_harmonization.py"), init_globals={"snakemake": snakemake})
        outputs = observed["step"]["outputs"]
        for role in ("vcf", "index", "evidence", "manifest"):
            assert outputs["support_" + role] == paths[role]
        assert outputs["vcf"] == output / "S1.harmonized.sv.vcf.gz"
        metadata = observed["complete"]["metadata"]
        assert metadata["minimum_caller_support"] == 2
        assert metadata["native_output_role_unchanged"] is True
        assert metadata["caller_count_is_confidence"] is False
    else:
        match = "rule output differs" if mode == "wrong_bundle" else "did not produce"
        with pytest.raises(InputValidationError, match=match):
            runpy.run_path(str(ROOT / "workflow/scripts/run_harmonization.py"), init_globals={"snakemake": snakemake})
        assert "complete" not in observed
        if mode == "wrong_bundle":
            assert "phase9" not in observed and "step" not in observed
    assert not output.exists() and not (root / "work").exists()
