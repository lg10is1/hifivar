"""Official Truvari 5.4 flags and distinct bench/refinement output contracts.

The process fixture is synthetic. It exercises the shared runner and native
filenames, not real Truvari matching or GIAB accuracy.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import pytest

from hifivar.benchmark import BenchmarkVariantClass
from hifivar.command import CommandRunner
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.reference import ReferenceGenome
from hifivar.truvari import (
    TruvariRequest,
    TruvariResultStatus,
    TruvariThresholds,
    TruvariWrapper,
    parse_truvari_refined_summary,
    parse_truvari_summary,
)


@pytest.fixture
def truvari_request(tmp_path: Path) -> TruvariRequest:
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    reference = ReferenceGenome.from_fasta(fasta, build="GRCh38")
    truth, query = tmp_path / "truth.vcf.gz", tmp_path / "query.vcf.gz"
    for path in (truth, query):
        path.write_bytes(b"synthetic input retained unchanged")
        Path(f"{path}.tbi").write_bytes(b"synthetic index")
    return TruvariRequest("HG002", reference, truth, query, tmp_path / "out")


def _with_policy(truvari_request: TruvariRequest, **values: object) -> TruvariRequest:
    return TruvariRequest(
        truvari_request.sample_id, truvari_request.reference, truvari_request.base_vcf,
        truvari_request.comparison_vcf, truvari_request.output_directory,
        thresholds=TruvariThresholds(**values),
    )


def _variant_summary() -> dict[str, int | float | None]:
    return {"TP-base": 9, "TP-comp": 8, "FP": 2, "FN": 1,
            "precision": 0.8, "recall": 0.9, "f1": 1.44 / 1.7,
            "base cnt": 10, "comp cnt": 10}


def test_default_command_remains_legacy_bench(truvari_request: TruvariRequest) -> None:
    command = TruvariWrapper().plan_command(truvari_request)
    assert command == (
        "truvari", "bench", "-b", str(truvari_request.base_vcf.absolute()),
        "-c", str(truvari_request.comparison_vcf.absolute()),
        "-f", str(truvari_request.reference.fasta.absolute()),
        "-o", str(truvari_request.output_directory.absolute()),
    )
    assert truvari_request.refinement_output_candidates == {}
    assert truvari_request.to_dict()["bench_summary_candidate"].endswith("summary.json")


@pytest.mark.parametrize("pick", ["single", "multi", "ac"])
def test_explicit_official_protocol_flags(truvari_request: TruvariRequest, pick: str) -> None:
    item = _with_policy(
        truvari_request, refdist=2000, chunksize=5000, pick=pick,
        pass_only=True, refine=True, sizemin=50, sizemax=100000,
    )
    command = TruvariWrapper().plan_command(item)
    for flag, value in (("--refdist", "2000"), ("--chunksize", "5000"), ("--pick", pick)):
        assert command[command.index(flag) + 1] == value
    assert "--passonly" in command and "--refine" in command
    policy = item.to_dict()["thresholds"]
    assert policy["pick"] == pick and policy["chunksize"] == 5000 and policy["refine"] is True
    candidates = item.refinement_output_candidates
    assert candidates["variant_summary"].name == "refine.variant_summary.json"
    assert candidates["region_summary"].name == "refine.region_summary.json"
    assert candidates["base_vcf"].name == "refine.base.vcf.gz"
    assert candidates["comparison_index"].name == "refine.comp.vcf.gz.tbi"


@pytest.mark.parametrize("field,value", [
    ("pick", "best"), ("pick", "AC"), ("pick", ""), ("pick", 1), ("pick", []),
    ("chunksize", 0), ("chunksize", -1), ("chunksize", True), ("chunksize", 1.5),
    ("chunksize", "5000"), ("refine", 1), ("refine", "true"), ("refine", None),
    ("pass_only", 0), ("pass_only", "false"), ("pass_only", None),
])
def test_invalid_protocol_values_fail_closed(field: str, value: object) -> None:
    with pytest.raises(InputValidationError):
        TruvariThresholds(**{field: value})


def test_bounds_are_checked_without_replacing_upstream_defaults() -> None:
    with pytest.raises(InputValidationError, match="sizemin"):
        TruvariThresholds(sizemin=100, sizemax=50)
    with pytest.raises(InputValidationError, match="chunksize"):
        TruvariThresholds(refdist=2000, chunksize=1999)
    assert TruvariThresholds(refdist=2000, chunksize=2000).chunksize == 2000
    assert TruvariThresholds(sizemin=50, sizemax=50).sizemin == 50
    assert TruvariThresholds().pick is None


def test_dry_run_has_no_tool_or_output_requirement(truvari_request: TruvariRequest) -> None:
    item = _with_policy(truvari_request, pick="ac", chunksize=5000, refdist=2000, refine=True)
    result = TruvariWrapper(executable="not-installed-truvari").run(item, dry_run=True)
    assert result.status is TruvariResultStatus.PLANNED
    assert "--refine" in result.command
    assert result.summary_path is None and result.refined_summary_path is None
    assert result.to_dict()["refinement_requested"] is True
    assert not truvari_request.output_directory.exists()


def test_refined_variant_summary_retains_distinct_tp_sides(tmp_path: Path) -> None:
    path = tmp_path / "refine.variant_summary.json"
    path.write_text(json.dumps(_variant_summary()), encoding="utf-8")
    metrics = {item.name: item for item in parse_truvari_refined_summary(
        path, variant_class=BenchmarkVariantClass.ASSEMBLY_SV,
    )}
    assert metrics["tp_base"].value == 9 and metrics["tp_call"].value == 8
    assert metrics["tp_call"].source_field == "TP-comp"
    assert metrics["f1"].variant_class is BenchmarkVariantClass.ASSEMBLY_SV


def test_null_rates_are_preserved_not_invented(tmp_path: Path) -> None:
    path = tmp_path / "refine.variant_summary.json"
    payload = {name: 0 for name in ("TP-base", "TP-comp", "FP", "FN", "base cnt", "comp cnt")}
    payload.update(precision=None, recall=None, f1=None)
    path.write_text(json.dumps(payload), encoding="utf-8")
    metrics = {item.name: item.value for item in parse_truvari_refined_summary(path)}
    assert metrics["precision"] is None and metrics["recall"] is None and metrics["f1"] is None


@pytest.mark.parametrize("name", ["summary.json", "refine.region_summary.json", "other.json"])
def test_region_or_bench_summary_is_not_refined_variant_metrics(tmp_path: Path, name: str) -> None:
    path = tmp_path / name
    path.write_text(json.dumps(_variant_summary()), encoding="utf-8")
    with pytest.raises(OutputValidationError, match="refine.variant_summary.json"):
        parse_truvari_refined_summary(path)


@pytest.mark.parametrize("change", [
    {"TP-base": True}, {"FP": -1}, {"FN": 1.0}, {"base cnt": 99},
    {"comp cnt": 99}, {"precision": float("nan")}, {"recall": 1.2}, {"f1": "0.9"},
])
def test_bad_refined_schema_is_rejected(tmp_path: Path, change: dict) -> None:
    path = tmp_path / "refine.variant_summary.json"
    payload = _variant_summary()
    payload.update(change)
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(OutputValidationError):
        parse_truvari_refined_summary(path)


_FAKE_PROCESS = '''
import gzip
import json
from pathlib import Path
import sys

args = sys.argv[1:]
if args == ["version"]:
    print("Truvari v5.4.0")
    raise SystemExit(0)
out = Path(args[args.index("-o") + 1])
out.mkdir()
(out / "native_argv.json").write_text(json.dumps(args), encoding="utf-8")
bench = {"TP-base": 2, "TP-comp": 2, "FP": 8, "FN": 8,
         "precision": .2, "recall": .2, "f1": .2}
(out / "summary.json").write_text(json.dumps(bench), encoding="utf-8")
if "--refine" in args:
    refined = {"TP-base": 9, "TP-comp": 8, "FP": 2, "FN": 1,
               "precision": .8, "recall": .9, "f1": 1.44 / 1.7,
               "base cnt": 10, "comp cnt": 10}
    (out / "refine.variant_summary.json").write_text(json.dumps(refined), encoding="utf-8")
    (out / "refine.region_summary.json").write_text(json.dumps({"TP": 1, "TN": 1}), encoding="utf-8")
    for name in ("refine.base.vcf.gz", "refine.comp.vcf.gz"):
        with gzip.open(out / name, "wt", encoding="utf-8") as stream:
            stream.write("##fileformat=VCFv4.3\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n")
        Path(str(out / name) + ".tbi").write_bytes(b"synthetic-index")
omit = OMIT
if omit is not None:
    (out / omit).unlink()
print("fake bench finished")
print("fake tool diagnostic", file=sys.stderr)
'''


class _ProcessRunner:
    """Only map a fake executable to Python; all execution uses CommandRunner."""

    def __init__(self, script: Path):
        self.script = script
        self.runner = CommandRunner()
        self.commands: list[tuple[str, ...]] = []

    def require_executable(self, executable: str) -> Path:
        return self.runner.require_executable(sys.executable)

    def run(self, command, **kwargs):
        self.commands.append(tuple(command))
        return self.runner.run(
            [sys.executable, "-I", "-B", str(self.script), *command[1:]], **kwargs,
        )


def _process_runner(tmp_path: Path, omit: str | None = None) -> _ProcessRunner:
    script = tmp_path / "fake_truvari.py"
    script.write_text(_FAKE_PROCESS.replace("OMIT", repr(omit)), encoding="utf-8")
    return _ProcessRunner(script)


def _source_hashes(truvari_request: TruvariRequest) -> dict[Path, str]:
    paths = (truvari_request.reference.fasta, truvari_request.reference.fai, truvari_request.base_vcf,
             Path(f"{truvari_request.base_vcf}.tbi"), truvari_request.comparison_vcf,
             Path(f"{truvari_request.comparison_vcf}.tbi"))
    return {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def test_fake_process_refinement_e2e_keeps_bench_and_refined_distinct(
    truvari_request: TruvariRequest, tmp_path: Path,
) -> None:
    item = _with_policy(truvari_request, pick="ac", pass_only=True, refdist=2000, chunksize=5000, refine=True)
    before = _source_hashes(item)
    runner = _process_runner(tmp_path)
    log = tmp_path / "logs" / "bench.log"
    result = TruvariWrapper(executable="fake-truvari", runner=runner).run(item, stderr_path=log)
    assert result.status is TruvariResultStatus.COMPLETED and result.version == "5.4.0"
    assert result.summary_path == item.output_directory / "summary.json"
    assert result.refined_summary_path == item.output_directory / "refine.variant_summary.json"
    assert result.refined_region_summary_path == item.output_directory / "refine.region_summary.json"
    assert result.refined_base_vcf_path == item.output_directory / "refine.base.vcf.gz"
    assert result.refined_comparison_vcf_path == item.output_directory / "refine.comp.vcf.gz"
    assert {m.name: m.value for m in parse_truvari_summary(result.summary_path)}["f1"] == .2
    assert {m.name: m.value for m in parse_truvari_refined_summary(result.refined_summary_path)}["precision"] == .8
    assert "--refine" in json.loads((item.output_directory / "native_argv.json").read_text())
    assert result.to_dict()["summary_kind"] == "bench"
    assert result.to_dict()["refined_summary_kind"] == "variant"
    assert result.to_dict()["refined_stratification_status"] == "NOT_RUN_NEEDS_REFINED_ASSIGNMENTS"
    assert "fake tool diagnostic" in log.read_text(encoding="utf-8")
    assert "fake bench finished" in log.with_name("bench.stdout.log").read_text(encoding="utf-8")
    assert _source_hashes(item) == before


@pytest.mark.parametrize("missing", [
    "refine.variant_summary.json", "refine.region_summary.json", "refine.base.vcf.gz",
    "refine.base.vcf.gz.tbi", "refine.comp.vcf.gz", "refine.comp.vcf.gz.tbi",
])
def test_exit_zero_without_required_refinement_artifact_fails(
    truvari_request: TruvariRequest, tmp_path: Path, missing: str,
) -> None:
    item = _with_policy(truvari_request, refine=True)
    before = _source_hashes(item)
    with pytest.raises(OutputValidationError):
        TruvariWrapper(runner=_process_runner(tmp_path, omit=missing)).run(item)
    assert (item.output_directory / "summary.json").is_file()
    assert _source_hashes(item) == before


def test_existing_output_is_not_overwritten(truvari_request: TruvariRequest) -> None:
    truvari_request.output_directory.mkdir()
    sentinel = truvari_request.output_directory / "owner.txt"
    sentinel.write_text("preserve", encoding="utf-8")
    with pytest.raises(OutputValidationError, match="exists"):
        _with_policy(truvari_request, refine=True)
    assert sentinel.read_text(encoding="utf-8") == "preserve"


def test_sizemax_minus_one_is_explicit_unlimited_not_default(truvari_request: TruvariRequest) -> None:
    item = _with_policy(truvari_request, sizemin=50, sizemax=-1)
    command = TruvariWrapper().plan_command(item)
    assert command[command.index("--sizemax") + 1] == "-1"
    assert TruvariThresholds().sizemax is None
    assert "--sizemax" not in TruvariWrapper().plan_command(truvari_request)
    with pytest.raises(InputValidationError, match="sizemax"):
        TruvariThresholds(sizemin=50, sizemax=-2)


@pytest.mark.parametrize("refined", [False, True])
@pytest.mark.parametrize("counts,rates", [
    ((0, 0, 0, 0), (None, None, None)),
    ((0, 0, 2, 0), (None, None, None)),
    ((0, 0, 0, 3), (None, None, None)),
    ((0, 0, 2, 3), (0.0, 0.0, None)),
    ((1, 0, 2, 1), (0.0, 0.5, 0.0)),
])
def test_native_null_and_asymmetric_count_semantics(tmp_path, refined, counts, rates) -> None:
    tp_base, tp_comp, fp, fn = counts
    payload = {"TP-base": tp_base, "TP-comp": tp_comp, "FP": fp, "FN": fn,
               "base cnt": tp_base + fn, "comp cnt": tp_comp + fp,
               "precision": rates[0], "recall": rates[1], "f1": rates[2],
               "gt_concordance": None}
    path = tmp_path / ("refine.variant_summary.json" if refined else "summary.json")
    path.write_text(json.dumps(payload), encoding="utf-8")
    parser = parse_truvari_refined_summary if refined else parse_truvari_summary
    result = {metric.name: metric.value for metric in parser(path)}
    assert (result["precision"], result["recall"], result["f1"]) == rates
    assert (result["tp_base"], result["tp_call"]) == (tp_base, tp_comp)


@pytest.mark.parametrize("refined", [False, True])
@pytest.mark.parametrize("change", [
    {"precision": None}, {"recall": None}, {"f1": None},
    {"precision": .1}, {"recall": .1}, {"f1": .1},
    {"precision": float("nan")}, {"recall": float("inf")},
    {"TP-base": -1}, {"FP": 1.0}, {"FN": True}, {"comp cnt": 999},
])
def test_native_summary_rejects_inconsistent_counts_rates_and_nulls(tmp_path, refined, change) -> None:
    payload = _variant_summary()
    payload.update(change)
    path = tmp_path / ("refine.variant_summary.json" if refined else "summary.json")
    path.write_text(json.dumps(payload), encoding="utf-8")
    parser = parse_truvari_refined_summary if refined else parse_truvari_summary
    with pytest.raises(OutputValidationError):
        parser(path)


@pytest.mark.parametrize("refined", [False, True])
def test_six_decimal_native_rates_accepted_without_replacing_values(tmp_path, refined) -> None:
    payload = _variant_summary()
    payload["f1"] = round(payload["f1"], 6)
    path = tmp_path / ("refine.variant_summary.json" if refined else "summary.json")
    path.write_text(json.dumps(payload), encoding="utf-8")
    parser = parse_truvari_refined_summary if refined else parse_truvari_summary
    assert {m.name: m.value for m in parser(path)}["f1"] == payload["f1"]


def test_legacy_abbreviated_rounded_bench_summary_remains_readable(tmp_path) -> None:
    payload = {"TP-base": 8, "TP-comp": 7, "FP": 2, "FN": 3,
               "precision": .777, "recall": .727, "f1": .75}
    path = tmp_path / "summary.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    result = {m.name: m.value for m in parse_truvari_summary(path)}
    assert (result["precision"], result["recall"], result["f1"]) == (.777, .727, .75)


@pytest.mark.parametrize("payload", [[], None, 1, {"TP-base": 0}])
def test_nonobject_or_incomplete_bench_summary_has_explicit_error(tmp_path, payload) -> None:
    path = tmp_path / "summary.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(OutputValidationError):
        parse_truvari_summary(path)


def test_bench_parser_cannot_read_refined_or_region_metrics(tmp_path) -> None:
    for name in ("refine.variant_summary.json", "refine.region_summary.json"):
        path = tmp_path / name
        path.write_text(json.dumps(_variant_summary()), encoding="utf-8")
        with pytest.raises(OutputValidationError, match="not a refined"):
            parse_truvari_summary(path)


def test_only_one_native_total_is_not_a_valid_legacy_summary(tmp_path) -> None:
    payload = _variant_summary()
    del payload["comp cnt"]
    path = tmp_path / "summary.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(OutputValidationError, match="both"):
        parse_truvari_summary(path)
