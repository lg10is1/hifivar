"""Synthetic process/manifest contracts, not real Truvari performance tests."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import struct
import sys
import zlib

import pytest
import yaml

from hifivar.cli import build_parser, main
from hifivar.command import CommandRunner
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.sv_benchmark import SourceImmutabilityError, load_sv_benchmark_job, run_sv_benchmark


def _bgzf(path: Path, payload: bytes) -> None:
    compressor = zlib.compressobj(level=6, wbits=-15)
    body = compressor.compress(payload) + compressor.flush()
    size = 18 + len(body) + 8
    header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", size - 1)
    path.write_bytes(header + body + struct.pack("<II", zlib.crc32(payload), len(payload)))


def _vcf(path: Path, sample: str = "HG002") -> None:
    text = ("##fileformat=VCFv4.3\n##reference=GRCh38\n##contig=<ID=chr1,length=1000>\n"
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
            "chr1\t101\tI1\tA\t" + "A" * 60 + "\t.\tPASS\tSVTYPE=INS;SVLEN=59\tGT\t0/1\n")
    _bgzf(path, text.encode())
    _bgzf(Path(f"{path}.tbi"), b"TBI\x01")


@pytest.fixture
def manifest(tmp_path: Path) -> Path:
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + "A" * 1000 + "\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t1000\t6\t1000\t1001\n", encoding="utf-8")
    (tmp_path / "regions.bed").write_text("chr1\t0\t1000\n", encoding="utf-8")
    for name in ("truth", "Q1", "Q2", "Q3", "Q4"):
        _vcf(tmp_path / f"{name}.vcf.gz")
    payload = {
        "schema_version": 1, "benchmark_id": "BATCH1", "sample_id": "HG002",
        "reference": {"fasta": "ref.fa", "build": "GRCh38"},
        "truth": {"vcf": "truth.vcf.gz", "confident_bed": "regions.bed", "version": "1.00", "source": "synthetic fixture", "reference_build": "GRCh38"},
        "tools": {"truvari_executable": "fake-truvari", "bcftools_executable": "fake-bcftools", "truvari_version": "5.4.0"},
        "thresholds": {"sizemin": 50, "refdist": 2000, "chunksize": 5000, "pick": "ac", "refine": True},
        "queries": [{"id": f"Q{i}", "vcf": f"Q{i}.vcf.gz", "caller": caller}
                    for i, caller in enumerate(("sawfish", "sniffles2", "harmonized", "pav"), start=1)],
    }
    path = tmp_path / "job.yaml"
    path.write_text(yaml.safe_dump(payload), encoding="utf-8")
    return path


def _change(path: Path, *location: str, value: object) -> None:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    parent = data
    for key in location[:-1]:
        parent = parent[key]
    parent[location[-1]] = value
    path.write_text(yaml.safe_dump(data), encoding="utf-8")


_PROCESS = r'''
import gzip
import json
from pathlib import Path
import sys

tool, *args = sys.argv[1:]
mode = MODE
if args in (["version"], ["--version"]):
    print("Truvari v5.4.0" if tool == "fake-truvari" else "bcftools 1.23.1")
    raise SystemExit(0)
if tool == "fake-bcftools":
    assert args[:4] == ["norm", "-N", "-c", "e"]
    assert "-Ou" in args and "-o" in args
    vcf = Path(args[-1])
    if mode == "bad_ref" and vcf.name == "Q2.vcf.gz":
        print("REF allele mismatch", file=sys.stderr)
        raise SystemExit(7)
    print("Reference check-only; no output VCF")
    raise SystemExit(0)
assert tool == "fake-truvari" and args[0] == "bench"
query = Path(args[args.index("-c") + 1])
if mode == "exit_fail" and query.name == "Q2.vcf.gz":
    print("deliberate synthetic caller-result failure", file=sys.stderr)
    raise SystemExit(4)
if mode == "mutate" and query.name == "Q2.vcf.gz":
    query.write_bytes(query.read_bytes() + b"changed")
out = Path(args[args.index("-o") + 1])
out.mkdir(parents=True)
bench = {"TP-base": 2, "TP-comp": 2, "FP": 8, "FN": 8,
         "base cnt": 10, "comp cnt": 10, "precision": .2, "recall": .2, "f1": .2}
if mode == "missing_denominator" and query.name == "Q2.vcf.gz":
    bench.pop("base cnt")
if not (mode == "missing_summary" and query.name == "Q2.vcf.gz"):
    (out / "summary.json").write_text(json.dumps(bench), encoding="utf-8")
header = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
for name in ("tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz"):
    if mode == "missing_assignment" and query.name == "Q2.vcf.gz" and name == "fp.vcf.gz":
        continue
    with gzip.open(out / name, "wt", encoding="utf-8") as f:
        f.write(header)
    Path(str(out / name) + ".tbi").write_bytes(b"synthetic-index")
if "--refine" in args:
    refined = {"TP-base": 9, "TP-comp": 8, "FP": 2, "FN": 1,
               "precision": .8, "recall": .9, "f1": 1.44 / 1.7,
               "base cnt": 10, "comp cnt": 10}
    (out / "refine.variant_summary.json").write_text(json.dumps(refined), encoding="utf-8")
    (out / "refine.region_summary.json").write_text(json.dumps({"TP": 1, "TN": 1}), encoding="utf-8")
    for name in ("refine.base.vcf.gz", "refine.comp.vcf.gz"):
        with gzip.open(out / name, "wt", encoding="utf-8") as f:
            f.write(header)
        Path(str(out / name) + ".tbi").write_bytes(b"synthetic-index")
print("synthetic benchmark complete")
print("retained tool stderr", file=sys.stderr)
'''


class ProcessRunner:
    """Map fake tool names to real child processes through CommandRunner."""

    def __init__(self, script: Path):
        self.script = script
        self.shared = CommandRunner()
        self.commands: list[tuple[str, ...]] = []

    def require_executable(self, executable: str) -> Path:
        return Path(sys.executable)

    def run(self, command, **kwargs):
        self.commands.append(tuple(command))
        return self.shared.run([sys.executable, "-I", "-B", str(self.script), *command], **kwargs)


def _runner(tmp_path: Path, mode: str = "success") -> ProcessRunner:
    script = tmp_path / "fake_tools.py"
    script.write_text(_PROCESS.replace("MODE", repr(mode)), encoding="utf-8")
    return ProcessRunner(script)


def _hashes(manifest: Path) -> dict[str, str]:
    return {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in load_sv_benchmark_job(manifest).source_paths}


def test_loader_resolves_paths_and_preserves_order(manifest: Path) -> None:
    job = load_sv_benchmark_job(manifest)
    assert job.reference.fasta == manifest.parent / "ref.fa"
    assert [q.id for q in job.queries] == ["Q1", "Q2", "Q3", "Q4"]
    assert [q.caller for q in job.queries] == ["sawfish", "sniffles2", "harmonized", "pav"]
    assert job.truth_version == "1.00" and job.thresholds.refine is True


@pytest.mark.parametrize("version", [True, False, "1", "1.0", 1.0, 2, None])
def test_schema_requires_integer_one(manifest: Path, version: object) -> None:
    _change(manifest, "schema_version", value=version)
    with pytest.raises(InputValidationError, match="integer 1"):
        load_sv_benchmark_job(manifest)


@pytest.mark.parametrize("location,value", [
    (("extra",), 1), (("reference", "extra"), 1),
    (("tools", "extra"), "x"), (("truth", "extra"), "x"),
    (("thresholds", "extra"), 1), (("sample_id",), ""),
    (("benchmark_id",), "../escape"), (("queries",), []),
    (("truth", "version"), "latest"), (("tools", "truvari_version"), "latest"),
    (("reference", "fasta"), "../ref.fa"), (("reference", "fasta"), "~/.ref.fa"),
    (("truth", "reference_build"), "GRCh37"), (("tools", "truvari_executable"), "truvari --flag"),
    (("thresholds", "pctseq"), float("nan")),
])
def test_invalid_manifest_fails_closed(manifest: Path, location: tuple, value: object) -> None:
    _change(manifest, *location, value=value)
    with pytest.raises(InputValidationError):
        load_sv_benchmark_job(manifest)


def test_yaml_duplicate_keys_are_rejected(manifest: Path) -> None:
    manifest.write_text(manifest.read_text() + "sample_id: OTHER\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="unique"):
        load_sv_benchmark_job(manifest)


@pytest.mark.parametrize("query", [
    {"id": "../Q", "vcf": "Q1.vcf.gz", "caller": "x"},
    {"id": "Q", "vcf": "Q1.vcf.gz", "caller": ""},
    {"id": "Q", "vcf": "Q1.vcf.gz", "caller": "x", "extra": 1},
    {"id": "Q", "vcf": "truth.vcf.gz", "caller": "x"},
])
def test_invalid_query_rejected(manifest: Path, query: dict) -> None:
    _change(manifest, "queries", value=[query])
    with pytest.raises(InputValidationError):
        load_sv_benchmark_job(manifest)


@pytest.mark.parametrize("ids,files", [(["Q", "q"], ["Q1.vcf.gz", "Q2.vcf.gz"]),
                                      (["Q1", "Q2"], ["Q1.vcf.gz", "Q1.vcf.gz"])])
def test_duplicate_query_id_or_source_rejected(manifest: Path, ids: list, files: list) -> None:
    _change(manifest, "queries", value=[{"id": i, "vcf": p, "caller": "x"} for i, p in zip(ids, files)])
    with pytest.raises(InputValidationError):
        load_sv_benchmark_job(manifest)


def test_dry_run_no_execution_no_output_and_no_sequence_pass(manifest: Path) -> None:
    class NeverRunner:
        def require_executable(self, *args, **kwargs):
            pytest.fail("dry-run must not locate or run native tools")
        def run(self, *args, **kwargs):
            pytest.fail("dry-run must not execute")
    before = _hashes(manifest)
    root = manifest.parent / "new" / "dry"
    report = run_sv_benchmark(manifest, root, dry_run=True, runner=NeverRunner())
    assert report["status"] == "PLANNED" and report["commands"] == []
    assert all(q["status"] == "PLANNED" and q["bench_metrics"] == [] for q in report["queries"])
    assert report["ref_sequence_check"] == "NOT_RUN_DRY_RUN"
    assert all(q["input_audit"]["sequence_check"] == "NOT_RUN_REQUIRES_BCFTOOLS" for q in report["queries"])
    assert not root.parent.exists() and _hashes(manifest) == before


def test_unknown_sample_dry_run_is_not_success(manifest: Path) -> None:
    _vcf(manifest.parent / "Q2.vcf.gz", "UnknownSampleName")
    report = run_sv_benchmark(manifest, manifest.parent / "dry", dry_run=True)
    assert report["status"] == "FAILED"
    assert [q["status"] for q in report["queries"]] == ["PLANNED", "FAILED", "PLANNED", "PLANNED"]
    assert not (manifest.parent / "dry").exists()


def test_real_fake_process_matrix_complete_and_metrics_separate(manifest: Path) -> None:
    before = _hashes(manifest)
    runner = _runner(manifest.parent)
    root = manifest.parent / "matrix"
    report = run_sv_benchmark(manifest, root, runner=runner)
    assert report["status"] == "COMPLETE" and report["schema_version"] == 1
    assert len(report["queries"]) == 4
    assert all(q["status"] == "COMPLETE" for q in report["queries"])
    for query in report["queries"]:
        assert {m["name"]: m["value"] for m in query["bench_metrics"]}["f1"] == .2
        assert {m["name"]: m["value"] for m in query["refined_metrics"]}["precision"] == .8
        assert query["result"]["summary_kind"] == "bench"
        assert query["refined_stratification_status"] == "NOT_RUN_NEEDS_REFINED_ASSIGNMENTS"
    norm = [command for command in runner.commands if command[1:2] == ("norm",)]
    assert len(norm) == 5 and all(command[2:5] == ("-N", "-c", "e") for command in norm)
    assert report["source_immutability"].startswith("PASS_FULL_SHA256")
    assert report["input_checksums_after"] == report["input_checksums"]
    assert _hashes(manifest) == before
    assert report["tools"]["observed_bcftools_version"] == "1.23.1"
    saved = json.loads((root / "benchmark_matrix.json").read_text(encoding="utf-8"))
    assert saved == report
    assert all(c["return_code"] == 0 and c["executed"] for c in report["commands"])
    assert all(Path(c["stdout_path"]).is_file() and Path(c["stderr_path"]).is_file() for c in report["commands"])
    assert any("retained tool stderr" in p.read_text() for p in (root / "logs").glob("*.stderr.log"))
    assert report["output_checksums"] and report["input_checksums"][str(manifest)]["sha256"] == before[str(manifest)]


@pytest.mark.parametrize("mode", ["exit_fail", "missing_summary", "missing_assignment", "missing_denominator", "bad_ref"])
def test_one_failure_isolated_and_remaining_queries_continue(manifest: Path, mode: str) -> None:
    before = _hashes(manifest)
    root = manifest.parent / "matrix"
    report = run_sv_benchmark(manifest, root, runner=_runner(manifest.parent, mode))
    assert report["status"] == "PARTIAL"
    assert [q["status"] for q in report["queries"]] == ["COMPLETE", "FAILED", "COMPLETE", "COMPLETE"]
    assert report["queries"][1]["error_type"]
    assert _hashes(manifest) == before
    if mode == "exit_fail":
        assert any(c["return_code"] == 4 for c in report["commands"])
    if mode == "bad_ref":
        assert report["queries"][1]["ref_sequence_check"] == "NOT_RUN"


def test_source_mutation_invalidates_matrix_and_stops_remaining(manifest: Path) -> None:
    runner = _runner(manifest.parent, "mutate")
    report = run_sv_benchmark(manifest, manifest.parent / "matrix", runner=runner)
    assert report["status"] == "FAILED" and report["source_immutability"] == "FAIL"
    assert report["metrics_valid_for_claims"] is False
    assert report["queries"][0]["status"] == "INVALIDATED_SOURCE_MUTATION"
    assert all(q["status"] == "NOT_RUN_SOURCE_MUTATION" for q in report["queries"][2:])
    assert not any("Q3.vcf.gz" in str(c) for c in runner.commands)


def test_version_mismatch_prevents_comparisons(manifest: Path) -> None:
    _change(manifest, "tools", "truvari_version", value="5.3.0")
    runner = _runner(manifest.parent)
    report = run_sv_benchmark(manifest, manifest.parent / "matrix", runner=runner)
    assert report["status"] == "FAILED" and report["error_type"] == "ToolVersionError"
    assert len(runner.commands) == 1
    assert all(q["status"] == "NOT_RUN" for q in report["queries"])


def test_fresh_output_required_and_owner_files_preserved(manifest: Path) -> None:
    root = manifest.parent / "matrix"
    root.mkdir()
    sentinel = root / "owner.txt"
    sentinel.write_text("preserve", encoding="utf-8")
    with pytest.raises(OutputValidationError, match="fresh"):
        run_sv_benchmark(manifest, root, dry_run=True)
    assert sentinel.read_text() == "preserve"


def test_output_traversal_and_input_alias_rejected(manifest: Path) -> None:
    with pytest.raises(InputValidationError, match="traversal"):
        run_sv_benchmark(manifest, manifest.parent / "x" / ".." / "new", dry_run=True)
    with pytest.raises(OutputValidationError):
        run_sv_benchmark(manifest, manifest.parent / "truth.vcf.gz", dry_run=True)


def test_cli_registration_and_global_dry_run(manifest: Path, capsys) -> None:
    args = build_parser().parse_args(["--dry-run", "benchmark", "sv", "--manifest", str(manifest), "--output-dir", str(manifest.parent / "dry")])
    assert args.command == "benchmark" and args.benchmark_command == "sv" and args.dry_run
    assert main(["--dry-run", "benchmark", "sv", "--manifest", str(manifest), "--output-dir", str(manifest.parent / "dry")]) == 0
    assert json.loads(capsys.readouterr().out)["status"] == "PLANNED"


def test_cli_partial_returns_nonzero_without_changing_existing_entries(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr("hifivar.sv_benchmark.run_sv_benchmark", lambda *a, **k: {"status": "PARTIAL", "queries": [{"status": "COMPLETE"}, {"status": "FAILED"}]})
    assert main(["benchmark", "sv", "--manifest", "unused.yaml", "--output-dir", str(tmp_path / "out")]) == 1
    assert build_parser().parse_args(["review", "prepare", "--output-dir", "out"]).command == "review"


def test_cli_benchmark_does_not_open_unowned_logging_path(manifest: Path) -> None:
    log = manifest.parent / "not-owned" / "outside.log"
    assert main(["--log-file", str(log), "--dry-run", "benchmark", "sv", "--manifest", str(manifest), "--output-dir", str(manifest.parent / "dry")]) == 0
    assert not log.parent.exists()


def test_manifest_change_between_parse_and_hash_rejected(manifest: Path, monkeypatch) -> None:
    original = load_sv_benchmark_job

    def changed(path):
        job = original(path)
        manifest.write_bytes(manifest.read_bytes() + b"\n")
        return job

    monkeypatch.setattr("hifivar.sv_benchmark.load_sv_benchmark_job", changed)
    root = manifest.parent / "matrix"
    with pytest.raises(SourceImmutabilityError, match="between parsing"):
        run_sv_benchmark(manifest, root, runner=_runner(manifest.parent))
    assert not root.exists()


def test_outputs_cannot_reenter_tool_installation(manifest: Path) -> None:
    executable = manifest.parent / "protected-tool" / "bin" / "tool"
    executable.parent.mkdir(parents=True)
    executable.write_text("synthetic executable", encoding="utf-8")

    class InstalledRunner:
        def require_executable(self, name):
            return executable

    root = executable.parent.parent / "analysis-output"
    with pytest.raises(OutputValidationError, match="inside a tool installation"):
        run_sv_benchmark(manifest, root, runner=InstalledRunner())
    assert not root.exists()


def test_nonrefined_protocol_does_not_invent_refined_metrics(manifest: Path) -> None:
    _change(manifest, "thresholds", "refine", value=False)
    report = run_sv_benchmark(manifest, manifest.parent / "matrix", runner=_runner(manifest.parent))
    assert report["status"] == "COMPLETE"
    assert all(q["refined_metrics"] == [] and q["refined_stratification_status"] == "NOT_REQUESTED"
               for q in report["queries"])
