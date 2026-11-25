"""SV-depth orchestration contracts. Synthetic tools are not accuracy evidence."""
from __future__ import annotations

import copy
import csv
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import struct
import sys
from types import SimpleNamespace
import zlib

import pytest
import yaml

from hifivar.command import CommandRunner
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.sv_depth_study import load_sv_depth_study, run_sv_depth_study


def _bgzf(data: bytes) -> bytes:
    c = zlib.compressobj(6, zlib.DEFLATED, -15)
    body = c.compress(data) + c.flush()
    size = 18 + len(body) + 8
    return (b"\x1f\x8b\x08\x04\0\0\0\0\0\xff\x06\0BC\x02\0" +
            struct.pack("<H", size - 1) + body + struct.pack("<II", zlib.crc32(data), len(data)))


_VCF = ('##fileformat=VCFv4.3\n##reference=GRCh38\n##contig=<ID=chr1,length=1000>\n'
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">\n'
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="length">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n'
        'chr1\t101\tI1\tA\t' + 'A' * 60 + '\t.\tPASS\tSVTYPE=INS;SVLEN=59\tGT\t0/1\n')


def _vcf(path):
    path.write_bytes(_bgzf(_VCF.encode()))
    Path(str(path) + ".tbi").write_bytes(_bgzf(b"TBI\x01"))


@pytest.fixture
def manifest(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + "A" * 1000 + "\n")
    Path(str(fasta) + ".fai").write_text("chr1\t1000\t6\t1000\t1001\n")
    truth = tmp_path / "truth.vcf.gz"
    _vcf(truth)
    (tmp_path / "confident.bed").write_text("chr1\t0\t1000\n")
    defaults = Path(__file__).parents[2] / "src/hifivar/resources/configs/default.yaml"
    source_sv = yaml.safe_load(defaults.read_text())["sv"]
    sv = copy.deepcopy({k: source_sv[k] for k in ("enabled", "overwrite", "sawfish", "sniffles2", "pbsv", "cutesv", "finalization")})
    for name in ("sawfish", "sniffles2", "pbsv", "cutesv"):
        sv[name]["enabled"] = name in {"sniffles2", "cutesv"}
        sv[name]["executable"] = "fake-" + name
    sv.update(enabled=True, overwrite=False)
    sv["finalization"] = {"bgzip_executable": "fake-bgzip", "tabix_executable": "fake-tabix"}
    units = []
    for nominal, fraction in [(5, .25), (20, 1)]:
        bam = tmp_path / f"{nominal}x.bam"
        bam.write_bytes(b"SYNTHETIC_BAM_" + str(nominal).encode())
        Path(str(bam) + ".bai").write_bytes(b"SYNTHETIC_INDEX")
        coverage = tmp_path / f"{nominal}x.coverage.json"
        coverage.write_text(json.dumps({"sample_id": "HG002", "alignment_bam": str(bam),
            "observed_coverage": nominal * 1.59, "reference_bases": 1000, "covered_bases": 1000}))
        units.append({"id": f"{nominal}x-r1", "nominal_depth": nominal, "replicate": 1,
                      "sampling_fraction": fraction, "bam": str(bam), "coverage_json": str(coverage)})
    payload = {"schema_version": 1, "study_id": "SYNTHETIC_DEPTH", "sample_id": "HG002",
        "reference": {"fasta": str(fasta), "build": "GRCh38"},
        "truth": {"vcf": str(truth), "confident_bed": str(tmp_path / "confident.bed"),
            "version": "synthetic1", "source": "synthetic test fixture", "reference_build": "GRCh38"},
        "tools": {"samtools_executable": "fake-samtools", "truvari_executable": "fake-truvari",
            "bcftools_executable": "fake-bcftools", "truvari_version": "5.4.0"},
        "thresholds": {"sizemin": 50, "sizemax": 50000, "refdist": 2000,
            "pctseq": .7, "pctsize": .7, "pass_only": True, "pick": "ac", "refine": True},
        "sv": sv, "design": {"sampling_unit": "pacbio_molecule", "nested": True,
            "biological_samples": 1, "full_input_nominal_depth": 20}, "units": units}
    path = tmp_path / "study.yaml"
    path.write_text(yaml.safe_dump(payload))
    return path


def _change(path, location, value):
    data = yaml.safe_load(path.read_text())
    parent = data
    for key in location[:-1]:
        parent = parent[key]
    parent[location[-1]] = value
    path.write_text(yaml.safe_dump(data))


_PROCESS = r'''
from pathlib import Path
import gzip, json, struct, sys, zlib
tool,*args=sys.argv[1:]
mode=MODE
tool=tool.removeprefix("fake-")
versions={"samtools":"samtools 1.24","bcftools":"bcftools 1.24","truvari":"Truvari v5.4.0","sniffles2":"Sniffles2 2.8.0","cutesv":"cuteSV 2.1.4","sawfish":"sawfish 2.2.1","pbsv":"pbsv 2.11.0","bgzip":"bgzip (htslib) 1.24","tabix":"tabix (htslib) 1.24"}
if args in (["version"],["--version"],["-v"]):
 print(versions[tool]);raise SystemExit(0)
def bgzf(payload):
 c=zlib.compressobj(6,zlib.DEFLATED,-15);body=c.compress(payload)+c.flush()
 size=18+len(body)+8
 return b"\x1f\x8b\x08\x04\0\0\0\0\0\xff\x06\0BC\x02\0"+struct.pack("<H",size-1)+body+struct.pack("<II",zlib.crc32(payload),len(payload))
vcf=VCF
def indexed(path):
 path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(bgzf(vcf.encode()))
 Path(str(path)+".tbi").write_bytes(bgzf(b"TBI\x01"))
if tool=="samtools":
 if args[0]=="quickcheck":raise SystemExit(0)
 if args[0]=="view" and "-H" in args:
  print("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:RG\tSM:HG002\tPL:PACBIO")
  raise SystemExit(0)
 if args[0]=="view" and "-c" in args:print("0");raise SystemExit(0)
 raise SystemExit(9)
if tool=="bcftools":
 assert args[:4]==["norm","-N","-c","e"];raise SystemExit(0)
if tool=="sniffles2":
 if mode=="fail_sniffles":
  print("deliberate test failure",file=sys.stderr);raise SystemExit(7)
 if mode=="mutate":
  p=Path(args[args.index("--input")+1]);p.write_bytes(p.read_bytes()+b"MUTATED")
 indexed(Path(args[args.index("--vcf")+1]));raise SystemExit(0)
if tool=="sawfish":
 out=Path(args[args.index("--output-dir")+1]);out.mkdir(parents=True,exist_ok=True)
 if args[0]=="joint-call":indexed(out/"genotyped.sv.vcf.gz")
 raise SystemExit(0)
if tool=="pbsv":
 out=Path(args[-1]);out.parent.mkdir(parents=True,exist_ok=True)
 out.write_bytes(b"SIG" if args[0]=="discover" else vcf.encode());raise SystemExit(0)
if tool=="cutesv":
 out=Path(args[2]);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(vcf);raise SystemExit(0)
if tool=="bgzip":sys.stdout.buffer.write(bgzf(Path(args[-1]).read_bytes()));raise SystemExit(0)
if tool=="tabix":Path(args[-1]+".tbi").write_bytes(bgzf(b"TBI\x01"));raise SystemExit(0)
assert tool=="truvari" and args[0]=="bench"
out=Path(args[args.index("-o")+1]);out.mkdir(parents=True)
values={"TP-base":9,"TP-comp":8,"FP":2,"FN":1,"base cnt":10,"comp cnt":10,"precision":.8,"recall":.9,"f1":1.44/1.7}
(out/"summary.json").write_text(json.dumps(values))
for name in ("tp-base.vcf.gz","tp-comp.vcf.gz","fp.vcf.gz","fn.vcf.gz"):indexed(out/name)
if "--refine" in args:
 (out/"refine.variant_summary.json").write_text(json.dumps(values))
 (out/"refine.region_summary.json").write_text(json.dumps({"TP":1,"TN":1}))
 for name in ("refine.base.vcf.gz","refine.comp.vcf.gz"):indexed(out/name)
print("synthetic comparator complete")
'''


class ProcessRunner:
    def __init__(self, script):
        self.script = script
        self.runner = CommandRunner()
        self.commands = []
        self.timeouts = []

    def require_executable(self, executable, **kwargs):
        return Path(sys.executable)

    def run(self, command, **kwargs):
        self.commands.append(tuple(command))
        self.timeouts.append(kwargs.get("timeout"))
        return self.runner.run([sys.executable, "-I", "-B", str(self.script), *command], **kwargs)


@pytest.fixture
def runtime(monkeypatch):
    monkeypatch.setattr("hifivar.sv_benchmark_selection._load_native_runtime", lambda: object())
    def derive(source, output, *, sample_id):
        output.mkdir(parents=True)
        query = output / "query.ins_del.vcf.gz"
        with gzip.open(source, "rt", encoding="utf-8") as handle:
            lines = [line for line in handle if line.startswith("#") or "SVTYPE=INS;" in line or "SVTYPE=DEL;" in line]
        query.write_bytes(_bgzf("".join(lines).encode()))
        shutil.copyfile(Path(str(source) + ".tbi"), Path(str(query) + ".tbi"))
        return {"status": "COMPLETE", "output_paths": {"vcf": str(query)}, "test_classifier": "synthetic_dependency"}
    monkeypatch.setattr("hifivar.sv_depth_study.derive_sv_benchmark_query", derive)


def _runner(path, mode="success"):
    script = path / ("fake_tools_" + mode + ".py")
    script.write_text(_PROCESS.replace("MODE", repr(mode)).replace("VCF", repr(_VCF + "chr1\t300\tINV1\tA\t<INV>\t.\tPASS\tSVTYPE=INV;SVLEN=70\tGT\t0/1\n")))
    return ProcessRunner(script)


def _hashes(paths):
    return {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def test_loader_nominal_observed_and_full_input_kind(manifest):
    job = load_sv_depth_study(manifest)
    assert job.caller_order == ("sniffles2", "cutesv")
    assert [(u.nominal_depth, u.observed_coverage) for u in job.units] == [(5, 7.95), (20, 31.8)]
    assert [u.replicate_kind for u in job.units] == ["NESTED_TECHNICAL_SUBSAMPLE", "FULL_INPUT_TECHNICAL_REPEAT"]
    assert job.design["biological_samples"] == 1


@pytest.mark.parametrize("location,value", [
    (("schema_version",), True), (("extra",), "forbidden"),
    (("truth", "reference_build"), "GRCh37"), (("truth", "version"), "latest"),
    (("tools", "truvari_version"), "5.3.0"), (("tools", "samtools_executable"), "samtools --bad"),
    (("thresholds", "sizemin"), 20), (("design", "biological_samples"), 3),
    (("design", "nested"), False), (("sv", "overwrite"), True),
    (("units", 0, "nominal_depth"), float("nan")), (("units", 0, "replicate"), True),
    (("units", 0, "id"), "../escape"), (("units", 0, "sampling_fraction"), 2),
    (("units", 1, "id"), "5x-r1"), (("units", 0, "sampling_fraction"), 1),
])
def test_invalid_manifest_fails_before_output(manifest, location, value):
    _change(manifest, location, value)
    with pytest.raises((InputValidationError, ValueError)):
        load_sv_depth_study(manifest)


def test_duplicate_yaml_keys_rejected(manifest):
    manifest.write_text(manifest.read_text() + "\nstudy_id: OTHER\n")
    with pytest.raises(InputValidationError, match="unique"):
        load_sv_depth_study(manifest)


@pytest.mark.parametrize("field,value", [
    ("sample_id", "OTHER"), ("alignment_bam", "elsewhere.bam"),
    ("observed_coverage", True), ("observed_coverage", float("inf")), ("reference_bases", 0),
])
def test_coverage_binding_and_positive_depth(manifest, field, value):
    data = yaml.safe_load(manifest.read_text())
    p = Path(data["units"][0]["coverage_json"])
    coverage = json.loads(p.read_text());coverage[field] = value;p.write_text(json.dumps(coverage))
    with pytest.raises(InputValidationError):
        load_sv_depth_study(manifest)


def test_dry_run_never_probes_executes_or_writes(manifest):
    class Never:
        def require_executable(self, *a, **k):pytest.fail("no probe")
        def run(self, *a, **k):pytest.fail("no process")
    out = manifest.parent / "new" / "plan"
    result = run_sv_depth_study(manifest, out, dry_run=True, unit_ids=["5x-r1"], callers=["sniffles2"], runner=Never())
    assert result["status"] == "PLANNED" and result["selected_run_ids"] == ["5x-r1/sniffles2"]
    assert len(result["runs"]) == 4 and not out.parent.exists()
    assert result["minimum_depth_claim"] == "NOT_ESTABLISHED"


@pytest.mark.parametrize("ids,callers", [([], None), (["unknown"], None), (None, []), (None, ["pav"]), (["5x-r1", "5x-r1"], None)])
def test_selection_requires_declared_unique_ids(manifest, ids, callers):
    with pytest.raises(InputValidationError):
        run_sv_depth_study(manifest, manifest.parent / "out", dry_run=True, unit_ids=ids, callers=callers)


def test_fake_process_full_matrix_and_separate_tp_semantics(manifest, runtime):
    job = load_sv_depth_study(manifest)
    before = _hashes([job.reference.fasta, job.truth["vcf"], *[u.bam for u in job.units]])
    runner = _runner(manifest.parent)
    out = manifest.parent / "out"
    result = run_sv_depth_study(manifest, out, runner=runner)
    assert result["status"] == result["selection_status"] == "COMPLETE"
    assert len(result["runs"]) == 4 and all(r["status"] == "COMPLETE" for r in result["runs"])
    for r in result["runs"]:
        assert r["metrics"]["tp_base"] == 9 and r["metrics"]["tp_call"] == 8
        assert r["metrics"]["precision"] == .8 and r["metrics"]["recall"] == .9
        assert r["summary_kind"] == "refined"
        assert Path(r["benchmark_report"]).is_file()
    assert before == _hashes([Path(p) for p in before])
    assert result["source_immutability"].startswith("PASS_FULL_SHA256")
    assert any(c[0] == "fake-bgzip" and "-c" in c for c in runner.commands)
    assert any(c[0] == "fake-cutesv" for c in runner.commands)
    assert any(c[0] == "fake-sniffles2" for c in runner.commands)
    with (out / "sv_depth_metrics.tsv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 4 and rows[0]["tp_base"] == "9" and rows[0]["tp_call"] == "8"


def test_all_four_real_wrapper_command_paths(manifest, runtime):
    data = yaml.safe_load(manifest.read_text())
    for name in ("sawfish", "pbsv"):
        data["sv"][name]["enabled"] = True
    manifest.write_text(yaml.safe_dump(data))
    result = run_sv_depth_study(manifest, manifest.parent / "out", unit_ids=["5x-r1"], runner=_runner(manifest.parent))
    assert result["selection_status"] == "COMPLETE" and result["status"] == "INCOMPLETE"
    assert [r["caller"] for r in result["runs"] if r["status"] == "COMPLETE"] == ["sawfish", "sniffles2", "pbsv", "cutesv"]


def test_bounded_execution_does_not_hash_unselected_bam(manifest, runtime, monkeypatch):
    import hifivar.sv_depth_study as module
    job = load_sv_depth_study(manifest)
    unselected = job.units[1].bam
    calls = []
    original = module.compute_sha256
    def checked(path):
        p = Path(path);calls.append(p)
        assert p != unselected, "unselected huge BAM must not be hashed"
        return original(p)
    monkeypatch.setattr(module, "compute_sha256", checked)
    result = run_sv_depth_study(manifest, manifest.parent / "out", unit_ids=["5x-r1"], callers=["sniffles2"], runner=_runner(manifest.parent))
    assert result["selection_status"] == "COMPLETE" and result["status"] == "INCOMPLETE"
    assert {r["id"] for r in result["runs"] if r["status"] == "COMPLETE"} == {"5x-r1/sniffles2"}
    assert calls.count(job.units[0].bam) == 2


def test_resume_verifies_complete_and_runs_only_pending_selection(manifest, runtime):
    out = manifest.parent / "out";runner = _runner(manifest.parent)
    first = run_sv_depth_study(manifest, out, unit_ids=["5x-r1"], callers=["sniffles2"], runner=runner)
    receipt = out / "units/5x-r1/sniffles2/unit_receipt.json"
    before = receipt.read_bytes()
    second_runner = _runner(manifest.parent, "success_again")
    second = run_sv_depth_study(manifest, out, resume=True, unit_ids=["5x-r1"], runner=second_runner)
    assert second["selection_status"] == "COMPLETE" and receipt.read_bytes() == before
    assert not any(c[0] == "fake-sniffles2" for c in second_runner.commands)
    assert any(c[0] == "fake-cutesv" for c in second_runner.commands)


def test_failed_caller_isolated_and_never_automatically_retried(manifest, runtime):
    out = manifest.parent / "out"
    first = run_sv_depth_study(manifest, out, unit_ids=["5x-r1"], runner=_runner(manifest.parent, "fail_sniffles"))
    assert first["status"] == "PARTIAL" and first["selection_status"] == "FAILED"
    assert [r["status"] for r in first["runs"][:2]] == ["FAILED", "COMPLETE"]
    failure = out / "units/5x-r1/sniffles2/failure.json";before = failure.read_bytes()
    runner = _runner(manifest.parent, "recovery_not_retry")
    second = run_sv_depth_study(manifest, out, resume=True, unit_ids=["5x-r1"], runner=runner)
    assert second["runs"][0]["status"] == "FAILED_PRESERVED_NO_AUTOMATIC_RETRY"
    assert failure.read_bytes() == before and runner.commands == []


def test_source_mutation_stops_remaining_and_invalidates_metrics(manifest, runtime):
    out = manifest.parent / "out"
    result = run_sv_depth_study(manifest, out, unit_ids=["5x-r1"], runner=_runner(manifest.parent, "mutate"))
    assert result["status"] == "FAILED" and result["source_immutability"] == "FAIL"
    assert result["metrics_valid_for_claims"] is False
    assert not (out / "units/5x-r1/cutesv").exists()
    with (out / "sv_depth_metrics.tsv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert all(not row["precision"] for row in rows)


def test_output_tampering_resume_is_failed_not_reused(manifest, runtime):
    out = manifest.parent / "out"
    run_sv_depth_study(manifest, out, unit_ids=["5x-r1"], callers=["sniffles2"], runner=_runner(manifest.parent))
    raw = out / "units/5x-r1/sniffles2/raw/HG002.sniffles2.sv.vcf.gz"
    raw.write_bytes(raw.read_bytes() + b"changed")
    runner = _runner(manifest.parent, "tamper_resume")
    with pytest.raises(OutputValidationError, match="output changed"):
        run_sv_depth_study(manifest, out, resume=True, unit_ids=["5x-r1"], callers=["sniffles2"], runner=runner)
    assert runner.commands == []


def test_existing_output_rejected_without_resume(manifest):
    out = manifest.parent / "owner";out.mkdir()
    (out / "keep.txt").write_text("preserve")
    with pytest.raises(OutputValidationError, match="fresh"):
        run_sv_depth_study(manifest, out, dry_run=True)
    assert (out / "keep.txt").read_text() == "preserve"


def test_resume_changed_manifest_rejected(manifest, runtime):
    out = manifest.parent / "out"
    run_sv_depth_study(manifest, out, unit_ids=["5x-r1"], callers=["sniffles2"], runner=_runner(manifest.parent))
    manifest.write_text(manifest.read_text() + "\n")
    with pytest.raises(OutputValidationError, match="identity"):
        run_sv_depth_study(manifest, out, resume=True, runner=_runner(manifest.parent, "changed"))


def test_optional_native_dependency_checked_before_caller_or_output(manifest, monkeypatch):
    def fail():
        raise InputValidationError("missing native classifier")
    monkeypatch.setattr("hifivar.sv_benchmark_selection._load_native_runtime", fail)
    out = manifest.parent / "out"
    with pytest.raises(InputValidationError, match="native classifier"):
        run_sv_depth_study(manifest, out)
    assert not out.exists()




def test_fraction_must_equal_nominal_ratio_not_just_increase(manifest):
    _change(manifest, ("units", 0, "sampling_fraction"), .3)
    with pytest.raises(InputValidationError, match="nominal_depth"):
        load_sv_depth_study(manifest)


def test_master_metrics_tamper_rejected_before_export_or_execution(manifest, runtime):
    out = manifest.parent / "out"
    run_sv_depth_study(manifest, out, unit_ids=["5x-r1"], callers=["sniffles2"], runner=_runner(manifest.parent))
    metrics_before = (out / "sv_depth_metrics.tsv").read_bytes()
    status = out / "study_report.json"
    data = json.loads(status.read_text())
    data["runs"][0]["metrics"]["precision"] = .1
    status.write_text(json.dumps(data))
    runner = _runner(manifest.parent, "master_tamper")
    with pytest.raises(OutputValidationError, match="sealed unit receipt"):
        run_sv_depth_study(manifest, out, resume=True, unit_ids=["5x-r1"], callers=["sniffles2"], runner=runner)
    assert runner.commands == [] and (out / "sv_depth_metrics.tsv").read_bytes() == metrics_before


def test_receipt_and_master_tamper_cannot_override_native_summary(manifest, runtime):
    out = manifest.parent / "out"
    run_sv_depth_study(manifest, out, unit_ids=["5x-r1"], callers=["sniffles2"], runner=_runner(manifest.parent))
    report_path = out / "study_report.json"
    receipt_path = out / "units/5x-r1/sniffles2/unit_receipt.json"
    report = json.loads(report_path.read_text())
    receipt = json.loads(receipt_path.read_text())
    for data in (report["runs"][0], receipt):
        data["metrics"]["precision"] = .1
    report_path.write_text(json.dumps(report));receipt_path.write_text(json.dumps(receipt))
    runner = _runner(manifest.parent, "both_tampered")
    with pytest.raises(OutputValidationError, match="canonical native"):
        run_sv_depth_study(manifest, out, resume=True, unit_ids=["5x-r1"], callers=["sniffles2"], runner=runner)
    assert runner.commands == []


def test_bounded_resume_checks_unselected_completed_outputs(manifest, runtime):
    out = manifest.parent / "out"
    result = run_sv_depth_study(manifest, out, runner=_runner(manifest.parent))
    assert result["status"] == "COMPLETE"
    unselected = out / "units/20x-r1/cutesv/raw/HG002.cutesv.sv.vcf.gz"
    unselected.write_bytes(unselected.read_bytes() + b"edited")
    runner = _runner(manifest.parent, "unselected_tamper")
    with pytest.raises(OutputValidationError, match="output changed"):
        run_sv_depth_study(manifest, out, resume=True, unit_ids=["5x-r1"], callers=["sniffles2"], runner=runner)
    assert runner.commands == []


def test_owned_lock_is_removed_on_terminal_and_foreign_lock_preserved(manifest, runtime):
    out = manifest.parent / "out"
    run_sv_depth_study(manifest, out, unit_ids=["5x-r1"], callers=["sniffles2"], runner=_runner(manifest.parent))
    lock = out / "study.lock";assert not lock.exists()
    lock.write_text("manual concurrent/stale marker")
    runner = _runner(manifest.parent, "lock")
    with pytest.raises(OutputValidationError, match="lock"):
        run_sv_depth_study(manifest, out, resume=True, unit_ids=["5x-r1"], callers=["cutesv"], runner=runner)
    assert lock.read_text() == "manual concurrent/stale marker" and runner.commands == []


def test_existing_incomplete_unit_is_preserved_not_retried(manifest, runtime):
    out = manifest.parent / "out"
    run_sv_depth_study(manifest, out, unit_ids=["5x-r1"], callers=["sniffles2"], runner=_runner(manifest.parent))
    partial = out / "units/5x-r1/cutesv";partial.mkdir()
    owner = partial / "keep.txt";owner.write_text("owner incomplete file")
    runner = _runner(manifest.parent, "incomplete")
    result = run_sv_depth_study(manifest, out, resume=True, unit_ids=["5x-r1"], callers=["cutesv"], runner=runner)
    assert result["runs"][1]["status"] == "FAILED_PRESERVED_NO_AUTOMATIC_RETRY"
    assert runner.commands == [] and owner.read_text() == "owner incomplete file"
    assert not (out / "study.lock").exists()


def test_manifest_change_between_parse_and_hash_fails_before_output(manifest, runtime, monkeypatch):
    import hifivar.sv_depth_study as module
    original = module.load_sv_depth_study
    def changed(path):
        result = original(path)
        manifest.write_text(manifest.read_text() + "\n")
        return result
    monkeypatch.setattr(module, "load_sv_depth_study", changed)
    out = manifest.parent / "out"
    with pytest.raises(module.SourceImmutabilityError, match="Manifest changed"):
        run_sv_depth_study(manifest, out, runner=_runner(manifest.parent))
    assert not out.exists()


@pytest.mark.parametrize("extra", [
    "chr1\t50\tSNV1\tA\tC\t.\tPASS\t.\tGT\t0/1\n",
    "chr1\t55\tSTAR1\tA\t*\t.\tPASS\t.\tGT\t0/1\n",
])
def test_original_mixed_truth_is_not_prefiltered(manifest, runtime, extra):
    job = load_sv_depth_study(manifest)
    job.truth["vcf"].write_bytes(_bgzf((_VCF + extra).encode()))
    original = job.truth["vcf"].read_bytes()
    runner = _runner(manifest.parent)
    result = run_sv_depth_study(manifest, manifest.parent / "out", unit_ids=["5x-r1"], callers=["sniffles2"], runner=runner)
    assert result["selection_status"] == "COMPLETE"
    assert result["truth_prefilter"] == "NONE_ORIGINAL_MIXED_TRUTH_PRESERVED"
    assert job.truth["vcf"].read_bytes() == original
    bench = [c for c in runner.commands if c[0] == "fake-truvari" and c[1] == "bench"]
    assert bench and all(c[c.index("-b") + 1] == str(job.truth["vcf"]) for c in bench)


def test_raw_all_type_caller_vcf_and_insdel_query_remain_separate(manifest, runtime):
    out = manifest.parent / "out"
    result = run_sv_depth_study(manifest, out, unit_ids=["5x-r1"], callers=["sniffles2"], runner=_runner(manifest.parent))
    assert result["selection_status"] == "COMPLETE"
    raw = out / "units/5x-r1/sniffles2/raw/HG002.sniffles2.sv.vcf.gz"
    query = out / "units/5x-r1/sniffles2/query/query.ins_del.vcf.gz"
    with gzip.open(raw, "rt", encoding="utf-8") as handle:raw_text = handle.read()
    with gzip.open(query, "rt", encoding="utf-8") as handle:query_text = handle.read()
    assert "SVTYPE=INV;" in raw_text and "SVTYPE=INV;" not in query_text
    assert "SVTYPE=INS;" in query_text and raw != query


def test_per_step_timeout_comes_from_declared_caller_resources(manifest, runtime):
    _change(manifest, ("sv", "sniffles2", "runtime_minutes"), 1)
    runner = _runner(manifest.parent)
    result = run_sv_depth_study(manifest, manifest.parent / "out", unit_ids=["5x-r1"], callers=["sniffles2"], runner=runner)
    assert result["selection_status"] == "COMPLETE"
    assert all(isinstance(v, (int, float)) and 0 < v <= 60 for v in runner.timeouts)
    probes = [timeout for command, timeout in zip(runner.commands, runner.timeouts) if len(command) > 1 and command[1] in {"--version", "version"}]
    assert probes and all(v <= 30 for v in probes)
    assert result["runs"][0]["per_command_timeout_seconds"] == 60


def test_resume_holds_lock_before_report_read_and_interleaved_resume(manifest, runtime, monkeypatch):
    import hifivar.sv_depth_study as module
    out = manifest.parent / "out"
    run_sv_depth_study(manifest, out, unit_ids=["5x-r1"], callers=["sniffles2"], runner=_runner(manifest.parent))
    original_read = Path.read_text
    observed = []
    def guarded_read(path, *args, **kwargs):
        if path == out / "study_report.json":
            assert (out / "study.lock").is_file(), "lock must precede mutable report read"
            observed.append("report_read_locked")
        return original_read(path, *args, **kwargs)
    monkeypatch.setattr(Path, "read_text", guarded_read)
    original = module._unchanged
    interleaved = []
    def interleave(snapshot, *, full=False):
        if full and not interleaved:
            interleaved.append("attempted")
            with pytest.raises(OutputValidationError, match="lock"):
                run_sv_depth_study(manifest, out, resume=True, unit_ids=["5x-r1"], callers=["cutesv"], runner=_runner(manifest.parent, "interleave_A"))
        return original(snapshot, full=full)
    monkeypatch.setattr(module, "_unchanged", interleave)
    result = run_sv_depth_study(manifest, out, resume=True, unit_ids=["5x-r1"], callers=["sniffles2"], runner=_runner(manifest.parent, "interleave_B"))
    assert observed and interleaved and result["selection_status"] == "COMPLETE"
    assert result["runs"][1]["status"] == "PENDING"
    assert not (out / "units/5x-r1/cutesv").exists()
    assert not (out / "study.lock").exists()


@pytest.mark.parametrize("location,value", [
    (("sv", "sniffles2", "threads"), 1.5),
    (("sv", "sniffles2", "threads"), True),
    (("sv", "sniffles2", "minimum_support"), True),
    (("sv", "cutesv", "minimum_support"), "10"),
    (("sv", "cutesv", "genotype"), "false"),
    (("sv", "sawfish", "disable_cnv"), "false"),
    (("sv", "cutesv", "diff_ratio_merging_ins"), True),
    (("sv", "cutesv", "minimum_suport"), 5),
    (("sv", "unconsumed_harmonization"), {"enabled": True}),
])
def test_manifest_rejects_silent_caller_parameter_coercion(manifest, location, value):
    _change(manifest, location, value)
    with pytest.raises(InputValidationError):
        load_sv_depth_study(manifest)


def test_coverage_change_between_parse_and_hash_invalidates_without_calling(manifest, runtime, monkeypatch):
    import hifivar.sv_depth_study as module
    original = module.load_sv_depth_study
    def changed(path):
        job = original(path)
        coverage = job.units[0].coverage_json
        payload = json.loads(coverage.read_text())
        payload["observed_coverage"] = 999.0
        coverage.write_text(json.dumps(payload))
        return job
    monkeypatch.setattr(module, "load_sv_depth_study", changed)
    runner = _runner(manifest.parent, "coverage_toctou")
    report = run_sv_depth_study(manifest, manifest.parent / "out", unit_ids=["5x-r1"], callers=["sniffles2"], runner=runner)
    assert report["status"] == "FAILED" and report["metrics_valid_for_claims"] is False
    assert report["source_immutability"] == "FAIL" and "Coverage changed" in report["error"]
    assert runner.commands == []


def test_fai_change_between_reference_parse_and_hash_rejected(manifest, runtime, monkeypatch):
    import hifivar.sv_depth_study as module
    original = module.load_sv_depth_study
    def changed(path):
        job = original(path)
        job.reference.fai.write_text(job.reference.fai.read_text() + "\n")
        return job
    monkeypatch.setattr(module, "load_sv_depth_study", changed)
    runner = _runner(manifest.parent, "fai_toctou")
    with pytest.raises(module.SourceImmutabilityError, match="FAI changed"):
        run_sv_depth_study(manifest, manifest.parent / "out", runner=runner)
    assert runner.commands == []


def test_actual_selection_membership_is_persisted_under_same_report(manifest, runtime):
    out = manifest.parent / "out"
    first = run_sv_depth_study(manifest, out, unit_ids=["5x-r1"], callers=["sniffles2"], runner=_runner(manifest.parent))
    assert first["selected_run_ids"] == ["5x-r1/sniffles2"]
    second = run_sv_depth_study(manifest, out, resume=True, unit_ids=["5x-r1"], callers=["cutesv"], runner=_runner(manifest.parent, "scope_resume"))
    assert second["selected_run_ids"] == ["5x-r1/cutesv"] and second["selection_status"] == "COMPLETE"
    stored = json.loads((out / "study_report.json").read_text())
    assert stored["selected_run_ids"] == second["selected_run_ids"]
