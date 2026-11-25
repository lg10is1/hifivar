"""Synthetic software contracts for explicit four-caller depth support analysis."""
from __future__ import annotations

import copy
import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys

import pytest
import yaml

from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.sv_depth_study import _snapshot, load_sv_depth_study, run_sv_depth_study
from hifivar.sv_depth_support import load_sv_depth_support_plan, run_sv_depth_support
from hifivar.sv_input_identity import POLICY, POLICY_SHA256

_spec = importlib.util.spec_from_file_location("_depth_producer_test_helpers", Path(__file__).with_name("test_sv_depth_study.py"))
_d = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_d)

_JASMINE = r'''
if tool=="jasmine":
 if mode=="fail_low" and "5x-r1" in " ".join(args):
  print("deliberate synthetic merge failure",file=sys.stderr);raise SystemExit(7)
 if mode=="mutate_source":
  p=Path(MUTATE);p.write_bytes(p.read_bytes()+b"MUTATED")
 values=dict(x.split("=",1) for x in args if "=" in x)
 paths=Path(values["file_list"]).read_text(encoding="utf-8").splitlines()
 sources=[]
 for path in paths:
  lines=Path(path).read_text(encoding="utf-8").splitlines()
  sources.append([x.split("\t") for x in lines if not x.startswith("#")])
 headers=[x for x in Path(paths[0]).read_text(encoding="utf-8").splitlines() if x.startswith("##") and not x.startswith("##hifivarSourceIdentity=")]
 headers+=['##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Sources">',
 '##INFO=<ID=SUPP,Number=1,Type=Integer,Description="Sources">',
 '##INFO=<ID=IDLIST,Number=.,Type=String,Description="IDs">']
 headers+=['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002']
 if mode=="empty_high":
  groups=[[(i,row)] for i,rows in enumerate(sources) for row in rows]
 else:
  groups=[[(i,rows[0]) for i,rows in enumerate(sources)],
          [(0,sources[0][1]),(1,sources[1][1])],[(2,sources[2][1])],[(3,sources[3][1])]]
 out=[]
 for n,group in enumerate(groups,1):
  row=group[0][1].copy()
  vector="".join("1" if any(i==j for i,_ in group) else "0" for j in range(4))
  row[2]="J"+str(n);row[7]+=";SUPP_VEC="+vector+";SUPP="+str(vector.count("1"))+";IDLIST="+",".join(r[2] for _,r in group)
  out.append("\t".join(row))
 Path(values["out_file"]).write_text("\n".join(headers+out)+"\n")
 raise SystemExit(0)
if tool=="bgzip":
 if args[0]=="-t":raise SystemExit(0)
 sys.stdout.buffer.write(bgzf(Path(args[-1]).read_bytes()));raise SystemExit(0)
if tool=="tabix":
 if args[0]=="-l":print("chr1");raise SystemExit(0)
 if args[0]=="-p":
  Path(args[-1]+".tbi").write_bytes(bgzf(b"TBI\x01"));raise SystemExit(0)
 contig,interval=args[1].split(":");start,end=map(int,interval.split("-"))
 with gzip.open(args[0],"rb") as handle:
  for line in handle:
   if not line.startswith(b"#"):
    parts=line.split(b"\t")
    if parts[0].decode()==contig and start<=int(parts[1])<=end:sys.stdout.buffer.write(line)
 raise SystemExit(0)
if tool=="truvari" and args[0]=="bench":
 source=Path(args[args.index("-c")+1])
 with gzip.open(source,"rt") as handle:empty=not any(not line.startswith("#") for line in handle)
 if empty:
  out=Path(args[args.index("-o")+1]);out.mkdir(parents=True)
  values={"TP-base":0,"TP-comp":0,"FP":0,"FN":10,"base cnt":10,"comp cnt":0,"precision":None,"recall":0.0,"f1":None}
  (out/"summary.json").write_text(json.dumps(values))
  for name in ("tp-base.vcf.gz","tp-comp.vcf.gz","fp.vcf.gz","fn.vcf.gz"):indexed(out/name)
  print("No regions to be refined. Exiting",file=sys.stderr)
  raise SystemExit(0)
'''


class _SupportRunner(_d.ProcessRunner):
    def __init__(self, script, root):
        super().__init__(script)
        self.tools = root / "tools/bin"
        self.tools.mkdir(parents=True, exist_ok=True)
        for name in ("jasmine","samtools","bcftools","truvari","sniffles2","cutesv","sawfish","pbsv","bgzip","tabix"):
            path = self.tools / ("fake-" + name)
            if not path.exists():path.write_text("SYNTHETIC_EXECUTABLE_IDENTITY_ONLY\n")
            path.chmod(0o700)

    def require_executable(self, executable, **kwargs):
        return self.tools / Path(executable).name

    def run(self, command, **kwargs):
        command = list(command)
        name = Path(command[0]).name
        if name.lower() in {"bash", "bash.exe"} and len(command) > 1 and Path(command[1]).name == "fake-jasmine":
            command = command[1:];name = "fake-jasmine"
        if name.startswith("fake-"):command[0] = name
        self.commands.append(tuple(command))
        self.timeouts.append(kwargs.get("timeout"))
        return self.runner.run([sys.executable, "-X", "utf8", "-I", "-B", str(self.script), *command], **kwargs)


def _runner(root, mode="success", mutation=None):
    code = _d._PROCESS.replace('"samtools":"samtools 1.24"', '"jasmine":"Jasmine 1.1.5","samtools":"samtools 1.24"')
    code = code.replace('if tool=="samtools":', _JASMINE.replace("Path(MUTATE)", "Path(" + repr(str(mutation)) + ")") + '\nif tool=="samtools":')
    code = code.replace("MODE", repr(mode)).replace("VCF", repr(_d._VCF + "chr1\t300\tINV1\tA\t<INV>\t.\tPASS\tSVTYPE=INV;SVLEN=70\tGT\t0/1\n"))
    script = root / ("support_fake_" + mode + ".py");script.write_text(code)
    return _SupportRunner(script, root)


@pytest.fixture
def completed(tmp_path, monkeypatch):
    manifest = _d.manifest.__wrapped__(tmp_path)
    data = yaml.safe_load(manifest.read_text())
    for caller in ("sawfish", "pbsv"):data["sv"][caller]["enabled"] = True
    manifest.write_text(yaml.safe_dump(data))
    _d.runtime.__wrapped__(monkeypatch)
    original = __import__("hifivar.sv_depth_study", fromlist=["derive_sv_benchmark_query"]).derive_sv_benchmark_query

    def derive(source, output, *, sample_id):
        value = original(source, output, sample_id=sample_id)
        query = Path(value["output_paths"]["vcf"])
        value.update(source_vcf=str(source), input_checksums=_snapshot((source, Path(str(source)+".tbi"))))
        value["output_paths"].update(index=str(query)+".tbi", manifest=str(output/"selection_manifest.json"))
        (output/"selection_manifest.json").write_text(json.dumps(value))
        return value
    monkeypatch.setattr("hifivar.sv_depth_study.derive_sv_benchmark_query", derive)
    monkeypatch.setattr("hifivar.sv_depth_support.derive_sv_benchmark_query", derive)
    output = tmp_path/"producer"
    result = run_sv_depth_study(manifest, output, runner=_runner(tmp_path))
    assert result["status"] == "COMPLETE" and len(result["runs"]) == 8
    policy = {"enabled": True, "policy_id": POLICY["policy_id"], "policy_version": POLICY["policy_version"],
              "policy_sha256": POLICY_SHA256, "mode": POLICY["mode"], "timeout_seconds": 300}
    plan = {"schema_version": 1, "study_manifest": str(manifest),
            "source_roots": {u["unit_id"]: {c: str(output) for c in ("sawfish","sniffles2","pbsv","cutesv")} for u in result["runs"]},
            "harmonization": {"jasmine_executable": "fake-jasmine", "bgzip_executable": "fake-bgzip", "tabix_executable": "fake-tabix",
                              "max_dist": 1000, "distance_type": "linear", "input_identity_policy": policy},
            "support_thresholds": [1,2,3,4]}
    path = tmp_path/"support.yaml";path.write_text(yaml.safe_dump(plan))
    return path, manifest, output, tmp_path


def _write(path, data):
    path.write_text(yaml.safe_dump(data))


def _hash_tree(root):
    return {p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in root.rglob("*") if p.is_file()}


def test_real_producer_receipts_then_actual_merge_selection_benchmark(completed):
    plan, manifest, producer, tmp = completed
    before = _hash_tree(producer)
    runner = _runner(tmp, "normal_support")
    result = run_sv_depth_support(plan, tmp/"response", runner=runner)
    assert result["status"] == result["selection_status"] == "COMPLETE"
    assert result["metrics_valid_for_claims"] is True
    assert result["caller_order"] == ["sawfish","sniffles2","pbsv","cutesv"]
    assert result["truth_prefilter"] == "NONE_ORIGINAL_MIXED_TRUTH_PRESERVED"
    assert len([c for c in runner.commands if c[0]=="fake-jasmine" and len(c)>2]) == 2
    assert len([c for c in runner.commands if c[:2]==("fake-truvari","bench")]) == 8
    benchmark = json.loads((tmp/"response/benchmark_manifest.json").read_text())
    original = yaml.safe_load(manifest.read_text())
    assert benchmark["truth"] == original["truth"]
    assert benchmark["thresholds"] == load_sv_depth_study(manifest).thresholds.to_dict()
    assert len(benchmark["queries"]) == 8
    for unit in result["units"]:
        assert [r["selection"]["selected_records"] for r in unit["thresholds"]] == [4,2,1,1]
        for row in unit["thresholds"]:
            assert row["status"] == "COMPLETE" and row["summary_kind"] == "refined"
            assert row["metrics"]["tp_base"] == 9 and row["metrics"]["tp_call"] == 8
            assert row["metrics"]["precision"] == .8 and row["metrics"]["recall"] == .9
            receipt = json.loads(Path(row["receipt_path"]).read_text())
            assert receipt["unit_id"] == unit["unit_id"] and receipt["support_ge"] == row["support_ge"]
            assert receipt["output_checksums"] and receipt["caller_order"] == result["caller_order"]
        assert unit["observed_coverage"] == unit["nominal_depth"] * 1.59
    with (tmp/"response/support_study_summary.tsv").open(newline="") as f:rows=list(csv.DictReader(f,delimiter="\t"))
    assert len(rows)==8 and rows[0]["tp_base"]=="9" and rows[0]["tp_call"]=="8"
    assert _hash_tree(producer) == before
    assert result["bam_integrity"].startswith("INHERITED_FULL_SHA_WITH_CURRENT_STAT")
    assert result["minimum_depth_claim"] == "NOT_ESTABLISHED"


def test_dry_run_reads_readiness_but_no_native_probe_or_output(completed):
    plan, _, producer, tmp = completed
    class Forbidden:
        def run(self,*a,**k):pytest.fail("no native run")
        def require_executable(self,*a,**k):pytest.fail("no executable probe")
    before=_hash_tree(producer)
    result=run_sv_depth_support(plan,tmp/"dry",dry_run=True,runner=Forbidden())
    assert result["selection_status"]=="PLANNED" and not (tmp/"dry").exists()
    assert _hash_tree(producer)==before


def test_subset_is_not_full_parent_completion(completed):
    plan,_,_,tmp=completed
    result=run_sv_depth_support(plan,tmp/"subset",unit_ids=["5x-r1"],runner=_runner(tmp,"subset"))
    assert result["selection_status"]=="COMPLETE" and result["status"]=="INCOMPLETE"
    assert result["metrics_valid_for_claims"] is False
    assert result["selected_unit_ids"]==["5x-r1"]


def test_declared_source_subset_is_not_full_parent(completed):
    plan,_,_,tmp=completed
    d=yaml.safe_load(plan.read_text());d["source_roots"].pop("20x-r1");_write(plan,d)
    result=run_sv_depth_support(plan,tmp/"declared",runner=_runner(tmp,"declared"))
    assert result["status"]=="INCOMPLETE" and result["selection_status"]=="COMPLETE"
    assert result["units"][1]["status"]=="NOT_DECLARED"


@pytest.mark.parametrize("case",["caller_failed","owner_mismatch","global_invalid","global_nonboolean","global_output_failed","active_lock","raw_unsealed","raw_tampered","wrong_depth_request","wrong_derivation_source","missing_source_selection"])
def test_unaccepted_source_never_merges_or_creates_accuracy(completed,case):
    plan,_,producer,tmp=completed
    master=producer/"study_report.json";report=json.loads(master.read_text())
    row=report["runs"][0];cell=producer/"units/5x-r1/sawfish";rp=cell/"unit_receipt.json";receipt=json.loads(rp.read_text())
    if case=="caller_failed":row["status"]="FAILED"
    elif case=="owner_mismatch":(producer/"study_owner.json").write_text("{}")
    elif case=="global_invalid":report["metrics_valid_for_claims"]=False
    elif case=="global_nonboolean":report["metrics_valid_for_claims"]="true"
    elif case=="global_output_failed":report["output_validation"]="FAIL"
    elif case=="active_lock":(producer/"study.lock").write_text("{}")
    elif case=="raw_unsealed":receipt["output_checksums"].pop(str(Path("raw")/"HG002.sawfish.sv.vcf.gz"))
    elif case=="raw_tampered":(cell/"raw/HG002.sawfish.sv.vcf.gz").write_bytes(b"TAMPER")
    elif case=="wrong_depth_request":
        row["caller_provenance"]["caller"]["request"]["alignment"]["path"]=str(tmp/"20x.bam")
        receipt["caller_provenance"]=copy.deepcopy(row["caller_provenance"])
    elif case=="wrong_derivation_source":
        row["derivation"]["source_vcf"]=str(tmp/"different.vcf.gz");receipt["derivation"]=copy.deepcopy(row["derivation"])
    elif case=="missing_source_selection":receipt["output_checksums"].pop(str(Path("query")/"selection_manifest.json"))
    master.write_text(json.dumps(report));rp.write_text(json.dumps(receipt))
    runner=_runner(tmp,"blocked_"+case)
    result=run_sv_depth_support(plan,tmp/"blocked",unit_ids=["5x-r1"],runner=runner)
    assert result["selection_status"]=="FAILED"
    if case=="active_lock":
        assert result["status"]=="FAILED"
        assert result["units"][0]["status"]=="INVALIDATED_SOURCE_MUTATION"
        assert result["source_immutability"]=="FAIL"
    else:
        assert result["units"][0]["status"]=="BLOCKED_REQUIRED_CALLERS"
    assert runner.commands==[]
    assert all(r["metrics"] is None for r in result["units"][0]["thresholds"])


def test_benchmark_only_recovery_is_not_complete_caller_receipt(completed):
    plan,_,producer,tmp=completed
    p=producer/"study_report.json";d=json.loads(p.read_text());d["runs"][0]["status"]="FAILED"
    d["runs"][0]["benchmark_only_recovery"]={"status":"COMPLETE","metrics":{"precision":1}}
    p.write_text(json.dumps(d))
    result=run_sv_depth_support(plan,tmp/"recovery",unit_ids=["5x-r1"],runner=_runner(tmp,"recovery"))
    assert result["selection_status"]=="FAILED" and result["units"][0]["status"]=="BLOCKED_REQUIRED_CALLERS"


def test_independent_merge_failure_preserves_other_unit(completed):
    plan,_,producer,tmp=completed;before=_hash_tree(producer)
    result=run_sv_depth_support(plan,tmp/"isolated",runner=_runner(tmp,"fail_low"))
    assert result["selection_status"]=="FAILED" and result["status"]=="PARTIAL"
    assert result["units"][0]["status"]=="FAILED_HARMONIZATION"
    assert result["units"][1]["status"]=="COMPLETE"
    assert all(r["metrics"] is None for r in result["units"][0]["thresholds"])
    assert _hash_tree(producer)==before


def test_shared_mutation_invalidates_all_metrics(completed):
    plan,manifest,_,tmp=completed
    truth=Path(yaml.safe_load(manifest.read_text())["truth"]["vcf"])
    result=run_sv_depth_support(plan,tmp/"mutated",runner=_runner(tmp,"mutate_source",truth))
    assert result["status"]==result["selection_status"]=="FAILED"
    assert result["source_immutability"]=="FAIL" and result["metrics_valid_for_claims"] is False
    assert all(r["metrics"] is None for u in result["units"] for r in u["thresholds"])

    assert all(u["status"] != "COMPLETE" for u in result["units"] if u["unit_id"] in result["selected_unit_ids"])


def test_empty_high_support_does_not_fake_refined_metrics(completed):
    plan,_,_,tmp=completed
    result=run_sv_depth_support(plan,tmp/"empty",unit_ids=["5x-r1"],runner=_runner(tmp,"empty_high"))
    unit=result["units"][0]
    assert unit["thresholds"][0]["status"]=="COMPLETE"
    for row in unit["thresholds"][1:]:
        assert row["selection"]["selected_records"]==0
        assert row["status"]=="FAILED_BENCHMARK" and row["metrics"] is None
        assert "summary_kind" not in row
    assert result["metrics_valid_for_claims"] is False


@pytest.mark.parametrize("key,value",[("support_thresholds",[1,2,3]),("support_thresholds",[True,2,3,4]),("unknown",1)])
def test_strict_plan_rejected_before_outputs(completed,key,value):
    plan,_,_,tmp=completed
    d=yaml.safe_load(plan.read_text());d[key]=value;_write(plan,d)
    with pytest.raises(InputValidationError):run_sv_depth_support(plan,tmp/"invalid")
    assert not (tmp/"invalid").exists()


def test_disabled_identity_rejected(completed):
    plan,_,_,tmp=completed
    d=yaml.safe_load(plan.read_text());d["harmonization"]["input_identity_policy"]={"enabled":False};_write(plan,d)
    with pytest.raises(InputValidationError):load_sv_depth_support_plan(plan)


@pytest.mark.parametrize("location",["inside_source","existing","input_parent"])
def test_output_alias_protection(completed,location):
    plan,_,producer,tmp=completed
    target=producer/"new_response" if location=="inside_source" else tmp if location=="existing" else tmp.parent
    with pytest.raises((InputValidationError,OutputValidationError)):run_sv_depth_support(plan,target)
    assert not (producer/"new_response").exists()


def test_fresh_output_beside_reference_is_allowed(completed):
    plan,_,_,tmp=completed
    result=run_sv_depth_support(plan,tmp/"fresh_beside_reference",dry_run=True)
    assert result["selection_status"]=="PLANNED"


def test_large_bam_not_freshly_hashed_by_postprocessing(completed,monkeypatch):
    plan,manifest,_,tmp=completed
    import hifivar.sv_depth_support as mod
    original=mod.compute_sha256
    bams={u.bam for u in load_sv_depth_study(manifest).units}
    def checked(path):
        assert Path(path) not in bams
        return original(path)
    monkeypatch.setattr(mod,"compute_sha256",checked)
    result=mod.run_sv_depth_support(plan,tmp/"no_bam_hash",unit_ids=["5x-r1"],runner=_runner(tmp,"nohash"))
    assert result["selection_status"]=="COMPLETE"


def test_no_resume_or_overwrite(completed):
    plan,_,_,tmp=completed
    out=tmp/"response_once"
    run_sv_depth_support(plan,out,unit_ids=["5x-r1"],runner=_runner(tmp,"once"))
    before=_hash_tree(out)
    with pytest.raises(InputValidationError):run_sv_depth_support(plan,out,unit_ids=["5x-r1"])
    assert _hash_tree(out)==before


def test_duplicate_json_control_keys_are_rejected(tmp_path):
    from hifivar.sv_depth_support import _ReadSet
    p = tmp_path / "duplicate.json";p.write_text('{"status":"COMPLETE","status":"FAILED"}')
    with pytest.raises(InputValidationError, match="Duplicate JSON"):
        _ReadSet().json(p)


def test_native_global_integrity_failure_stays_failed(completed, monkeypatch):
    import hifivar.sv_depth_support as mod
    plan, _, _, tmp = completed;original = mod.run_sv_benchmark
    def failed(*args, **kwargs):
        value = original(*args, **kwargs)
        value.update(status="FAILED", source_immutability="FAIL", metrics_valid_for_claims=False)
        (Path(args[1]) / "benchmark_matrix.json").write_text(json.dumps(value))
        return value
    monkeypatch.setattr(mod, "run_sv_benchmark", failed)
    result = mod.run_sv_depth_support(plan, tmp/"native_global_fail", unit_ids=["5x-r1"], runner=_runner(tmp,"globalfail"))
    assert result["status"] == result["selection_status"] == "FAILED"
    assert result["source_immutability"] == "FAIL"
    assert all(r["metrics"] is None for u in result["units"] for r in u["thresholds"])


def test_native_summary_change_after_metric_read_is_detected(completed, monkeypatch):
    import os
    import hifivar.sv_depth_support as mod
    plan, _, _, tmp = completed;original = mod._metrics;changed = []
    def mutate(path, **kwargs):
        metrics = original(path, **kwargs)
        if not changed:
            p = Path(path);before = p.stat();raw = p.read_bytes();changed_raw = raw.replace(b'"TP-base": 9', b'"TP-base": 8', 1)
            assert changed_raw != raw and len(changed_raw) == len(raw);p.write_bytes(changed_raw)
            os.utime(p, ns=(before.st_atime_ns,before.st_mtime_ns));changed.append(str(p))
        return metrics
    monkeypatch.setattr(mod, "_metrics", mutate)
    result = mod.run_sv_depth_support(plan, tmp/"summary_mutation", unit_ids=["5x-r1"], runner=_runner(tmp,"summarychange"))
    assert changed and result["status"] == "FAILED"
    assert result["selection_metrics_valid"] is False
    assert all(r["metrics"] is None for u in result["units"] for r in u["thresholds"])

    assert all(u["status"] != "COMPLETE" for u in result["units"] if u["unit_id"] in result["selected_unit_ids"])


def test_static_admission_reason_does_not_echo_source_text(completed):
    plan,_,producer,tmp=completed
    p=producer/"study_owner.json"
    p.write_text('{"schema_version":1,"study_id":"PRIVATE_SOURCE_MARKER_74fc","manifest_sha256":"bad"}')
    result=run_sv_depth_support(plan,tmp/"safe_reason",unit_ids=["5x-r1"],runner=_runner(tmp,"safe"))
    row=result["units"][0]
    assert row["error_code"] == "SUPPORT_CONTRACT_REJECTED" and "owner" in row["error_reason"]
    assert "PRIVATE_SOURCE_MARKER_74fc" not in json.dumps(result)


def test_output_oserror_does_not_claim_source_mutation(completed, monkeypatch):
    import hifivar.sv_depth_support as mod
    plan, _, producer, tmp = completed;before = _hash_tree(producer)
    def unavailable(*args, **kwargs):
        raise OSError("PRIVATE_IO_MARKER_85cb")
    monkeypatch.setattr(mod,"run_sv_benchmark",unavailable)
    result=mod.run_sv_depth_support(plan,tmp/"output_io",unit_ids=["5x-r1"],runner=_runner(tmp,"outputio"))
    assert result["status"] == result["selection_status"] == "FAILED"
    assert result["source_immutability"].startswith("PASS_")
    assert all(r["metrics"] is None for u in result["units"] for r in u["thresholds"])
    assert _hash_tree(producer) == before
    assert "PRIVATE_IO_MARKER_85cb" not in json.dumps(result)

    assert all(u["status"] != "COMPLETE" for u in result["units"] if u["unit_id"] in result["selected_unit_ids"])


def test_native_returned_output_seal_is_authoritative(completed, monkeypatch):
    import os
    import hifivar.sv_depth_support as mod
    plan, _, _, tmp = completed;original = mod.run_sv_benchmark
    def changed(*args, **kwargs):
        result=original(*args,**kwargs)
        p=next((Path(args[1])/"queries").glob("*/refine.variant_summary.json"))
        old=p.stat();raw=p.read_bytes();new=raw.replace(b'"TP-base": 9',b'"TP-base": 8',1)
        assert new != raw and len(new)==len(raw)
        p.write_bytes(new);os.utime(p,ns=(old.st_atime_ns,old.st_mtime_ns))
        return result
    monkeypatch.setattr(mod,"run_sv_benchmark",changed)
    result=mod.run_sv_depth_support(plan,tmp/"native_seal",unit_ids=["5x-r1"],runner=_runner(tmp,"seal"))
    assert result["status"] == "FAILED" and result["output_validation"] == "FAIL"
    assert result["error_code"] == "SUPPORT_CONTRACT_REJECTED"
    assert all(r["metrics"] is None for u in result["units"] for r in u["thresholds"])


def test_caller_version_drift_across_ready_units_blocks_merge(completed):
    plan,_,producer,tmp=completed
    master=producer/"study_report.json";value=json.loads(master.read_text())
    row=next(r for r in value["runs"] if r["id"]=="20x-r1/sniffles2")
    row["caller_provenance"]["caller"]["tool_version"]="99.1"
    receipt=producer/"units/20x-r1/sniffles2/unit_receipt.json";saved=json.loads(receipt.read_text())
    saved["caller_provenance"]=copy.deepcopy(row["caller_provenance"])
    master.write_text(json.dumps(value));receipt.write_text(json.dumps(saved))
    runner=_runner(tmp,"drift")
    result=run_sv_depth_support(plan,tmp/"version_drift",runner=runner)
    assert result["selection_status"] == "FAILED" and runner.commands == []
    assert result["caller_versions"]["sniffles2"] is None
    assert result["caller_version_observations"]["sniffles2"] == ["2.8.0","99.1"]
    assert all(u["status"]=="BLOCKED_CALLER_VERSION_MISMATCH" for u in result["units"])
    assert all(r["metrics"] is None for u in result["units"] for r in u["thresholds"])
