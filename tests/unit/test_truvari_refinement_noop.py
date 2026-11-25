"""Strict no-op state/provenance contracts; synthetic fixtures are not accuracy data.

Most cases inject a labelled partition verifier because core Windows has no
pysam. The separate native-reader test uses actual BGZF/TBI when pysam exists.
"""
from __future__ import annotations

import copy
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from hifivar import truvari_refinement as noop
from hifivar.command import CommandResult
from hifivar.exceptions import OutputValidationError
from hifivar.reference import ReferenceGenome
from hifivar.truvari import TruvariRequest, TruvariThresholds, TruvariWrapper, TruvariResult, TruvariResultStatus
from hifivar.validation import compute_sha256
from hifivar.sv_depth_study import _verify_receipt, _output_checksums

POLICY = TruvariThresholds(refdist=2000, pctseq=.7, pctsize=.7, sizemin=50,
                          sizemax=50000, pass_only=True, pick="ac", chunksize=5000, refine=True)


def _write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def _text_vcf(count, sample="HG002"):
    header = ('##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n'
              '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + sample + '\n')
    return header + "".join(f"chr1\t{100+i}\tV{i}\tA\tAT\t.\tPASS\tSVTYPE=INS\tGT\t0/1\n" for i in range(count))


def _params(request):
    result = {"base": str(request.base_vcf.absolute()), "comp": str(request.comparison_vcf.absolute()),
              "reference": str(request.reference.fasta.absolute()), "output": str(request.output_directory.absolute()),
              "includebed": str(request.confident_regions.absolute()) if request.confident_regions else None,
              "bSample": request.sample_id, "cSample": request.sample_id, "skip_gt": False,
              "passonly": request.thresholds.pass_only, "sizefilt": 30}
    for key, value in request.thresholds.to_dict().items():
        if key not in {"pass_only", "refine"} and value is not None:
            result[key] = value
    return result


def _log(request, params):
    refined = {"benchdir": str(request.output_directory.absolute()),
               "reference": str(request.reference.fasta.absolute()),
               "regions": str(request.output_directory.absolute() / "candidate.refine.bed"),
               "use_original_vcfs": False, "coords": "R"}
    return ("[INFO] Truvari v5.4.0\n[INFO] Params:\n" + json.dumps(params) +
            "\n[INFO] Finished bench\n[INFO] Params:\n" + json.dumps(refined) +
            "\n[INFO] 0 regions to be refined\n[CRITICAL] No regions to be refined. Exiting\n")


def _emit(request, log):
    out = request.output_directory
    out.mkdir(parents=True, exist_ok=False)
    params = _params(request)
    summary = {"TP-base": 3, "TP-comp": 2, "FP": 0, "FN": 0,
               "base cnt": 3, "comp cnt": 2, "precision": 1., "recall": 1., "f1": 1.}
    _write_json(out / "params.json", params)
    _write_json(out / "summary.json", summary)
    (out / "candidate.refine.bed").write_text("chr1\t0\t1000\n" * 5)
    (out / "refine.log.txt").write_bytes(b"")
    for name, count in (("tp-base", 3), ("tp-comp", 2), ("fp", 0), ("fn", 0)):
        with gzip.open(out / (name + ".vcf.gz"), "wt") as stream:
            stream.write(_text_vcf(count))
        (out / (name + ".vcf.gz.tbi")).write_bytes(b"UNIT_PARTITION_READER_TEST_DOUBLE")
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(_log(request, params), encoding="utf-8")


def _partition_test_double(vcf, index, sample_id):
    """Tests proof state only; this is intentionally not a native TBI validator."""
    if index.read_bytes() != b"UNIT_PARTITION_READER_TEST_DOUBLE":
        raise OutputValidationError("Synthetic index verifier rejected corruption")
    with gzip.open(vcf, "rt") as stream:
        text = stream.read()
    assert text.split("#CHROM", 1)[1].splitlines()[0].endswith(sample_id)
    records = [line for line in text.splitlines() if line and not line.startswith("#")]
    return {"records": len(records), "indexed_records": len(records),
            "record_stream_sha256": hashlib.sha256(("\n".join(records) + ("\n" if records else "")).encode()).hexdigest()}


@pytest.fixture
def bundle(tmp_path, monkeypatch):
    fasta = tmp_path / "ref.fa";fasta.write_text(">chr1\n" + "A"*1000 + "\n")
    Path(str(fasta)+".fai").write_text("chr1\t1000\t6\t1000\t1001\n")
    truth = tmp_path / "truth.vcf.gz";truth.write_bytes(b"synthetic protected truth")
    Path(str(truth)+".tbi").write_bytes(b"synthetic protected index")
    unit = tmp_path / "unit";query = unit / "query/query.ins_del.vcf.gz";query.parent.mkdir(parents=True)
    query.write_bytes(b"synthetic protected query");Path(str(query)+".tbi").write_bytes(b"synthetic protected index")
    bed = tmp_path/"regions.bed";bed.write_text("chr1\t0\t1000\n")
    request = TruvariRequest("HG002", ReferenceGenome.from_fasta(fasta, build="GRCh38"),
                            truth, query, unit/"benchmark/queries/sniffles2",
                            confident_regions=bed, thresholds=POLICY)
    log = unit/"benchmark/logs/current.stderr.log"
    _emit(request, log)
    monkeypatch.setattr(noop, "_partition_counts", _partition_test_double)
    return SimpleNamespace(request=request, log=log, unit=unit)


def _proof(bundle, **updates):
    kwargs = dict(version="5.4.0", returncode=0, executed=True, stderr_path=bundle.log, output_absent_before=True)
    kwargs.update(updates)
    return noop.verify_noop_refinement(bundle.request.to_dict(), TruvariWrapper().plan_command(bundle.request), **kwargs)


def _result(bundle):
    proof = _proof(bundle);summary = bundle.request.output_directory/"summary.json"
    result = TruvariResult(bundle.request, TruvariResultStatus.COMPLETED,
        TruvariWrapper().plan_command(bundle.request), "5.4.0", 0., summary,
        refinement_status="NOT_NEEDED", refinement_proof=proof,
        effective_summary_path=summary, effective_summary_sha256=compute_sha256(summary))
    return result.to_dict()


def test_verified_positive_tp_noop_is_stable_readonly_and_keeps_bench(bundle):
    before = {str(p): p.read_bytes() for p in bundle.unit.rglob("*") if p.is_file()}
    one, two = _proof(bundle), _proof(bundle)
    assert one == two and one["refinement_status"] == "NOT_NEEDED"
    assert one["partition_counts"]["tp-base"]["records"] == 3
    assert one["partition_counts"]["tp-comp"]["records"] == 2
    assert one["effective_summary"]["kind"] == "bench"
    assert before == {str(p): p.read_bytes() for p in bundle.unit.rglob("*") if p.is_file()}
    result = _result(bundle)
    assert result["refinement_requested"] is True
    assert result["refined_summary_path"] is None
    assert noop.verify_noop_result(result) == one


@pytest.mark.parametrize("change", [
    {"returncode": 1}, {"returncode": False}, {"executed": False}, {"version": "5.3.0"},
    {"version": "5.4.1"}, {"output_absent_before": False}, {"stderr_path": None},
])
def test_current_execution_evidence_is_mandatory(bundle, change):
    with pytest.raises(OutputValidationError):
        _proof(bundle, **change)


@pytest.mark.parametrize("field,value", [
    ("base", "/other/base.vcf.gz"), ("comp", "/other/query.vcf.gz"),
    ("output", "/other/output"), ("reference", "/other/ref.fa"), ("includebed", None),
    ("bSample", "OTHER"), ("cSample", "OTHER"), ("skip_gt", True),
    ("refdist", 500), ("pctseq", .5), ("pctsize", .5), ("passonly", False), ("pick", "single"),
])
def test_native_parameter_binding_is_required(bundle, field, value):
    p=bundle.request.output_directory/"params.json";data=json.loads(p.read_text());data[field]=value
    _write_json(p,data)
    with pytest.raises(OutputValidationError):
        _proof(bundle)


@pytest.mark.parametrize("name", [
    "refine.variant_summary.json", "refine.region_summary.json",
    "refine.base.vcf.gz", "refine.base.vcf.gz.tbi",
    "refine.comp.vcf.gz", "refine.comp.vcf.gz.tbi",
    "refine.base.vcf.gz.partial", "phab.output.vcf.gz", "phab_partial",
])
def test_any_partial_refine_or_phab_artifact_is_failure(bundle, name):
    (bundle.request.output_directory/name).write_bytes(b"partial")
    with pytest.raises(OutputValidationError):
        _proof(bundle)


@pytest.mark.parametrize("mode", ["missing_log", "one_marker", "wrong_hook", "wrong_logged_params", "wrong_version", "traceback"])
def test_both_current_native_log_context_and_zero_region_markers_required(bundle, mode):
    text=bundle.log.read_text()
    if mode=="missing_log":bundle.log.unlink()
    else:
        if mode=="one_marker":text=text.replace("[INFO] 0 regions to be refined\n","")
        elif mode=="wrong_hook":text=text.replace('"use_original_vcfs": false','"use_original_vcfs": true')
        elif mode=="wrong_logged_params":text=text.replace('"bSample": "HG002"','"bSample": "OTHER"')
        elif mode=="wrong_version":text=text.replace("Truvari v5.4.0","Truvari v5.3.0")
        else:text+="Traceback: synthetic failure\n"
        bundle.log.write_text(text)
    with pytest.raises((OutputValidationError, FileNotFoundError)):
        _proof(bundle)


@pytest.mark.parametrize("mode", ["fp", "fn", "empty_both", "missing_totals", "wrong_rate", "boolean_count"])
def test_native_summary_must_be_positive_tp_perfect_with_exact_totals_and_rates(bundle, mode):
    path=bundle.request.output_directory/"summary.json";s=json.loads(path.read_text())
    if mode=="fp":
        s.update(FP=1, **{"comp cnt":3, "precision":2/3,"f1":.8})
    elif mode=="fn":
        s.update(FN=1, **{"base cnt":4, "recall":.75,"f1":6/7})
    elif mode=="empty_both":
        s={k:0 for k in ("TP-base","TP-comp","FP","FN","base cnt","comp cnt")}
        s.update(precision=None,recall=None,f1=None)
    elif mode=="missing_totals":s.pop("base cnt")
    elif mode=="wrong_rate":s["precision"]=.5
    else:s["TP-base"]=True
    _write_json(path,s)
    with pytest.raises(OutputValidationError):
        _proof(bundle)


def test_partition_counts_and_index_validity_cannot_be_replaced_by_summary(bundle):
    out=bundle.request.output_directory
    with gzip.open(out/"tp-base.vcf.gz","wt") as f:f.write(_text_vcf(2))
    with pytest.raises(OutputValidationError,match="partition count"):
        _proof(bundle)
    with gzip.open(out/"tp-base.vcf.gz","wt") as f:f.write(_text_vcf(3))
    (out/"fp.vcf.gz.tbi").write_bytes(b"corrupt")
    with pytest.raises(OutputValidationError,match="index"):
        _proof(bundle)


@pytest.mark.parametrize("field,value", [
    ("status","planned"), ("status","FORGED"),
    ("refinement_status","COMPLETE"), ("refinement_requested",False),
    ("effective_summary_kind","refined"), ("effective_summary_sha256","0"*64),
    ("refined_summary_path","fake.refine.json"), ("refinement_proof",None),
])
def test_stored_flags_do_not_grant_noop_success(bundle,field,value):
    result=_result(bundle);result[field]=value
    with pytest.raises(OutputValidationError):
        noop.verify_noop_result(result)


def test_proof_hash_and_summary_tampering_is_rejected(bundle):
    result=_result(bundle)
    proof=result["refinement_proof"];name=next(iter(proof["proof_files"]))
    proof["proof_files"][name]["sha256"]="0"*64
    with pytest.raises(OutputValidationError,match="proof"):
        noop.verify_noop_result(result)
    result=_result(bundle)
    path=bundle.request.output_directory/"summary.json";s=json.loads(path.read_text());s["TP-base"]=2;s["base cnt"]=2;_write_json(path,s)
    with pytest.raises(OutputValidationError,match="partition count"):
        noop.verify_noop_result(result)


def _seal_study_bundle(bundle):
    req=bundle.request;result=_result(bundle);root=bundle.unit
    protocol={"schema_version":1,"benchmark_id":"fixture","sample_id":"HG002",
        "reference":{"fasta":str(req.reference.fasta),"build":"GRCh38"},
        "truth":{"vcf":str(req.base_vcf),"confident_bed":str(req.confident_regions),"version":"synthetic","source":"unit fixture","reference_build":"GRCh38"},
        "tools":{"truvari_executable":"truvari","bcftools_executable":"bcftools","truvari_version":"5.4.0"},
        "thresholds":POLICY.to_dict(),"queries":[{"id":"sniffles2","vcf":str(req.comparison_vcf),"caller":"sniffles2"}]}
    _write_json(root/"benchmark_manifest.json",protocol)
    matrix={"status":"COMPLETE","commands":[{"command":result["command"],"return_code":0,"executed":True,"stderr_path":str(bundle.log)}],
        "queries":[{"id":"sniffles2","status":"COMPLETE","command":result["command"],"result":result,
        "refinement_status":"NOT_NEEDED","refinement_proof":result["refinement_proof"],"refined_metrics":[]}]}
    _write_json(root/"benchmark/benchmark_matrix.json",matrix)
    entry={"id":"u/sniffles2","unit_id":"u","caller":"sniffles2","status":"COMPLETE","manifest_sha256":"a"*64,
        "nominal_depth":5.,"observed_coverage":7.95,"replicate":1,"sampling_fraction":.25,
        "replicate_kind":"NESTED_TECHNICAL_SUBSAMPLE","coverage_reference_bases":1000,"coverage_sha256":"b"*64,
        "metrics":{"tp_base":3,"tp_call":2,"fp":0,"fn":0,"precision":1.,"recall":1.,"f1":1.},
        "summary_kind":"bench","benchmark_report":str(root/"benchmark/benchmark_matrix.json"),
        "caller_provenance":{"fixture":True},"derivation":{"fixture":True},
        "precision_formula":"TP_COMP/(TP_COMP+FP)","recall_formula":"TP_BASE/(TP_BASE+FN)"}
    for k in ("refinement_requested","refinement_status","refinement_proof","effective_summary_kind","effective_summary_path","effective_summary_sha256"):
        entry[k]=result[k]
    receipt={**entry,"output_checksums":_output_checksums(root)}
    _write_json(root/"unit_receipt.json",receipt)
    return entry,receipt


def test_study_receipt_revalidates_noop_without_fabricating_refined_files(bundle):
    entry,_=_seal_study_bundle(bundle)
    assert _verify_receipt(bundle.unit,entry,"a"*64)["summary_kind"]=="bench"
    assert not (bundle.request.output_directory/"refine.variant_summary.json").exists()


@pytest.mark.parametrize("mode",["master_flag","joint_flag","joint_summary","missing_proof","partial_after_seal"])
def test_resume_cannot_be_tricked_by_flags_receipt_or_native_summary(bundle,mode):
    entry,receipt=_seal_study_bundle(bundle)
    if mode=="master_flag":entry["refinement_status"]="COMPLETE"
    elif mode=="joint_flag":
        entry["refinement_requested"]=receipt["refinement_requested"]=False
    elif mode=="missing_proof":
        entry["refinement_proof"]=receipt["refinement_proof"]=None
    elif mode=="partial_after_seal":
        (bundle.request.output_directory/"phab_partial").write_bytes(b"x")
    else:
        path=bundle.request.output_directory/"summary.json";s=json.loads(path.read_text());s["TP-base"]=2;s["base cnt"]=2;_write_json(path,s)
        receipt["output_checksums"]=_output_checksums(bundle.unit)
        entry["metrics"]=dict(entry["metrics"],tp_base=2);receipt["metrics"]=entry["metrics"]
    _write_json(bundle.unit/"unit_receipt.json",receipt)
    with pytest.raises(OutputValidationError):
        _verify_receipt(bundle.unit,entry,"a"*64)


def test_wrapper_returns_truthful_not_needed_with_no_result_rewrite(bundle,tmp_path):
    # Fresh output owned by a labelled synthetic runner, not a native process.
    req=TruvariRequest(bundle.request.sample_id,bundle.request.reference,bundle.request.base_vcf,
        bundle.request.comparison_vcf,tmp_path/"fresh",confident_regions=bundle.request.confident_regions,thresholds=POLICY)
    log=tmp_path/"fresh.stderr.log"
    class Runner:
        def require_executable(self,exe):return Path(exe)
        def run(self,command,**kw):
            if command[1] in {"version","--version"}:
                return SimpleNamespace(stdout="Truvari v5.4.0",stderr="",returncode=0)
            _emit(req,log)
            return SimpleNamespace(args=tuple(command),returncode=0,executed=True,duration_seconds=.01,stderr_path=log)
    r=TruvariWrapper(runner=Runner()).run(req)
    assert r.refinement_status=="NOT_NEEDED" and r.refined_summary_path is None
    assert r.to_dict()["effective_summary_kind"]=="bench"
    assert json.loads((req.output_directory/"summary.json").read_text())["TP-comp"]==2


def test_native_pysam_partition_reader_uses_actual_bgzf_and_index(tmp_path):
    pysam=pytest.importorskip("pysam")
    plain=tmp_path/"source.vcf";plain.write_text(_text_vcf(3),encoding="utf-8")
    vcf=tmp_path/"source.vcf.gz"
    pysam.tabix_compress(str(plain),str(vcf),force=False)
    pysam.tabix_index(str(vcf),preset="vcf",force=False)
    index=Path(str(vcf)+".tbi")
    result=noop._partition_counts(vcf,index,"HG002")
    assert result["records"]==result["indexed_records"]==3
    index.write_bytes(b"invalid native index")
    with pytest.raises(OutputValidationError):
        noop._partition_counts(vcf,index,"HG002")


@pytest.mark.parametrize("requested,bad_status",[
    (True,"NOT_REQUESTED"),(True,"PLANNED"),(True,"FORGED_STATE"),
    (False,"COMPLETE"),(False,"PLANNED"),(False,"NOT_NEEDED")])
def test_explicit_normal_refinement_state_must_match_protocol(bundle,requested,bad_status):
    entry,receipt=_normal_study_bundle(bundle,requested)
    entry["refinement_status"]=receipt["refinement_status"]=bad_status
    _write_json(bundle.unit/"unit_receipt.json",receipt)
    with pytest.raises(OutputValidationError,match="state conflicts"):
        _verify_receipt(bundle.unit,entry,"a"*64)


@pytest.mark.parametrize("requested",[False,True])
def test_historical_missing_normal_refinement_state_remains_compatible(bundle,requested):
    entry,receipt=_normal_study_bundle(bundle,requested)
    for data in (entry,receipt):
        for k in ("refinement_requested","refinement_status","refinement_proof",
                  "effective_summary_kind","effective_summary_path","effective_summary_sha256"):
            data.pop(k,None)
    _write_json(bundle.unit/"unit_receipt.json",receipt)
    assert _verify_receipt(bundle.unit,entry,"a"*64)["summary_kind"]==("refined" if requested else "bench")


def _normal_study_bundle(bundle,requested):
    entry,receipt=_seal_study_bundle(bundle)
    root=bundle.unit;control=root/"benchmark_manifest.json"
    protocol=json.loads(control.read_text());protocol["thresholds"]["refine"]=requested;_write_json(control,protocol)
    kind="refined" if requested else "bench"
    summary=bundle.request.output_directory/("refine.variant_summary.json" if requested else "summary.json")
    if requested:summary.write_bytes((bundle.request.output_directory/"summary.json").read_bytes())
    for data in (entry,receipt):
        data.update(summary_kind=kind,refinement_requested=requested,
            refinement_status="COMPLETE" if requested else "NOT_REQUESTED",refinement_proof=None,
            effective_summary_kind=kind,effective_summary_path=str(summary),effective_summary_sha256=compute_sha256(summary))
    receipt["output_checksums"]=_output_checksums(root)
    _write_json(root/"unit_receipt.json",receipt)
    return entry,receipt


@pytest.mark.parametrize("mode",[
    "row_failed","result_planned","record_nonzero","record_unexecuted",
    "record_other_log","record_missing","record_duplicate","row_other_command"])
def test_noop_resume_requires_completed_matching_current_invocation(bundle,mode):
    entry,receipt=_seal_study_bundle(bundle)
    path=bundle.unit/"benchmark/benchmark_matrix.json";matrix=json.loads(path.read_text())
    row=matrix["queries"][0];record=matrix["commands"][0]
    if mode=="row_failed":row["status"]="FAILED"
    elif mode=="result_planned":row["result"]["status"]="planned"
    elif mode=="record_nonzero":record["return_code"]=1
    elif mode=="record_unexecuted":record["executed"]=False
    elif mode=="record_other_log":record["stderr_path"]="different.stderr.log"
    elif mode=="record_missing":matrix["commands"]=[]
    elif mode=="record_duplicate":matrix["commands"].append(dict(record))
    else:row["command"]=["other","bench"]
    _write_json(path,matrix);receipt["output_checksums"]=_output_checksums(bundle.unit)
    _write_json(bundle.unit/"unit_receipt.json",receipt)
    with pytest.raises(OutputValidationError):
        _verify_receipt(bundle.unit,entry,"a"*64)
