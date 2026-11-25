from __future__ import annotations
import html
import json
from pathlib import Path
import pytest
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.report import FinalRunReport, FinalStatus, ReportArtifact, ToolRecord, TrackReport


def _report(tmp_path: Path, states=(FinalStatus.COMPLETE, FinalStatus.DISABLED)) -> FinalRunReport:
    artifact=ReportArtifact(tmp_path/"S1.small.vcf.gz","small-vcf","vcf","S1",selected_for_bundle=True)
    tracks=(TrackReport("small","3","small",states[0],True,(artifact,),("S1",),qc={"status":"PASS"}),
            TrackReport("assembly","7","assembly",states[1],False))
    return FinalRunReport("RUN-1","abc123",{"build":"GRCh38","sha256":"refsum"},
        ({"sample_id":"S1"},),tracks,{"api_token":"secret","paths":{"outdir":"results"}},
        (ToolRecord("DeepVariant","1.10.0","run_deepvariant","apptainer","PASS"),),
        provenance={"commands_recorded":True,"service_secret":"hidden-provenance"})


def test_final_report_serialization_and_offline_human_reports(tmp_path: Path) -> None:
    report=_report(tmp_path); root=tmp_path/"reports"
    report.write_json(root/"report.json"); report.write_yaml(root/"report.yaml")
    report.write_markdown(root/"report.md"); report.write_html(root/"report.html")
    payload=json.loads((root/"report.json").read_text(encoding="utf-8"))
    assert payload["schema_version"]=="1.0" and payload["status"]=="COMPLETE"
    assert payload["config"]["api_token"]=="***"
    assert payload["provenance"]["service_secret"]=="***"
    assert payload["interpretation_policy"]["clinical_interpretation_performed"] is False
    markdown=(root/"report.md").read_text(encoding="utf-8")
    assert "Small variants" in markdown and "Benchmark" in markdown and "clinical interpretation" in markdown
    html=(root/"report.html").read_text(encoding="utf-8")
    assert "<!doctype html>" in html and "http://" not in html and "https://" not in html


@pytest.mark.parametrize("states,expected",[
    ((FinalStatus.PARTIAL,FinalStatus.DISABLED),FinalStatus.PARTIAL),
    ((FinalStatus.NOT_RUN,FinalStatus.DISABLED),FinalStatus.NOT_RUN),
    ((FinalStatus.FAILED,FinalStatus.DISABLED),FinalStatus.FAILED),
])
def test_final_report_partial_and_failed_states(tmp_path: Path,states,expected) -> None:
    assert _report(tmp_path,states).status is expected


def test_track_enabled_status_contract_and_report_overwrite(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError): TrackReport("x","1","small",FinalStatus.DISABLED,True)
    report=_report(tmp_path); path=tmp_path/"report.json"; report.write_json(path)
    with pytest.raises(OutputValidationError): report.write_json(path)

def test_report_unicode_paths_and_html_escaping(tmp_path: Path) -> None:
    item=_report(tmp_path/"路径")
    path=tmp_path/"报告"/"结果.html"; item.write_html(path)
    assert path.is_file() and "HiFiVar run report" in path.read_text(encoding="utf-8")

def test_complete_plus_not_run_aggregates_to_partial(tmp_path: Path) -> None:
    complete=TrackReport("small","3","small",FinalStatus.COMPLETE,True)
    not_run=TrackReport("benchmark","13","benchmark",FinalStatus.NOT_RUN,True)
    item=FinalRunReport("R",None,{"build":"GRCh38"},(),(complete,not_run),{})
    assert item.status is FinalStatus.PARTIAL


@pytest.mark.parametrize("status,category,phase,message,warnings", [
    (FinalStatus.COMPLETE, "review", "10",
     "IGV execution is technical evidence only; human review remains NOT_REVIEWED.", ()),
    (FinalStatus.PARTIAL, "small", "3",
     "Expected unit receipts or source materialization are incomplete.",
     ("Missing sample HG002 receipt; no completion inferred.",)),
    (FinalStatus.FAILED, "sv", "4",
     "Native caller failed with exit 9; retained output is not accepted.",
     ("Inspect the recorded stderr before an authorized retry.",)),
])
def test_html_preserves_track_status_explanations(status, category, phase, message, warnings) -> None:
    track = TrackReport("selected_track", phase, category, status, True,
                        message=message, warnings=warnings)
    report = FinalRunReport("R", None, {"build": "GRCh38"}, (), (track,), {})
    rendered = report._html()
    detail = rendered.split("<div class='track-details'>", 1)[1].split("</div>", 1)[0]
    assert "<h3>selected_track</h3>" in detail
    assert html.escape(message) in detail
    assert all(html.escape(warning) in detail for warning in warnings)
    assert f"class='status {status.value.lower()}'>{status.value}</td>" in rendered
    if category == "review":
        assert "NOT_REVIEWED" in detail
    assert report.status is status


def test_html_track_details_escape_text_and_keep_categories_associated() -> None:
    review_name = '复核 <script>track()</script> & "A"'
    review_message = '人工未判定 <img src=x onerror="review()"> & 中文说明'
    review_warning = '保留原始 VCF <tag> & "review-warning"'
    cohort_name = '批量结果 <cohort> & "B"'
    cohort_message = '缺失样本不得当作参考基因型 <missing> & 中文'
    cohort_warning = 'cohort-only <script>warning()</script> & "quoted"'
    tracks = (
        TrackReport(review_name, '10<&"\'>', "review", FinalStatus.COMPLETE, True,
                    message=review_message, warnings=(review_warning,)),
        TrackReport(cohort_name, "12", "cohort", FinalStatus.PARTIAL, True,
                    message=cohort_message, warnings=(cohort_warning,)),
    )
    report = FinalRunReport("R", None, {"build": "GRCh38"}, (), tracks, {})
    before_dict, before_markdown, before_status = report.to_dict(), report._markdown(), report.status
    rendered = report._html()
    details = [part.split("</div>", 1)[0] for part in rendered.split("<div class='track-details'>")[1:]]
    assert len(details) == 2
    for detail, own, other in (
        (details[0], (review_name, review_message, review_warning), (cohort_name, cohort_message, cohort_warning)),
        (details[1], (cohort_name, cohort_message, cohort_warning), (review_name, review_message, review_warning)),
    ):
        assert all(html.escape(value) in detail for value in own)
        assert all(html.escape(value) not in detail for value in other)
    assert "<script>" not in rendered and "<img " not in rendered and "<missing>" not in rendered
    assert html.escape(tracks[0].phase) in rendered
    assert report.to_dict() == before_dict
    assert report._markdown() == before_markdown
    assert report.status is before_status


def test_html_absent_track_details_do_not_create_empty_placeholders() -> None:
    tracks = (
        TrackReport("none", "3", "small", FinalStatus.COMPLETE, True),
        TrackReport("empty", "4", "sv", FinalStatus.PARTIAL, True,
                    message="", warnings=("",)),
    )
    report = FinalRunReport("R", None, {"build": "GRCh38"}, (), tracks, {})
    rendered = report._html()
    assert "track-details" not in rendered
    assert "<h3>" not in rendered and "<li></li>" not in rendered
    assert "<table>" in rendered and "COMPLETE" in rendered and "PARTIAL" in rendered
