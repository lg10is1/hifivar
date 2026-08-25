from __future__ import annotations
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
