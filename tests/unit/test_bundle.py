from __future__ import annotations
import json
from pathlib import Path, PurePosixPath
import pytest
from hifivar.bundle import BundleItem, ReproducibilityRecord, create_release_bundle, selected_report_items
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.report import FinalRunReport, FinalStatus, ReportArtifact, TrackReport


def report(tmp_path: Path) -> FinalRunReport:
    vcf=tmp_path/"S1.small.vcf.gz"; vcf.write_bytes(b"selected-vcf")
    bam=tmp_path/"S1.bam"; bam.write_bytes(b"large-alignment")
    track=TrackReport("small","3","small",FinalStatus.COMPLETE,True,
        (ReportArtifact(vcf,"small-vcf","vcf","S1",sha256="configured",selected_for_bundle=True),
         ReportArtifact(bam,"alignment","bam","S1",selected_for_bundle=True)),("S1",))
    return FinalRunReport("R1","deadbeef",{"build":"GRCh38"},({"sample_id":"S1"},),(track,),{"password":"do-not-store"})


def test_release_bundle_selection_large_pointer_and_reproducibility(tmp_path: Path) -> None:
    item_report=report(tmp_path); sheet=tmp_path/"samples.tsv"; sheet.write_text("sample_id\nS1\n")
    repro=ReproducibilityRecord({"hifivar":"0.0.1.dev0"},(("tool","--token","s3cr3t"),),
        {"API_TOKEN":"hidden","python":"3.12"},{"build":"GRCh38"},("s3cr3t",),sheet)
    result=create_release_bundle(item_report,tmp_path/"bundle",items=selected_report_items(item_report),reproducibility=repro)
    assert len(result.copied)==1 and result.copied[0].name.endswith(".vcf.gz")
    assert len(result.pointers)==1 and not any(path.suffix==".bam" for path in result.copied)
    manifest=json.loads(result.manifest.read_text(encoding="utf-8")); assert manifest["large_primary_data_included"] is False
    config=(result.root/"configs"/"effective_config.yaml").read_text(); assert "do-not-store" not in config and "'***'" in config
    commands=(result.root/"provenance"/"commands.json").read_text(); assert "s3cr3t" not in commands and "***" in commands
    environment=(result.root/"provenance"/"environment.json").read_text(); assert "hidden" not in environment
    assert (result.root/"provenance"/"sample_sheet.tsv").is_file()


def test_bundle_rejects_escape_duplicate_and_nonempty_destination(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError): BundleItem(tmp_path/"x",PurePosixPath("../x"),"x")
    source=tmp_path/"x.tsv"; source.write_text("x")
    item=BundleItem(source,PurePosixPath("results/x.tsv"),"table")
    with pytest.raises(InputValidationError): create_release_bundle(report(tmp_path),tmp_path/"bundle",items=(item,item))
    occupied=tmp_path/"occupied"; occupied.mkdir(); (occupied/"user.txt").write_text("keep")
    with pytest.raises(OutputValidationError): create_release_bundle(report(tmp_path),occupied)
    with pytest.raises(OutputValidationError): create_release_bundle(report(tmp_path),occupied,overwrite=True)
    assert (occupied/"user.txt").read_text()=="keep"
