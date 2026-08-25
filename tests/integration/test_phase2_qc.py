"""Phase 1 to Phase 2.1 lightweight input-QC handoff integration test."""

from __future__ import annotations

import json
from pathlib import Path

import yaml

from hifivar.config import load_config
from hifivar.context import AnalysisContext
from hifivar.qc import QCStatus, run_input_qc


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = (
    PROJECT_ROOT / "src" / "hifivar" / "resources" / "configs" / "default.yaml"
)
STANDARD_PRESET = (
    PROJECT_ROOT
    / "src"
    / "hifivar"
    / "resources"
    / "configs"
    / "presets"
    / "standard.yaml"
)


def test_phase2_lightweight_qc_end_to_end(tmp_path: Path) -> None:
    """Build a mixed context and write matching Unicode JSON/YAML QC reports."""
    project = tmp_path / "科研项目" / "输入数据"
    project.mkdir(parents=True)
    fasta = project / "参考.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    fastq = project / "样本一.fastq"
    fastq.write_text("@read1\nACGT\n+\nIIII\n", encoding="utf-8")
    bam = project / "样本二.bam"
    bam.write_bytes(b"alignment-placeholder")
    sheet = project / "samples.tsv"
    sheet.write_text(
        "sample_id\tinput\n"
        "S1\t样本一.fastq\n"
        "S2\t样本二.bam\n",
        encoding="utf-8",
    )
    user_config = project / "analysis.yaml"
    user_config.write_text(
        yaml.safe_dump(
            {
                "reference": {"fasta": "参考.fa", "build": "GRCh38"},
                "samples": {"sheet": "samples.tsv"},
            },
            allow_unicode=True,
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    effective = load_config(DEFAULT_CONFIG, STANDARD_PRESET, user_config)
    context = AnalysisContext.from_config(effective)
    report = run_input_qc(context)
    json_path = report.write_json(tmp_path / "结果" / "input-qc.json")
    yaml_path = report.write_yaml(tmp_path / "结果" / "input-qc.yaml")
    json_payload = json.loads(json_path.read_text(encoding="utf-8"))
    yaml_payload = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))

    assert context.sample_ids == ("S1", "S2")
    assert [result.status for result in report.sample_results] == [
        QCStatus.PASS,
        QCStatus.WARN,
    ]
    assert report.overall_status is QCStatus.WARN
    assert report.sample_results[1].issues[0].code == "ALIGNMENT_INDEX_MISSING"
    assert report.get_metric("reference_build").value == "GRCh38"
    assert json_payload == yaml_payload == report.to_dict()
    assert "输入数据" in json_path.read_text(encoding="utf-8")
    assert not Path(f"{bam}.bai").exists()
