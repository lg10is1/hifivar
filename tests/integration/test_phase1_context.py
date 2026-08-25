"""Phase 1 seal test from effective config through portable manifests."""

from __future__ import annotations

import json
from pathlib import Path

import yaml

from hifivar.config import load_config
from hifivar.context import AnalysisContext
from hifivar.manifest import RunManifest


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


def test_phase1_config_context_manifest_end_to_end(tmp_path: Path) -> None:
    """Seal Phase 1 with tiny Unicode-path reference and FASTQ inputs."""
    project = tmp_path / "科研项目" / "输入数据"
    project.mkdir(parents=True)
    fasta = project / "参考.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    fastq = project / "样本.fastq"
    fastq.write_text("@read1\nACGT\n+\nIIII\n", encoding="utf-8")
    sheet = project / "samples.tsv"
    sheet.write_text("sample_id\tinput\nS1\t样本.fastq\n", encoding="utf-8")
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
    manifest = RunManifest.from_context(context)
    json_path = manifest.write_json(tmp_path / "结果" / "run-manifest.json")
    yaml_path = manifest.write_yaml(tmp_path / "结果" / "run-manifest.yaml")

    json_payload = json.loads(json_path.read_text(encoding="utf-8"))
    yaml_payload = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
    assert context.reference.fasta == fasta
    assert context.sample_ids == ("S1",)
    assert context.samples[0].sample.input.files == (fastq,)
    assert json_payload == yaml_payload == manifest.to_dict()
    assert json_payload["inputs"][0]["sha256"] is None
    assert json_payload["reference"]["sha256"] is None
    assert Path(json_payload["reference"]["fasta"]).is_absolute()
    assert Path(json_payload["source_sample_sheet"]).is_absolute()
