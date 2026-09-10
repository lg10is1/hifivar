from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest
import yaml

from hifivar.config import load_config, validate_config
from hifivar.exceptions import ConfigurationError


ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "src" / "hifivar" / "resources" / "configs" / "default.yaml"


def _enabled() -> dict[str, object]:
    return {
        "enabled": True,
        "study_id": "HG002-low-depth",
        "sample_id": "HG002",
        "input_mode": "aligned_bam",
        "source_bam": "aligned/HG002.bam",
        "source_fastqs": [],
        "baseline_coverage": 30,
        "target_coverages": [5, 10, 20, 30],
        "replicates": 3,
        "base_seed": 1729,
        "profiles": [
            {"name": "baseline", "settings": {}},
            {"name": "candidate-a", "settings": {"small.model_type": "PACBIO"}},
        ],
        "thresholds": [
            {"variant_class": "SNP", "metric": "f1", "minimum": 0.98},
            {"variant_class": "INDEL", "metric": "f1", "minimum": 0.90},
        ],
        "minimum_successful_replicates": 3,
        "minimum_pass_fraction": 1.0,
        "samtools_executable": "samtools",
        "threads": 4,
        "memory_mb": 8000,
        "runtime_minutes": 120,
        "overwrite": False,
    }


def test_default_low_depth_branch_is_disabled_and_valid() -> None:
    loaded = load_config(DEFAULT)
    section = loaded["low_depth"]
    assert section["enabled"] is False
    assert section["thresholds"] == []


def test_enabled_low_depth_config_is_valid_and_path_expands(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump({"low_depth": _enabled()}), encoding="utf-8")
    loaded = load_config(DEFAULT, user_config=config_path)
    assert loaded["low_depth"]["study_id"] == "HG002-low-depth"
    assert len(loaded["low_depth"]["profiles"]) == 2


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("target_coverages", [10, 5], "target_coverages"),
        ("target_coverages", [5, 31], "target_coverages"),
        ("minimum_successful_replicates", 4, "cannot exceed"),
        ("thresholds", [], "incomplete"),
        ("minimum_pass_fraction", 0, "minimum_pass_fraction"),
    ],
)
def test_invalid_low_depth_config_is_rejected(field: str, value: object, message: str) -> None:
    section = _enabled()
    section[field] = value
    with pytest.raises(ConfigurationError, match=message):
        validate_config({"low_depth": section})


def test_duplicate_profiles_and_thresholds_are_rejected() -> None:
    section = _enabled()
    section["profiles"] = [deepcopy(section["profiles"][0]), deepcopy(section["profiles"][0])]
    with pytest.raises(ConfigurationError, match="unique"):
        validate_config({"low_depth": section})


def test_study_identifier_cannot_escape_output_root() -> None:
    section = _enabled()
    section["study_id"] = "../outside"
    with pytest.raises(ConfigurationError, match="study_id"):
        validate_config({"low_depth": section})

    section = _enabled()
    section["thresholds"] = [deepcopy(section["thresholds"][0]), deepcopy(section["thresholds"][0])]
    with pytest.raises(ConfigurationError, match="unique"):
        validate_config({"low_depth": section})


def test_raw_fastq_mode_requires_unique_sources_and_not_a_bam() -> None:
    section = _enabled()
    section.update(
        {
            "input_mode": "raw_fastq",
            "source_bam": None,
            "source_fastqs": ["reads/movie1.fastq.gz", "reads/movie2.fastq.gz"],
        }
    )
    validate_config({"low_depth": section})
    section["source_fastqs"] = []
    with pytest.raises(ConfigurationError, match="source_fastqs"):
        validate_config({"low_depth": section})
    section["source_fastqs"] = ["reads.fastq.gz", "reads.fastq.gz"]
    with pytest.raises(ConfigurationError, match="unique"):
        validate_config({"low_depth": section})


def test_end_to_end_requires_raw_inputs_reference_and_truth() -> None:
    section = _enabled()
    section.update({
        "input_mode": "raw_fastq", "source_bam": None,
        "source_fastqs": ["reads.fastq.gz"],
        "profiles": [{"name": "baseline", "settings": {}}],
        "end_to_end": {
            "enabled": True, "truth_vcf": "truth.vcf.gz",
            "truth_version": "GIAB-4.2.1", "truth_source": "GIAB",
            "confident_bed": "confident.bed",
            "confident_bed_version": "GIAB-4.2.1",
            "metrics_compression": "gzip",
        },
    })
    validate_config({
        "reference": {"fasta": "ref.fa", "build": "GRCh38"},
        "low_depth": section,
    })
    del section["end_to_end"]["truth_version"]
    with pytest.raises(ConfigurationError, match="incomplete"):
        validate_config({
            "reference": {"fasta": "ref.fa", "build": "GRCh38"},
            "low_depth": section,
        })


def test_end_to_end_rejects_unvalidated_optimization_profiles() -> None:
    section = _enabled()
    section.update({
        "input_mode": "raw_fastq", "source_bam": None,
        "source_fastqs": ["reads.fastq.gz"],
        "end_to_end": {
            "enabled": True, "truth_vcf": "truth.vcf.gz",
            "truth_version": "GIAB-4.2.1", "truth_source": "GIAB",
            "confident_bed": "confident.bed",
            "confident_bed_version": "GIAB-4.2.1",
        },
    })
    with pytest.raises(ConfigurationError, match="exactly one profile"):
        validate_config({
            "reference": {"fasta": "ref.fa", "build": "GRCh38"},
            "low_depth": section,
        })
