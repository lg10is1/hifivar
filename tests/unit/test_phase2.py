"""Unit tests for Phase 2 settings and orchestration boundary checks."""

from __future__ import annotations

from pathlib import Path

import pytest

from hifivar.alignment import (
    AlignmentOutputFormat,
    AlignmentResources,
    AlignmentTool,
)
from hifivar.alignment_postprocess import AlignmentIndexFormat
from hifivar.exceptions import ConfigurationError, InputValidationError
from hifivar.phase2 import Phase2Settings, run_phase2


def alignment_config(**overrides: object) -> dict[str, object]:
    section: dict[str, object] = {
        "tool": "pbmm2",
        "output_format": "bam",
        "threads": 16,
        "memory_mb": 64000,
        "runtime_minutes": 720,
        "overwrite": False,
        "index_threads": 4,
        "bam_index_format": "auto",
    }
    section.update(overrides)
    return {"alignment": section}


def test_phase2_settings_load_complete_alignment_config() -> None:
    settings = Phase2Settings.from_config(alignment_config())
    assert settings.tool is AlignmentTool.PBMM2
    assert settings.output_format is AlignmentOutputFormat.BAM
    assert settings.resources == AlignmentResources(16, 64000, 720)
    assert settings.index_threads == 4
    assert settings.bam_index_format is None
    assert settings.to_dict()["bam_index_format"] == "auto"


@pytest.mark.parametrize(
    "configured,expected",
    (("bai", AlignmentIndexFormat.BAI), ("csi", AlignmentIndexFormat.CSI)),
)
def test_phase2_settings_accept_explicit_bam_index_format(
    configured: str,
    expected: AlignmentIndexFormat,
) -> None:
    settings = Phase2Settings.from_config(
        alignment_config(bam_index_format=configured)
    )
    assert settings.bam_index_format is expected


def test_phase2_execution_rejects_unimplemented_minimap2_and_cram() -> None:
    with pytest.raises(ConfigurationError, match="pbmm2 only"):
        Phase2Settings(
            AlignmentTool.MINIMAP2,
            AlignmentOutputFormat.BAM,
            AlignmentResources(),
            False,
            1,
        )
    with pytest.raises(ConfigurationError, match="output_format=bam"):
        Phase2Settings(
            AlignmentTool.PBMM2,
            AlignmentOutputFormat.CRAM,
            AlignmentResources(),
            False,
            1,
        )


@pytest.mark.parametrize(
    "config",
    (
        {},
        {"alignment": {}},
        alignment_config(index_threads=0),
        alignment_config(bam_index_format="bad"),
    ),
)
def test_phase2_settings_reject_incomplete_or_invalid_config(
    config: dict[str, object],
) -> None:
    with pytest.raises((ConfigurationError, InputValidationError)):
        Phase2Settings.from_config(config)


def test_run_phase2_rejects_wrong_context(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="AnalysisContext"):
        run_phase2(object(), tmp_path)  # type: ignore[arg-type]
