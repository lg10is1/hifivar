import pytest

from hifivar.config import validate_config
from hifivar.exceptions import ConfigurationError


def test_cohort_config_accepts_independent_tracks() -> None:
    validate_config({"cohort": {
        "enabled": True,
        "cohort_id": "C1",
        "input_manifest": "cohort.tsv",
        "small_variants": {"enabled": True, "glnexus_executable": "glnexus_cli", "bcftools_executable": "bcftools", "preset": "DeepVariantWGS", "threads": 8, "memory_gb": 32, "runtime_minutes": 100},
        "sv": {"enabled": False, "threads": 1, "memory_mb": 1000, "runtime_minutes": 10},
        "tr": {"enabled": True, "threads": 1, "memory_mb": 1000, "runtime_minutes": 10},
    }})


@pytest.mark.parametrize("config, message", [
    ({"cohort": {"enabled": True, "cohort_id": "C1", "input_manifest": None}}, "requires"),
    ({"cohort": {"enabled": False, "small_variants": {"threads": 0}}}, "positive"),
    ({"cohort": {"enabled": True, "cohort_id": "C1", "input_manifest": "cohort.tsv", "small_variants": {"enabled": True}}}, "explicit positive"),
    ({"cohort": {"enabled": False, "sv": {"enabled": "yes"}}}, "boolean"),
    ({"cohort": {"enabled": False, "tr": {"mystery": 1}}}, "Unknown"),
])
def test_cohort_config_rejects_incomplete_unknown_or_wrong_types(config, message) -> None:
    with pytest.raises(ConfigurationError, match=message):
        validate_config(config)


def test_disabled_small_cohort_allows_unset_memory() -> None:
    validate_config(
        {
            "cohort": {
                "enabled": False,
                "small_variants": {"enabled": False, "memory_gb": None},
            }
        }
    )
