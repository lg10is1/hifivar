"""Tests for the shared HiFiVar YAML configuration system."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest
import yaml

from hifivar.config import (
    HiFiVarConfig,
    deep_merge,
    load_config,
    load_yaml,
    merge_configs,
    validate_config,
    write_effective_config,
)
from hifivar.exceptions import ConfigurationError


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_RESOURCE_DIRECTORY = (
    PROJECT_ROOT / "src" / "hifivar" / "resources" / "configs"
)
DEFAULT_CONFIG = CONFIG_RESOURCE_DIRECTORY / "default.yaml"
PRESET_DIRECTORY = CONFIG_RESOURCE_DIRECTORY / "presets"
PRESET_NAMES = ("fast", "standard", "comprehensive", "cohort", "trio")


def complete_config() -> dict[str, object]:
    """Return an independent, schema-complete configuration mapping."""
    return {
        "project": {"name": "hifivar"},
        "reference": {"fasta": None, "build": None},
        "samples": {"sheet": None},
        "runtime": {"threads": 1, "tmpdir": None},
        "paths": {"workdir": None, "outdir": None},
        "logging": {"level": "INFO", "file": None},
        "workflow": {"preset": "standard"},
        "alignment": {
            "tool": "pbmm2",
            "output_format": "bam",
            "threads": 8,
            "memory_mb": 32000,
            "runtime_minutes": 1440,
            "overwrite": False,
            "pbmm2_preset": "CCS",
            "pbmm2_log_level": "INFO",
            "index_threads": 4,
            "bam_index_format": "auto",
        },
    }


def write_yaml(path: Path, data: object) -> None:
    """Write a UTF-8 YAML fixture."""
    path.write_text(
        yaml.safe_dump(data, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )


def test_load_yaml_reads_basic_mapping(tmp_path: Path) -> None:
    """A valid YAML mapping should load without conversion surprises."""
    config_path = tmp_path / "config.yaml"
    config_path.write_text("runtime:\n  threads: 8\n", encoding="utf-8")

    loaded = load_yaml(config_path)

    assert loaded == {"runtime": {"threads": 8}}


def test_load_yaml_treats_empty_document_as_empty_mapping(tmp_path: Path) -> None:
    """The behavior of an empty YAML document should be explicit."""
    config_path = tmp_path / "empty.yaml"
    config_path.write_text("", encoding="utf-8")

    assert load_yaml(config_path) == {}


def test_deep_merge_preserves_unmodified_nested_values() -> None:
    """Overriding one runtime value should not discard its siblings."""
    default = {"runtime": {"threads": 1, "tmpdir": None}}
    override = {"runtime": {"threads": 8}}

    merged = deep_merge(default, override)

    assert merged == {"runtime": {"threads": 8, "tmpdir": None}}


def test_load_config_precedence_is_user_then_preset_then_default(
    tmp_path: Path,
) -> None:
    """Later layers should deterministically override earlier layers."""
    default_path = tmp_path / "default.yaml"
    preset_path = tmp_path / "preset.yaml"
    user_path = tmp_path / "user.yaml"
    write_yaml(default_path, complete_config())
    write_yaml(
        preset_path,
        {"runtime": {"threads": 8}, "workflow": {"preset": "fast"}},
    )
    write_yaml(user_path, {"runtime": {"threads": 32}})

    config = load_config(default_path, preset_path, user_path)

    assert config["runtime"]["threads"] == 32  # type: ignore[index]
    assert config["runtime"]["tmpdir"] is None  # type: ignore[index]
    assert config["workflow"]["preset"] == "fast"  # type: ignore[index]


def test_lists_are_replaced_instead_of_appended() -> None:
    """Lists are atomic override values rather than mergeable collections."""
    merged = deep_merge({"items": ["a", "b"]}, {"items": ["c"]})

    assert merged["items"] == ["c"]


def test_merge_configs_applies_layers_left_to_right() -> None:
    """The generic merge helper should use its documented priority order."""
    merged = merge_configs(
        {"runtime": {"threads": 1, "tmpdir": None}},
        {"runtime": {"threads": 8}},
        {"runtime": {"threads": 64}},
    )

    assert merged == {"runtime": {"threads": 64, "tmpdir": None}}


def test_missing_file_raises_configuration_error(tmp_path: Path) -> None:
    """Filesystem details should be wrapped in a HiFiVar exception."""
    missing_path = tmp_path / "missing.yaml"

    with pytest.raises(ConfigurationError, match="missing.yaml"):
        load_yaml(missing_path)


def test_yaml_syntax_error_raises_configuration_error(tmp_path: Path) -> None:
    """PyYAML parser failures should identify the source file."""
    config_path = tmp_path / "invalid.yaml"
    config_path.write_text("runtime:\n  threads: [\n", encoding="utf-8")

    with pytest.raises(ConfigurationError, match="invalid.yaml"):
        load_yaml(config_path)


def test_non_mapping_yaml_root_is_rejected(tmp_path: Path) -> None:
    """A list cannot serve as the root configuration value."""
    config_path = tmp_path / "list.yaml"
    config_path.write_text("- a\n- b\n", encoding="utf-8")

    with pytest.raises(ConfigurationError, match="YAML root"):
        load_yaml(config_path)


def test_unknown_top_level_key_is_rejected() -> None:
    """Top-level typos must not be silently ignored."""
    with pytest.raises(ConfigurationError, match="unknown_section"):
        validate_config({"unknown_section": {"abc": 1}})


def test_unknown_nested_key_is_rejected() -> None:
    """Nested typos must identify their fully qualified key."""
    with pytest.raises(ConfigurationError, match=r"runtime\.threeds"):
        validate_config({"runtime": {"threeds": 32}})


def test_positive_integer_threads_is_valid() -> None:
    """A positive integer thread count should pass validation."""
    validate_config({"runtime": {"threads": 64}})


@pytest.mark.parametrize(
    "reference",
    (
        {"fasta": None, "build": None},
        {"fasta": "/references/genome.fa", "build": "GRCh38"},
        {"build": "my_species_v1"},
    ),
)
def test_reference_config_accepts_optional_string_metadata(
    reference: dict[str, object],
) -> None:
    """Reference config remains optional until a model is constructed."""
    validate_config({"reference": reference})


@pytest.mark.parametrize(
    "reference,keyword",
    (
        ({"fasta": 42}, "reference.fasta"),
        ({"build": ""}, "reference.build"),
        ({"build": 38}, "reference.build"),
    ),
)
def test_reference_config_rejects_invalid_metadata_types(
    reference: dict[str, object],
    keyword: str,
) -> None:
    """Reference paths/build labels use explicit nullable string types."""
    with pytest.raises(ConfigurationError, match=keyword):
        validate_config({"reference": reference})


@pytest.mark.parametrize("sheet", (None, "/project/samples.tsv"))
def test_samples_sheet_accepts_nullable_string_path(sheet: str | None) -> None:
    """The Phase 1.4 sample-sheet entry follows existing config path types."""
    validate_config({"samples": {"sheet": sheet}})


def test_samples_sheet_unknown_nested_key_is_rejected() -> None:
    """Unknown sample configuration cannot be silently ignored."""
    with pytest.raises(ConfigurationError, match=r"samples\.spreadsheet"):
        validate_config({"samples": {"spreadsheet": "samples.tsv"}})


@pytest.mark.parametrize("sheet", (42, True, ["samples.tsv"]))
def test_samples_sheet_rejects_wrong_type(sheet: object) -> None:
    """Sample sheet settings remain nullable scalar path strings."""
    with pytest.raises(ConfigurationError, match=r"samples\.sheet"):
        validate_config({"samples": {"sheet": sheet}})


@pytest.mark.parametrize("threads", (0, -5, "64", True))
def test_invalid_threads_are_rejected(threads: object) -> None:
    """Zero, negatives, strings, and booleans are not thread counts."""
    with pytest.raises(ConfigurationError, match="positive integer"):
        validate_config({"runtime": {"threads": threads}})


def test_alignment_config_accepts_phase2_settings() -> None:
    validate_config(
        {
            "alignment": {
                "tool": "pbmm2",
                "output_format": "bam",
                "threads": 24,
                "memory_mb": 64000,
                "runtime_minutes": 720,
                "overwrite": False,
                "pbmm2_preset": "HIFI",
                "pbmm2_log_level": "DEBUG",
                "index_threads": 4,
                "bam_index_format": "csi",
            }
        }
    )


@pytest.mark.parametrize(
    "key,value",
    (
        ("tool", "bwa"),
        ("output_format", "sam"),
        ("threads", 0),
        ("memory_mb", True),
        ("runtime_minutes", -1),
        ("overwrite", "false"),
        ("pbmm2_preset", "SUBREAD"),
        ("pbmm2_log_level", "VERBOSE"),
        ("index_threads", 0),
        ("bam_index_format", "crai"),
    ),
)
def test_alignment_config_rejects_invalid_settings(
    key: str,
    value: object,
) -> None:
    with pytest.raises(ConfigurationError, match=f"alignment.{key}"):
        validate_config({"alignment": {key: value}})


@pytest.mark.parametrize(
    "level",
    ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL", "debug"),
)
def test_supported_logging_levels_are_valid(level: str) -> None:
    """Config validation should reuse logging's accepted level parser."""
    validate_config({"logging": {"level": level}})


def test_invalid_logging_level_is_rejected() -> None:
    """Logging typos should surface as configuration errors."""
    with pytest.raises(ConfigurationError, match="INVALID_LEVEL"):
        validate_config({"logging": {"level": "INVALID_LEVEL"}})


def test_utf8_yaml_is_loaded_without_data_loss(tmp_path: Path) -> None:
    """Project names and future metadata may contain Unicode text."""
    config_path = tmp_path / "unicode.yaml"
    config_path.write_text(
        'project:\n  name: "HiFiVar 测试"\n',
        encoding="utf-8",
    )

    loaded = load_yaml(config_path)

    assert loaded["project"]["name"] == "HiFiVar 测试"  # type: ignore[index]


def test_effective_config_can_be_written_and_reloaded(tmp_path: Path) -> None:
    """Effective YAML should preserve validated configuration values."""
    config = load_config(
        DEFAULT_CONFIG,
        PRESET_DIRECTORY / "standard.yaml",
    )
    output_path = tmp_path / "nested" / "effective.yaml"

    write_effective_config(config, output_path)
    reloaded = load_yaml(output_path)

    assert output_path.is_file()
    assert reloaded == config.to_dict()
    validate_config(reloaded, require_complete=True)


def test_deep_merge_does_not_mutate_inputs() -> None:
    """Merge operations should not alter either caller-owned mapping."""
    base = {"runtime": {"threads": 1}, "items": ["a"]}
    override = {"runtime": {"threads": 8}, "items": ["b"]}
    original_base = deepcopy(base)
    original_override = deepcopy(override)

    deep_merge(base, override)

    assert base == original_base
    assert override == original_override


def test_repeated_loads_do_not_share_nested_state() -> None:
    """Mutating one result must not contaminate a subsequent load."""
    first = load_config(DEFAULT_CONFIG)
    first_runtime = first["runtime"]
    assert isinstance(first_runtime, dict)
    first_runtime["threads"] = 100

    second = load_config(DEFAULT_CONFIG)

    assert second["runtime"]["threads"] == 1  # type: ignore[index]


@pytest.mark.parametrize("preset_name", PRESET_NAMES)
def test_repository_presets_parse_and_validate(preset_name: str) -> None:
    """Every shipped preset should be a valid minimal config fragment."""
    preset = load_yaml(PRESET_DIRECTORY / f"{preset_name}.yaml")

    validate_config(preset)

    assert preset["workflow"]["preset"] == preset_name  # type: ignore[index]


def test_config_records_source_provenance(tmp_path: Path) -> None:
    """Loaded configs should retain simple default/preset/user provenance."""
    user_path = tmp_path / "user.yaml"
    write_yaml(user_path, {"runtime": {"threads": 2}})

    config = load_config(
        DEFAULT_CONFIG,
        PRESET_DIRECTORY / "fast.yaml",
        user_path,
    )

    assert config.sources == {
        "default": DEFAULT_CONFIG,
        "preset": PRESET_DIRECTORY / "fast.yaml",
        "user": user_path,
    }


def test_user_paths_expand_tilde_without_existence_check(tmp_path: Path) -> None:
    """A leading tilde should expand without resolving or requiring a path."""
    default_path = tmp_path / "default.yaml"
    data = complete_config()
    data["paths"]["outdir"] = "~/future-output"  # type: ignore[index]
    write_yaml(default_path, data)

    config = load_config(default_path)

    assert config["paths"]["outdir"] == str(  # type: ignore[index]
        Path("~/future-output").expanduser()
    )


def test_loading_does_not_modify_user_yaml(tmp_path: Path) -> None:
    """User input is read-only unless effective output is explicitly requested."""
    user_path = tmp_path / "user.yaml"
    original_text = "runtime:\n  threads: 7\n"
    user_path.write_text(original_text, encoding="utf-8")

    load_config(DEFAULT_CONFIG, user_config=user_path)

    assert user_path.read_text(encoding="utf-8") == original_text


def test_incomplete_effective_config_is_rejected(tmp_path: Path) -> None:
    """The final merged config must contain all current foundation fields."""
    default_path = tmp_path / "incomplete.yaml"
    write_yaml(default_path, {"runtime": {"threads": 1}})

    with pytest.raises(ConfigurationError, match="Missing required"):
        load_config(default_path)


def test_config_repr_does_not_dump_configuration_values() -> None:
    """Default representations should not expose future secret-like values."""
    config = HiFiVarConfig(
        {"project": {"name": "private-value"}},
        {"default": Path("default.yaml")},
    )

    assert "private-value" not in repr(config)


@pytest.mark.parametrize(
    "small,pattern",
    (
        ({"execution_mode": "docker", "deepvariant_image": None}, "required"),
        ({"execution_mode": "apptainer"}, "required"),
        (
            {"execution_mode": "native", "deepvariant_image": "unexpected.sif"},
            "must be null",
        ),
        ({"model_type": "WGS"}, "model_type"),
        ({"threads": 0}, "threads"),
        ({"max_concurrent_samples": 0}, "max_concurrent_samples"),
        ({"overwrite": "yes"}, "overwrite"),
    ),
)
def test_small_variant_config_rejects_invalid_values(
    small: dict[str, object],
    pattern: str,
) -> None:
    with pytest.raises(ConfigurationError, match=pattern):
        validate_config({"small": small})


def test_small_variant_container_config_is_valid() -> None:
    validate_config(
        {
            "small": {
                "enabled": True,
                "execution_mode": "docker",
                "deepvariant_executable": "run_deepvariant",
                "deepvariant_image": "google/deepvariant:1.10.0",
                "model_type": "PACBIO",
                "threads": 16,
                "max_concurrent_samples": 1,
                "memory_mb": 64000,
                "runtime_minutes": 2880,
                "overwrite": False,
            }
        }
    )


def test_phase4_sv_config_accepts_strict_nested_callers() -> None:
    validate_config(
        {
            "sv": {
                "enabled": True,
                "overwrite": False,
                "sawfish": {"enabled": True, "executable": "sawfish", "threads": 16, "disable_cnv": False},
                "sniffles2": {"enabled": True, "executable": "sniffles", "minimum_support": None, "minimum_sv_length": 50},
                "pbsv": {"enabled": True, "executable": "pbsv", "threads": 8},
                "cutesv": {"enabled": True, "executable": "cuteSV", "diff_ratio_merging_ins": 0.9, "genotype": True},
                "finalization": {"bgzip_executable": "bgzip", "tabix_executable": "tabix"},
            }
        }
    )


@pytest.mark.parametrize(
    "sv,pattern",
    (
        ({"enabled": "yes"}, r"sv\.enabled"),
        ({"sawfish": {"threeds": 8}}, r"sv\.sawfish\.threeds"),
        ({"sniffles2": {"minimum_sv_length": 0}}, r"minimum_sv_length"),
        ({"cutesv": {"diff_ratio_merging_ins": 1.5}}, r"diff_ratio_merging_ins"),
        ({"cutesv": {"genotype": "yes"}}, r"genotype"),
        ({"finalization": {"bgzip_executable": ""}}, r"bgzip_executable"),
    ),
)
def test_phase4_sv_config_rejects_invalid_or_unknown_values(sv, pattern) -> None:
    with pytest.raises(ConfigurationError, match=pattern):
        validate_config({"sv": sv})


def test_phase5_tr_config_accepts_strict_settings() -> None:
    validate_config({
        "tr": {
            "enabled": True,
            "catalog": "/reference/trgt.bed",
            "catalog_reference_build": "GRCh38",
            "executable": "trgt",
            "bcftools_executable": "bcftools",
            "samtools_executable": "samtools",
            "threads": 8,
            "memory_mb": 16000,
            "runtime_minutes": 720,
            "preset": "wgs",
            "karyotype": "auto",
            "overwrite": False,
        }
    })


@pytest.mark.parametrize(
    "tr,pattern",
    [
        ({"enabled": "yes"}, r"tr\.enabled"),
        ({"catalog": 42}, r"tr\.catalog"),
        ({"preset": "unknown"}, r"tr\.preset"),
        ({"karyotype": "unknown"}, r"tr\.karyotype"),
        ({"threeds": 8}, r"tr\.threeds"),
    ],
)
def test_phase5_tr_config_rejects_invalid_or_unknown_values(tr, pattern) -> None:
    with pytest.raises(ConfigurationError, match=pattern):
        validate_config({"tr": tr})
