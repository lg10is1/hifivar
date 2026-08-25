"""Policy tests for the Phase 2-4 external-tool runtime specification."""

from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[2]
MATRIX = ROOT / "workflow" / "envs" / "external-tools.yaml"
EXPECTED = {
    "pbmm2": "1.17.0",
    "deepvariant": "1.10.0",
    "sawfish": "2.2.1",
    "sniffles2": "2.8.0",
    "pbsv": "2.11.0",
    "cutesv": "2.1.4",
}


def _matrix() -> dict[str, object]:
    loaded = yaml.safe_load(MATRIX.read_text(encoding="utf-8"))
    assert isinstance(loaded, dict)
    return loaded


def test_external_tool_matrix_is_complete_and_exactly_versioned() -> None:
    document = _matrix()
    assert document["schema_version"] == 1
    tools = document["tools"]
    assert isinstance(tools, dict)
    assert set(tools) == set(EXPECTED)

    serialized = MATRIX.read_text(encoding="utf-8").lower()
    assert "latest" not in serialized
    for name, version in EXPECTED.items():
        spec = tools[name]
        assert spec["target_version"] == version
        assert spec["version_constraint"] == f"=={version}"
        for field in (
            "installation_source",
            "execution_backend",
            "executable_or_container",
            "wrapper_module",
            "version_command",
            "wrapper_compatibility_assumption",
            "config_section",
            "mock_verification_status",
            "official_cli_verification_status",
            "linux_real_verification_status",
        ):
            assert spec[field]


def test_native_conda_specs_pin_the_selected_tool() -> None:
    tools = _matrix()["tools"]
    for name in ("pbmm2", "sawfish", "sniffles2", "pbsv", "cutesv"):
        spec = tools[name]
        assert spec["execution_backend"] == "native-conda"
        environment = ROOT / spec["environment_spec"]
        loaded = yaml.safe_load(environment.read_text(encoding="utf-8"))
        dependencies = loaded["dependencies"]
        package = "sniffles" if name == "sniffles2" else name
        assert f"{package}={EXPECTED[name]}" in dependencies
        assert loaded["channels"][-1] == "nodefaults"


def test_deepvariant_production_backend_is_pinned_apptainer() -> None:
    deepvariant = _matrix()["tools"]["deepvariant"]
    assert deepvariant["execution_backend"] == "apptainer"
    assert deepvariant["executable_or_container"] == (
        "docker://google/deepvariant:1.10.0"
    )
    assert deepvariant["optional_backends"]["docker"]["status"] == "OPTIONAL"
    assert deepvariant["optional_backends"]["native"]["status"] == (
        "OPTIONAL_UNPACKAGED"
    )
    assert deepvariant["immutable_digest"] == "VERSION_PENDING_LINUX_VERIFICATION"
