"""Static and package-level checks for public distribution definitions."""

from __future__ import annotations

import importlib.util
import tomllib
from pathlib import Path

import yaml
from hifivar import __version__


ROOT = Path(__file__).resolve().parents[2]


def _text(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def _load_public_audit_module():
    script = ROOT / "scripts" / "audit_public_release.py"
    spec = importlib.util.spec_from_file_location("audit_public_release", script)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_public_license_and_python_metadata_are_consistent() -> None:
    license_text = _text("LICENSE")
    metadata = tomllib.loads(_text("pyproject.toml"))["project"]
    citation = yaml.safe_load(_text("CITATION.cff"))

    assert "Apache License" in license_text
    assert "Version 2.0, January 2004" in license_text
    assert metadata["license"] == "Apache-2.0"
    assert metadata["license-files"] == ["LICENSE"]
    assert "snakemake>=8,<10" in metadata["optional-dependencies"]["workflow"]
    assert citation["version"] == __version__


def test_conda_core_definitions_are_bounded_and_non_monolithic() -> None:
    environment = yaml.safe_load(_text("environment.yml"))
    recipe = _text("conda-recipe/meta.yaml")
    installation = _text("docs/installation.md")
    dependencies = [str(item) for item in environment["dependencies"] if isinstance(item, str)]

    assert environment["name"] == "hifivar"
    assert "python=3.12" in dependencies
    assert "snakemake-minimal>=8,<10" in dependencies
    assert f'version = "{__version__}"' in recipe
    assert "noarch: python" in recipe
    assert "license: Apache-2.0" in recipe
    assert "deepvariant" not in recipe.lower()
    assert "annovar" not in recipe.lower()
    assert "latest" not in recipe.lower()
    build_command = "conda build -c conda-forge -c bioconda conda-recipe"
    assert build_command in recipe
    assert build_command in installation


def test_docker_core_is_non_root_and_does_not_bundle_callers() -> None:
    dockerfile = _text("Dockerfile")
    dockerignore = _text(".dockerignore")

    assert "FROM python:3.12-slim-bookworm" in dockerfile
    assert "USER 10001:10001" in dockerfile
    assert 'ENTRYPOINT ["hifivar"]' in dockerfile
    assert "hifivar --version" in dockerfile
    assert "latest" not in dockerfile.lower()
    assert "00.hifivar" not in dockerfile
    assert "*.sif" in dockerignore


def test_apptainer_core_has_build_and_test_contract() -> None:
    definition = _text("containers/hifivar.def")

    assert definition.startswith("Bootstrap: docker\nFrom: python:3.12-slim-bookworm")
    assert "%files" in definition
    assert "%post" in definition
    assert "%runscript" in definition
    assert "%test" in definition
    assert "hifivar config validate" in definition
    assert "latest" not in definition.lower()


def test_minimal_example_is_safe_and_caller_disabled() -> None:
    config = yaml.safe_load(_text("examples/minimal/config.yaml"))
    sample_lines = _text("examples/minimal/samples.tsv").splitlines()

    assert config["reference"]["build"] == "GRCh38"
    assert config["samples"]["sheet"] == "/work/samples.tsv"
    for section in (
        "small",
        "sv",
        "tr",
        "phasing",
        "assembly",
        "assembly_sv",
        "review",
        "annotation",
        "cohort",
        "benchmark",
    ):
        assert config[section]["enabled"] is False
    assert sample_lines[0] == "sample_id\tinput\tinput_type\tsex"
    assert sample_lines[1].endswith("\tbam\tmale")


def test_public_docs_disclose_execution_and_container_boundaries() -> None:
    readme = _text("README.md")
    readme_zh = _text("README.zh-CN.md")
    quickstart = _text("docs/quickstart.md")
    quickstart_zh = _text("docs/zh_CN/quickstart.md")
    installation_zh = _text("docs/zh_CN/installation.md")
    containers = _text("docs/containers.md")

    assert "0.1.0rc4" in readme
    assert "0.1.0rc4" in readme_zh
    assert __version__ in readme
    assert __version__ in readme_zh
    assert "[简体中文](README.zh-CN.md)" in readme
    # The tagged installation guide must agree with the documented launcher.
    assert "[`hifivar run`](docs/indexed_bam_run.md)" in readme
    assert "docs/raw_input_run.md" in readme
    assert "releases/tag/v0.1.0-rc5" in readme
    assert "15 October payload" in _text("docs/rc5_candidate.md")
    candidate_guide = _text("docs/rc5_candidate.md")
    assert "different manifests and validation scopes" in candidate_guide
    assert "CRAM remains deferred at the main launcher" in readme
    assert "docs/capabilities.md" in readme and "docs/capabilities.md" in readme_zh
    assert "README.zh-CN.md" in _text("MANIFEST.in")
    capabilities = _text("docs/capabilities.md")
    assert "GENERATED" in capabilities and "NOT_REVIEWED" in capabilities
    assert "separate" in capabilities and "report" in capabilities
    assert "docs/raw_input_run.md" in readme_zh
    raw_guide = _text("docs/raw_input_run.md")
    assert "RAW_INPUT_PLAN_ONLY" in raw_guide
    assert "not a single raw-input" in raw_guide
    assert "no FASTQ conversion" in raw_guide
    assert "consume an existing indexed" in readme
    assert "BAM/CRAM" in readme
    assert "FASTQ boundary" in quickstart
    assert "`/data/...` and" in quickstart
    assert "`/work/...` path" in quickstart
    assert "RAW_INPUT_PLAN_ONLY" in quickstart_zh
    assert "run --run-dir ./runs/sample01 --resume" in quickstart_zh
    assert "0.1.0rc4" in installation_zh
    assert __version__ in installation_zh
    assert "does not contain PAV, DeepVariant" in containers


def test_ci_separates_python310_core_from_workflow_and_forces_utf8() -> None:
    ci = _text(".github/workflows/ci.yml")

    assert 'PYTHONUTF8: "1"' in ci
    assert "actions/checkout@v7" in ci
    assert "actions/setup-python@v7" in ci
    assert "if: matrix.python == '3.10'" in ci
    assert 'python -m pip install ".[dev]"' in ci
    assert "python -m pytest tests/unit -p no:cacheprovider" in ci
    assert "if: matrix.python != '3.10'" in ci
    assert 'python -m pip install ".[dev,workflow]"' in ci


def test_distribution_ci_uses_current_node24_actions() -> None:
    distribution = _text(".github/workflows/distribution.yml")

    assert distribution.count("actions/checkout@v7") == 2
    assert "actions/setup-python@v7" in distribution
    assert "actions/checkout@v4" not in distribution
    assert "actions/setup-python@v5" not in distribution


def test_public_audit_ignores_named_local_virtual_environments(tmp_path: Path) -> None:
    audit_module = _load_public_audit_module()
    local_environment = tmp_path / ".venv-validation" / "Lib" / "site-packages"
    local_environment.mkdir(parents=True)
    private_cache_text = "maintainer" + "@" + "example.org\n" + "C:" + "\\Users\\My\\local-only-cache\n"
    (local_environment / "third_party_metadata.txt").write_text(
        private_cache_text,
        encoding="utf-8",
    )
    (tmp_path / "README.md").write_text("Public project documentation.\n", encoding="utf-8")

    assert audit_module.audit(tmp_path, max_bytes=1024 * 1024) == []


def test_publication_visuals_follow_shared_nature_profile() -> None:
    rules = _text("PROJECT_VISUAL_RULES.md")
    style = yaml.safe_load(_text("project_visuals/visual_style.yaml"))
    plotting = _text("project_visuals/scripts/plot_hg002_benchmarks.py")

    assert style["profile"] == "nature_family"
    assert style["publication_geometry"]["single_column_width_mm"] == 89
    assert style["publication_geometry"]["double_column_width_mm"] == 183
    assert style["typography"]["minimum_pt"] == 5
    assert style["typography"]["maximum_pt"] == 7
    assert style["lines"]["minimum_pt"] >= 0.5
    assert style["accessibility"]["forbid_red_green_only_contrast"] is True
    assert "invented_error_bars" in style["prohibited"]
    assert "truncated_bar_axis" in style["prohibited"]
    assert "Nature-family production target" in rules
    assert "#0072B2" in plotting
    assert "#009E73" in plotting
    assert "'pdf.fonttype': 42" in plotting
    assert ".grid(" not in plotting
