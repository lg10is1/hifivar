"""Static and package-level checks for public distribution definitions."""

from __future__ import annotations

import tomllib
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[2]


def _text(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def test_public_license_and_python_metadata_are_consistent() -> None:
    license_text = _text("LICENSE")
    metadata = tomllib.loads(_text("pyproject.toml"))["project"]

    assert "Apache License" in license_text
    assert "Version 2.0, January 2004" in license_text
    assert metadata["license"] == "Apache-2.0"
    assert metadata["license-files"] == ["LICENSE"]
    assert "snakemake>=8,<10" in metadata["optional-dependencies"]["workflow"]


def test_conda_core_definitions_are_bounded_and_non_monolithic() -> None:
    environment = yaml.safe_load(_text("environment.yml"))
    recipe = _text("conda-recipe/meta.yaml")
    installation = _text("docs/installation.md")
    dependencies = [str(item) for item in environment["dependencies"] if isinstance(item, str)]

    assert environment["name"] == "hifivar"
    assert "python=3.12" in dependencies
    assert "snakemake-minimal>=8,<10" in dependencies
    assert "version = \"0.1.0rc1\"" in recipe
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
    quickstart = _text("docs/quickstart.md")
    containers = _text("docs/containers.md")

    assert "does not yet expose one unified `hifivar run`" in readme
    assert "consume an existing indexed" in readme
    assert "BAM/CRAM" in readme
    assert "FASTQ boundary" in quickstart
    assert "`/data/...` and" in quickstart
    assert "`/work/...` path" in quickstart
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
