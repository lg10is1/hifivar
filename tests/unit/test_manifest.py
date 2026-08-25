"""Tests for Phase 1.4 portable run provenance manifests."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import pytest
import yaml

from hifivar import __version__
from hifivar.context import AnalysisContext
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.manifest import MANIFEST_SCHEMA_VERSION, RunManifest
from hifivar.reference import ReferenceGenome
from hifivar.sample import InputDataset, InputType, Sample
from hifivar.sample_sheet import SampleRecord, SampleSheet


def write_reference(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    fasta = root / "reference.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    return fasta


def write_fastq(root: Path, name: str = "sample.fastq") -> Path:
    root.mkdir(parents=True, exist_ok=True)
    fastq = root / name
    fastq.write_text("@read1\nACGT\n+\nIIII\n", encoding="utf-8")
    return fastq


def make_context(
    tmp_path: Path,
    *,
    config: dict[str, object] | None = None,
    with_sheet: bool = False,
) -> AnalysisContext:
    reference = ReferenceGenome.from_fasta(write_reference(tmp_path), build="GRCh38")
    sample = Sample("S1", InputDataset.from_files((write_fastq(tmp_path),)))
    effective = config or {
        "reference": {"fasta": str(reference.fasta.absolute()), "build": "GRCh38"},
        "samples": {"sheet": None},
        "runtime": {"threads": 2},
    }
    if not with_sheet:
        return AnalysisContext.from_sample(reference, sample, effective)
    sheet_path = tmp_path / "samples.tsv"
    sheet_path.write_text("sample_id\tinput\nS1\tsample.fastq\n", encoding="utf-8")
    sheet = SampleSheet(sheet_path, (SampleRecord(sample),))
    return AnalysisContext.from_sample_sheet(reference, sheet, effective)


def assert_standard_types(value: object) -> None:
    if isinstance(value, dict):
        assert all(isinstance(key, str) for key in value)
        for item in value.values():
            assert_standard_types(item)
    elif isinstance(value, list):
        for item in value:
            assert_standard_types(item)
    else:
        assert value is None or isinstance(value, (str, int, float, bool))


def test_manifest_has_schema_version(tmp_path: Path) -> None:
    manifest = RunManifest.from_context(make_context(tmp_path))
    assert manifest.schema_version == MANIFEST_SCHEMA_VERSION == "1.0"


def test_manifest_records_hifivar_version(tmp_path: Path) -> None:
    assert RunManifest.from_context(make_context(tmp_path)).hifivar_version == __version__


def test_manifest_timestamp_is_utc_iso8601(tmp_path: Path) -> None:
    timestamp = RunManifest.from_context(make_context(tmp_path)).created_at
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}T.+Z", timestamp)


def test_reference_paths_are_absolute_without_mutating_context(tmp_path: Path) -> None:
    context = make_context(tmp_path)
    original = context.reference.fasta
    manifest = RunManifest.from_context(context)
    assert Path(manifest.reference["fasta"]).is_absolute()
    assert Path(manifest.reference["fai"]).is_absolute()
    assert context.reference.fasta == original


def test_sample_input_paths_are_absolute(tmp_path: Path) -> None:
    manifest = RunManifest.from_context(make_context(tmp_path))
    input_metadata = manifest.samples[0]["input"]
    assert isinstance(input_metadata, dict)
    assert Path(input_metadata["files"][0]).is_absolute()


def test_inputs_have_sample_type_path_and_size(tmp_path: Path) -> None:
    context = make_context(tmp_path)
    manifest = RunManifest.from_context(context)
    input_record = manifest.inputs[0]
    assert input_record["sample_id"] == "S1"
    assert input_record["input_type"] == "fastq"
    assert Path(input_record["path"]).is_absolute()
    assert input_record["size_bytes"] == context.samples[0].sample.input.files[0].stat().st_size


def test_input_checksum_is_disabled_by_default(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def unexpected_checksum(path: object) -> str:
        raise AssertionError(f"unexpected checksum for {path}")

    monkeypatch.setattr("hifivar.manifest.validation.compute_sha256", unexpected_checksum)
    manifest = RunManifest.from_context(make_context(tmp_path))
    assert manifest.inputs[0]["sha256"] is None


def test_input_checksum_can_be_enabled_explicitly(tmp_path: Path) -> None:
    context = make_context(tmp_path)
    manifest = RunManifest.from_context(context, compute_input_checksums=True)
    fastq = context.samples[0].sample.input.files[0]
    assert manifest.inputs[0]["sha256"] == hashlib.sha256(fastq.read_bytes()).hexdigest()


def test_reference_checksum_is_not_forced(tmp_path: Path) -> None:
    manifest = RunManifest.from_context(make_context(tmp_path))
    assert manifest.reference["sha256"] is None


def test_effective_config_is_recorded(tmp_path: Path) -> None:
    manifest = RunManifest.from_context(make_context(tmp_path))
    assert manifest.config["runtime"] == {"threads": 2}


def test_sensitive_config_keys_are_redacted_recursively(tmp_path: Path) -> None:
    config = {
        "credentials": {
            "token": "a",
            "PASSWORD": "b",
            "nested": [{"Api_Key": "c"}, {"secret": "d"}],
        }
    }
    manifest = RunManifest.from_context(make_context(tmp_path, config=config))
    credentials = manifest.config["credentials"]
    assert credentials == {
        "token": "***",
        "PASSWORD": "***",
        "nested": [{"Api_Key": "***"}, {"secret": "***"}],
    }


def test_non_secret_substring_key_is_not_redacted(tmp_path: Path) -> None:
    config = {"metadata": {"secretary": "Alice", "token_count": 3}}
    manifest = RunManifest.from_context(make_context(tmp_path, config=config))
    assert manifest.config == config


def test_composite_sensitive_keys_are_redacted_consistently(tmp_path: Path) -> None:
    config = {
        "credentials": {
            "access_token": "a",
            "client-secret": "b",
            "service_api_key": "c",
            "token_count": 3,
            "secretary": "Alice",
        }
    }
    manifest = RunManifest.from_context(make_context(tmp_path, config=config))
    assert manifest.config["credentials"] == {
        "access_token": "***",
        "client-secret": "***",
        "service_api_key": "***",
        "token_count": 3,
        "secretary": "Alice",
    }


def test_manifest_payload_contains_only_standard_types(tmp_path: Path) -> None:
    payload = RunManifest.from_context(make_context(tmp_path)).to_dict()
    assert_standard_types(payload)


def test_path_and_enum_values_are_standardized(tmp_path: Path) -> None:
    config = {"future": {"path": Path("relative/path"), "type": InputType.FASTQ}}
    manifest = RunManifest.from_context(make_context(tmp_path, config=config))
    assert manifest.config == {
        "future": {"path": str(Path("relative/path")), "type": "fastq"}
    }


def test_to_dict_returns_independent_data(tmp_path: Path) -> None:
    manifest = RunManifest.from_context(make_context(tmp_path))
    payload = manifest.to_dict()
    payload["reference"]["build"] = "changed"  # type: ignore[index]
    assert manifest.reference["build"] == "GRCh38"


def test_json_serialization_round_trip(tmp_path: Path) -> None:
    manifest = RunManifest.from_context(make_context(tmp_path))
    path = manifest.write_json(tmp_path / "provenance" / "run.json")
    assert json.loads(path.read_text(encoding="utf-8")) == manifest.to_dict()


def test_yaml_serialization_round_trip(tmp_path: Path) -> None:
    manifest = RunManifest.from_context(make_context(tmp_path))
    path = manifest.write_yaml(tmp_path / "provenance" / "run.yaml")
    assert yaml.safe_load(path.read_text(encoding="utf-8")) == manifest.to_dict()


def test_unicode_paths_remain_readable_in_json_and_yaml(tmp_path: Path) -> None:
    context = make_context(tmp_path / "科研数据")
    manifest = RunManifest.from_context(context)
    json_path = manifest.write_json(tmp_path / "结果" / "运行.json")
    yaml_path = manifest.write_yaml(tmp_path / "结果" / "运行.yaml")
    assert "科研数据" in json_path.read_text(encoding="utf-8")
    assert "科研数据" in yaml_path.read_text(encoding="utf-8")


@pytest.mark.parametrize("method,suffix", (("write_json", "json"), ("write_yaml", "yaml")))
def test_existing_manifest_is_not_overwritten_by_default(
    tmp_path: Path, method: str, suffix: str
) -> None:
    manifest = RunManifest.from_context(make_context(tmp_path))
    path = tmp_path / f"run.{suffix}"
    path.write_text("original", encoding="utf-8")
    with pytest.raises(OutputValidationError, match="already exists"):
        getattr(manifest, method)(path)
    assert path.read_text(encoding="utf-8") == "original"


@pytest.mark.parametrize("method,suffix", (("write_json", "json"), ("write_yaml", "yaml")))
def test_overwrite_true_replaces_existing_manifest(
    tmp_path: Path, method: str, suffix: str
) -> None:
    manifest = RunManifest.from_context(make_context(tmp_path))
    path = tmp_path / f"run.{suffix}"
    path.write_text("original", encoding="utf-8")
    getattr(manifest, method)(path, overwrite=True)
    assert path.read_text(encoding="utf-8") != "original"


def test_failed_atomic_replace_removes_owned_temp_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest = RunManifest.from_context(make_context(tmp_path))
    destination = tmp_path / "run.json"

    def fail_replace(source: object, target: object) -> None:
        raise OSError("replace failed")

    monkeypatch.setattr("hifivar.serialization.os.replace", fail_replace)
    with pytest.raises(OutputValidationError, match="replace failed"):
        manifest.write_json(destination)
    assert not destination.exists()
    assert list(tmp_path.glob(".run.json.*.tmp")) == []


def test_json_reader_returns_manifest_data_object(tmp_path: Path) -> None:
    original = RunManifest.from_context(make_context(tmp_path))
    path = original.write_json(tmp_path / "run.json")
    assert RunManifest.from_json(path).to_dict() == original.to_dict()


def test_yaml_reader_returns_manifest_data_object(tmp_path: Path) -> None:
    original = RunManifest.from_context(make_context(tmp_path))
    path = original.write_yaml(tmp_path / "run.yaml")
    assert RunManifest.from_yaml(path).to_dict() == original.to_dict()


def test_reader_rejects_missing_required_key(tmp_path: Path) -> None:
    path = tmp_path / "bad.json"
    path.write_text('{"schema_version": "1.0"}', encoding="utf-8")
    with pytest.raises(InputValidationError, match="missing required"):
        RunManifest.from_json(path)


def test_single_sample_manifest_has_null_sheet_source(tmp_path: Path) -> None:
    manifest = RunManifest.from_context(make_context(tmp_path))
    assert manifest.source_sample_sheet is None


def test_sheet_source_is_absolute(tmp_path: Path) -> None:
    manifest = RunManifest.from_context(make_context(tmp_path, with_sheet=True))
    assert manifest.source_sample_sheet is not None
    assert Path(manifest.source_sample_sheet).is_absolute()


def test_reference_compatibility_status_is_preserved(tmp_path: Path) -> None:
    manifest = RunManifest.from_context(make_context(tmp_path))
    assert manifest.reference_compatibility[0]["status"] == "not_applicable"  # type: ignore[index]


def test_from_context_rejects_wrong_object() -> None:
    with pytest.raises(InputValidationError, match="AnalysisContext"):
        RunManifest.from_context(object())  # type: ignore[arg-type]


def test_unsupported_config_value_type_is_rejected(tmp_path: Path) -> None:
    config = {"future": {"value": object()}}
    with pytest.raises(InputValidationError, match="unsupported type"):
        RunManifest.from_context(make_context(tmp_path, config=config))
