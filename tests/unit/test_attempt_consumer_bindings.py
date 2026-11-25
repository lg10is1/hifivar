"""Synthetic source/stage linkage tests, not a Python security sandbox.

Positive NATIVE tests execute a real isolated Python child through CommandRunner.
The runtime fingerprint and small owned control documents are unit fixtures;
these tests are not native biological-tool or scientific-accuracy validation.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest
import yaml

from hifivar import launcher, package_resources
import hifivar.attempt_receipts as receipt_module
from hifivar.attempt_receipts import (
    ATTEMPT_CONTEXT_ENV, artifact_identity, attempt_step, create_attempt_context,
    read_step_receipt, receipt_name,
)
from hifivar.command import CommandRunner
from hifivar.exceptions import InputValidationError, OutputValidationError


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


def native_outputs(paths):
    for path in paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    env = {key: value for key, value in os.environ.items() if key.upper() in {"SYSTEMROOT", "WINDIR"}}
    CommandRunner().run(
        [sys.executable, "-X", "utf8", "-I", "-B", "-c",
         "import json,sys;from pathlib import Path;"
         "[(Path(p).write_bytes(('synthetic:'+str(i)).encode())) for i,p in enumerate(json.loads(sys.argv[1]))]",
         json.dumps([str(path) for path in paths.values()])], inherit_env=False, env=env, timeout=20,
    )


def fixture_attempt(tmp_path, monkeypatch, *, track, generated=True, plain=False, second_state=None):
    root, external = tmp_path / "run", tmp_path / "external"
    root.mkdir()
    external.mkdir()
    attempt = root / "attempts/attempt_test"
    attempt.mkdir(parents=True)
    annotation = track.startswith("annotation.")
    group, category = ("annotation", "small") if annotation else ("cohort", "sv")
    samples = ("S1", "S2") if second_state else ("S1",)
    config = {
        "workflow": {"downstream_source_plan": str(root / "downstream_plan.yaml") if generated else None},
        "small": {"enabled": True},
        "annotation": {"enabled": annotation, "vep_enabled": True, "overlap_enabled": True},
        "cohort": {"enabled": not annotation, "sv": {"enabled": True}},
    }
    input_manifest = (attempt / "downstream" if generated else root) / (group + ".inputs.tsv")
    config[group]["input_manifest"] = None if generated else str(input_manifest)
    controls = ["effective.yaml", "samples.tsv"]
    sheet = "sample_id\tinput\n" + "".join(sample + "\t/external/" + sample + ".bam\n" for sample in samples)
    (root / "samples.tsv").write_text(sheet, encoding="utf-8")
    (root / "effective.yaml").write_text(yaml.safe_dump(config), encoding="utf-8")
    runtime = {"hifivar_version": "synthetic-unit-runtime", "payload_sha256": "a" * 64}
    monkeypatch.setattr(launcher, "runtime_identity", lambda _: runtime)
    monkeypatch.setattr(package_resources, "installed_workflow_root", lambda: root)
    if generated:
        (root / "downstream_plan.yaml").write_text("schema_version: 1\n", encoding="utf-8")
        controls.append("downstream_plan.yaml")
    source = (root / "results/generated.vcf.gz") if generated else external / ("input.vcf" if plain else "input.vcf.gz")
    index = Path(str(source) + ".tbi") if not plain or generated else None
    source_pair = {"vcf": source, "index": index}
    if not generated:
        source.write_bytes(b"synthetic-source")
        if index is not None:
            index.write_bytes(b"synthetic-index")

    def manifest_text():
        if annotation:
            return ("sample\tvariant_category\tsource_vcf\tsource_tool\tsource_variant_ids\n"
                    + f"S1\tsmall\t{source}\tunit-source\tv1\n")
        rows = f"S1\tsv\tCALLED\t{source}\t{index}\tunit-source\t1\tGRCh38\t\n"
        if second_state:
            rows += f"S2\tsv\t{second_state}\t\t\t\t\tGRCh38\t\n"
        return "sample\ttrack\tstate\tsource_path\tindex_path\tsource_tool\tsource_version\treference_build\tcatalog_id\n" + rows

    if not generated:
        input_manifest.write_text(manifest_text(), encoding="utf-8")
        controls.append(input_manifest.name)
    identity = {"schema": 2, "run_dir": str(root), "config": config, "sample_sheet": sheet,
                "runtime": runtime, "inputs": [artifact_identity(path) for path in (source, index) if path is not None] if not generated else []}
    write_json(root / "launch_manifest.json", {"identity": identity, "control_sha256": {
        name: artifact_identity(root / name)["sha256"] for name in controls}})
    write_json(attempt / "status.json", {"status": "PRODUCER_DAG_RUNNING" if generated else "DAG_RUNNING", "dry_run": False})
    context = create_attempt_context(root, attempt, dry_run=False)
    monkeypatch.setenv(ATTEMPT_CONTEXT_ENV, str(context))
    if generated:
        producer_outputs = {"vcf": source, "index": index, "gvcf": root / "results/generated.g.vcf.gz",
                            "gvcf_index": root / "results/generated.g.vcf.gz.tbi"}
        with attempt_step("small", ("S1",), outputs=producer_outputs) as receipt:
            native_outputs(producer_outputs)
            receipt.complete(outputs=producer_outputs, tool_versions={"unit-process": "1"}, validation={"validated": True})
        input_manifest.parent.mkdir()
        input_manifest.write_text(manifest_text(), encoding="utf-8")
        effective = dict(config)
        effective[group] = {**config[group], "input_manifest": str(input_manifest)}
        effective_path = input_manifest.parent / "effective.yaml"
        effective_path.write_text(yaml.safe_dump(effective), encoding="utf-8")
        producer_path = attempt / "receipts" / receipt_name("small", ("S1",))
        materialization = {
            "schema_version": 1, "kind": "hifivar_downstream_materialization", "attempt_id": attempt.name,
            "context_sha256": artifact_identity(context)["sha256"],
            "plan_sha256": artifact_identity(root / "downstream_plan.yaml")["sha256"],
            "controls": {p.name: artifact_identity(p) for p in (input_manifest, effective_path)},
            "source_receipts": [{**artifact_identity(producer_path), "track": "small", "samples": ["S1"], "output_role": "vcf"}],
            "source_validation": {group: [{"sample_id": "S1", "category" if annotation else "track": category, "source_origin": "GENERATED"}]},
        }
        write_json(input_manifest.parent / "materialization.json", materialization)
        write_json(attempt / "status.json", {"status": "DOWNSTREAM_DAG_RUNNING", "dry_run": False,
            "downstream_materialization_sha256": artifact_identity(input_manifest.parent / "materialization.json")["sha256"]})
    outputs = ({"tsv": root / "results/annotation/result.tsv"} if annotation else
               {"sites_tsv": root / "results/cohort/sites.tsv", "sample_matrix_tsv": root / "results/cohort/matrix.tsv",
                "provenance": root / "results/cohort/provenance.json"})
    return SimpleNamespace(root=root, attempt=attempt, track=track, samples=samples, outputs=outputs,
        unit="small" if annotation else None, manifest=input_manifest,
        inputs={"S1": source_pair, **({"S2": {}} if second_state else {})},
        source=source, index=index, generated=generated)


def step(fixture, **overrides):
    options = {"outputs": fixture.outputs, "unit_id": fixture.unit, "inputs": fixture.inputs,
               "input_manifest": fixture.manifest}
    options.update(overrides)
    return attempt_step(fixture.track, fixture.samples, **options)


def finish(fixture, receipt):
    native = fixture.track == "annotation.vep"
    if native:
        native_outputs(fixture.outputs)
    else:
        for path in fixture.outputs.values():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"synthetic-in-process-table")
    receipt.complete(outputs=fixture.outputs, tool_versions={"unit-process": "1"}, validation={"validated": True},
                     execution_kind="NATIVE" if native else "IN_PROCESS")


@pytest.mark.parametrize("track", ["annotation.vep", "cohort.sv"])
@pytest.mark.parametrize("generated", [False, True])
def test_correct_native_and_inprocess_sources_bind_to_their_stage(tmp_path, monkeypatch, track, generated):
    fixture = fixture_attempt(tmp_path, monkeypatch, track=track, generated=generated)
    with step(fixture) as receipt:
        finish(fixture, receipt)
    value = read_step_receipt(fixture.root, fixture.attempt.name, track, fixture.samples, unit_id=fixture.unit)
    binding = value["consumer_binding"]
    assert binding["source_inputs"]["S1"]["vcf"] == artifact_identity(fixture.source)
    assert binding["mode"] == ("MATERIALIZED_PLAN" if generated else "EXTERNAL_READY")
    assert binding["execution_stage"] == ("DOWNSTREAM_DAG_RUNNING" if generated else "DAG_RUNNING")
    assert value["scientific_accuracy"] == "NOT_ASSESSED"


def test_legacy_plain_annotation_has_explicit_no_index_witness(tmp_path, monkeypatch):
    fixture = fixture_attempt(tmp_path, monkeypatch, track="annotation.vep", generated=False, plain=True)
    with step(fixture) as receipt:
        finish(fixture, receipt)
    value = read_step_receipt(fixture.root, fixture.attempt.name, fixture.track, fixture.samples, unit_id=fixture.unit)
    assert value["consumer_binding"]["source_inputs"]["S1"]["index"] is None
    assert not Path(str(fixture.source) + ".tbi").exists()


@pytest.mark.parametrize("track", ["annotation.vep", "cohort.sv"])
@pytest.mark.parametrize("stage", ["DAG_RUNNING", "PRODUCER_DAG_RUNNING", "DOWNSTREAM_PREPARATION"])
def test_consumer_cannot_execute_early_in_a_generated_plan(tmp_path, monkeypatch, track, stage):
    fixture = fixture_attempt(tmp_path, monkeypatch, track=track)
    status = json.loads((fixture.attempt / "status.json").read_text(encoding="utf-8"))
    status["status"] = stage
    write_json(fixture.attempt / "status.json", status)
    with pytest.raises(InputValidationError):
        with step(fixture):
            pytest.fail("Early-stage consumer reached execution")


@pytest.mark.parametrize("track", ["annotation.vep", "cohort.sv"])
@pytest.mark.parametrize("change", ["missing_inputs", "wrong_vcf", "wrong_index", "wrong_manifest"])
def test_actual_request_source_set_must_match_not_just_metadata(tmp_path, monkeypatch, track, change):
    fixture = fixture_attempt(tmp_path, monkeypatch, track=track)
    options = {}
    if change == "missing_inputs":
        options["inputs"] = None
    elif change == "wrong_manifest":
        other = fixture.root / "other.tsv"
        other.write_bytes(fixture.manifest.read_bytes())
        options["input_manifest"] = other
    else:
        other = fixture.root / "other.vcf.gz"
        other.write_bytes(fixture.source.read_bytes())
        options["inputs"] = {"S1": {**fixture.inputs["S1"], "vcf" if change == "wrong_vcf" else "index": other}}
    with pytest.raises(InputValidationError):
        with step(fixture, **options):
            pytest.fail("Wrong native source reached execution")


@pytest.mark.parametrize("track", ["annotation.vep", "cohort.sv"])
@pytest.mark.parametrize("change", ["marker", "materialization_and_marker", "effective", "manifest", "source", "index"])
def test_change_during_execution_cannot_seal_a_consumer_receipt(tmp_path, monkeypatch, track, change):
    fixture = fixture_attempt(tmp_path, monkeypatch, track=track)
    with pytest.raises((InputValidationError, OutputValidationError)):
        with step(fixture) as receipt:
            if change == "marker":
                status = json.loads((fixture.attempt / "status.json").read_text(encoding="utf-8"))
                status["downstream_materialization_sha256"] = "0" * 64
                write_json(fixture.attempt / "status.json", status)
            elif change == "materialization_and_marker":
                materialization = fixture.manifest.parent / "materialization.json"
                document = json.loads(materialization.read_text(encoding="utf-8"))
                document["changed_after_entry"] = True
                write_json(materialization, document)
                status = json.loads((fixture.attempt / "status.json").read_text(encoding="utf-8"))
                status["downstream_materialization_sha256"] = artifact_identity(materialization)["sha256"]
                write_json(fixture.attempt / "status.json", status)
            else:
                path = {"effective": fixture.manifest.parent / "effective.yaml", "manifest": fixture.manifest,
                        "source": fixture.source, "index": fixture.index}[change]
                with path.open("ab") as handle:
                    handle.write(b"changed")
            finish(fixture, receipt)
    assert not (fixture.attempt / "receipts" / receipt_name(track, fixture.samples, fixture.unit)).exists()


@pytest.mark.parametrize("track", ["annotation.vep", "cohort.sv"])
@pytest.mark.parametrize("change", ["missing_binding", "old_stage", "wrong_input_hash", "marker"])
def test_reader_does_not_upgrade_unbound_or_tampered_consumers(tmp_path, monkeypatch, track, change):
    fixture = fixture_attempt(tmp_path, monkeypatch, track=track)
    with step(fixture) as receipt:
        finish(fixture, receipt)
    path = fixture.attempt / "receipts" / receipt_name(track, fixture.samples, fixture.unit)
    record = json.loads(path.read_text(encoding="utf-8"))
    if change == "missing_binding":
        record.pop("consumer_binding")
    elif change == "old_stage":
        record["consumer_binding"]["execution_stage"] = "PRODUCER_DAG_RUNNING"
    elif change == "wrong_input_hash":
        record["consumer_binding"]["source_inputs"]["S1"]["vcf"]["sha256"] = "0" * 64
    else:
        status = json.loads((fixture.attempt / "status.json").read_text(encoding="utf-8"))
        status["downstream_materialization_sha256"] = "0" * 64
        write_json(fixture.attempt / "status.json", status)
    write_json(path, record)
    with pytest.raises((InputValidationError, OutputValidationError)):
        read_step_receipt(fixture.root, fixture.attempt.name, track, fixture.samples, unit_id=fixture.unit, verify_outputs=False)


def test_external_source_changed_before_entry_is_not_adopted(tmp_path, monkeypatch):
    fixture = fixture_attempt(tmp_path, monkeypatch, track="cohort.sv", generated=False)
    fixture.source.write_bytes(b"different-bytes")
    with pytest.raises(InputValidationError, match="frozen launch identity"):
        with step(fixture):
            pytest.fail("Changed external source reached execution")


def test_noncallable_sample_requires_explicit_empty_inputs_not_homref(tmp_path, monkeypatch):
    fixture = fixture_attempt(tmp_path, monkeypatch, track="cohort.sv", second_state="MISSING_INPUT")
    with step(fixture) as receipt:
        finish(fixture, receipt)
    value = read_step_receipt(fixture.root, fixture.attempt.name, fixture.track, fixture.samples)
    assert value["consumer_binding"]["source_inputs"]["S2"] == {}


def test_no_calls_without_real_vcf_is_still_rejected(tmp_path, monkeypatch):
    fixture = fixture_attempt(tmp_path, monkeypatch, track="cohort.sv", second_state="NO_CALLS")
    with pytest.raises(InputValidationError, match="NO_CALLS"):
        with step(fixture):
            pytest.fail("Missing NO_CALLS artifact cannot become hom-ref")


@pytest.mark.parametrize("track", ["annotation.vep", "cohort.sv"])
def test_unused_producer_gvcf_bytes_are_not_rehashed_as_consumed_sources(tmp_path, monkeypatch, track):
    fixture = fixture_attempt(tmp_path, monkeypatch, track=track)
    unused = {fixture.root / "results/generated.g.vcf.gz", fixture.root / "results/generated.g.vcf.gz.tbi"}
    for path in unused:
        path.write_bytes(b"changed-unconsumed-output-not-certified-by-this-consumer")
    original = receipt_module.artifact_identity
    checked = []

    def selected_only(path):
        assert Path(path) not in unused, "Unconsumed large gVCF was unnecessarily scanned"
        checked.append(Path(path))
        return original(path)

    monkeypatch.setattr(receipt_module, "artifact_identity", selected_only)
    with step(fixture) as receipt:
        finish(fixture, receipt)
    value = read_step_receipt(fixture.root, fixture.attempt.name, track, fixture.samples, unit_id=fixture.unit)
    assert fixture.source in checked and fixture.index in checked
    assert len(value["consumer_binding"]["producer_receipts"]) == 1
    assert value["consumer_binding"]["source_validation_scope"] == "ACTUAL_SELECTED_SOURCE_ROLES_ONLY_OTHER_OUTPUT_BYTES_NOT_REHASHED"


@pytest.mark.parametrize("track", ["annotation.vep", "cohort.sv"])
@pytest.mark.parametrize("role", ["vcf", "index"])
def test_reading_consumer_rechecks_selected_source_bytes_even_without_output_verification(tmp_path, monkeypatch, track, role):
    fixture = fixture_attempt(tmp_path, monkeypatch, track=track)
    with step(fixture) as receipt:
        finish(fixture, receipt)
    fixture.inputs["S1"][role].write_bytes(b"changed-actually-consumed-source")
    with pytest.raises((InputValidationError, OutputValidationError)):
        read_step_receipt(fixture.root, fixture.attempt.name, track, fixture.samples,
                          unit_id=fixture.unit, verify_outputs=False)
