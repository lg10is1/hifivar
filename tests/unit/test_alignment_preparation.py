"""Real subprocess witnesses over synthetic bytes, not native BAM validation."""
from copy import deepcopy
import csv
import io
import json
from pathlib import Path
import sys

import pytest
import yaml

from hifivar import alignment_preparation as prep, launcher, package_resources
from hifivar.attempt_receipts import artifact_identity, create_attempt_context, read_step_receipt, _context
from hifivar.command import CommandRunner
from hifivar.command_observation import observe_commands
from hifivar.exceptions import InputValidationError, OutputValidationError


def put(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


@pytest.fixture
def case(tmp_path, monkeypatch):
    source = tmp_path / "reads.fastq"
    source.write_text("@m/1/ccs\nACGT\n+\nIIII\n")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n")
    root = tmp_path / "run"
    root.mkdir()
    workflow = tmp_path / "workflow"
    workflow.mkdir()
    (workflow / "Snakefile").write_text("rule all:\n    input: []\n")
    monkeypatch.setattr(package_resources, "installed_workflow_root", lambda: workflow)
    runtime = launcher.runtime_identity(workflow)
    config = {"alignment": {"enabled": True, "tool": "pbmm2", "pbmm2_executable": sys.executable},
              "reference": {"fasta": str(ref)}, "samples": {"sheet": str(root / "samples.tsv")}}
    text = "sample_id\tinput\tinput_type\tsex\nS1\t" + str(source) + "\tfastq\tfemale\n"
    (root / "samples.tsv").write_text(text)
    (root / "effective.yaml").write_text(yaml.safe_dump(config))
    manifest = {"identity": {"schema": 2, "run_dir": str(root), "config": config, "sample_sheet": text,
                            "runtime": runtime, "inputs": [artifact_identity(source), artifact_identity(ref)]},
                "control_sha256": {name: artifact_identity(root / name)["sha256"] for name in ("samples.tsv", "effective.yaml")}}
    put(root / "launch_manifest.json", manifest)
    attempt = root / "attempts" / "attempt_one"
    put(attempt / "status.json", {"status": "ALIGNMENT_PREPARATION", "dry_run": False})
    bam = root / "results/alignment/S1/attempt_native/S1.bam"
    bam.parent.mkdir(parents=True)
    bam.write_bytes(b"synthetic-bam-not-native")
    index = bam.with_suffix(".bam.bai")
    index.write_bytes(b"synthetic-index-not-native")
    return root, attempt, {"sample_id": "S1", "kind": "FASTQ", "inputs": [str(source)]}, bam, index


def observed_validation(bam, *, mapper=None):
    """Run harmless subprocesses; command suffixes model expected contracts."""
    runner = CommandRunner()
    if mapper:
        runner.run([sys.executable, "-c", "pass", "align", *mapper])
    steps = []
    for args in (("quickcheck", str(bam)), ("view", "-H", str(bam)), ("view", "-X", str(bam))):
        command = [sys.executable, "-c", "pass", *args]
        runner.run(command)
        steps.append({"command": command, "returncode": 0})
    return {"sample_id": "S1", "quickcheck": "PASS", "header_identity": "PASS", "indexed_access": "PASS",
            "samtools_version": "synthetic-process-contract", "steps": steps}


def publish_new(case):
    root, attempt, route, bam, index = case
    with observe_commands() as commands:
        validation = observed_validation(bam, mapper=[str(root.parent / "ref.fa"), *route["inputs"]])
        return prep.publish_preparation(root, attempt, route, bam, index, validation, commands,
                                        action="ALIGN_NEW", execution={"tool_version": "synthetic-not-a-mapper"})


def calling(case):
    root, attempt, route, bam, index = case
    rows = list(csv.DictReader(io.StringIO(json.loads((root / "launch_manifest.json").read_text())["identity"]["sample_sheet"]), delimiter="\t"))
    rows[0].update(input=str(bam), input_type="bam")
    with (root / "calling.samples.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    config = yaml.safe_load((root / "effective.yaml").read_text())
    config["samples"]["sheet"] = str(root / "calling.samples.tsv")
    (root / "calling.yaml").write_text(yaml.safe_dump(config))


def test_new_alignment_context_and_report_reader(case):
    root, attempt, route, bam, index = case
    record = publish_new(case)
    assert record["execution_kind"] == "NATIVE"
    assert record["producer_attempt_id"] == attempt.name
    calling(case)
    binding = prep.seal_handoff(root, attempt)
    assert binding["declared_sample_order"] == ["S1"]
    context = create_attempt_context(root, attempt, dry_run=False)
    assert read_step_receipt(root, attempt.name, "alignment", ("S1",)) == record
    assert json.loads(context.read_text())["alignment_handoff"] == artifact_identity(attempt / "alignment_handoff.json")


def test_unobserved_results_cannot_certify_alignment(case):
    root, attempt, route, bam, index = case
    with observe_commands() as commands, pytest.raises(OutputValidationError, match="observed"):
        prep.publish_preparation(root, attempt, route, bam, index, {}, commands,
                                 action="ALIGN_NEW", execution={"tool_version": "fake"})
    assert not prep.preparation_path(attempt, "S1").exists()


def test_old_observation_cannot_certify_alignment(case):
    root, attempt, route, bam, index = case
    with observe_commands() as old:
        check = observed_validation(bam)
    with observe_commands(), pytest.raises(OutputValidationError, match="current"):
        prep.publish_preparation(root, attempt, route, bam, index, check, old, action="ALIGN_NEW")


def test_validation_only_cannot_be_new_alignment(case):
    root, attempt, route, bam, index = case
    with observe_commands() as commands:
        check = observed_validation(bam)
        with pytest.raises(OutputValidationError, match="mapper"):
            prep.publish_preparation(root, attempt, route, bam, index, check, commands,
                                     action="ALIGN_NEW", execution={"tool_version": "fake"})


@pytest.mark.parametrize("field,value", [("sample_id", "S2"), ("quickcheck", "FAIL"), ("header_identity", "FAIL"), ("indexed_access", "FAIL"), ("steps", [])])
def test_invalid_or_unobserved_bam_validation(case, field, value):
    root, attempt, route, bam, index = case
    with observe_commands() as commands:
        check = observed_validation(bam)
        check[field] = value
        with pytest.raises(OutputValidationError):
            prep.publish_preparation(root, attempt, route, bam, index, check, commands, action="ALIGN_NEW")


@pytest.mark.parametrize("which", ["calling.yaml", "calling.samples.tsv", "receipt", "handoff"])
def test_calling_controls_and_receipts_bound_to_context(case, which):
    root, attempt, route, bam, index = case
    publish_new(case)
    calling(case)
    prep.seal_handoff(root, attempt)
    context = create_attempt_context(root, attempt, dry_run=False)
    path = (prep.preparation_path(attempt, "S1") if which == "receipt" else
            attempt / "alignment_handoff.json" if which == "handoff" else root / which)
    path.write_text(path.read_text() + "\n ")
    with pytest.raises(InputValidationError):
        _context(context, live=False)


@pytest.mark.parametrize("which", ["bam", "index"])
def test_selected_artifact_mutation_rejected(case, which):
    root, attempt, route, bam, index = case
    publish_new(case)
    (bam if which == "bam" else index).write_bytes(b"changed")
    with pytest.raises(OutputValidationError, match="changed"):
        prep.read_preparation(root, attempt, "S1")


@pytest.mark.parametrize("change", ["order", "metadata", "config", "missing_receipt"])
def test_handoff_rejects_unbound_selection(case, change):
    root, attempt, route, bam, index = case
    if change != "missing_receipt":
        publish_new(case)
    calling(case)
    if change in {"order", "metadata"}:
        sheet = root / "calling.samples.tsv"
        sheet.write_text(sheet.read_text().replace("S1", "S2") if change == "order" else sheet.read_text().replace("female", "male"))
    if change == "config":
        path = root / "calling.yaml"
        value = yaml.safe_load(path.read_text())
        value["alignment"]["tool"] = "minimap2"
        path.write_text(yaml.safe_dump(value))
    with pytest.raises(InputValidationError):
        prep.seal_handoff(root, attempt)


def test_prior_reuse_preserves_original_producer_and_validation_only(case):
    root, first, route, bam, index = case
    original = publish_new(case)
    completed = root / "results/alignment/S1/completed.json"
    put(completed, {"sample_id": "S1", "inputs": route["inputs"], "bam": str(bam),
                    "artifacts": list(original["selected_artifacts"].values()), "producer_attempt_id": first.name,
                    "producer_receipt": artifact_identity(prep.preparation_path(first, "S1"))})
    before = artifact_identity(completed)
    second = root / "attempts/attempt_two"
    put(second / "status.json", {"status": "ALIGNMENT_PREPARATION", "dry_run": False})
    with observe_commands() as commands:
        check = observed_validation(bam)
        record = prep.publish_preparation(root, second, route, bam, index, check, commands,
                                          action="VALIDATE_PRIOR_ALIGNMENT_REUSE", prior=completed)
    assert record["producer_attempt_id"] == first.name
    assert record["execution_kind"] == "VALIDATION_ONLY"
    assert set(record["outputs"]) == {"validation"}
    assert record["selected_artifacts"] == original["selected_artifacts"]
    assert artifact_identity(completed) == before


def test_legacy_reuse_does_not_retrofit_producer_identity(case):
    root, attempt, route, bam, index = case
    prior = root / "results/alignment/S1/completed.json"
    put(prior, {"sample_id": "S1", "inputs": route["inputs"], "bam": str(bam), "artifacts": [artifact_identity(bam), artifact_identity(index)]})
    with observe_commands() as commands:
        record = prep.publish_preparation(root, attempt, route, bam, index, observed_validation(bam), commands,
                                          action="VALIDATE_PRIOR_ALIGNMENT_REUSE", prior=prior)
    assert record["producer_attempt_id"] is None
    assert record["metadata"]["producer_origin"] == "NOT_PROVEN_LEGACY_OR_EXTERNAL"


def test_external_reuse_requires_frozen_bam_and_index(case):
    root, attempt, route, bam, index = case
    path = root / "launch_manifest.json"
    manifest = json.loads(path.read_text())
    route.update(kind="ALIGNED_BAM", inputs=[str(bam)])
    text = manifest["identity"]["sample_sheet"].replace(str(root.parent / "reads.fastq"), str(bam)).replace("fastq", "bam")
    manifest["identity"].update(sample_sheet=text, inputs=[artifact_identity(bam), artifact_identity(index)])
    (root / "samples.tsv").write_text(text)
    manifest["control_sha256"]["samples.tsv"] = artifact_identity(root / "samples.tsv")["sha256"]
    put(path, manifest)
    with observe_commands() as commands:
        record = prep.publish_preparation(root, attempt, route, bam, index, observed_validation(bam), commands,
                                          action="VALIDATE_EXTERNAL_REUSE")
    assert record["execution_kind"] == "VALIDATION_ONLY"
    assert record["producer_attempt_id"] is None
    assert set(record["outputs"]) == {"validation"}


def test_legacy_context_without_preparation_does_not_infer_alignment(case):
    root, attempt, *_ = case
    create_attempt_context(root, attempt, dry_run=False)
    assert read_step_receipt(root, attempt.name, "alignment", ("S1",)) is None
