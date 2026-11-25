"""Report genuine tiny receipt/API execution, not fabricated native success."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path

import pytest
import yaml

from hifivar import launcher
from hifivar.attempt_receipts import (
    ATTEMPT_CONTEXT_ENV, artifact_identity, attempt_step, create_attempt_context,
    read_step_receipt, receipt_name, required_output_roles,
)
from hifivar.command import CommandRunner
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.workflow_report import build_launcher_report


def now():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False), encoding="utf-8")


def fingerprint(root):
    return {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in root.rglob("*") if path.is_file()}


def payload(result):
    return json.loads(Path(result["output_paths"]["json"]).read_text(encoding="utf-8"))


def tracks(result):
    return {item["name"]: item for item in payload(result)["tracks"]}


@pytest.fixture
def make_run(tmp_path, monkeypatch):
    def make(config=None, *, samples=("S2", "S1"), annotation_rows=(), source_plan=None,
             cohort_rows=None, dry_run=False):
        root = tmp_path / "run"
        root.mkdir()
        config = deepcopy(config if config is not None else {"small": {"enabled": True}})
        config.setdefault("reference", {"build": "GRCh38", "fasta": "/external/reference.fa"})
        controls = {}
        declared_inputs = {}

        def source_pair(sample, category):
            source = root / "work" / "declared_sources" / f"{sample}.{category}.vcf.gz"
            index = Path(str(source) + ".tbi")
            source.parent.mkdir(parents=True, exist_ok=True)
            source.write_bytes(b"synthetic-input-bytes-not-a-native-VCF-test")
            index.write_bytes(b"synthetic-index-bytes-not-a-native-index-test")
            for item in (source, index):
                declared_inputs[str(item)] = artifact_identity(item)
            return source, index

        if config.get("sv", {}).get("harmonization", {}).get("enabled"):
            config['sv']['harmonization']['input_vcfs'] = {
                sample: {'pav': str(source_pair(sample, 'pav.assembly.sv')[0])} for sample in samples}
        if annotation_rows:
            path = root / "annotation.inputs.tsv"
            config["annotation"]["input_manifest"] = str(path)
            with path.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                writer.writerow(("sample", "variant_category", "source_vcf", "source_tool", "source_variant_ids"))
                for sample, category in annotation_rows:
                    source, _ = source_pair(sample, category)
                    writer.writerow((sample, category, str(source), "declared-source", "v1"))
            controls[path.name] = artifact_identity(path)["sha256"]
        if config.get("cohort", {}).get("enabled") and source_plan is None:
            path = root / "cohort.inputs.tsv"
            config["cohort"]["input_manifest"] = str(path)
            with path.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                writer.writerow(("sample", "track", "state", "source_path", "index_path", "source_tool",
                                 "source_version", "reference_build", "catalog_id"))
                rows = cohort_rows
                if rows is None:
                    rows = [(sample, track, "CALLED")
                            for track in ("small_variants", "sv", "tr")
                            if config["cohort"].get(track, {}).get("enabled")
                            for sample in samples]
                for sample, track, state in rows:
                    source, index = source_pair(sample, track)
                    writer.writerow((sample, track, state, str(source),
                                      str(index), "declared-source", "1.0", "GRCh38",
                                     "sha256:" + "c" * 64 if track == "tr" else ""))
            controls[path.name] = artifact_identity(path)["sha256"]
        if source_plan is not None:
            path = root / "downstream_plan.yaml"
            config.setdefault("workflow", {})["downstream_source_plan"] = str(path)
            path.write_text(yaml.safe_dump(source_plan), encoding="utf-8")
            controls[path.name] = artifact_identity(path)["sha256"]
        sheet = "sample_id\tinput\n" + "".join(f"{sample}\t/external/{sample}.bam\n" for sample in samples)
        (root / "effective.yaml").write_text(yaml.safe_dump(config), encoding="utf-8")
        (root / "samples.tsv").write_text(sheet, encoding="utf-8")
        for name in ("effective.yaml", "samples.tsv"):
            controls[name] = artifact_identity(root / name)["sha256"]
        runtime = {"payload_sha256": "a" * 64, "hifivar_version": "synthetic-receipt-test"}
        monkeypatch.setattr(launcher, "runtime_identity", lambda _: runtime)
        import hifivar.package_resources
        monkeypatch.setattr(hifivar.package_resources, "installed_workflow_root", lambda: root)
        manifest = {"identity": {"schema": 2, "kind": "hifivar_two_stage_launcher", "run_dir": str(root),
                                 "config": config, "sample_sheet": sheet, "runtime": runtime,
                                 "inputs": list(declared_inputs.values())},
                    "control_sha256": controls}
        write_json(root / "launch_manifest.json", manifest)
        attempt = root / "attempts" / "attempt_a"
        write_json(attempt / "status.json", {"started": now(), "status": "PRODUCER_DAG_RUNNING" if source_plan is not None else "DAG_RUNNING", "dry_run": dry_run,
                                            "returncode": None, "finished": None})
        context = create_attempt_context(root, attempt, dry_run=dry_run)
        monkeypatch.setenv(ATTEMPT_CONTEXT_ENV, str(context))
        return root, attempt, config, samples
    return make


def consumer_inputs(run, track, samples, unit_id):
    """Explicit synthetic actual inputs, not a self-reported metadata claim."""
    if track == 'sv.harmonization':
        root, _, config, _ = run
        sample = samples[0]
        source = config['sv']['harmonization']['input_vcfs'][sample]['pav']
        return {'inputs': {'pav': {'vcf': Path(source), 'index': Path(source + '.tbi')}},
                'input_manifest': root / 'effective.yaml'}
    if not track.startswith(("annotation.", "cohort.")):
        return {}
    root, attempt, config, _ = run
    group = track.partition(".")[0]
    planned = bool(config.get("workflow", {}).get("downstream_source_plan"))
    owner = attempt / "downstream"
    value = owner / (group + ".inputs.tsv") if planned else config[group].get("input_manifest")
    if not value or not Path(value).is_file():
        return {"inputs": {sample: {} for sample in samples},
                "input_manifest": root / "missing-consumer-manifest.tsv"}
    path = Path(value)
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    inputs = {sample: {} for sample in samples}
    category = unit_id if group == "annotation" else track.partition(".")[2]
    for row in rows:
        if row["sample"] not in samples or row.get("variant_category" if group == "annotation" else "track") != category:
            continue
        source = row.get("source_vcf" if group == "annotation" else "source_path")
        if not source:
            continue
        index = row.get("index_path")
        if group == "annotation" and not planned and source.lower().endswith(".vcf.gz"):
            index = source + ".tbi"
        if group == "annotation" and planned:
            plan = yaml.safe_load((root / "downstream_plan.yaml").read_text(encoding="utf-8"))
            selected = next(item for item in plan["annotation"]
                            if item["sample_id"] == row["sample"] and item["track"] == category)
            producer_track = {"small.deepvariant": "small", "tr.trgt": "tr", "phasing.hiphase": "phasing"}.get(selected["producer"], selected["producer"])
            producer = read_step_receipt(root, attempt.name, producer_track, (row["sample"],))
            index_role = "gvcf_index" if selected["output_role"] == "gvcf" else "index"
            index = producer["outputs"][index_role]["path"]
        inputs[row["sample"]] = {"vcf": Path(source), "index": Path(index) if index else None}
    return {"inputs": inputs, "input_manifest": path}


def complete_step(run, track, samples, *, unit_id=None, metadata=None, tool_versions=None):
    root, _, _, _ = run
    filename = receipt_name(track, samples, unit_id).replace(".json", ".data")
    outputs = {role: root / "results" / (filename + "." + role)
               for role in sorted(required_output_roles(track) | {"test_output"})}
    if "vcf" in outputs:
        outputs["vcf"] = root / "results" / (filename + ".vcf.gz")
        if "index" in outputs:
            outputs["index"] = Path(str(outputs["vcf"]) + ".tbi")
    if "gvcf" in outputs:
        outputs["gvcf"] = root / "results" / (filename + ".g.vcf.gz")
        outputs["gvcf_index"] = Path(str(outputs["gvcf"]) + ".tbi")
    output = outputs["test_output"]
    output.parent.mkdir(parents=True, exist_ok=True)
    environment = {}
    system_root = os.environ.get("SystemRoot") or os.environ.get("SYSTEMROOT")
    if os.name == "nt" and system_root:
        environment["SystemRoot"] = system_root
    with attempt_step(track, samples, outputs=outputs, unit_id=unit_id,
                      **consumer_inputs(run, track, samples, unit_id)) as receipt:
        CommandRunner().run([
            sys.executable, "-I", "-B", "-c",
            "import pathlib,sys;[pathlib.Path(p).write_bytes(b'controlled-synthetic-step-output') for p in sys.argv[1:]]",
            *[str(path) for path in outputs.values()],
        ], env=environment, inherit_env=False)
        receipt.complete(
            outputs=outputs, tool_versions=tool_versions or {"synthetic-python-process": "1"},
            validation={"validated": True, "scope": "SYNTHETIC_TEST_BYTES_NOT_BIOINFORMATICS"},
            metadata=metadata or {"tool_version_sources": {"synthetic-python-process": "EXPLICIT_TEST_PROCESS_VERSION"}},
        )
    return output


def materialize_fixture(run, *, wrong_source=False):
    """Small exact-schema handoff fixture; no real caller/annotation is claimed."""
    root, attempt, config, _ = run
    source_plan = yaml.safe_load((root / "downstream_plan.yaml").read_text(encoding="utf-8"))
    owner = attempt / "downstream"
    owner.mkdir()
    effective = deepcopy(config)
    links = []
    controls = {}
    for group in ("annotation", "cohort"):
        rows = source_plan.get(group, [])
        if not rows:
            continue
        target = owner / (group + ".inputs.tsv")
        effective[group]["input_manifest"] = str(target)
        fields = (("sample", "variant_category", "source_vcf", "source_tool", "source_variant_ids")
                  if group == "annotation" else
                  ("sample", "track", "state", "source_path", "index_path", "source_tool", "source_version", "reference_build", "catalog_id"))
        with target.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(fields)
            for row in rows:
                track = {"small.deepvariant": "small", "tr.trgt": "tr", "phasing.hiphase": "phasing"}.get(row["producer"], row["producer"])
                receipt = read_step_receipt(root, attempt.name, track, (row["sample_id"],))
                role = row["output_role"]
                source = "/external/unrelated.vcf.gz" if wrong_source else receipt["outputs"][role]["path"]
                index = receipt["outputs"]["gvcf_index" if role == "gvcf" else "index"]["path"]
                tool = "deepvariant"
                source_receipt = attempt / "receipts" / receipt_name(track, (row["sample_id"],))
                links.append({**artifact_identity(source_receipt), "producer_attempt_id": attempt.name,
                              "track": track, "samples": [row["sample_id"]], "output_role": role,
                              "source_tool": tool, "source_version": receipt["tool_versions"][tool]})
                if group == "annotation":
                    writer.writerow((row["sample_id"], row["track"], source, tool, ";".join(row.get("variant_ids", []))))
                else:
                    writer.writerow((row["sample_id"], row["track"], "CALLED", source, index, tool,
                                     receipt["tool_versions"][tool], row["reference_build"], row.get("catalog_id") or ""))
        controls[target.name] = artifact_identity(target)
    for group in ("alignment", "small", "sv", "tr", "phasing", "assembly", "assembly_sv"):
        effective.setdefault(group, {})["enabled"] = False
    effective["workflow"]["downstream_source_plan"] = None
    target = owner / "effective.yaml"
    target.write_text(yaml.safe_dump(effective), encoding="utf-8")
    controls[target.name] = artifact_identity(target)
    source_validation = {group: [{"sample_id": row["sample_id"],
                                 "category" if group == "annotation" else "track": row["track"],
                                 "source_origin": "GENERATED"} for row in source_plan.get(group, [])]
                         for group in ("annotation", "cohort")}
    marker = {"schema_version": 1, "kind": "hifivar_downstream_materialization", "attempt_id": attempt.name,
              "context_sha256": artifact_identity(attempt / "execution_context.json")["sha256"],
              "plan_sha256": artifact_identity(root / "downstream_plan.yaml")["sha256"],
              "controls": controls, "source_receipts": links, "source_validation": source_validation}
    write_json(owner / "materialization.json", marker)
    status = json.loads((attempt / "status.json").read_text(encoding="utf-8"))
    status.update(status="DOWNSTREAM_DAG_RUNNING",
                  downstream_materialization_sha256=artifact_identity(owner / "materialization.json")["sha256"])
    write_json(attempt / "status.json", status)
    return owner, marker


def finish(run, state="WORKFLOW_SUCCEEDED", returncode=0):
    _, attempt, _, _ = run
    status = json.loads((attempt / "status.json").read_text(encoding="utf-8"))
    status.update(status=state, returncode=returncode, finished=now(), input_immutability="PASS",
                  snakemake_version="test-scheduler-record")
    write_json(attempt / "status.json", status)


def test_all_explicit_sample_receipts_close_module_and_expose_verified_artifacts(make_run, tmp_path):
    run = make_run()
    root, _, _, samples = run
    outputs = [complete_step(run, "small", (sample,)) for sample in samples]
    finish(run)
    before = fingerprint(root)
    result = build_launcher_report(root, tmp_path / "report")
    small = tracks(result)["small"]
    assert result["status"] == small["status"] == "COMPLETE"
    assert small["qc"]["expected_units"] == small["qc"]["validated_units"] == 2
    assert all(row["producer_attempt_id"] == "attempt_a" for row in small["qc"]["units"])
    assert all(row["provenance_class"] == "THIS_ATTEMPT_EXECUTED_NOT_REUSED" for row in small["qc"]["units"])
    assert all(row["commands"][0]["returncode"] == 0 for row in small["qc"]["units"])
    output_records = [item for item in small["artifacts"] if item["role"] == "test_output"]
    assert {item["path"] for item in output_records} == {str(path) for path in outputs}
    assert all(item["sha256"] == hashlib.sha256(Path(item["path"]).read_bytes()).hexdigest() for item in output_records)
    data = payload(result)
    assert data["provenance"]["scope"] == "LAUNCHER_ATTEMPT_RECEIPTS_AND_OUTPUT_SHA256"
    assert data["provenance"]["source_input_checksum_verification"] == "NOT_RUN_METADATA_ONLY"
    assert data["provenance"]["scientific_accuracy"] == "NOT_ASSESSED"
    assert data["tools"][-1]["validation_status"] == "EXPLICIT_TEST_PROCESS_VERSION"
    assert fingerprint(root) == before


def test_missing_one_sample_receipt_keeps_module_partial(make_run, tmp_path):
    run = make_run()
    complete_step(run, "small", ("S1",))
    finish(run)
    result = build_launcher_report(run[0], tmp_path / "partial-report")
    small = tracks(result)["small"]
    assert small["status"] == result["status"] == "PARTIAL"
    assert [unit["status"] for unit in small["qc"]["units"]] == ["MISSING_RECEIPT", "COMPLETE"]
    assert small["qc"]["validated_units"] == 1


def test_context_without_receipts_does_not_upgrade_scope(make_run, tmp_path):
    run = make_run()
    finish(run)
    result = build_launcher_report(run[0], tmp_path / "missing-report")
    assert result["status"] == "PARTIAL"
    assert payload(result)["provenance"]["scope"] == "LAUNCHER_ATTEMPT_METADATA_ONLY"
    assert payload(result)["provenance"]["validated_receipt_units"] == 0


def test_selected_failed_scheduler_remains_failed_despite_all_module_receipts(make_run, tmp_path):
    run = make_run()
    for sample in run[3]:
        complete_step(run, "small", (sample,))
    finish(run, state="FAILED", returncode=7)
    result = build_launcher_report(run[0], tmp_path / "failed-report")
    assert result["status"] == "FAILED"
    assert tracks(result)["launcher_execution"]["status"] == "FAILED"
    assert tracks(result)["small"]["status"] == "COMPLETE"


def test_dry_run_never_reads_native_receipts(make_run, tmp_path, monkeypatch):
    run = make_run(dry_run=True)
    finish(run, state="DRY_RUN_PASS")
    import hifivar.attempt_receipts
    monkeypatch.setattr(hifivar.attempt_receipts, "read_step_receipt", lambda *a, **k: pytest.fail("dry-run receipt probe"))
    result = build_launcher_report(run[0], tmp_path / "dry-report")
    assert result["status"] == tracks(result)["small"]["status"] == "NOT_RUN"


def test_disabled_tracks_do_not_probe_receipts(make_run, tmp_path, monkeypatch):
    run = make_run({"small": {"enabled": False}})
    finish(run)
    import hifivar.attempt_receipts
    monkeypatch.setattr(hifivar.attempt_receipts, "read_step_receipt", lambda *a, **k: pytest.fail("disabled receipt probe"))
    result = build_launcher_report(run[0], tmp_path / "disabled-report")
    assert tracks(result)["small"]["status"] == "DISABLED"


def test_output_mutation_is_rejected_before_report_destination_creation(make_run, tmp_path):
    run = make_run()
    output = complete_step(run, "small", ("S1",))
    finish(run)
    output.write_bytes(b"not-the-recorded-output")
    with pytest.raises(OutputValidationError, match="bytes differ"):
        build_launcher_report(run[0], tmp_path / "must-not-exist")
    assert not (tmp_path / "must-not-exist").exists()


@pytest.mark.parametrize("change", ["missing", "mutated"])
def test_failed_attempt_delivers_failed_artifact_diagnosis_without_borrowing_success(make_run, tmp_path, change):
    run = make_run()
    output = complete_step(run, "small", ("S1",))
    complete_step(run, "small", ("S2",))
    finish(run, state="FAILED", returncode=9)
    if change == "missing":
        preserved = output.with_suffix(output.suffix + ".preserved")
        output.rename(preserved)
    else:
        output.write_bytes(b"different-from-receipt")
    before = fingerprint(run[0])
    result = build_launcher_report(run[0], tmp_path / "failed-artifact-report")
    small = tracks(result)["small"]
    assert result["status"] == small["status"] == "FAILED"
    assert small["qc"]["artifact_validation"] == "FAIL_OUTPUT_VALIDATION"
    assert small["qc"]["expected_units"] == 2
    assert small["qc"]["validated_units"] == 1
    assert {row["status"] for row in small["qc"]["units"]} == {"COMPLETE", "FAILED_OUTPUT_VALIDATION"}
    assert all(item["path"] != str(output) for item in small["artifacts"])
    assert set(result["output_paths"]) == {"json", "yaml", "md", "html"}
    assert payload(result)["provenance"]["scientific_accuracy"] == "NOT_ASSESSED"
    assert fingerprint(run[0]) == before


def test_failed_attempt_does_not_relax_frozen_control_validation(make_run, tmp_path):
    run = make_run()
    complete_step(run, "small", ("S1",))
    finish(run, state="FAILED", returncode=9)
    with (run[0] / "effective.yaml").open("a", encoding="utf-8") as stream:
        stream.write("# changed control\n")
    with pytest.raises(InputValidationError, match="LAUNCHER_CONTROL_CHANGED"):
        build_launcher_report(run[0], tmp_path / "unsafe-control-report")
    assert not (tmp_path / "unsafe-control-report").exists()


def test_receipt_cannot_be_reused_as_new_attempt_completion(make_run, tmp_path, monkeypatch):
    run = make_run()
    root, attempt, _, samples = run
    complete_step(run, "small", ("S1",))
    finish(run)
    second = root / "attempts" / "attempt_b"
    write_json(second / "status.json", {"started": now(), "status": "DAG_RUNNING", "dry_run": False})
    create_attempt_context(root, second, dry_run=False)
    status = json.loads((second / "status.json").read_text(encoding="utf-8"))
    status.update(status="WORKFLOW_SUCCEEDED", returncode=0, finished=now(), input_immutability="PASS")
    write_json(second / "status.json", status)
    result = build_launcher_report(root, tmp_path / "new-attempt-report", "attempt_b")
    assert tracks(result)["small"]["qc"]["validated_units"] == 0
    assert result["status"] == "PARTIAL"
    receipt = attempt / "receipts" / receipt_name("small", ("S1",))
    target = second / "receipts" / receipt.name
    target.parent.mkdir()
    target.write_bytes(receipt.read_bytes())
    with pytest.raises(InputValidationError, match="selected native attempt"):
        build_launcher_report(root, tmp_path / "no-relabel", "attempt_b")


@pytest.mark.parametrize("track", ["tr", "phasing", "assembly", "sv.sawfish", "sv.sniffles2", "sv.pbsv",
                                    "sv.cutesv", "sv.harmonization", "assembly_sv.pav", "assembly_sv.svim_asm"])
def test_expected_native_modules_are_enumerated_explicitly(make_run, tmp_path, track):
    section, _, child = track.partition(".")
    config = {section: {"enabled": True}}
    if child:
        config[section][child] = {"enabled": True}
    run = make_run(config)
    for sample in run[3]:
        complete_step(run, track, (sample,))
    finish(run)
    result = build_launcher_report(run[0], tmp_path / "module-report")
    assert tracks(result)[track]["status"] == "COMPLETE"


def test_annotation_requires_every_enabled_backend_and_category(make_run, tmp_path):
    config = {"annotation": {"enabled": True, "annovar_enabled": True, "vep_enabled": True, "overlap_enabled": True}}
    run = make_run(config, annotation_rows=[("S1", "small"), ("S1", "sv"), ("S2", "tr")])
    for sample, category in (("S1", "small"), ("S1", "sv"), ("S2", "tr")):
        for backend in ("annovar", "vep"):
            complete_step(run, "annotation." + backend, (sample,), unit_id=category)
        if category != "small":
            complete_step(run, "annotation.overlap", (sample,), unit_id=category)
    finish(run)
    result = build_launcher_report(run[0], tmp_path / "annotation-report")
    annotation = tracks(result)["annotation"]
    assert annotation["status"] == "COMPLETE"
    assert annotation["qc"]["expected_units"] == annotation["qc"]["validated_units"] == 8
    assert {unit["unit_id"] for unit in annotation["qc"]["units"]} == {"small", "sv", "tr"}
    assert "annotation.inputs.tsv" in payload(result)["provenance"]["control_sha256_verified"]


def test_annotation_partial_backend_does_not_invent_success(make_run, tmp_path):
    config = {"annotation": {"enabled": True, "annovar_enabled": True, "vep_enabled": True}}
    run = make_run(config, annotation_rows=[("S1", "sv")])
    complete_step(run, "annotation.annovar", ("S1",), unit_id="sv",
                  metadata={"tool_version_sources": {"synthetic-python-process": "USER_DECLARED_NOT_NATIVE_PROBED"}})
    finish(run, state="FAILED", returncode=1)
    result = build_launcher_report(run[0], tmp_path / "annotation-partial")
    assert tracks(result)["annotation"]["status"] == "PARTIAL"
    assert result["status"] == "FAILED"
    units = tracks(result)["annotation"]["qc"]["units"]
    assert [(unit["track"], unit["status"]) for unit in units] == [
        ("annotation.annovar", "COMPLETE"), ("annotation.vep", "MISSING_RECEIPT")]
    assert payload(result)["tools"][-1]["validation_status"] == "USER_DECLARED_NOT_NATIVE_PROBED"


def test_generated_annotation_plan_remains_partial_until_materialization_verified(make_run, tmp_path):
    config = {"annotation": {"enabled": True, "input_manifest": None, "vep_enabled": True},
              "small": {"enabled": True}}
    source_plan = {"schema_version": 1, "annotation": [{"sample_id": "S1", "track": "small", "origin": "GENERATED",
                                                       "producer": "small.deepvariant", "output_role": "vcf"}]}
    run = make_run(config, source_plan=source_plan)
    # A source selector and even unrelated valid consumer bytes do not prove a
    # verified, attempt-bound READY-manifest handoff.
    with pytest.raises(InputValidationError):
        complete_step(run, "annotation.vep", ("S1",), unit_id="small")
    finish(run)
    result = build_launcher_report(run[0], tmp_path / "generated-pending")
    annotation = tracks(result)["annotation"]
    assert annotation["status"] == "PARTIAL"
    assert annotation["qc"]["pending_requirements"] == ["SOURCE_PLAN_MATERIALIZATION_NOT_VERIFIED"]
    assert annotation["qc"]["source_rows"][0]["origin"] == "GENERATED"
    assert "downstream_plan.yaml" in payload(result)["provenance"]["control_sha256_verified"]


def generated_run(make_run, *, cohort=False):
    config = {"annotation": {"enabled": not cohort, "input_manifest": None, "vep_enabled": not cohort},
              "cohort": {"enabled": cohort, "input_manifest": None, "small_variants": {"enabled": cohort}},
              "small": {"enabled": True}}
    group = "cohort" if cohort else "annotation"
    source_plan = {"schema_version": 1, group: [{"sample_id": "S1", "track": "small_variants" if cohort else "small",
                                              "origin": "GENERATED", "producer": "small.deepvariant",
                                              "output_role": "gvcf" if cohort else "vcf", "reference_build": "GRCh38"}]}
    run = make_run(config, samples=("S1",), source_plan=source_plan)
    complete_step(run, "small", ("S1",), tool_versions={"deepvariant": "SYNTHETIC_TEST_VERSION"},
                  metadata={"reference_build": "GRCh38", "tool_version_sources": {"deepvariant": "TEST_DECLARATION_NOT_REAL_DEEPVARIANT"}})
    return run


@pytest.mark.parametrize("cohort", [False, True])
def test_valid_materialization_and_actual_receipts_close_generated_consumer(make_run, tmp_path, cohort):
    run = generated_run(make_run, cohort=cohort)
    materialize_fixture(run)
    if cohort:
        complete_step(run, "cohort.small_variants", ("S1",))
    else:
        complete_step(run, "annotation.vep", ("S1",), unit_id="small")
    finish(run)
    before = fingerprint(run[0])
    result = build_launcher_report(run[0], tmp_path / "materialized-report")
    assert tracks(result)["cohort" if cohort else "annotation"]["status"] == "COMPLETE"
    assert result["status"] == "COMPLETE"
    assert payload(result)["provenance"]["scientific_accuracy"] == "NOT_ASSESSED"
    assert fingerprint(run[0]) == before


@pytest.mark.parametrize("cohort", [False, True])
def test_failed_generated_attempt_preserves_diagnosis_when_producer_source_is_missing(make_run, tmp_path, cohort):
    run = generated_run(make_run, cohort=cohort)
    materialize_fixture(run)
    if cohort:
        complete_step(run, "cohort.small_variants", ("S1",))
    else:
        complete_step(run, "annotation.vep", ("S1",), unit_id="small")
    producer = read_step_receipt(run[0], run[1].name, "small", ("S1",))
    source = Path(producer["outputs"]["gvcf" if cohort else "vcf"]["path"])
    source.rename(source.with_suffix(source.suffix + ".preserved"))
    finish(run, state="FAILED", returncode=9)
    before = fingerprint(run[0])
    result = build_launcher_report(run[0], tmp_path / "failed-generated-report")
    observed = tracks(result)
    downstream = observed["cohort" if cohort else "annotation"]
    assert result["status"] == observed["small"]["status"] == downstream["status"] == "FAILED"
    assert downstream["qc"]["pending_requirements"] == ["SOURCE_ARTIFACT_VALIDATION_FAILED"]
    assert downstream["qc"]["validated_units"] == 0
    assert downstream["qc"]["source_rows"] == []
    assert all(Path(path).is_file() for path in result["output_paths"].values())
    assert fingerprint(run[0]) == before


@pytest.mark.parametrize("change", ["attempt_id", "context_sha256", "plan_sha256"])
def test_materialization_identity_mismatch_fails_before_output(make_run, tmp_path, change):
    run = generated_run(make_run)
    owner, marker = materialize_fixture(run)
    marker[change] = "wrong"
    write_json(owner / "materialization.json", marker)
    finish(run)
    with pytest.raises(InputValidationError, match="bound to"):
        build_launcher_report(run[0], tmp_path / "no-wrong-handoff")
    assert not (tmp_path / "no-wrong-handoff").exists()


def test_materialized_control_mutation_is_rejected(make_run, tmp_path):
    run = generated_run(make_run)
    owner, _ = materialize_fixture(run)
    with (owner / "annotation.inputs.tsv").open("a", encoding="utf-8") as stream:
        stream.write("changed\n")
    finish(run)
    with pytest.raises(InputValidationError, match="control identity"):
        build_launcher_report(run[0], tmp_path / "no-control-change")


@pytest.mark.parametrize("change", ["missing", "sha256", "producer_attempt_id", "source_version"])
def test_materialization_source_receipt_linkage_is_not_assumed(make_run, tmp_path, change):
    run = generated_run(make_run)
    owner, marker = materialize_fixture(run)
    if change == "missing":
        marker["source_receipts"] = []
    else:
        marker["source_receipts"][0][change] = "wrong"
    write_json(owner / "materialization.json", marker)
    finish(run)
    with pytest.raises(InputValidationError):
        build_launcher_report(run[0], tmp_path / "no-invalid-link")


def test_hash_valid_but_wrong_materialized_source_is_rejected(make_run, tmp_path):
    run = generated_run(make_run)
    materialize_fixture(run, wrong_source=True)
    finish(run)
    with pytest.raises(InputValidationError, match="annotation source differs"):
        build_launcher_report(run[0], tmp_path / "no-wrong-source")


@pytest.mark.parametrize("value", [["sv"], {"unexpected": "value"}, None])
def test_wrong_typed_planned_category_is_explicit_validation_error(make_run, tmp_path, value):
    config = {"annotation": {"enabled": True, "input_manifest": None, "vep_enabled": True}}
    run = make_run(config, source_plan={"schema_version": 1, "annotation": [{"sample_id": "S1", "track": value,
                                                                          "origin": "GENERATED"}]})
    finish(run)
    with pytest.raises(InputValidationError):
        build_launcher_report(run[0], tmp_path / "no-invalid-plan-category")


def test_cohort_complete_requires_all_selected_tracks_and_declared_sample_order(make_run, tmp_path):
    config = {"cohort": {"enabled": True, "small_variants": {"enabled": True},
                         "sv": {"enabled": True}, "tr": {"enabled": True}}}
    run = make_run(config)
    for track in ("small_variants", "sv", "tr"):
        complete_step(run, "cohort." + track, run[3])
    finish(run)
    result = build_launcher_report(run[0], tmp_path / "cohort-report")
    cohort = tracks(result)["cohort"]
    assert cohort["status"] == "COMPLETE" and cohort["qc"]["validated_units"] == 3
    assert all(unit["samples"] == ["S2", "S1"] for unit in cohort["qc"]["units"])
    assert "cohort.inputs.tsv" in payload(result)["provenance"]["control_sha256_verified"]


def test_cohort_missing_track_is_partial_not_all_vcf_unification(make_run, tmp_path):
    config = {"cohort": {"enabled": True, "sv": {"enabled": True}, "tr": {"enabled": True}}}
    run = make_run(config)
    complete_step(run, "cohort.sv", run[3])
    finish(run)
    result = build_launcher_report(run[0], tmp_path / "cohort-partial")
    assert tracks(result)["cohort"]["status"] == "PARTIAL"
    assert [unit["track"] for unit in tracks(result)["cohort"]["qc"]["units"]] == ["cohort.sv", "cohort.tr"]


def test_additional_frozen_manifest_change_is_rejected_before_receipt_lookup(make_run, tmp_path):
    run = make_run({"annotation": {"enabled": True, "vep_enabled": True}}, annotation_rows=[("S1", "sv")])
    finish(run)
    with (run[0] / "annotation.inputs.tsv").open("a", encoding="utf-8") as stream:
        stream.write("changed\n")
    with pytest.raises(InputValidationError, match="LAUNCHER_CONTROL_CHANGED"):
        build_launcher_report(run[0], tmp_path / "no-changed-control")
    assert not (tmp_path / "no-changed-control").exists()


def test_unfrozen_annotation_source_manifest_is_not_trusted(make_run, tmp_path):
    config = {"annotation": {"enabled": True, "vep_enabled": True, "input_manifest": "/external/unfrozen.tsv"}}
    run = make_run(config)
    finish(run)
    with pytest.raises(InputValidationError, match="owned frozen"):
        build_launcher_report(run[0], tmp_path / "no-external-manifest")


def test_context_identity_change_is_rejected_even_without_receipts(make_run, tmp_path):
    run = make_run()
    finish(run)
    context = run[1] / "execution_context.json"
    record = json.loads(context.read_text(encoding="utf-8"))
    record["attempt_id"] = "different_attempt"
    write_json(context, record)
    with pytest.raises(InputValidationError, match="context"):
        build_launcher_report(run[0], tmp_path / "no-wrong-context")


@pytest.mark.parametrize("name", ["../escape.tsv", "/absolute.tsv", "a\\b.tsv", "C:relative.tsv"])
def test_unsafe_extra_control_names_fail_closed(make_run, tmp_path, name):
    run = make_run()
    finish(run)
    path = run[0] / "launch_manifest.json"
    value = json.loads(path.read_text(encoding="utf-8"))
    value["control_sha256"][name] = "a" * 64
    write_json(path, value)
    with pytest.raises(InputValidationError, match="Unsafe launcher control"):
        build_launcher_report(run[0], tmp_path / "no-escape")


def test_ready_harmonization_is_not_hidden_by_disabled_read_calling(make_run, tmp_path):
    run = make_run({'sv': {'enabled': False, 'harmonization': {'enabled': True}}})
    for sample in run[3]:
        complete_step(run, 'sv.harmonization', (sample,))
    finish(run)
    result = build_launcher_report(run[0], tmp_path / 'ready-sv-report')
    assert tracks(result)['sv.harmonization']['status'] == 'COMPLETE'
    assert tracks(result)['sv.sawfish']['status'] == 'DISABLED'
