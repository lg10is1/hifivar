"""Adapt an explicitly selected launcher attempt into a local final report.

Legacy launcher metadata remains conservative. New attempts can additionally
prove bounded step execution through explicitly selected, validated receipts
and streamed output hashes. Neither mode infers success from result filenames
or turns engineering completion into benchmark accuracy.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import re
from collections import Counter
from datetime import datetime
from pathlib import Path

import yaml

from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.report import FinalRunReport, FinalStatus, ReportArtifact, ToolRecord, TrackReport


_ATTEMPT_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}")
_ACTIVE_STATES = {
    "PREFLIGHT_RUNNING", "ALIGNMENT_PREPARATION", "DAG_RUNNING",
    "PRODUCER_DAG_RUNNING", "DOWNSTREAM_PREPARATION", "DOWNSTREAM_DAG_RUNNING",
    "REVIEW_PREPARATION", "REVIEW_DAG_RUNNING", "BENCHMARK_PREPARATION", "BENCHMARK_DAG_RUNNING",
}
_PLANNED_STATES = {"DRY_RUN_PASS", "RAW_INPUT_PLAN_ONLY", "GENERATED_SOURCE_PLAN_ONLY", "REVIEW_SOURCE_PLAN_ONLY", "BENCHMARK_SOURCE_PLAN_ONLY"}
_KNOWN_STATES = _ACTIVE_STATES | _PLANNED_STATES | {"WORKFLOW_SUCCEEDED", "FAILED", "INTERRUPTED_LOCK_RETAINED"}
_HASH = re.compile(r"[0-9a-f]{64}\Z")
_CONTROL_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_./-]{0,511}\Z")
_SCIENTIFIC_BOUNDARY = (
    "Execution status is separate from scientific accuracy. Accuracy is NOT_ASSESSED; "
    "this adapter does not run a matched-truth benchmark, repeat calling or a full raw-read audit, "
    "or substitute caller support/manual review for truth."
)


def build_launcher_report(
    run_dir: Path,
    output_dir: Path,
    attempt_id: str | None = None,
) -> dict[str, object]:
    """Write fresh JSON/YAML/Markdown/HTML reports from launcher metadata.

    ``run_dir`` must be an owned schema-2 launcher directory. ``output_dir`` must
    not exist and must be outside it. Redirected paths and ``..`` are refused.
    With multiple attempts an explicit immediate-child ``attempt_id`` is
    mandatory: folder names and mtimes never choose which attempt is reported.
    No-attempt/missing-status cases are reported as NOT_RUN, not success.

    Every declared frozen control is checked. Legacy attempts without a linked
    execution context remain NOT_RUN/PARTIAL per module. For new contexts, only
    explicitly expected step receipts are read; their declared output hashes
    are streamed and verified. For an already FAILED scheduler attempt, an
    unavailable or changed artifact is delivered as a failed validation unit,
    not as module completion. Invalid controls/receipt identity still fail
    closed. No result-directory discovery or reuse inference
    occurs. Receipt-bound selected source roles may be revalidated on readback;
    this is not a full raw-input audit or proof of scientific accuracy. Private
    paths remain local; the adapter never publishes any report.
    """
    root = _safe_path(run_dir, "Run directory", must_exist=True)
    destination = _safe_path(output_dir, "Report directory", must_exist=False)
    if destination.is_relative_to(root) or root.is_relative_to(destination):
        raise OutputValidationError("Report output must be independent of the original run tree.")
    if destination.exists():
        raise OutputValidationError(f"Report output already exists: '{destination}'.")
    manifest_path, config_path, sheet_path = (root / name for name in ("launch_manifest.json", "effective.yaml", "samples.tsv"))
    manifest = _mapping_file(manifest_path, yaml_input=False)
    identity = manifest.get("identity")
    if not isinstance(identity, dict) or identity.get("schema") != 2 or identity.get("kind") != "hifivar_two_stage_launcher":
        raise InputValidationError("Unsupported or missing launcher manifest identity/schema.")
    if identity.get("run_dir") != str(root):
        raise InputValidationError("LAUNCHER_RUN_IDENTITY_MISMATCH: manifest does not identify this run directory.")
    config = _mapping_file(config_path, yaml_input=True)
    if identity.get("config") != config:
        raise InputValidationError("LAUNCHER_CONFIG_MISMATCH: effective config differs from the manifest.")
    control_records = _verify_controls(root, manifest.get("control_sha256"))
    try:
        sheet = sheet_path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as error:
        raise InputValidationError(f"Unable to read launcher sample sheet: {error}") from error
    if identity.get("sample_sheet") != sheet:
        raise InputValidationError("LAUNCHER_SAMPLE_SHEET_MISMATCH: sample-sheet bytes differ from identity.")
    samples = _samples(sheet)
    attempt, status = _select_attempt(root, attempt_id)
    execution, execution_status, execution_warnings = _execution_status(status)
    warnings = [_SCIENTIFIC_BOUNDARY, *execution_warnings]
    if attempt is None:
        warnings.append("No launcher attempt is recorded; enabled modules are NOT_RUN.")
    context = _linked_context(root, attempt, status, identity, manifest_path)
    can_read_receipts = context is not None and status is not None and status.get("dry_run") is False
    artifacts = [
        ReportArtifact(manifest_path, "launcher-manifest", "json"),
        ReportArtifact(config_path, "effective-config", "yaml"),
        ReportArtifact(sheet_path, "declared-sample-sheet", "tsv"),
    ]
    if attempt is not None and status is not None:
        artifacts.append(ReportArtifact(attempt / "status.json", "selected-attempt-status", "json"))
    sample_ids = tuple(str(sample["sample_id"]) for sample in samples)
    tracks = [TrackReport(
        "launcher_execution", "launcher", "sample", execution, True,
        tuple(artifacts), sample_ids,
        qc={"execution_status": execution_status, "accuracy_status": "NOT_ASSESSED"},
        message=f"Selected attempt execution: {execution_status}. {_SCIENTIFIC_BOUNDARY}",
    )]
    validated_receipts: list[dict] = []
    for name, phase, category, enabled in _configured_tracks(config):
        module_status = FinalStatus.DISABLED
        if enabled:
            module_status = FinalStatus.PARTIAL if execution in {FinalStatus.COMPLETE, FinalStatus.PARTIAL, FinalStatus.FAILED} and status and not status.get("dry_run", False) else FinalStatus.NOT_RUN
        module_artifacts: tuple[ReportArtifact, ...] = ()
        qc = {"artifact_validation": "NOT_RUN", "attempt_linked_result_record": "NOT_AVAILABLE"} if enabled else {}
        message = "Enabled, but no attempt-linked per-track receipt is available." if enabled else "Not enabled in the frozen effective config."
        if enabled and can_read_receipts:
            failed_source_artifact = False
            try:
                units, pending, source_rows = _expected_units(root, attempt, config, control_records, name, sample_ids)
            except OutputValidationError:
                if execution is not FinalStatus.FAILED:
                    raise
                # Only a terminal failed attempt may produce this diagnostic.
                # Never turn unavailable producer bytes into a READY handoff.
                units, source_rows = [], []
                pending = ["SOURCE_ARTIFACT_VALIDATION_FAILED"]
                failed_source_artifact = True
            module_artifacts, unit_records, actual_receipts = _read_units(
                root, attempt.name, units, diagnose_failed_artifacts=execution is FinalStatus.FAILED,
            )
            validated_receipts.extend(actual_receipts)
            complete = bool(units) and not pending and all(item["status"] == "COMPLETE" for item in unit_records)
            failed_artifacts = failed_source_artifact or any(item["status"] == "FAILED_OUTPUT_VALIDATION" for item in unit_records)
            module_status = FinalStatus.FAILED if failed_artifacts else FinalStatus.COMPLETE if complete else FinalStatus.PARTIAL
            qc = {
                "artifact_validation": "FAIL_OUTPUT_VALIDATION" if failed_artifacts else "PASS_OUTPUT_SHA256" if actual_receipts else "NOT_RUN",
                "attempt_linked_result_record": "ALL_EXPECTED_VERIFIED" if complete else "INCOMPLETE",
                "expected_units": len(units), "validated_units": len(actual_receipts),
                "units": unit_records, "pending_requirements": pending,
                "source_rows": source_rows, "accuracy_status": "NOT_ASSESSED",
                "source_input_revalidation": ("EXPLICIT_BND_BUNDLE_AND_SELECTED_RECEIPT_BINDINGS_NOT_RAW_DATA_AUDIT"
                    if any(row.get("bnd_bundle_manifest") for row in source_rows) else "NOT_RUN_BY_REPORT"),
            }
            message = (
                "The selected failed attempt has unavailable or changed recorded artifacts; no completion is inferred for failed units."
                if failed_artifacts else
                "All explicit units have validated attempt-linked receipts and unchanged output bytes; alignment VALIDATION_ONLY means verified reuse, not a new alignment. Scientific accuracy is NOT_ASSESSED."
                if complete else "Expected unit receipts or source materialization are incomplete; no completion was inferred from existing result files."
            )
            if name == "alignment":
                qc["alignment_actions"] = [item["metadata"]["action"] for item in actual_receipts]
                qc["source_input_revalidation"] = "SELECTED_BAM_INDEX_SHA256_NOT_ORIGINAL_FASTQ_REHASH"
            if name == "benchmark":
                qc["evidence_kind"] = "CONTEXT_SPECIFIC_COMPARATOR_OUTPUTS_NOT_UNIVERSAL_ACCURACY"
                qc["truth_scope"] = source_rows
                qc["comparison_statuses"] = [item["metadata"].get("benchmark_statuses", {}) for item in actual_receipts]
                qc["accuracy_status"] = "EXPLICIT_BENCHMARK_ARTIFACTS_AVAILABLE_NOT_UNIVERSAL_ACCURACY" if complete else "NOT_ASSESSED"
                message += " Benchmark completion describes comparator execution; inspect comparison_statuses and metric artifacts for PASS/PARTIAL and their exact truth/region policy."
            if name == "review":
                qc["human_review_status"] = "NOT_REVIEWED"
                qc["evidence_kind"] = "TECHNICAL_IGV_EXECUTION_NOT_HUMAN_OR_TRUTH_VALIDATION"
                qc["screenshot_artifacts"] = sum(role.startswith("screenshot_") for item in actual_receipts for role in item["outputs"])
                qc["source_input_revalidation"] = "BOUND_RECORDS_VCF_INDEX_AND_ACTUAL_BAM_PROBES_NOT_FULL_RAW_READ_AUDIT"
                message += " IGV execution is technical evidence only; human review remains NOT_REVIEWED."
        tracks.append(TrackReport(
            name, phase, category, module_status, enabled, artifacts=module_artifacts,
            samples=((config["benchmark"]["sample_id"],) if name=="benchmark" else sample_ids) if enabled else (), qc=qc, message=message,
        ))
    if config.get("assembly_sv", {}).get("input_manifest"):
        reuse_status = FinalStatus.NOT_RUN
        reuse_artifacts = ()
        reuse_qc = {"execution_kind": "VALIDATION_ONLY", "new_hifiasm_execution": False,
                    "accuracy_status": "NOT_ASSESSED", "input_revalidation": "NOT_RUN"}
        if can_read_receipts:
            from hifivar.assembly_reuse import HANDOFF_NAME, read_ready_handoff
            try:
                preparation, ready_items = read_ready_handoff(root, attempt, context.get("ready_assembly_handoff"), verify_inputs=True)
                reuse_status = FinalStatus.COMPLETE
                reuse_artifacts = (ReportArtifact(attempt / HANDOFF_NAME, "ready-assembly-validation-not-native-calling", "json"),)
                reuse_qc.update(input_revalidation="FULL_SELECTED_FASTA_MANIFEST_SHA256_AND_STRUCTURE",
                                declared_sample_order=preparation["declared_sample_order"],
                                assemblies=[item.to_dict() for item in ready_items])
            except (InputValidationError, OutputValidationError) as error:
                if execution is not FinalStatus.FAILED:
                    raise
                reuse_status = FinalStatus.FAILED
                reuse_qc.update(input_revalidation="FAILED", diagnostic=str(error))
        tracks.append(TrackReport("assembly_input_reuse", "8", "assembly", reuse_status, True,
                                  artifacts=reuse_artifacts, samples=sample_ids, qc=reuse_qc,
                                  message="Validation-only reuse of user-declared assembly inputs; not a new hifiasm execution or phasing accuracy result."))
    if validated_receipts:
        warnings.append(
            "Only explicitly selected this-attempt receipts and their output byte hashes were verified. "
            "Receipt-bound source/control checks retain their recorded scope; whole raw-input identity, biological accuracy and human decisions are not re-established by reporting."
        )
    else:
        warnings.append(
            "No validated per-track receipt is available for this selected attempt. "
            "Enabled module completion and output integrity are NOT independently verified; "
            "old results are not evidence for this selected attempt."
        )
    runtime = identity.get("runtime") if isinstance(identity.get("runtime"), dict) else {}
    reference = dict(config.get("reference", {}))
    reference["checksum_verification"] = "NOT_RUN_METADATA_ONLY"
    inputs = identity.get("inputs", [])
    if isinstance(inputs, list):
        for item in inputs:
            if isinstance(item, dict) and item.get("path") == reference.get("fasta"):
                reference["recorded_sha256"] = item.get("sha256")
    tools: tuple[ToolRecord, ...] = ()
    if status and isinstance(status.get("snakemake_version"), str):
        tools = (ToolRecord("snakemake", status["snakemake_version"], validation_status="RECORDED_NOT_REPROBED"),)
    tools += _receipt_tools(validated_receipts)
    report = FinalRunReport(
        run_id=root.name, git_sha=None, reference=reference, samples=tuple(samples), tracks=tuple(tracks),
        config=config, tools=tools, warnings=tuple(warnings),
        project_version=str(runtime.get("hifivar_version", "NOT_RECORDED")),
        provenance={
            "adapter": "hifivar.workflow_report", "scope": "LAUNCHER_ATTEMPT_RECEIPTS_AND_OUTPUT_SHA256" if validated_receipts else "LAUNCHER_ATTEMPT_METADATA_ONLY",
            "run_dir": str(root), "attempt_id": attempt.name if attempt else None,
            "execution_status": execution_status, "selected_attempt_status": status,
            "runtime_identity": runtime, "control_sha256_verified": control_records,
            "source_input_checksum_verification": "NOT_RUN_METADATA_ONLY",
            "scientific_accuracy": "NOT_ASSESSED", "publication": "LOCAL_ONLY_NOT_PUBLISHED",
            "git_sha": "NOT_RECORDED_BY_LAUNCHER", "original_run_modified": False,
            "validated_receipt_units": len(validated_receipts),
            "output_checksum_verification": "PASS_EXPLICIT_RECEIPT_ARTIFACTS_ONLY" if validated_receipts else "NOT_RUN",
            "receipt_reuse": "NOT_IMPLEMENTED_NO_OLD_ATTEMPT_INFERENCE",
        },
    )
    destination.mkdir(parents=True, exist_ok=False)
    paths = {}
    for suffix, writer in (("json", report.write_json), ("yaml", report.write_yaml), ("md", report.write_markdown), ("html", report.write_html)):
        target = destination / f"final_report.{suffix}"
        writer(target)
        paths[suffix] = str(target)
    return {
        "status": report.status.value, "execution_status": execution_status,
        "scientific_accuracy": "NOT_ASSESSED", "attempt_id": attempt.name if attempt else None,
        "run_dir": str(root), "output_paths": paths,
        "output_publication": "LOCAL_ONLY_NOT_PUBLISHED",
    }


def _verify_controls(root: Path, value: object) -> dict[str, str]:
    if not isinstance(value, dict) or not {"effective.yaml", "samples.tsv"}.issubset(value):
        raise InputValidationError("Launcher manifest lacks required control_sha256 records.")
    result = {}
    for name, expected in value.items():
        if not isinstance(name, str) or not _CONTROL_NAME.fullmatch(name) or ".." in Path(name).parts:
            raise InputValidationError("Unsafe launcher control path.")
        if not isinstance(expected, str) or not _HASH.fullmatch(expected):
            raise InputValidationError("Malformed launcher control SHA256.")
        path = root / name
        digest = _metadata_digest(path)
        if digest != expected:
            raise InputValidationError(f"LAUNCHER_CONTROL_CHANGED: '{name}' checksum differs.")
        result[name] = digest
    return result


def _metadata_digest(path: Path) -> str:
    _safe_file(path)
    before = path.stat()
    data = path.read_bytes()
    after = path.stat()
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (
        after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns
    ):
        raise InputValidationError("Launcher metadata changed while being read.")
    return hashlib.sha256(data).hexdigest()


def _linked_context(root: Path, attempt: Path | None, status: dict | None,
                    identity: dict, manifest_path: Path) -> dict | None:
    if attempt is None:
        return None
    path = attempt / "execution_context.json"
    if not path.exists() and not path.is_symlink():
        return None
    context = _mapping_file(path, yaml_input=False)
    if (type(context.get("schema_version")) is not int or context["schema_version"] != 1
            or context.get("kind") != "hifivar_attempt_execution_context"
            or context.get("run_dir") != str(root) or context.get("attempt_id") != attempt.name
            or context.get("runtime") != identity.get("runtime")
            or context.get("manifest_sha256") != _metadata_digest(manifest_path)
            or not isinstance(context.get("dry_run"), bool)):
        raise InputValidationError("Execution context does not match the selected launcher attempt.")
    if status is not None and context["dry_run"] != status.get("dry_run"):
        raise InputValidationError("Execution context and attempt dry_run states differ.")
    _timestamp(context.get("created"), "context created")
    return context


def _owned_control(root: Path, value: object, controls: dict[str, str]) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise InputValidationError("Downstream manifest requires an explicit frozen control path.")
    path = Path(value)
    if ".." in path.parts:
        raise InputValidationError("Downstream control must not contain '..'.")
    if not path.is_absolute():
        path = root / path
    if not path.is_relative_to(root):
        raise InputValidationError("Downstream manifest must be an owned frozen launcher control.")
    name = path.relative_to(root).as_posix()
    if name not in controls:
        raise InputValidationError("Downstream manifest is not declared in frozen control_sha256.")
    _safe_file(path)
    return path


def _planned_source_rows(root: Path, config: dict, controls: dict[str, str],
                         group: str, sample_ids: tuple[str, ...]) -> list[dict]:
    section = config.get("workflow", {})
    if not isinstance(section, dict):
        raise InputValidationError("Frozen workflow config must be a mapping.")
    value = section.get("downstream_source_plan")
    if value is None:
        return []
    path = _owned_control(root, value, controls)
    document = _mapping_file(path, yaml_input=True)
    if (type(document.get("schema_version")) is not int or document["schema_version"] != 1
            or set(document).difference({"schema_version", "annotation", "cohort"})):
        raise InputValidationError("Unsupported frozen downstream source-plan schema.")
    rows = document.get(group, [])
    if not isinstance(rows, list) or len(rows) > 4096:
        raise InputValidationError("Frozen downstream source-plan rows must be bounded lists.")
    result, seen = [], set()
    allowed = {"small", "sv", "tr"} if group == "annotation" else {"small_variants", "sv", "tr"}
    for row in rows:
        if (not isinstance(row, dict) or not isinstance(row.get("sample_id"), str) or row["sample_id"] not in sample_ids
                or not isinstance(row.get("track"), str) or row["track"] not in allowed
                or not isinstance(row.get("origin"), str) or row["origin"] not in {"GENERATED", "EXTERNAL"}):
            raise InputValidationError("Frozen downstream source-plan row has invalid sample/category/origin.")
        key = row["sample_id"], row["track"]
        if key in seen:
            raise InputValidationError("Frozen source plan contains duplicate sample/category rows.")
        seen.add(key)
        result.append({"sample_id": key[0], "track": key[1], "origin": row["origin"],
                       "producer": row.get("producer"), "state": row.get("state", "CALLED"),
                       "source_plan_sha256": controls[path.relative_to(root).as_posix()]})
    return result


def _materialized_controls(root: Path, attempt: Path, config: dict,
                            controls: dict[str, str]) -> tuple[Path, dict, dict[str, str]] | None:
    """Verify the fixed attempt-owned handoff, never discover nearby files."""
    from hifivar.attempt_receipts import read_step_receipt, receipt_name
    owner = attempt / "downstream"
    marker = owner / "materialization.json"
    if not marker.exists() and not marker.is_symlink():
        return None
    record = _mapping_file(marker, yaml_input=False)
    plan_path = _owned_control(root, config["workflow"]["downstream_source_plan"], controls)
    if (type(record.get("schema_version")) is not int or record["schema_version"] != 1
            or record.get("kind") != "hifivar_downstream_materialization"
            or record.get("attempt_id") != attempt.name
            or record.get("context_sha256") != _metadata_digest(attempt / "execution_context.json")
            or record.get("plan_sha256") != controls[plan_path.relative_to(root).as_posix()]):
        raise InputValidationError("Downstream materialization is not bound to the selected context and source plan.")
    declared = record.get("controls")
    required = {"effective.yaml"} | {
        group + ".inputs.tsv" for group in ("annotation", "cohort") if config.get(group, {}).get("enabled") is True
    }
    if not isinstance(declared, dict) or set(declared) != required:
        raise InputValidationError("Materialization lacks the exact enabled downstream control set.")
    verified = {}
    for name, item in declared.items():
        target = owner / name
        if (not isinstance(item, dict) or item.get("path") != str(target)
                or type(item.get("bytes")) is not int or item["bytes"] <= 0
                or not isinstance(item.get("sha256"), str) or not _HASH.fullmatch(item["sha256"])
                or _metadata_digest(target) != item["sha256"] or target.stat().st_size != item["bytes"]):
            raise InputValidationError("Materialized downstream control identity changed.")
        verified[name] = item["sha256"]
    effective = _mapping_file(owner / "effective.yaml", yaml_input=True)
    for group in ("annotation", "cohort"):
        expected = dict(config.get(group, {}))
        if expected.get("enabled"):
            expected["input_manifest"] = str(owner / (group + ".inputs.tsv"))
        if effective.get(group, {}) != expected:
            raise InputValidationError("Materialized downstream settings differ from the frozen source intent.")
    if effective.get("reference", {}) != config.get("reference", {}):
        raise InputValidationError("Materialized reference differs from the frozen launcher reference.")
    if effective.get("paths", {}) != config.get("paths", {}) or effective.get("runtime", {}) != config.get("runtime", {}):
        raise InputValidationError("Materialized output/resource settings differ from the frozen launcher.")
    for group in ("alignment", "small", "sv", "tr", "phasing", "assembly", "assembly_sv"):
        section = effective.get(group, {})
        if not isinstance(section, dict) or section.get("enabled", False) is not False:
            raise InputValidationError("Materialized downstream DAG must not enable producer reruns.")
    integration = effective.get("sv", {}).get("harmonization", {})
    if not isinstance(integration, dict) or integration.get("enabled", False) is not False:
        raise InputValidationError("Materialized downstream DAG must not enable independent harmonization producer reruns.")
    from hifivar.breakends import BndPolicy
    from hifivar.sv_input_identity import SourceIdentityPolicy
    if (BndPolicy.from_config(integration.get("bnd_representation", {})) is not None
            or SourceIdentityPolicy.from_config(integration.get("input_identity", {})) is not None):
        raise InputValidationError("Materialized downstream DAG must disable producer-only representation policies.")
    if not isinstance(effective.get("workflow", {}), dict) or effective.get("workflow", {}).get("downstream_source_plan") is not None:
        raise InputValidationError("Materialized downstream configuration cannot recursively launch a source plan.")
    plan = _mapping_file(plan_path, yaml_input=True)
    expected_links = Counter()
    for group in ("annotation", "cohort"):
        rows = plan.get(group, [])
        if not isinstance(rows, list) or len(rows) > 4096:
            raise InputValidationError("Invalid materialized source-plan rows.")
        for row in rows:
            if not isinstance(row, dict):
                raise InputValidationError("Invalid materialized source-plan row.")
            if row.get("origin") == "GENERATED":
                producer = row.get("producer")
                if not isinstance(producer, str):
                    raise InputValidationError("Invalid generated producer in the frozen source plan.")
                track = {"small.deepvariant": "small", "phasing.hiphase": "phasing", "tr.trgt": "tr"}.get(producer, producer)
                if not isinstance(track, str) or not isinstance(row.get("sample_id"), str) or row.get("output_role") not in {"vcf", "gvcf", "bnd_vcf"}:
                    raise InputValidationError("Invalid generated producer linkage in the frozen source plan.")
                if row["output_role"] == "bnd_vcf" and (group != "annotation" or row.get("track") != "sv" or track != "sv.harmonization"):
                    raise InputValidationError("Generated BND role is only an explicit harmonization SV annotation source.")
                expected_links[(track, (row["sample_id"],), row["output_role"])] += 1
    links = record.get("source_receipts")
    if not isinstance(links, list) or len(links) > 8192:
        raise InputValidationError("Materialization source receipts must be a bounded list.")
    observed_links = Counter()
    producer_records = {}
    producer_sources = {}
    for link in links:
        if not isinstance(link, dict) or not isinstance(link.get("track"), str) or not isinstance(link.get("samples"), list):
            raise InputValidationError("Invalid materialization source receipt link.")
        samples = link["samples"]
        if not samples or any(not isinstance(sample, str) for sample in samples):
            raise InputValidationError("Invalid materialization source receipt samples.")
        if not isinstance(link.get("output_role"), str):
            raise InputValidationError("Invalid materialization producer output role.")
        key = link["track"], tuple(samples), link["output_role"]
        if key not in expected_links or link.get("producer_attempt_id") != attempt.name:
            raise InputValidationError("Materialization links an unexpected or reused producer.")
        receipt_path = attempt / "receipts" / receipt_name(link["track"], samples)
        if (link.get("path") != str(receipt_path) or link.get("sha256") != _metadata_digest(receipt_path)
                or type(link.get("bytes")) is not int or link["bytes"] != receipt_path.stat().st_size):
            raise InputValidationError("Materialization source receipt hash/path differs.")
        producer_key = key[:2]
        if producer_key not in producer_records:
            producer_records[producer_key] = read_step_receipt(root, attempt.name, link["track"], samples, verify_outputs=True)
        source = producer_records[producer_key]
        if source is None or link.get("output_role") not in source["outputs"]:
            raise InputValidationError("Materialized source has no accepted selected-attempt producer output.")
        tool = {"small": "deepvariant", "tr": "trgt", "phasing": "hiphase", "sv.harmonization": "jasmine"}.get(link["track"], link["track"].split(".")[-1])
        version_keys = [name for name in source["tool_versions"] if name.casefold().replace("-", "_") == tool]
        index_role = {"gvcf": "gvcf_index", "bnd_vcf": "bnd_index"}.get(key[2], "index")
        if link.get("index_role", index_role) != index_role:
            raise InputValidationError("Materialized producer index selector differs from its output role.")
        bundle_manifest = None
        if key[2] == "bnd_vcf":
            if link.get("bnd_bundle_manifest_role") != "bnd_manifest":
                raise InputValidationError("Materialized BND source lacks its explicit bundle selector.")
            from hifivar.downstream_bnd import generated_bnd_roles
            from hifivar.reference import ReferenceGenome
            ref = config["reference"]
            genome = ReferenceGenome.from_fasta(Path(ref["fasta"]), build=ref["build"])
            bound_roles = generated_bnd_roles(source, genome, samples[0])
            bundle_manifest = bound_roles["bnd_bundle_manifest"]["path"]
        elif link.get("bnd_bundle_manifest_role") is not None:
            raise InputValidationError("Native generated source cannot attach a BND bundle selector.")
        if (len(version_keys) != 1 or index_role not in source["outputs"]
                or link.get("source_tool") != tool
                or link.get("source_version") != source["tool_versions"][version_keys[0]]):
            raise InputValidationError("Materialized producer version or index role differs from its receipt.")
        producer_sources[key] = {
            "path": source["outputs"][key[2]]["path"], "index_path": source["outputs"][index_role]["path"],
            "source_tool": tool, "source_version": link["source_version"], "metadata": source.get("metadata", {}),
            "bnd_bundle_manifest": bundle_manifest,
        }
        observed_links[key] += 1
    if observed_links != expected_links:
        raise InputValidationError("Materialization does not account for every explicit generated source row.")
    _validate_materialized_rows(owner, plan, producer_sources)
    return owner, effective, verified


def _validate_materialized_rows(owner: Path, plan: dict, producers: dict) -> None:
    for group in ("annotation", "cohort"):
        expected = plan.get(group, [])
        if not expected:
            continue
        target = owner / (group + ".inputs.tsv")
        with target.open("r", encoding="utf-8", newline="") as stream:
            observed = list(csv.DictReader(stream, delimiter="\t"))
        if len(observed) != len(expected):
            raise InputValidationError("Materialized row count differs from the frozen source plan.")
        category_field = "variant_category" if group == "annotation" else "track"
        actual = {}
        for row in observed:
            if None in row or any(value is None for value in row.values()):
                raise InputValidationError("Malformed materialized manifest row.")
            key = row.get("sample"), row.get(category_field)
            if key in actual:
                raise InputValidationError("Duplicate materialized source row.")
            actual[key] = row
        for row in expected:
            entry = actual.get((row.get("sample_id"), row.get("track")))
            if entry is None:
                raise InputValidationError("Materialized sample/category differs from the source plan.")
            source = row
            if row.get("origin") == "GENERATED":
                producer = row.get("producer")
                track = {"small.deepvariant": "small", "phasing.hiphase": "phasing", "tr.trgt": "tr"}.get(producer, producer)
                source = producers[(track, (row["sample_id"],), row["output_role"])]
                metadata = source["metadata"]
                if (not isinstance(metadata, dict) or metadata.get("reference_build") != row.get("reference_build")
                        or row.get("track") == "tr" and metadata.get("catalog_id") != row.get("catalog_id")):
                    raise InputValidationError("Generated source reference/catalog identity differs from its frozen plan.")
            if group == "annotation":
                if (entry.get("source_vcf") != str(source.get("path")) or entry.get("source_tool") != source.get("source_tool")
                        or entry.get("source_variant_ids") != ";".join(row.get("variant_ids", []))
                        or (entry.get("bnd_bundle_manifest") or "") != (str(source.get("bnd_bundle_manifest")) if source.get("bnd_bundle_manifest") else "")):
                    raise InputValidationError("Materialized annotation source differs from its explicit producer/plan.")
            else:
                allowed_states = {"CALLED", "NO_CALLS"} if row.get("origin") == "GENERATED" else {row.get("state", "CALLED")}
                if (entry.get("source_path") != (str(source["path"]) if source.get("path") else "")
                        or entry.get("index_path") != (str(source["index_path"]) if source.get("index_path") else "")
                        or entry.get("source_tool") != (source.get("source_tool") or "")
                        or entry.get("source_version") != (source.get("source_version") or "")
                        or entry.get("state") not in allowed_states
                        or entry.get("reference_build") != row.get("reference_build")
                        or entry.get("catalog_id") != (row.get("catalog_id") or "")):
                    raise InputValidationError("Materialized cohort source/state differs from its explicit producer/plan.")


def _annotation_sources(root: Path, attempt: Path, config: dict, controls: dict[str, str],
                        sample_ids: tuple[str, ...]) -> tuple[list[dict], list[str]]:
    from hifivar.annotation import read_annotation_input_rows
    section = config["annotation"]
    planned = _planned_source_rows(root, config, controls, "annotation", sample_ids)
    value = section.get("input_manifest")
    owner = root
    if planned:
        if value is not None:
            raise InputValidationError("Conflicting ready annotation manifest and source-plan rows.")
        materialized = _materialized_controls(root, attempt, config, controls)
        if materialized is None:
            return planned, ["SOURCE_PLAN_MATERIALIZATION_NOT_VERIFIED"]
        owner, effective, controls = materialized
        value = effective["annotation"].get("input_manifest")
    if value is None:
        return [], ["NO_FROZEN_ANNOTATION_SOURCE_ROWS"]
    path = _owned_control(owner, value, controls)
    rows = read_annotation_input_rows(path)
    result = []
    for row in rows:
        if row["sample"].strip() not in sample_ids:
            raise InputValidationError("Frozen annotation manifest contains an undeclared sample.")
        key = row["sample"].strip(), row["variant_category"].strip().lower()
        origin = next((item["origin"] for item in planned if (item["sample_id"], item["track"]) == key), "READY_EXTERNAL_MANIFEST")
        result.append({"sample_id": key[0], "track": key[1], "origin": origin, "source_tool": row["source_tool"],
                       "source_vcf": row["source_vcf"], "source_variant_ids": row["source_variant_ids"],
                       "manifest_sha256": controls[path.relative_to(owner).as_posix()],
                       "source_materialization": "VERIFIED_ATTEMPT_HANDOFF" if planned else "READY_EXTERNAL_FROZEN_CONTROL"})
        if row.get("bnd_bundle_manifest"):
            result[-1]["bnd_bundle_manifest"] = row["bnd_bundle_manifest"]
            result[-1]["source_representation_role"] = "EXPLICIT_VERIFIED_BND_DERIVED_WITH_NATIVE_LEDGER"
    if planned and {(item["sample_id"], item["track"]) for item in result} != {(item["sample_id"], item["track"]) for item in planned}:
        raise InputValidationError("Materialized annotation rows differ from the frozen source plan.")
    return result, []


def _cohort_sources(root: Path, attempt: Path, config: dict, controls: dict[str, str],
                    sample_ids: tuple[str, ...]) -> tuple[list[dict], list[str]]:
    from hifivar.cohort import SampleCallState
    planned = _planned_source_rows(root, config, controls, "cohort", sample_ids)
    value = config["cohort"].get("input_manifest")
    owner = root
    if planned:
        if value is not None:
            raise InputValidationError("Conflicting ready cohort manifest and source-plan rows.")
        materialized = _materialized_controls(root, attempt, config, controls)
        if materialized is None:
            return planned, ["SOURCE_PLAN_MATERIALIZATION_NOT_VERIFIED"]
        owner, effective, controls = materialized
        value = effective["cohort"].get("input_manifest")
    if value is None:
        return [], ["NO_FROZEN_COHORT_SOURCE_ROWS"]
    path = _owned_control(owner, value, controls)
    fields = ("sample", "track", "state", "source_path", "index_path", "source_tool",
              "source_version", "reference_build", "catalog_id")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != fields:
            raise InputValidationError("Frozen cohort manifest has invalid columns.")
        rows = []
        seen = set()
        for row in reader:
            if (set(row) != set(fields) or any(value is None for value in row.values())
                    or row["sample"] not in sample_ids or row["track"] not in {"small_variants", "sv", "tr"}
                    or row["state"] not in {item.value for item in SampleCallState}):
                raise InputValidationError("Invalid frozen cohort manifest row.")
            key = row["sample"], row["track"]
            if key in seen:
                raise InputValidationError("Duplicate frozen cohort sample/track row.")
            seen.add(key)
            rows.append(dict(row))
            if len(rows) > 4096:
                raise InputValidationError("Frozen cohort manifest exceeds 4096 rows.")
    for track in ("small_variants", "sv", "tr"):
        section = config["cohort"].get(track, {})
        if not isinstance(section, dict) or not isinstance(section.get("enabled", False), bool):
            raise InputValidationError("Cohort track enabled flags must be booleans.")
        if section.get("enabled") is True:
            if {row["sample"] for row in rows if row["track"] == track} != set(sample_ids):
                raise InputValidationError("Frozen cohort manifest omits a declared sample in an enabled track.")
    if planned and {(item["sample"], item["track"]) for item in rows} != {(item["sample_id"], item["track"]) for item in planned}:
        raise InputValidationError("Materialized cohort rows differ from the frozen source plan.")
    return rows, []


def _expected_units(root: Path, attempt: Path, config: dict, controls: dict[str, str], name: str,
                     sample_ids: tuple[str, ...]) -> tuple[list[tuple[str, tuple[str, ...], str | None]], list[str], list[dict]]:
    if name == 'benchmark':
        if not (attempt/'benchmark/materialization.json').exists():
            return [], ['BENCHMARK_SOURCE_MATERIALIZATION_NOT_AVAILABLE'], []
        from hifivar.benchmark_handoff import validate_materialization
        _,_,sources = validate_materialization(root,attempt)
        return [('benchmark',(config['benchmark']['sample_id'],),None)], [], [
            {'sample_id':config['benchmark']['sample_id'],'track':track,'origin':row['origin'],
             'source_version_scope':row['source_version_scope'],'truth_context':row['truth_metadata']} for track,row in sources.items()]
    if name == 'review':
        if not config.get('review',{}).get('source_plan'):
            return [], ['LEGACY_REVIEW_HAS_NO_SOURCE_PLAN_RECEIPT_BRIDGE'], []
        path = attempt/'review/materialization.json'
        if not path.exists():
            return [], ['REVIEW_SOURCE_MATERIALIZATION_NOT_AVAILABLE'], []
        from hifivar.review_handoff import validate_review_materialization
        record,_,resolved,selected = validate_review_materialization(root,attempt)
        rows = [{'sample_id':b['sample_id'],'track':b['track'],'origin':b['origin'],
                 'candidate_review_ids':b['candidate_review_ids'],
                 'source_version_scope':b['version_source'],'human_review':'NOT_REVIEWED'} for b in resolved.source_bindings]
        return [('review',selected,None)], [], rows
    if name == "annotation":
        rows, pending = _annotation_sources(root, attempt, config, controls, sample_ids)
        units = []
        for row in rows:
            applicable = []
            for backend in ("annovar", "vep", "overlap"):
                enabled = config[name].get(backend + "_enabled", False)
                if not isinstance(enabled, bool):
                    raise InputValidationError("Annotation backend enabled flags must be booleans.")
                if enabled and (backend != "overlap" or row["track"] in {"sv", "tr"}):
                    applicable.append(backend)
                    units.append(("annotation." + backend, (row["sample_id"],), row["track"]))
            if not applicable:
                pending.append("NO_APPLICABLE_BACKEND_FOR_" + row["track"])
        if config[name].get("functional_enabled"):
            pending.append("FUNCTIONAL_BACKEND_NOT_VALIDATED_BY_STEP_RECEIPTS")
        return units, pending, rows
    if name == "cohort":
        rows, pending = _cohort_sources(root, attempt, config, controls, sample_ids)
        units = []
        for track in ("small_variants", "sv", "tr"):
            section = config[name].get(track, {})
            if not isinstance(section, dict) or not isinstance(section.get("enabled", False), bool):
                raise InputValidationError("Cohort track enabled flags must be booleans.")
            if section.get("enabled", False):
                units.append(("cohort." + track, sample_ids, None))
        if not units:
            pending.append("NO_ENABLED_COHORT_TRACK")
        return units, pending, rows
    supported = {"alignment", "small", "tr", "phasing", "assembly"} | {
        "sv." + caller for caller in ("sawfish", "sniffles2", "pbsv", "cutesv", "harmonization")
    } | {"assembly_sv." + caller for caller in ("pav", "svim_asm")}
    if name not in supported:
        return [], ["NO_SUPPORTED_ATTEMPT_RECEIPT_CONTRACT"], []
    return [(name, (sample,), None) for sample in sample_ids], [], []


def _read_units(root: Path, attempt_id: str, units: list[tuple[str, tuple[str, ...], str | None]], *,
                diagnose_failed_artifacts: bool = False) -> tuple[tuple[ReportArtifact, ...], list[dict], list[dict]]:
    from hifivar.attempt_receipts import read_step_receipt, receipt_name
    artifacts, summaries, actual = [], [], []
    for track, samples, unit_id in units:
        summary = {"track": track, "samples": list(samples), "unit_id": unit_id,
                   "status": "MISSING_RECEIPT", "producer_attempt_id": None}
        try:
            record = read_step_receipt(root, attempt_id, track, samples, verify_outputs=True, unit_id=unit_id)
        except OutputValidationError:
            if not diagnose_failed_artifacts:
                raise
            # Preserve the failure without copying external diagnostics, trusting
            # a changed artifact, or adding a missing file as a valid deliverable.
            summary.update(status="FAILED_OUTPUT_VALIDATION", error_type="OutputValidationError",
                           diagnostic="An explicitly recorded artifact is missing or differs from its recorded identity.")
            summaries.append(summary)
            continue
        if record is not None:
            metadata = record.get("metadata", {})
            if not isinstance(metadata, dict):
                raise InputValidationError("Receipt metadata must be a mapping.")
            receipt_path = root / "attempts" / attempt_id / "receipts" / receipt_name(track, samples, unit_id)
            receipt_sha = _metadata_digest(receipt_path)
            if _mapping_file(receipt_path, yaml_input=False) != record:
                raise InputValidationError("Receipt changed after output verification.")
            summary.update(status="COMPLETE", producer_attempt_id=record["producer_attempt_id"],
                           provenance_class=record["provenance_class"], execution_kind=record["execution_kind"],
                           commands=record["commands"], tool_versions=record["tool_versions"],
                           validation=record["validation"], metadata=metadata,
                           receipt_path=str(receipt_path), receipt_sha256=receipt_sha,
                           output_roles=list(record["outputs"]))
            if track == "alignment":
                summary.update(alignment_action=record["action"],
                               selected_input_artifacts=record["selected_artifacts"])
            artifacts.append(ReportArtifact(receipt_path, "execution-receipt", "json", samples[0] if len(samples) == 1 else None, receipt_sha))
            for role, item in record["outputs"].items():
                artifacts.append(ReportArtifact(Path(item["path"]), role, "receipt-declared-output",
                                                samples[0] if len(samples) == 1 else None, item["sha256"]))
            actual.append(record)
        summaries.append(summary)
    return tuple(artifacts), summaries, actual


def _receipt_tools(records: list[dict]) -> tuple[ToolRecord, ...]:
    result, seen = [], set()
    for record in records:
        sources = record.get("metadata", {}).get("tool_version_sources", {})
        if not isinstance(sources, dict):
            raise InputValidationError("Receipt tool version sources must be a mapping.")
        for name, version in record["tool_versions"].items():
            source = sources.get(name, "RECEIPT_RECORDED_NOT_REPROBED_BY_REPORT")
            if not isinstance(source, str) or not source.strip():
                raise InputValidationError("Receipt tool version source must be nonempty text.")
            key = name, version, record["track"], source
            if key not in seen:
                result.append(ToolRecord(name, version, backend=record["track"], validation_status=source))
                seen.add(key)
    return tuple(result)


def _safe_path(value: Path, label: str, *, must_exist: bool) -> Path:
    if not isinstance(value, (str, Path)) or (isinstance(value, str) and not value.strip()):
        raise InputValidationError(f"{label} must be a nonempty path.")
    path = Path(value).expanduser()
    if ".." in path.parts:
        raise InputValidationError(f"{label} must not contain '..'.")
    absolute = path.absolute()
    for current in (absolute, *absolute.parents):
        if current.is_symlink() or (hasattr(current, "is_junction") and current.is_junction()):
            raise InputValidationError(f"Redirected {label.lower()} is forbidden: '{current}'.")
    if absolute == Path(absolute.anchor):
        raise InputValidationError(f"{label} must not be a filesystem root.")
    if must_exist and not absolute.is_dir():
        raise InputValidationError(f"{label} is not an existing directory: '{absolute}'.")
    return absolute


def _safe_file(path: Path) -> None:
    _safe_path(path.parent, "Metadata parent", must_exist=True)
    if path.is_symlink() or (hasattr(path, "is_junction") and path.is_junction()):
        raise InputValidationError(f"Redirected metadata is forbidden: '{path}'.")
    if not path.is_file() or path.stat().st_size == 0:
        raise InputValidationError(f"Required launcher metadata missing or empty: '{path}'.")
    if path.stat().st_size > 16 * 1024 * 1024:
        raise InputValidationError(f"Launcher metadata exceeds the bounded 16 MiB limit: '{path}'.")


def _mapping_file(path: Path, *, yaml_input: bool) -> dict:
    _safe_file(path)
    try:
        text = path.read_text(encoding="utf-8")
        value = yaml.safe_load(text) if yaml_input else json.loads(text)
    except (OSError, UnicodeError, ValueError, yaml.YAMLError) as error:
        raise InputValidationError(f"Invalid launcher metadata '{path}': {error}") from error
    if not isinstance(value, dict):
        raise InputValidationError(f"Launcher metadata must be a mapping: '{path}'.")
    return value


def _samples(sheet: str) -> list[dict[str, object]]:
    reader = csv.DictReader(io.StringIO(sheet), delimiter="\t")
    if not reader.fieldnames or "sample_id" not in reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
        raise InputValidationError("Launcher sample sheet has invalid column names.")
    result, seen = [], set()
    for row in reader:
        sample = row.get("sample_id")
        if not isinstance(sample, str) or not sample or any(char.isspace() for char in sample) or sample in seen or None in row or any(value is None for value in row.values()):
            raise InputValidationError("Launcher sample sheet has missing/duplicate sample identities or malformed rows.")
        seen.add(sample)
        result.append(dict(row))
    if not result:
        raise InputValidationError("Launcher sample sheet contains no samples.")
    return result


def _select_attempt(root: Path, attempt_id: str | None) -> tuple[Path | None, dict | None]:
    attempts = root / "attempts"
    if attempt_id is not None and (not isinstance(attempt_id, str) or _ATTEMPT_ID.fullmatch(attempt_id) is None):
        raise InputValidationError("attempt_id must identify one immediate attempt-directory name.")
    if attempts.is_symlink() or (hasattr(attempts, "is_junction") and attempts.is_junction()):
        raise InputValidationError("Redirected attempts directory is forbidden.")
    if not attempts.exists():
        if attempt_id is not None:
            raise InputValidationError(f"Requested attempt does not exist: '{attempt_id}'.")
        return None, None
    _safe_path(attempts, "Attempts directory", must_exist=True)
    entries = list(attempts.iterdir())
    for path in entries:
        if path.is_symlink() or (hasattr(path, "is_junction") and path.is_junction()):
            raise InputValidationError(f"Redirected attempt entry is forbidden: '{path}'.")
    choices = [path for path in entries if path.is_dir()]
    if attempt_id is None:
        if len(choices) > 1:
            raise InputValidationError("MULTIPLE_ATTEMPTS: supply an explicit attempt_id; names/mtimes are not selection evidence.")
        if not choices:
            return None, None
        chosen = choices[0]
    else:
        chosen = attempts / attempt_id
        if not chosen.is_dir():
            raise InputValidationError(f"Requested attempt does not exist: '{attempt_id}'.")
    _safe_path(chosen, "Selected attempt", must_exist=True)
    status_file = chosen / "status.json"
    if not status_file.exists() and not status_file.is_symlink():
        return chosen, None
    status = _mapping_file(status_file, yaml_input=False)
    started = _timestamp(status.get("started"), "started")
    if status.get("finished") is not None and _timestamp(status["finished"], "finished") < started:
        raise InputValidationError("Attempt finished timestamp precedes started timestamp.")
    return chosen, status


def _timestamp(value: object, field: str) -> datetime:
    if not isinstance(value, str):
        raise InputValidationError(f"Attempt requires a traceable ISO8601 {field} timestamp.")
    try:
        timestamp = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as error:
        raise InputValidationError(f"Invalid attempt {field} timestamp.") from error
    if timestamp.tzinfo is None or timestamp.utcoffset() is None:
        raise InputValidationError(f"Attempt {field} timestamp requires an explicit timezone.")
    return timestamp


def _execution_status(status: dict | None) -> tuple[FinalStatus, str, list[str]]:
    if status is None:
        return FinalStatus.NOT_RUN, "NO_STATUS_RECORD", []
    state = status.get("status")
    if not isinstance(state, str) or state not in _KNOWN_STATES:
        raise InputValidationError(f"Unsupported launcher attempt status: {state!r}.")
    if not isinstance(status.get("dry_run"), bool):
        raise InputValidationError("Attempt status requires a boolean dry_run field.")
    code = status.get("returncode")
    if code is not None and (not isinstance(code, int) or isinstance(code, bool)):
        raise InputValidationError("Attempt returncode must be an integer, not a boolean/string.")
    if state == "FAILED" or (code is not None and code != 0):
        return FinalStatus.FAILED, str(state), ["Selected attempt failed; pre-existing successful outputs do not override its failure."]
    if state in _PLANNED_STATES:
        if not status["dry_run"]:
            raise InputValidationError("Planning-only status conflicts with dry_run=false.")
        return FinalStatus.NOT_RUN, str(state), ["Dry-run/planning success is not analysis execution."]
    if state == "WORKFLOW_SUCCEEDED":
        if status["dry_run"] or code != 0 or status.get("finished") is None or status.get("input_immutability") != "PASS":
            raise InputValidationError("WORKFLOW_SUCCEEDED lacks a coherent finalized execution/immutability record.")
        return FinalStatus.COMPLETE, state, []
    return FinalStatus.PARTIAL, str(state), ["Selected attempt is running, interrupted or not finalized; no process liveness check was performed."]


def _configured_tracks(config: dict) -> list[tuple[str, str, str, bool]]:
    definitions = [
        ("alignment", "2", "sample"), ("small", "3", "small"),
        ("tr", "5", "tr"), ("phasing", "6", "phasing"),
        ("assembly", "7", "assembly"), ("review", "10", "review"),
        ("annotation", "11", "annotation"), ("cohort", "12", "cohort"),
        ("benchmark", "13", "benchmark"), ("low_depth", "15", "benchmark"),
    ]
    def enabled(section: object) -> bool:
        if not isinstance(section, dict) or not isinstance(section.get("enabled", False), bool):
            raise InputValidationError("Effective config track sections require boolean enabled values.")
        return section.get("enabled", False)
    tracks = [(name, phase, category, enabled(config.get(name, {}))) for name, phase, category in definitions]
    for group, callers, phase in (("sv", ("sawfish", "sniffles2", "pbsv", "cutesv", "harmonization"), "4"), ("assembly_sv", ("pav", "svim_asm"), "8")):
        section = config.get(group, {})
        group_enabled = enabled(section)
        for caller in callers:
            caller_enabled = enabled(section.get(caller, {}))
            tracks.append((f"{group}.{caller}", "9" if caller == "harmonization" else phase, "sv", caller_enabled and (group_enabled or caller == "harmonization")))
    return tracks


__all__ = ["build_launcher_report"]
