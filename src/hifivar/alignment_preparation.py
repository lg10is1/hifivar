"""Attempt-bound alignment preparation; reuse certifies validation, not alignment.

These receipts precede the calling DAG. They never turn external or previously
produced BAMs into fresh caller outputs, and do not establish biological truth.
"""
from __future__ import annotations

import csv
import io
from copy import deepcopy
from pathlib import Path

import yaml

from hifivar.attempt_receipts import (
    _allowed, _commands, _path, _read, _write_new, artifact_identity, receipt_name,
)
from hifivar.command_observation import is_current_observation
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.serialization import utc_now_iso8601


def _state(root: Path, attempt: Path, *, live: bool) -> dict:
    root, attempt = _path(root), _path(attempt, root=root)
    if attempt.parent != root / "attempts":
        raise InputValidationError("Alignment preparation requires an immediate launcher attempt.")
    manifest = _read(root / "launch_manifest.json")
    identity = manifest.get("identity", {})
    if identity.get("schema") != 2 or identity.get("run_dir") != str(root):
        raise InputValidationError("Alignment preparation launcher identity mismatch.")
    controls = manifest.get("control_sha256", {})
    if not {"effective.yaml", "samples.tsv"}.issubset(controls):
        raise InputValidationError("Alignment preparation lacks frozen controls.")
    for name, digest in controls.items():
        if not isinstance(name, str) or Path(name).is_absolute() or ".." in Path(name).parts:
            raise InputValidationError("Unsafe alignment preparation control.")
        if artifact_identity(_path(root / name, root=root))["sha256"] != digest:
            raise InputValidationError("Alignment preparation control changed.")
    if live:
        status = _read(attempt / "status.json")
        if status.get("status") != "ALIGNMENT_PREPARATION" or status.get("dry_run") is not False:
            raise InputValidationError("Alignment preparation has no active native attempt.")
        from hifivar.launcher import runtime_identity
        from hifivar.package_resources import installed_workflow_root
        if runtime_identity(installed_workflow_root()) != identity.get("runtime"):
            raise InputValidationError("Alignment preparation runtime differs from frozen source.")
    return identity


def preparation_path(attempt: Path, sample: str) -> Path:
    """Return one deterministic receipt, never discover results by globbing."""
    return attempt / "receipts" / receipt_name("alignment", (sample,))


def publish_preparation(root: Path, attempt: Path, route: dict, bam: Path,
                        index: Path, validation: dict, commands, *, action: str,
                        execution: dict | None = None, prior: Path | None = None) -> dict:
    """Seal one observed alignment or this-attempt validation of reused bytes."""
    identity = _state(root, attempt, live=True)
    sample = route["sample_id"]
    _allowed(identity, "alignment", (sample,))
    if not is_current_observation(commands):
        raise OutputValidationError("Alignment preparation requires current CommandRunner observation.")
    observed = _commands(commands, "NATIVE")
    if action not in {"ALIGN_NEW", "VALIDATE_EXTERNAL_REUSE", "VALIDATE_PRIOR_ALIGNMENT_REUSE"}:
        raise InputValidationError("Unknown alignment preparation action.")
    if (validation.get("sample_id") != sample or validation.get("quickcheck") != "PASS"
            or validation.get("header_identity") != "PASS" or validation.get("indexed_access") != "PASS"):
        raise OutputValidationError("Alignment preparation requires successful indexed BAM validation.")
    steps = validation.get("steps", [])
    if (not steps or any(step.get("command") not in [item["args"] for item in observed] for step in steps)
            or not all(any(token in step.get("command", []) for step in steps) for token in ("quickcheck", "-H", "-X"))):
        raise OutputValidationError("BAM validation commands were not observed in this preparation.")
    bam, index = _path(bam), _path(index)
    if bam == index or bam.samefile(index):
        raise InputValidationError("Alignment BAM and index must be distinct artifacts.")
    selected = {"bam": artifact_identity(bam), "index": artifact_identity(index)}
    frozen = {item["path"]: item for item in identity["inputs"]}
    rows = list(csv.DictReader(io.StringIO(identity["sample_sheet"]), delimiter="\t"))
    declared = next(row for row in rows if row["sample_id"] == sample)
    if route["inputs"] != declared["input"].split(";"):
        raise InputValidationError("Alignment preparation source paths differ from declared sample.")
    if any(path not in frozen for path in route["inputs"]):
        raise InputValidationError("Alignment preparation source is not frozen.")
    source_inputs = [frozen[path] for path in route["inputs"]]
    prior_identity = artifact_identity(prior) if prior is not None else None
    producer = None
    if action == "ALIGN_NEW":
        owner = root / "results" / "alignment" / sample
        if not bam.is_relative_to(owner) or bam.parent.parent != owner or not bam.parent.name.startswith("attempt_"):
            raise InputValidationError("New alignment must own a fresh isolated output directory.")
        tool = identity["config"]["alignment"]["tool"]
        executable = identity["config"]["alignment"][tool + "_executable"]
        if not any(item["args"][0] == executable and all(path in item["args"] for path in route["inputs"])
                   and identity["config"]["reference"]["fasta"] in item["args"]
                   and ("align" in item["args"] if tool == "pbmm2" else "-a" in item["args"])
                   for item in observed):
            raise OutputValidationError("New alignment lacks an observed configured mapper/source/reference command.")
        if not execution or not execution.get("tool_version"):
            raise OutputValidationError("New alignment lacks its detected tool version.")
        if prior is not None:
            raise InputValidationError("New alignment cannot claim a prior completion.")
        producer = attempt.name
    elif action == "VALIDATE_EXTERNAL_REUSE":
        if prior is not None or route["kind"] != "ALIGNED_BAM" or route["inputs"] != [str(bam)]:
            raise InputValidationError("External BAM reuse differs from its declared source.")
        if any(frozen.get(item["path"]) != item for item in selected.values()):
            raise InputValidationError("External reused BAM/index differs from frozen input bytes.")
    else:
        expected = root / "results" / "alignment" / sample / "completed.json"
        if prior != expected:
            raise InputValidationError("Prior alignment reuse requires its exact completion record.")
        saved = _read(prior)
        if (saved.get("sample_id") != sample or saved.get("inputs") != route["inputs"]
                or saved.get("bam") != str(bam) or saved.get("artifacts") != list(selected.values())):
            raise InputValidationError("Prior alignment completion no longer seals selected BAM/index.")
        # Legacy records have no producer receipt. Do not retrofit an origin.
        producer = saved.get("producer_attempt_id")
        if producer is not None:
            origin = _path(root / "attempts" / producer, root=root)
            if origin.parent != root / "attempts" or origin == attempt or artifact_identity(preparation_path(origin, sample)) != saved.get("producer_receipt"):
                raise InputValidationError("Prior alignment origin receipt has changed.")
            old = read_preparation(root, origin, sample)
            if old is None or old.get("action") != "ALIGN_NEW" or old.get("selected_artifacts") != selected:
                raise InputValidationError("Prior alignment origin does not certify the selected bytes.")
    validation_path = attempt / "alignment_validation" / (sample + ".json")
    _write_new(validation_path, validation)
    outputs = {"validation": artifact_identity(validation_path)}
    if action == "ALIGN_NEW":
        outputs.update(selected)
    record = {
        "schema_version": 1, "kind": "hifivar_alignment_preparation_receipt", "status": "COMPLETE",
        "track": "alignment", "samples": [sample], "unit_id": None,
        "attempt_id": attempt.name, "producer_attempt_id": producer,
        "manifest_sha256": artifact_identity(root / "launch_manifest.json")["sha256"],
        "runtime": identity["runtime"], "action": action,
        "execution_kind": "NATIVE" if action == "ALIGN_NEW" else "VALIDATION_ONLY",
        "provenance_class": "THIS_ATTEMPT_ALIGNED" if action == "ALIGN_NEW" else "THIS_ATTEMPT_VALIDATED_REUSE_NOT_NEW_ALIGNMENT",
        "source_inputs": source_inputs, "selected_artifacts": selected, "prior_completion": prior_identity,
        "outputs": outputs, "commands": observed,
        "tool_versions": {"samtools": validation["samtools_version"], **({"mapper": execution["tool_version"]} if execution else {})},
        "validation": validation, "metadata": {"action": action, "producer_origin": producer or "NOT_PROVEN_LEGACY_OR_EXTERNAL"},
        "scientific_accuracy": "NOT_ASSESSED", "created": utc_now_iso8601(),
    }
    _write_new(preparation_path(attempt, sample), record)
    return record


def read_preparation(root: Path, attempt: Path, sample: str, *, verify_outputs: bool = True) -> dict | None:
    """Verify a selected preparation; validation-only reuse remains explicit."""
    target = preparation_path(attempt, sample)
    if not target.exists():
        return None
    identity = _state(root, attempt, live=False)
    _allowed(identity, "alignment", (sample,))
    record = _read(target)
    expected = {"schema_version": 1, "kind": "hifivar_alignment_preparation_receipt", "status": "COMPLETE",
                "track": "alignment", "samples": [sample], "unit_id": None, "attempt_id": attempt.name,
                "manifest_sha256": artifact_identity(root / "launch_manifest.json")["sha256"],
                "runtime": identity["runtime"], "scientific_accuracy": "NOT_ASSESSED"}
    if any(record.get(key) != value for key, value in expected.items()):
        raise InputValidationError("Alignment preparation receipt/attempt mismatch.")
    action = record.get("action")
    fresh = action == "ALIGN_NEW"
    if action not in {"ALIGN_NEW", "VALIDATE_EXTERNAL_REUSE", "VALIDATE_PRIOR_ALIGNMENT_REUSE"}:
        raise InputValidationError("Invalid recorded alignment action.")
    if (record.get("execution_kind") != ("NATIVE" if fresh else "VALIDATION_ONLY")
            or record.get("provenance_class") != ("THIS_ATTEMPT_ALIGNED" if fresh else "THIS_ATTEMPT_VALIDATED_REUSE_NOT_NEW_ALIGNMENT")
            or (fresh and record.get("producer_attempt_id") != attempt.name)):
        raise InputValidationError("Alignment reuse cannot certify new alignment execution.")
    _commands(record.get("commands"), "NATIVE", reading=True)
    outputs, selected = record.get("outputs", {}), record.get("selected_artifacts", {})
    if set(selected) != {"bam", "index"} or set(outputs) != ({"validation", "bam", "index"} if fresh else {"validation"}):
        raise InputValidationError("Alignment preparation artifact roles are incomplete.")
    if verify_outputs:
        items = [*outputs.values(), *selected.values()]
        if record.get("prior_completion") is not None:
            items.append(record["prior_completion"])
        for item in items:
            if artifact_identity(Path(item["path"])) != item:
                raise OutputValidationError("Alignment preparation artifact changed.")
    return record


def seal_handoff(root: Path, attempt: Path) -> dict:
    """Freeze exactly the declared ordered BAMs and calling controls for a DAG."""
    identity = _state(root, attempt, live=False)
    declared = list(csv.DictReader(io.StringIO(identity["sample_sheet"]), delimiter="\t"))
    with (root / "calling.samples.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    expected_config = deepcopy(identity["config"])
    expected_config["samples"]["sheet"] = str(root / "calling.samples.tsv")
    if yaml.safe_load((root / "calling.yaml").read_text(encoding="utf-8")) != expected_config:
        raise InputValidationError("Calling config is not the frozen config with only its sample sheet replaced.")
    if [row["sample_id"] for row in rows] != [row["sample_id"] for row in declared]:
        raise InputValidationError("Calling sample order/set differs from the declared cohort.")
    receipts, selected = [], []
    for original, row in zip(declared, rows, strict=True):
        record = read_preparation(root, attempt, row["sample_id"])
        if record is None or row["input_type"] != "bam" or row["input"] != record["selected_artifacts"]["bam"]["path"]:
            raise InputValidationError("Calling sample has no exact attempt-linked prepared BAM.")
        if any(row[key] != value for key, value in original.items() if key not in {"input", "input_type"}):
            raise InputValidationError("Calling sample metadata differs from its declared source.")
        receipts.append(artifact_identity(preparation_path(attempt, row["sample_id"])))
        selected.extend(record["selected_artifacts"].values())
    value = {"schema_version": 1, "kind": "hifivar_alignment_calling_handoff", "attempt_id": attempt.name,
             "declared_sample_order": [row["sample_id"] for row in declared], "receipts": receipts,
             "selected_artifacts": selected, "controls": [artifact_identity(root / name) for name in ("calling.yaml", "calling.samples.tsv")]}
    _write_new(attempt / "alignment_handoff.json", value)
    return value


def verify_handoff(root: Path, attempt: Path, binding: dict, *, verify_bams: bool = False) -> dict:
    """Cheap controls per caller; full selected BAM hashing at launcher boundaries."""
    if artifact_identity(attempt / "alignment_handoff.json") != binding:
        raise InputValidationError("Alignment calling handoff changed.")
    value = _read(attempt / "alignment_handoff.json")
    if value.get("kind") != "hifivar_alignment_calling_handoff" or value.get("attempt_id") != attempt.name:
        raise InputValidationError("Alignment calling handoff attempt mismatch.")
    expected_controls = {str(root / name) for name in ("calling.yaml", "calling.samples.tsv")}
    if {item["path"] for item in value.get("controls", [])} != expected_controls:
        raise InputValidationError("Alignment calling handoff lacks exact calling controls.")
    for item in [*value["controls"], *value["receipts"], *(value["selected_artifacts"] if verify_bams else [])]:
        if artifact_identity(Path(item["path"])) != item:
            raise InputValidationError("Alignment calling handoff control, receipt or selected artifact changed.")
    return value
