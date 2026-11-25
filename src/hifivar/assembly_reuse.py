"""Attempt-bound validation of external assemblies; never a native caller receipt."""
from __future__ import annotations

import csv
import io
import os
from pathlib import Path
from typing import Sequence

from hifivar.assembly import AssemblyRole, HaplotypeAssemblyArtifact
from hifivar.exceptions import InputValidationError
from hifivar.ready_assemblies import ReadyAssemblyArtifact, load_ready_assembly_manifest

HANDOFF_NAME = "ready_assembly_input_validation.json"
_KIND = "hifivar_ready_assembly_input_validation"
_SCOPE = "THIS_ATTEMPT_VALIDATED_REUSE_NOT_NEW_ASSEMBLY"


def _frozen(root: Path, attempt: Path) -> tuple[dict, dict, tuple[str, ...]]:
    from hifivar.attempt_receipts import _path, _read, artifact_identity
    root, attempt = _path(root), _path(attempt, root=root)
    if attempt.parent != root / "attempts":
        raise InputValidationError("Ready assembly validation needs an owned attempt directory.")
    manifest = _read(root / "launch_manifest.json")
    identity = manifest.get("identity", {})
    config = identity.get("config", {})
    section = config.get("assembly_sv", {})
    if (identity.get("schema") != 2 or identity.get("run_dir") != str(root)
            or not section.get("enabled") or not section.get("input_manifest")
            or config.get("assembly", {}).get("enabled")):
        raise InputValidationError("Frozen run does not select explicit ready assembly reuse.")
    samples = tuple(row["sample_id"] for row in csv.DictReader(io.StringIO(identity["sample_sheet"]), delimiter="\t"))
    if not samples or len(set(samples)) != len(samples):
        raise InputValidationError("Frozen ready assembly sample IDs are missing or repeated.")
    controls = manifest.get("control_sha256", {})
    if not {"effective.yaml", "samples.tsv"} <= set(controls):
        raise InputValidationError("Ready assembly validation lacks frozen launcher controls.")
    for name, digest in controls.items():
        target = _path(root / name, root=root)
        if artifact_identity(target)["sha256"] != digest:
            raise InputValidationError("Ready assembly launcher control changed.")
    return manifest, identity, samples


def _input_identity(item: ReadyAssemblyArtifact) -> dict:
    return {"path": str(item.path), "bytes": item.file_size, "sha256": item.sha256}


def _artifact(row: dict) -> ReadyAssemblyArtifact:
    """Refuse relabelled/invented source metadata in the preparation record."""
    try:
        item = ReadyAssemblyArtifact(
            row["sample_id"], AssemblyRole(row["role"]), Path(row["path"]), row["sha256"],
            row["file_size"], row["sequence_count"], row["base_count"],
            row["source_tool"], row["source_version"], Path(row["manifest_path"]), row["manifest_sha256"],
        )
    except (KeyError, TypeError, ValueError) as error:
        raise InputValidationError("Ready assembly handoff contains a malformed artifact.") from error
    if item.to_dict() != row:
        raise InputValidationError("Ready assembly handoff origin/role metadata differs from its contract.")
    return item


def prepare_ready_handoff(root: Path, attempt: Path) -> Path:
    """Validate original manifest/FASTAs before the execution context is sealed."""
    from hifivar.attempt_receipts import _write_new, artifact_identity
    _, identity, samples = _frozen(root, attempt)
    ready = load_ready_assembly_manifest(Path(identity["config"]["assembly_sv"]["input_manifest"]), expected_samples=samples)
    frozen = identity.get("inputs", [])
    selected = [artifact_identity(ready.path), *(_input_identity(item) for item in ready.assemblies)]
    if any(item not in frozen for item in selected):
        raise InputValidationError("Ready assembly bytes/path differ from frozen launch inputs.")
    record = {"schema_version": 1, "kind": _KIND, "status": "VALIDATED_INPUTS",
              "execution_kind": "VALIDATION_ONLY", "provenance_class": _SCOPE,
              "run_dir": str(root), "attempt_id": attempt.name,
              "manifest_sha256": artifact_identity(root / "launch_manifest.json")["sha256"],
              "samples": list(samples), "input_manifest": selected[0],
              "declared_sample_order": list(ready.declared_sample_order),
              "assemblies": [item.to_dict() for item in ready.assemblies],
              "native_commands": [], "new_hifiasm_execution": False, "phase_accuracy": "NOT_ASSESSED"}
    path = attempt / HANDOFF_NAME
    _write_new(path, record)
    return path


def read_ready_handoff(root: Path, attempt: Path, binding: dict, *,
                       verify_inputs: bool = False, sample_id: str | None = None) -> tuple[dict, tuple[ReadyAssemblyArtifact, ...]]:
    """Verify exact attempt/control binding; optionally rehash selected FASTAs."""
    from hifivar.attempt_receipts import _read, artifact_identity
    _, identity, samples = _frozen(root, attempt)
    path = attempt / HANDOFF_NAME
    if artifact_identity(path) != binding:
        raise InputValidationError("Ready assembly validation record differs from the sealed context.")
    record = _read(path)
    if type(record.get("schema_version")) is not int or record.get("new_hifiasm_execution") is not False:
        raise InputValidationError("Ready input validation has invalid schema or execution semantics.")
    expected = {"schema_version": 1, "kind": _KIND, "status": "VALIDATED_INPUTS",
                "execution_kind": "VALIDATION_ONLY", "provenance_class": _SCOPE,
                "run_dir": str(root), "attempt_id": attempt.name,
                "manifest_sha256": artifact_identity(root / "launch_manifest.json")["sha256"],
                "samples": list(samples), "native_commands": [], "new_hifiasm_execution": False,
                "phase_accuracy": "NOT_ASSESSED"}
    if any(record.get(key) != value for key, value in expected.items()):
        raise InputValidationError("Ready input reuse is not linked to this frozen attempt.")
    selected_manifest = artifact_identity(Path(identity["config"]["assembly_sv"]["input_manifest"]))
    if record.get("input_manifest") != selected_manifest or selected_manifest not in identity["inputs"]:
        raise InputValidationError("Ready assembly manifest differs from the frozen input identity.")
    rows = record.get("assemblies")
    if not isinstance(rows, list) or not rows or any(not isinstance(row, dict) for row in rows):
        raise InputValidationError("Ready input validation has no explicit assembly rows.")
    items = tuple(_artifact(row) for row in rows)
    if set(item.sample_id for item in items) != set(samples):
        raise InputValidationError("Ready assembly handoff sample set differs from the frozen run.")
    if (len({(item.sample_id, item.role) for item in items}) != len(items)
            or len({item.path for item in items}) != len(items)):
        raise InputValidationError("Ready assembly handoff duplicates roles or input paths.")
    try:
        physical = [(item.path.stat().st_dev, item.path.stat().st_ino) for item in items]
    except OSError as error:
        raise InputValidationError("Ready assembly input is missing or unreadable.") from error
    if len(set(physical)) != len(physical):
        raise InputValidationError("Ready assembly handoff inputs now alias physical files.")
    order = tuple(dict.fromkeys(item.sample_id for item in items))
    if record.get("declared_sample_order") != list(order):
        raise InputValidationError("Ready assembly declared sample order changed.")
    for item in items:
        if (_input_identity(item) not in identity["inputs"] or str(item.manifest_path) != selected_manifest["path"]
                or item.manifest_sha256 != selected_manifest["sha256"]):
            raise InputValidationError("Ready assembly row is not bound to frozen input bytes/manifest.")
    selected = tuple(item for item in items if sample_id is None or item.sample_id == sample_id)
    if not selected:
        raise InputValidationError("Selected sample is not present in the ready assembly validation.")
    if verify_inputs:
        from hifivar.ready_assemblies import protect_ready_assembly_inputs
        protect_ready_assembly_inputs(selected)
    return record, tuple(sorted(selected, key=lambda item: (samples.index(item.sample_id), item.role.value)))


def workflow_assemblies(sample_id: str, paths: Sequence[Path], config: dict) -> tuple:
    """Resolve actual script inputs without converting ready FASTAs to generated ones."""
    from hifivar.attempt_receipts import ATTEMPT_CONTEXT_ENV, _context, read_step_receipt
    context_value = os.environ.get(ATTEMPT_CONTEXT_ENV)
    context_path = Path(context_value) if context_value else None
    section = config.get("assembly_sv", {})
    if section.get("input_manifest"):
        if context_path is None:
            from hifivar.context import AnalysisContext
            ready = load_ready_assembly_manifest(Path(section["input_manifest"]),
                                                 expected_samples=AnalysisContext.from_config(config).sample_ids)
            items = ready.for_sample(sample_id)
        else:
            context, _, root = _context(context_path, live=True)
            _, items = read_ready_handoff(root, context_path.parent, context.get("ready_assembly_handoff"),
                                          verify_inputs=True, sample_id=sample_id)
        if tuple(Path(p).resolve() for p in paths) != tuple(item.path for item in items):
            raise InputValidationError("Actual ready caller input paths/roles differ from the frozen declaration.")
        return items
    if len(paths) != 2:
        raise InputValidationError("Generated assembly caller requires explicit h1/h2 FASTAs.")
    version, details, reference_independent = None, {}, None
    if context_path is not None:
        context, _, root = _context(context_path, live=True)
        producer = read_step_receipt(root, context["attempt_id"], "assembly", (sample_id,))
        if producer is None:
            raise InputValidationError("Generated assembly source has no current-attempt assembly receipt.")
        version = producer["tool_versions"].get("hifiasm")
        details = {entry["role"]: entry for entry in producer.get("metadata", {}).get("assembly_artifacts", [])}
        reference_independent = producer.get("metadata", {}).get("reference_independent")
    return tuple(HaplotypeAssemblyArtifact(sample_id, role, Path(path),
                 Path(details[role.value]["source_gfa"]) if details.get(role.value, {}).get("source_gfa") else None,
                 version, Path(path).stat().st_size, reference_independent)
                 for role, path in zip((AssemblyRole.HAPLOTYPE1, AssemblyRole.HAPLOTYPE2), paths))


def ready_work_directory(work_root: Path, caller: str, sample_id: str) -> str:
    """Keep explicit resume attempts out of previous caller working files."""
    from hifivar.attempt_receipts import ATTEMPT_CONTEXT_ENV, _context
    from hifivar.sample import validate_sample_id
    validate_sample_id(sample_id)
    if caller not in {"pav", "svim_asm"}:
        raise InputValidationError("Unsupported ready assembly caller work directory.")
    base = Path(work_root) / caller / sample_id
    context_value = os.environ.get(ATTEMPT_CONTEXT_ENV)
    if context_value:
        context, _, root = _context(Path(context_value), live=False)
        if Path(work_root).absolute() != root / "work":
            raise InputValidationError("Ready assembly work root differs from its owned attempt.")
        base = base / context["attempt_id"]
    return str(base)


__all__ = ["HANDOFF_NAME", "prepare_ready_handoff", "read_ready_handoff", "workflow_assemblies", "ready_work_directory"]
