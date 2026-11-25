"""Offline review preparation and explicit, independently stored decisions."""
from __future__ import annotations

import csv
import io
import json
import re
from collections import Counter
from copy import deepcopy
from datetime import datetime
from pathlib import Path

from hifivar import __version__
from hifivar.alignment_postprocess import find_alignment_index
from hifivar.config import HiFiVarConfig
from hifivar.context import AnalysisContext
from hifivar.exceptions import ConfigurationError, InputValidationError, OutputValidationError
from hifivar.igv import IgvWrapper
from hifivar.junction_consumers import revalidate_junction_record
from hifivar.logging_utils import get_logger
from hifivar.review import ReviewEvidence, ReviewManifest, ReviewResult, ReviewStatus, read_review_selection
from hifivar.sample import InputType
from hifivar.serialization import redact_sensitive_data, utc_now_iso8601, write_json_atomic, write_yaml_atomic
from hifivar.validation import compute_sha256, validate_file

_LOGGER = get_logger(__name__)
_DECISION_FIELDS = ("review_id", "sample", "variant_id", "status", "reviewer", "reviewed_at", "comment")
_SEALED_FILES = ("review.igv.batch", "review_manifest.json", "review_manifest.yaml",
                 "review_manifest.tsv", "effective_config.yaml")


def _linked(path: Path) -> bool:
    """Detect symbolic links and Windows junctions where supported by pathlib."""
    return path.is_symlink() or getattr(path, "is_junction", lambda: False)()


def _output_path(value: str | Path) -> Path:
    """Reject occupied or symlink-routed destinations before any write."""
    output = Path(value).expanduser().absolute()
    if output == Path(output.anchor) or output == Path.home():
        raise OutputValidationError("Review output must be a new, specific project subdirectory.")
    if any(_linked(part) for part in (output, *output.parents)):
        raise OutputValidationError("Review output cannot use a symlink directory.")
    if output.exists():
        raise OutputValidationError(f"Review output already exists: '{output}'. Use a new directory.")
    return output


def _reserve(output: Path) -> None:
    """Own a fresh leaf directory; retain partial evidence on write failure."""
    _output_path(output)
    try:
        output.mkdir(parents=True, exist_ok=False)
    except OSError as error:
        raise OutputValidationError(f"Unable to create review directory '{output}': {error}") from error


def _write_text(path: Path, content: str) -> None:
    """Write exclusively inside an owned directory; never replace old files."""
    try:
        with path.open("xt", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
    except OSError as error:
        raise OutputValidationError(f"Unable to write review artifact '{path}': {error}") from error


def _rows(path: Path, fields: tuple[str, ...] | None = None) -> list[dict[str, str]]:
    """Read explicitly selected review tables with strict header and row shape."""
    validate_file(path, require_nonempty=True)
    try:
        with path.open("rt", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, [])
            if not header or len(header) != len(set(header)) or any(not name for name in header):
                raise InputValidationError("Review TSV has empty or duplicate column names.")
            if fields is not None and set(header) != set(fields):
                raise InputValidationError(f"Decision TSV requires exactly these columns: {fields!r}.")
            result = []
            for line, values in enumerate(reader, 2):
                if not values:
                    continue
                if len(values) != len(header):
                    raise InputValidationError(f"Review TSV line {line} has an incorrect number of fields.")
                if any(any(ord(character) < 32 for character in value) for value in values):
                    raise InputValidationError(f"Review TSV line {line} contains a control character.")
                result.append(dict(zip(header, values, strict=True)))
    except (OSError, UnicodeError, csv.Error) as error:
        raise InputValidationError(f"Unable to read review TSV '{path}': {error}") from error
    return result


def _tsv(rows: list[dict], fields: tuple[str, ...]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows({key: row.get(key, "") for key in fields} for row in rows)
    return stream.getvalue()


def _metadata(paths: list[Path]) -> dict[str, dict[str, int]]:
    """Record stat metadata only, explicitly not full raw-data hashes."""
    result = {}
    for path in paths:
        try:
            stat = path.stat()
        except OSError as error:
            raise InputValidationError(f"Unable to inspect review input '{path}': {error}") from error
        result[str(path.absolute())] = {"bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    return result


def prepare_review(config: HiFiVarConfig, output_directory: str | Path, *, dry_run: bool = False) -> dict:
    """Prepare a batch and pending manifest without invoking IGV or callers."""
    output = _output_path(output_directory)
    section = config["review"]
    selection_value = section.get("selection_file")
    if not isinstance(selection_value, str) or not selection_value.strip():
        raise ConfigurationError("review prepare requires review.selection_file in YAML.")
    selection = Path(selection_value).expanduser()
    if not selection.is_absolute():
        user = config.sources.get("user")
        if user is None:
            raise ConfigurationError("Relative review.selection_file requires a user YAML source.")
        selection = user.absolute().parent / selection
    selection = selection.absolute()
    _rows(selection)
    selection_sha = compute_sha256(selection)
    context = AnalysisContext.from_config(config)
    if context.reference.build is None:
        raise ConfigurationError("review prepare requires an explicit reference.build.")
    alignments = {record.sample.sample_id: record.sample.input.files[0].absolute()
                  for record in context.samples if record.sample.input.input_type in (InputType.BAM, InputType.CRAM)}
    targets = read_review_selection(selection, alignments=alignments,
                                    reference_fasta=context.reference.fasta.absolute(),
                                    output_directory=output, default_flank_bp=section["flank_bp"], reference_build=context.reference.build)
    for target in targets:
        if target._junction_record is None and target.end > context.reference.get_contig(target.contig).length:
            raise InputValidationError(f"Review target '{target.review_id}' extends beyond its reference contig.")
        if target._junction_record is None and target.mate_contig and target.mate_position > context.reference.get_contig(target.mate_contig).length:
            raise InputValidationError(f"BND mate for '{target.review_id}' extends beyond its reference contig.")
    executable = section["igv_executable"]
    if ("/" in executable or "\\" in executable) and not Path(executable).is_absolute():
        user = config.sources.get("user")
        if user is None:
            raise ConfigurationError("Relative IGV executable requires a user YAML source.")
        executable = str((user.absolute().parent / executable).absolute())
    plan = IgvWrapper(executable=executable).plan(targets, batch_path=output / "review.igv.batch",
                                                snapshot_directory=output / "screenshots")
    manifest = ReviewManifest(tuple(ReviewResult(target, ReviewEvidence(
        target, plan.batch_path, plan.args, target.screenshot_paths, generated=False)) for target in targets))
    inputs = [context.reference.fasta, context.reference.fai, selection]
    if context.source_sample_sheet:
        inputs.append(context.source_sample_sheet)
    for target in targets:
        inputs.extend((target.alignment_path, target.source_vcf))
        if target.bnd_bundle_manifest:
            inputs.append(target.bnd_bundle_manifest)
        index = find_alignment_index(target.alignment_path)
        if index is not None:
            inputs.append(index)
        if str(target.source_vcf).endswith(".vcf.gz"):
            inputs.append(Path(f"{target.source_vcf}.tbi"))
        if target.trgt_visualization_path:
            inputs.append(target.trgt_visualization_path)
    before = _metadata(inputs)
    if compute_sha256(selection) != selection_sha:
        raise InputValidationError("Review selection changed during preparation.")
    result = {"status": "REVIEW_PLAN" if dry_run else "REVIEW_PREPARED",
              "output_directory": str(output), "targets": len(targets),
              "expected_screenshots": len(plan.screenshots), "screenshots_generated": False,
              "manual_decisions_completed": 0, "igv_command": list(plan.args)}
    for target in targets:
        if target._junction_record:
            revalidate_junction_record(target._junction_record)
    if dry_run:
        return result
    _reserve(output)
    (output / "screenshots").mkdir()
    _write_text(plan.batch_path, plan.batch_text)
    manifest.write_json(output / "review_manifest.json")
    manifest.write_yaml(output / "review_manifest.yaml")
    manifest.write_tsv(output / "review_manifest.tsv")
    write_yaml_atomic(redact_sensitive_data(config.to_dict()), output / "effective_config.yaml")
    pending = [{"review_id": target.review_id, "sample": target.sample_id,
                "variant_id": target.variant_id, "status": ReviewStatus.NOT_REVIEWED.value,
                "reviewer": "", "reviewed_at": "", "comment": ""} for target in targets]
    _write_text(output / "manual_decisions.tsv", _tsv(pending, _DECISION_FIELDS))
    for target in targets:
        if target._junction_record:
            revalidate_junction_record(target._junction_record)
    if _metadata(inputs) != before or compute_sha256(selection) != selection_sha:
        raise InputValidationError("Review input metadata changed; bundle has no completion receipt.")
    receipt = {"schema_version": 1, "kind": "hifivar_review_preparation",
               "hifivar_version": __version__, "created_at": utc_now_iso8601(),
               "result": result, "reference_build": context.reference.build,
               "selection_sha256": selection_sha, "source_metadata": before,
               "raw_input_hash_policy": "STAT_ONLY_NOT_FULL_CONTENT_HASH",
               "sealed_artifacts": {name: compute_sha256(output / name) for name in _SEALED_FILES},
               "editable_template": "manual_decisions.tsv",
               "scientific_policy": manifest.to_dict()["scientific_policy"],
               "validation_scope": "Lightweight Phase 10 inputs; not BAM header/record or selected VCF-record identity validation."}
    # Written last: absent receipt means incomplete preparation, not success.
    write_json_atomic(receipt, output / "review_bundle.json")
    _LOGGER.info("Review prepared: %d targets, %d planned screenshots; no human decisions", len(targets), len(plan.screenshots))
    return result


def _load_bundle(directory: Path) -> tuple[dict, dict]:
    """Check the sealed review metadata; do not evaluate embedded commands."""
    try:
        if _linked(directory) or _linked(directory / "review_bundle.json"):
            raise InputValidationError("Review bundle cannot be a symlink.")
        receipt = json.loads(validate_file(directory / "review_bundle.json").read_text(encoding="utf-8"))
        if type(receipt.get("schema_version")) is not int or receipt.get("schema_version") != 1 or receipt.get("kind") not in {"hifivar_review_preparation", "hifivar_review_snapshots"}:
            raise InputValidationError("Unsupported review bundle identity/schema.")
        seals = receipt.get("sealed_artifacts")
        if not isinstance(seals, dict) or set(seals) != set(_SEALED_FILES):
            raise InputValidationError("Review bundle sealed file list is invalid.")
        for name, expected in seals.items():
            path = directory / name
            if _linked(path) or compute_sha256(path) != expected:
                raise InputValidationError(f"Review bundle artifact changed: {name}.")
        manifest = json.loads((directory / "review_manifest.json").read_text(encoding="utf-8"))
        rows = manifest.get("results")
        if not isinstance(rows, list):
            raise InputValidationError("Review bundle results must be a list.")
        ids = [row["review_id"] for row in rows]
        if len(ids) != len(set(ids)) or any(row["status"] != "NOT_REVIEWED" for row in rows):
            raise InputValidationError("Prepared review IDs/statuses are invalid.")
        if any(row.get("status_is_truth") is not False for row in rows):
            raise InputValidationError("Review bundle must retain non-truth semantics.")
        if manifest.get("scientific_policy", {}).get("manual_status_is_truth") is not False:
            raise InputValidationError("Review bundle scientific policy is invalid.")
        if receipt['kind'] == 'hifivar_review_snapshots':
            result = directory / 'snapshot_result.json'
            if _linked(result) or compute_sha256(result) != receipt.get('snapshot_result_sha256'):
                raise InputValidationError('Snapshot execution record changed.')
            shots = receipt.get('snapshot_sha256')
            expected = {str(Path(path).relative_to(directory).as_posix()) for row in rows for path in row['screenshots']}
            if not isinstance(shots,dict) or not shots or set(shots) != expected:
                raise InputValidationError('Snapshot receipt does not seal the exact expected screenshot set.')
            for name, expected_sha in shots.items():
                if re.fullmatch(r'screenshots/[A-Za-z0-9_.-]+\.png',name) is None:
                    raise InputValidationError('Invalid sealed screenshot relative path.')
                path = directory / name
                if _linked(path) or _linked(path.parent) or compute_sha256(path) != expected_sha:
                    raise InputValidationError('Sealed snapshot image changed.')
            sessions = receipt.get('session_sha256')
            expected_sessions = {str(Path(name).with_suffix('.session.xml').as_posix()) for name in shots}
            if (receipt.get('session_validation_policy') != 'REFERENCE_LOCUS_AND_SOURCE_RESOURCES'
                    or not isinstance(sessions,dict) or set(sessions) != expected_sessions):
                raise InputValidationError('Snapshot receipt lacks exact per-image session bindings.')
            for name,expected_sha in sessions.items():
                path=directory/name
                if _linked(path) or compute_sha256(path)!=expected_sha:
                    raise InputValidationError('Sealed snapshot session changed.')
    except (OSError, UnicodeError, ValueError, AttributeError, KeyError, TypeError) as error:
        raise InputValidationError(f"Unable to read prepared review bundle: {error}") from error
    return receipt, manifest


def import_review_decisions(bundle_directory: str | Path, decisions_file: str | Path,
                            output_directory: str | Path, *, dry_run: bool = False) -> dict:
    """Import only user-declared decisions; preserve the original review bundle."""
    output = _output_path(output_directory)
    bundle = Path(bundle_directory).expanduser().absolute()
    if output.resolve().is_relative_to(bundle.resolve()):
        raise OutputValidationError("Imported decisions must use a separate directory outside the prepared bundle.")
    decisions = Path(decisions_file).expanduser().absolute()
    receipt, original = _load_bundle(bundle)
    decision_sha = compute_sha256(decisions)
    prepared = {row["review_id"]: row for row in original["results"]}
    supplied = {}
    for row in _rows(decisions, _DECISION_FIELDS):
        identifier = row["review_id"]
        if identifier not in prepared or identifier in supplied:
            raise InputValidationError(f"Unknown or duplicate review_id: '{identifier}'.")
        source = prepared[identifier]
        if row["sample"] != source["sample"] or row["variant_id"] != source["variant_id"]:
            raise InputValidationError(f"Sample/variant identity mismatch for '{identifier}'.")
        try:
            status = ReviewStatus(row["status"])
        except ValueError as error:
            raise InputValidationError(f"Unsupported manual-review status for '{identifier}'.") from error
        if status is not ReviewStatus.NOT_REVIEWED:
            if any(not row[key].strip() for key in ("reviewer", "reviewed_at", "comment")):
                raise InputValidationError(f"Decision for '{identifier}' requires reviewer, reviewed_at and comment.")
            try:
                timestamp = datetime.fromisoformat(row["reviewed_at"].replace("Z", "+00:00"))
                if timestamp.utcoffset() is None:
                    raise ValueError("timezone is required")
            except ValueError as error:
                raise InputValidationError(f"Decision for '{identifier}' requires a timezone-aware ISO timestamp.") from error
        elif any(row[key] for key in ("reviewer", "reviewed_at", "comment")):
            raise InputValidationError("NOT_REVIEWED rows must leave reviewer, date and comment empty.")
        supplied[identifier] = row
    updated = deepcopy(original)
    for row in updated["results"]:
        declared = supplied.get(row["review_id"])
        row.update({"reviewer": "", "reviewed_at": "", "comment": ""})
        if declared:
            row.update({key: declared[key] for key in ("status", "reviewer", "reviewed_at", "comment")})
            row["notes"] = declared["comment"]
        row["reviewer_identity_verified"] = False
    counts = dict(Counter(row["status"] for row in updated["results"]))
    result = {"status": "DECISIONS_PLAN" if dry_run else "DECISIONS_IMPORTED",
              "output_directory": str(output), "total_targets": len(prepared),
              "reviewed_targets": sum(value for key, value in counts.items() if key != "NOT_REVIEWED"),
              "status_counts": counts, "missing_decision_ids": sorted(set(prepared) - set(supplied)),
              "reviewer_identity_verified": False, "status_is_truth": False}
    if compute_sha256(decisions) != decision_sha:
        raise InputValidationError("Decision TSV changed during import.")
    if dry_run:
        return result
    _reserve(output)
    updated["decision_provenance"] = {
        "imported_at": utc_now_iso8601(), "hifivar_version": __version__,
        "prepared_bundle": str(bundle), "prepared_receipt_sha256": compute_sha256(bundle / "review_bundle.json"),
        "decision_tsv": str(decisions), "decision_tsv_sha256": decision_sha,
        "reviewer_identity_verified": False, "source": "USER_SUPPLIED_TSV",
        "raw_input_revalidation": "NOT_PERFORMED_DURING_DECISION_IMPORT"}
    write_json_atomic(updated, output / "review_manifest.adjudicated.json")
    write_yaml_atomic(updated, output / "review_manifest.adjudicated.yaml")
    _write_text(output / "review_decisions.adjudicated.tsv", _tsv(updated["results"], _DECISION_FIELDS))
    _load_bundle(bundle)
    if compute_sha256(decisions) != decision_sha:
        raise InputValidationError("Decision TSV changed; imported directory has no completion receipt.")
    write_json_atomic({"schema_version": 1, "kind": "hifivar_review_decision_import", "result": result,
                       "prepared_seals": receipt["sealed_artifacts"],
                       "provenance": updated["decision_provenance"]}, output / "decision_import.json")
    _LOGGER.info("Imported %d user-declared decisions of %d targets; not independently authenticated", result["reviewed_targets"], len(prepared))
    return result


__all__ = ["prepare_review", "import_review_decisions"]
