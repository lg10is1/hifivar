"""Attempt-linked execution receipts, never inferred from result filenames.

Receipts certify bounded engineering facts: observed CommandRunner exits,
explicit wrapper validation and output byte identities. They do not establish
biological accuracy, call confidence, or independent human review.
"""

from __future__ import annotations

import hashlib
import csv
import io
import json
import math
import os
import re
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Mapping, Sequence

import yaml

from hifivar.command_observation import CommandObservation, ObservedCommand, is_current_observation, observe_commands
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.serialization import redact_sensitive_data, standardize_data, utc_now_iso8601
from hifivar.sample import validate_sample_id


ATTEMPT_CONTEXT_ENV = "HIFIVAR_ATTEMPT_CONTEXT"
_TOKEN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}")
_HASH = re.compile(r"[0-9a-f]{64}")
_LIMIT = 1024 * 1024


def required_output_roles(track: str) -> frozenset[str]:
    """Minimum native wrapper output contract, distinct from data accuracy."""
    explicit = {
        "small": {"vcf", "index", "gvcf", "gvcf_index"},
        "tr": {"vcf", "index", "spanning_bam", "spanning_bam_index"},
        "phasing": {"vcf", "index"},
        "assembly": {"primary_fasta", "haplotype1_fasta", "haplotype2_fasta"},
        "sv.harmonization": {"vcf", "index", "evidence", "provenance"},
        "annotation.annovar": {"tsv", "vcf"},
        "annotation.vep": {"tsv"}, "annotation.overlap": {"tsv"},
        "cohort.small_variants": {"bcf", "bcf_index", "vcf", "index", "provenance"},
        "cohort.sv": {"sites_tsv", "sample_matrix_tsv", "provenance"},
        "cohort.tr": {"loci_tsv", "sample_matrix_tsv", "provenance"},
        "review": {"manifest_json", "manifest_yaml", "manifest_tsv", "batch", "provenance"},
        "benchmark": {"manifest", "metrics"},
    }
    if track in explicit:
        return frozenset(explicit[track])
    if track in {"sv." + name for name in ("sawfish", "sniffles2", "pbsv", "cutesv")} | {
        "assembly_sv.pav", "assembly_sv.svim_asm"
    }:
        return frozenset({"vcf", "index"})
    raise InputValidationError("Receipt track has no implemented output contract.")


def _path(value: Path | str, *, root: Path | None = None) -> Path:
    if not isinstance(value, (str, Path)):
        raise InputValidationError("Attempt paths must be explicit text or pathlib paths.")
    path = Path(value)
    if not path.is_absolute() or ".." in path.parts:
        raise InputValidationError("Attempt paths must be absolute and contain no '..'.")
    for part in (path, *path.parents):
        if part.is_symlink() or (hasattr(part, "is_junction") and part.is_junction()):
            raise InputValidationError("Redirected attempt paths are forbidden.")
    if path == Path(path.anchor) or (root is not None and not path.is_relative_to(root)):
        raise InputValidationError("Attempt path escapes its owned run directory.")
    return path


def artifact_identity(path: Path) -> dict:
    """Stream-hash a regular file and reject changes during the read."""
    path = _path(path)
    if not path.is_file():
        raise OutputValidationError(f"Receipt artifact is not a regular file: '{path}'.")
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    after = path.stat()
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (
        after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns
    ):
        raise OutputValidationError("Receipt artifact changed during checksum verification.")
    return {"path": str(path), "bytes": after.st_size, "sha256": digest.hexdigest()}


def _read(path: Path) -> dict:
    path = _path(path)
    if not path.is_file() or not 0 < path.stat().st_size <= _LIMIT:
        raise InputValidationError("Attempt metadata is missing, empty or exceeds 1 MiB.")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError) as error:
        raise InputValidationError("Invalid attempt receipt metadata.") from error
    if not isinstance(value, dict):
        raise InputValidationError("Attempt metadata must be a mapping.")
    return value


def _write_new(path: Path, value: dict) -> None:
    """Publish a complete sibling file atomically, without replacing any file."""
    path = _path(path)
    if path.exists():
        raise OutputValidationError(f"Attempt metadata already exists: '{path}'.")
    payload = (json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")
    if len(payload) > _LIMIT:
        raise OutputValidationError("Attempt metadata exceeds 1 MiB.")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=".receipt-", suffix=".partial", dir=path.parent)
    temporary_path = Path(temporary)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary_path, path)
        except FileExistsError as error:
            raise OutputValidationError("Refusing to overwrite an existing attempt receipt.") from error
    finally:
        temporary_path.unlink(missing_ok=True)


def _samples(value: Sequence[str]) -> tuple[str, ...]:
    if isinstance(value, (str, bytes)):
        raise InputValidationError("Receipt samples must be an explicit sequence.")
    samples = tuple(value)
    if not samples or any(not isinstance(item, str) for item in samples) or len(set(samples)) != len(samples):
        raise InputValidationError("Receipt sample identities are empty, repeated or unsafe.")
    for sample in samples:
        validate_sample_id(sample)
    return samples


def receipt_name(track: str, samples: Sequence[str], unit_id: str | None = None) -> str:
    """Map an explicit unit to a deterministic filename, without globbing."""
    if not isinstance(track, str) or not _TOKEN.fullmatch(track):
        raise InputValidationError("Receipt track must be a safe nonempty identifier.")
    if unit_id is not None and (not isinstance(unit_id, str) or not _TOKEN.fullmatch(unit_id)):
        raise InputValidationError("Receipt unit must be an explicit safe identifier.")
    digest = hashlib.sha256(json.dumps([track, _samples(samples), unit_id], separators=(",", ":")).encode()).hexdigest()
    return f"{track}-{digest}.json"


def create_attempt_context(root: Path, attempt: Path, *, dry_run: bool) -> Path:
    """Seal a new context only after the launcher has frozen its manifest."""
    root, attempt = _path(root), _path(attempt, root=root)
    if attempt.parent != root / "attempts" or not _TOKEN.fullmatch(attempt.name):
        raise InputValidationError("Attempt context must belong to an immediate attempt directory.")
    if not isinstance(dry_run, bool):
        raise InputValidationError("Attempt dry_run must be boolean.")
    manifest_path = root / "launch_manifest.json"
    manifest = _read(manifest_path)
    identity = manifest.get("identity", {})
    if identity.get("schema") != 2 or identity.get("run_dir") != str(root):
        raise InputValidationError("Attempt context does not match the launcher manifest.")
    path = attempt / "execution_context.json"
    record = {
        "schema_version": 1, "kind": "hifivar_attempt_execution_context",
        "run_dir": str(root), "attempt_id": attempt.name, "dry_run": dry_run,
        "manifest_sha256": artifact_identity(manifest_path)["sha256"],
        "runtime": identity.get("runtime"), "created": utc_now_iso8601(),
    }
    handoff = attempt / "alignment_handoff.json"
    if handoff.exists() and identity.get("config", {}).get("alignment", {}).get("enabled"):
        from hifivar.alignment_preparation import verify_handoff
        binding = artifact_identity(handoff)
        verify_handoff(root, attempt, binding)
        record["alignment_handoff"] = binding
    if identity.get("config", {}).get("assembly_sv", {}).get("input_manifest"):
        from hifivar.assembly_reuse import HANDOFF_NAME, read_ready_handoff
        binding = artifact_identity(attempt / HANDOFF_NAME)
        read_ready_handoff(root, attempt, binding)
        record["ready_assembly_handoff"] = binding
    _write_new(path, record)
    return path


def _context(path: Path, *, live: bool) -> tuple[dict, dict, Path]:
    path = _path(path)
    context = _read(path)
    if context.get("schema_version") != 1 or context.get("kind") != "hifivar_attempt_execution_context":
        raise InputValidationError("Unsupported attempt execution context.")
    root = _path(context.get("run_dir", ""))
    attempt_id = context.get("attempt_id")
    if not isinstance(attempt_id, str) or not _TOKEN.fullmatch(attempt_id):
        raise InputValidationError("Invalid context attempt identity.")
    attempt = root / "attempts" / attempt_id
    if path != attempt / "execution_context.json" or not isinstance(context.get("dry_run"), bool):
        raise InputValidationError("Attempt context path/identity mismatch.")
    manifest_path = root / "launch_manifest.json"
    if artifact_identity(manifest_path)["sha256"] != context.get("manifest_sha256"):
        raise InputValidationError("Attempt launcher manifest has changed.")
    manifest = _read(manifest_path)
    identity = manifest.get("identity", {})
    if identity.get("schema") != 2 or identity.get("run_dir") != str(root) or identity.get("runtime") != context.get("runtime"):
        raise InputValidationError("Attempt runtime/manifest identity mismatch.")
    controls = manifest.get("control_sha256")
    if not isinstance(controls, dict) or not {"effective.yaml", "samples.tsv"}.issubset(controls):
        raise InputValidationError("Attempt manifest lacks frozen config/sample controls.")
    for name, digest in controls.items():
        if not isinstance(name, str) or not name or Path(name).is_absolute() or ".." in Path(name).parts:
            raise InputValidationError("Unsafe launcher control path.")
        if artifact_identity(_path(root / name, root=root))["sha256"] != digest:
            raise InputValidationError("Frozen launcher control changed.")
    if context.get("alignment_handoff") is not None:
        from hifivar.alignment_preparation import verify_handoff
        verify_handoff(root, attempt, context["alignment_handoff"])
    if identity.get("config", {}).get("assembly_sv", {}).get("input_manifest") or context.get("ready_assembly_handoff") is not None:
        from hifivar.assembly_reuse import read_ready_handoff
        read_ready_handoff(root, attempt, context.get("ready_assembly_handoff"))
    if live:
        if identity.get("config", {}).get("alignment", {}).get("enabled") and context.get("alignment_handoff") is None:
            raise InputValidationError("Calling after alignment requires an attempt-bound prepared BAM handoff.")
        status = _read(attempt / "status.json")
        if status.get("status") not in {"DAG_RUNNING", "PRODUCER_DAG_RUNNING", "DOWNSTREAM_DAG_RUNNING", "REVIEW_DAG_RUNNING", "BENCHMARK_DAG_RUNNING"} or status.get("dry_run") is not False:
            raise InputValidationError("No active native launcher attempt owns this receipt.")
        from hifivar.launcher import runtime_identity
        from hifivar.package_resources import installed_workflow_root
        actual = runtime_identity(installed_workflow_root())
        if actual != context.get("runtime"):
            raise InputValidationError("Receipt producer runtime differs from the frozen launcher runtime.")
    return context, identity, root


def _allowed(identity: dict, track: str, samples: tuple[str, ...]) -> None:
    import csv
    import io
    rows = csv.DictReader(io.StringIO(identity.get("sample_sheet", "")), delimiter="\t")
    declared = tuple(row.get("sample_id") for row in rows)
    if not declared or any(item not in declared for item in samples):
        raise InputValidationError("Receipt sample is not declared in the frozen run.")
    config = identity.get("config", {})
    known = {"small", "tr", "phasing", "assembly", "alignment", "review", "benchmark"} | {
        "sv." + name for name in ("sawfish", "sniffles2", "pbsv", "cutesv", "harmonization")
    } | {"assembly_sv." + name for name in ("pav", "svim_asm")} | {
        "annotation." + name for name in ("annovar", "vep", "overlap")
    } | {"cohort." + name for name in ("small_variants", "sv", "tr")}
    if track not in known:
        raise InputValidationError("Unsupported receipt track.")
    parts = track.split(".")
    section = config.get(parts[0], {})
    if not isinstance(section, dict) or (section.get("enabled") is not True and track != "sv.harmonization"):
        raise InputValidationError("Receipt track is not enabled in the frozen run.")
    if len(parts) == 2:
        if parts[0] == "annotation":
            sub_enabled = section.get(parts[1] + "_enabled") is True
        else:
            sub = section.get(parts[1], {})
            sub_enabled = isinstance(sub, dict) and sub.get("enabled") is True
        if not sub_enabled:
            raise InputValidationError("Receipt subtrack is not enabled in the frozen run.")
    elif len(parts) != 1:
        raise InputValidationError("Unsupported receipt track.")
    if parts[0] == "cohort" and samples != declared:
        raise InputValidationError("Cohort receipt must preserve all declared samples and their order.")


def _consumer(track: str) -> bool:
    return track in {"review", "benchmark"} or track.startswith(("annotation.", "cohort."))


def _stage(context_path: Path, identity: dict, track: str, *, live: bool) -> str:
    generated = bool(identity.get("config", {}).get("workflow", {}).get("downstream_source_plan"))
    expected = ("DOWNSTREAM_DAG_RUNNING" if _consumer(track) else "PRODUCER_DAG_RUNNING") if generated else "DAG_RUNNING"
    if track == "review":
        expected = "REVIEW_DAG_RUNNING"
    if track == "benchmark":
        expected = "BENCHMARK_DAG_RUNNING"
    if live and _read(context_path.parent / "status.json").get("status") != expected:
        raise InputValidationError("Receipt belongs to a different producer/consumer execution stage.")
    return expected


def _source_paths(inputs: Mapping[str, Mapping[str, Path | None]] | None, samples: tuple[str, ...]) -> dict:
    if not isinstance(inputs, Mapping) or set(inputs) != set(samples):
        raise InputValidationError("Consumer inputs must explicitly name exactly every selected sample.")
    result = {}
    for sample in samples:
        roles = inputs[sample]
        if not isinstance(roles, Mapping) or set(roles) - {"vcf", "index", "bnd_bundle_manifest"}:
            raise InputValidationError("Consumer source roles must be the explicit vcf/index pair.")
        if not roles:
            result[sample] = {}
            continue
        if roles.get("vcf") is None:
            raise InputValidationError("A callable consumer source needs a real VCF path.")
        result[sample] = {
            "vcf": _path(Path(roles["vcf"]).absolute()),
            "index": _path(Path(roles["index"]).absolute()) if roles.get("index") is not None else None,
        }
        if roles.get('bnd_bundle_manifest') is not None:
            result[sample]['bnd_bundle_manifest'] = _path(Path(roles['bnd_bundle_manifest']).absolute())
    return result


def _manifest_sources(path: Path, track: str, samples: tuple[str, ...], unit_id: str | None) -> dict:
    if not path.is_file() or not 0 < path.stat().st_size <= _LIMIT:
        raise InputValidationError("Consumer source manifest is missing or exceeds 1 MiB.")
    annotation = track.startswith("annotation.")
    fields = (("sample", "variant_category", "source_vcf", "source_tool", "source_variant_ids") if annotation else
              ("sample", "track", "state", "source_path", "index_path", "source_tool", "source_version", "reference_build", "catalog_id"))
    try:
        reader = csv.DictReader(io.StringIO(path.read_text(encoding="utf-8-sig")), delimiter="\t")
        actual_fields = tuple(reader.fieldnames or ())
        if actual_fields != fields and not (annotation and actual_fields == (*fields, 'bnd_bundle_manifest')):
            raise InputValidationError("Bound consumer source manifest has unexpected columns/order.")
        rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise InputValidationError("Cannot parse bound consumer source manifest.") from error
    seen, selected = set(), {}
    category = unit_id if annotation else track.split(".", 1)[1]
    for row in rows:
        if None in row or any(not isinstance(value, str) for value in row.values()):
            raise InputValidationError("Malformed bound consumer source row.")
        row = {key: value.strip() for key, value in row.items()}
        key = (row["sample"], row["variant_category" if annotation else "track"])
        if key in seen:
            raise InputValidationError("Duplicate bound consumer sample/category row.")
        seen.add(key)
        if key[0] not in samples or key[1] != category:
            continue
        source = row["source_vcf" if annotation else "source_path"]
        if annotation:
            if not source:
                raise InputValidationError("Annotation source path is missing.")
            selected[key[0]] = {"vcf": _path(source), "index": _path(source + ".tbi") if source.lower().endswith(".vcf.gz") else None}
            if row.get('bnd_bundle_manifest'):
                selected[key[0]]['bnd_bundle_manifest'] = _path(row['bnd_bundle_manifest'])
        elif row["state"] in {"CALLED", "NO_CALLS"}:
            if not source or not row["index_path"]:
                raise InputValidationError("CALLED/NO_CALLS cohort rows require real VCF/index sources.")
            selected[key[0]] = {"vcf": _path(source), "index": _path(row["index_path"])}
        elif row["state"] in {"NOT_RUN", "FAILED", "DISABLED", "MISSING_INPUT", "NOT_OBSERVED"}:
            if source or row["index_path"]:
                raise InputValidationError("Noncallable cohort state cannot attach a future source.")
            selected[key[0]] = {}
        else:
            raise InputValidationError("Unknown bound cohort source state.")
    if set(selected) != set(samples):
        raise InputValidationError("Bound consumer selection does not contain exactly the requested samples/unit.")
    return {sample: selected[sample] for sample in samples}


def _recorded_artifact(value: object) -> dict:
    if (not isinstance(value, dict) or type(value.get("bytes")) is not int or value["bytes"] < 0
            or not isinstance(value.get("sha256"), str) or not _HASH.fullmatch(value["sha256"])):
        raise InputValidationError("Bound source artifact identity is malformed.")
    return {"path": str(_path(value.get("path", ""))), "bytes": value["bytes"], "sha256": value["sha256"]}


def _harmonization_binding(context_path, context, identity, root, samples, inputs, input_manifest, *, live):
    """Bind each actual canonical source to frozen input or this-attempt producer."""
    from hifivar.harmonization_inputs import harmonization_sources
    if len(samples) != 1 or not isinstance(inputs, Mapping):
        raise InputValidationError('Harmonization requires one sample and explicit caller source pairs.')
    expected_manifest = root / 'effective.yaml'
    if input_manifest is None or _path(Path(input_manifest).absolute()) != expected_manifest:
        raise InputValidationError('Harmonization source binding requires the frozen effective configuration.')
    selected = harmonization_sources(identity['config'], samples[0], root)
    if set(inputs) != {item['caller'] for item in selected}:
        raise InputValidationError('Harmonization caller membership differs from frozen source selection.')
    frozen = {item['path']: _recorded_artifact(item) for item in identity.get('inputs', ())}
    from hifivar.assembly_source_roles import source_haplotype_binding
    witnesses, origins, producers, haplotype_roles = {}, {}, [], {}
    for item in selected:
        caller = item['caller']
        actual = inputs[caller]
        if not isinstance(actual, Mapping) or set(actual) != {'vcf', 'index'}:
            raise InputValidationError('Harmonization requires exact VCF/index source roles.')
        if any(value is None or _path(Path(value).absolute()) != item[role] for role, value in actual.items()):
            raise InputValidationError('Actual harmonization sources differ from frozen selected paths.')
        witness = {role: artifact_identity(item[role]) for role in ('vcf', 'index')}
        if item['origin'] == 'EXTERNAL':
            if any(frozen.get(info['path']) != info for info in witness.values()):
                raise InputValidationError('External harmonization source bytes differ from launch identity.')
        else:
            producer = read_step_receipt(root, context['attempt_id'], item['track'], samples, verify_outputs=False)
            if producer is None or any(producer['outputs'].get(role) != witness[role] for role in witness):
                raise InputValidationError('Generated harmonization source lacks a matching current-attempt producer receipt.')
            receipt_path = context_path.parent / 'receipts' / receipt_name(item['track'], samples)
            producers.append({'caller': caller, 'track': item['track'], 'receipt': artifact_identity(receipt_path)})
        haplotype_roles[caller] = source_haplotype_binding(samples[0], caller, external=item['origin']=='EXTERNAL', context_path=context_path)
        witnesses[caller], origins[caller] = witness, item['origin']
    return {'schema_version': 1, 'mode': 'PER_SAMPLE_MULTISOURCE_HARMONIZATION',
            'execution_stage': _stage(context_path, identity, 'sv.harmonization', live=live),
            'input_manifest': artifact_identity(expected_manifest), 'source_inputs': witnesses,
            'source_origins': origins, 'source_order': [item['caller'] for item in selected],
            'selection': {'track': 'sv.harmonization', 'samples': list(samples)},
            'producer_receipts': producers, 'source_haplotype_roles': haplotype_roles,
            'scientific_accuracy': 'NOT_ASSESSED'}


def _consumer_binding(context_path: Path, context: dict, identity: dict, root: Path, track: str,
                      samples: tuple[str, ...], unit_id: str | None, *, inputs, input_manifest,
                      live: bool) -> dict | None:
    """Bind actual request sources to sealed launcher controls, not caller metadata.

    This checks non-malicious integration mistakes and persisted byte identity;
    it is not a security sandbox against arbitrary code in this Python process.
    """
    if track == "sv.harmonization":
        return _harmonization_binding(context_path, context, identity, root, samples, inputs, input_manifest, live=live)
    if track == "benchmark":
        from hifivar.benchmark_handoff import receipt_binding
        return receipt_binding(context_path, context, identity, root, samples, inputs, input_manifest, live=live)
    if track == "review":
        from hifivar.review_handoff import review_receipt_binding
        return review_receipt_binding(context_path, context, identity, root, samples, inputs, input_manifest, live=live)
    if not _consumer(track):
        return None
    actual = _source_paths(inputs, samples)
    group = track.split(".", 1)[0]
    manifest = _read(root / "launch_manifest.json")
    frozen_controls = dict(manifest["control_sha256"])
    generated = bool(identity.get("config", {}).get("workflow", {}).get("downstream_source_plan"))
    producer_artifacts, selected_producers, origins, materialization_identity = [], [], {}, None
    if generated:
        status = _read(context_path.parent / "status.json")
        materialization = _path(context_path.parent / "downstream/materialization.json", root=root)
        materialization_identity = artifact_identity(materialization)
        if materialization_identity["sha256"] != status.get("downstream_materialization_sha256"):
            raise InputValidationError("Downstream materialization does not match the owned stage marker.")
        bound = _read(materialization)
        if (bound.get("schema_version") != 1 or bound.get("kind") != "hifivar_downstream_materialization"
                or bound.get("attempt_id") != context["attempt_id"]
                or bound.get("context_sha256") != artifact_identity(context_path)["sha256"]
                or bound.get("plan_sha256") != artifact_identity(root / "downstream_plan.yaml")["sha256"]):
            raise InputValidationError("Downstream materialization context/plan identity differs.")
        controls = bound.get("controls")
        if not isinstance(controls, dict) or not {"effective.yaml", group + ".inputs.tsv"}.issubset(controls):
            raise InputValidationError("Downstream materialization lacks its consumer controls.")
        for name, info in controls.items():
            if name not in {"effective.yaml", "annotation.inputs.tsv", "cohort.inputs.tsv"}:
                raise InputValidationError("Unexpected materialized consumer control.")
            expected_path = context_path.parent / "downstream" / name
            if _recorded_artifact(info) != artifact_identity(expected_path):
                raise InputValidationError("Materialized consumer control changed.")
        effective = context_path.parent / "downstream/effective.yaml"
        expected_manifest = context_path.parent / "downstream" / (group + ".inputs.tsv")
        # Validate links independently; consumer manifests must not substitute a
        # different index for the generated producer's selected output pair.
        links = bound.get("source_receipts", ())
        if not isinstance(links, (list, tuple)):
            raise InputValidationError("Generated producer links must be an explicit sequence.")
        for link in links:
            if not isinstance(link, dict) or _consumer(str(link.get("track", ""))):
                raise InputValidationError("Invalid generated producer receipt link.")
            expected_link = context_path.parent / "receipts" / receipt_name(link.get("track", ""), link.get("samples", ()))
            if _recorded_artifact(link) != artifact_identity(expected_link):
                raise InputValidationError("Selected generated producer receipt changed.")
            source_receipt = read_step_receipt(root, context["attempt_id"], link["track"], link["samples"], verify_outputs=False)
            if source_receipt is None:
                raise InputValidationError("Selected generated producer receipt is absent.")
            role = link.get("output_role")
            index_role = {"gvcf": "gvcf_index", "bnd_vcf": "bnd_index"}.get(role, "index")
            if link.get("index_role", index_role) != index_role:
                raise InputValidationError("Generated source link changed its selected index role.")
            if role not in source_receipt["outputs"] or index_role not in source_receipt["outputs"]:
                raise InputValidationError("Generated producer lacks its selected VCF/index roles.")
            item = {"receipt": artifact_identity(expected_link), "samples": source_receipt["samples"],
                    "vcf": source_receipt["outputs"][role], "index": source_receipt["outputs"][index_role]}
            if role == "bnd_vcf":
                if link["track"] != "sv.harmonization" or link.get("bnd_bundle_manifest_role") != "bnd_manifest":
                    raise InputValidationError("Generated BND source link lacks its explicit harmonization bundle role.")
                from hifivar.downstream_bnd import generated_bnd_roles
                from hifivar.reference import ReferenceGenome
                ref = identity["config"]["reference"]
                genome = ReferenceGenome.from_fasta(Path(ref["fasta"]), build=ref["build"])
                bound_roles = generated_bnd_roles(source_receipt, genome, source_receipt["samples"][0])
                item["bnd_bundle_manifest"] = bound_roles["bnd_bundle_manifest"]
            elif link.get("bnd_bundle_manifest_role") is not None:
                raise InputValidationError("Native generated source cannot attach a derived BND role.")
            producer_artifacts.append(item)
        source_validation = bound.get("source_validation")
        if not isinstance(source_validation, dict) or not isinstance(source_validation.get(group), list):
            raise InputValidationError("Materialized source selection audit is missing.")
        audits = source_validation[group]
        for row in audits:
            if not isinstance(row, dict):
                raise InputValidationError("Materialized source selection audit row is malformed.")
            key = (row.get("sample_id"), row.get("category" if group == "annotation" else "track"))
            if key in origins or row.get("source_origin") not in {"GENERATED", "EXTERNAL"}:
                raise InputValidationError("Materialized source origin is missing or repeated.")
            origins[key] = row["source_origin"]
    else:
        effective = root / "effective.yaml"
        expected_manifest = root / (group + ".inputs.tsv")
        if frozen_controls.get(expected_manifest.name) != artifact_identity(expected_manifest)["sha256"]:
            raise InputValidationError("External consumer manifest is not a frozen launcher control.")
    effective_identity = artifact_identity(effective)
    if not 0 < effective_identity["bytes"] <= _LIMIT:
        raise InputValidationError("Consumer effective configuration exceeds its bounded size.")
    try:
        effective_config = yaml.safe_load(effective.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as error:
        raise InputValidationError("Cannot parse bound consumer effective configuration.") from error
    if (not isinstance(effective_config, dict) or not isinstance(effective_config.get(group), dict)
            or effective_config[group].get("input_manifest") != str(expected_manifest)):
        raise InputValidationError("Consumer manifest does not match its effective configuration.")
    if input_manifest is None or _path(Path(input_manifest).absolute()) != expected_manifest:
        raise InputValidationError("Actual consumer manifest differs from the bound normalized manifest.")
    selected = _manifest_sources(expected_manifest, track, samples, unit_id)
    if actual != selected:
        raise InputValidationError("Actual consumer VCF/index sources differ from the bound manifest selection.")
    frozen_inputs = {}
    for item in identity.get("inputs", ()):
        info = _recorded_artifact(item)
        if info["path"] in frozen_inputs and frozen_inputs[info["path"]] != info:
            raise InputValidationError("Conflicting frozen external input identities.")
        frozen_inputs[info["path"]] = info
    witnesses = {}
    category = unit_id if group == "annotation" else track.split(".", 1)[1]
    for sample, roles in selected.items():
        witnesses[sample] = {role: artifact_identity(path) if path is not None else None for role, path in roles.items()}
        if not roles:
            continue
        origin = origins.get((sample, category)) if generated else "EXTERNAL"
        if origin == "GENERATED":
            pairs = [item for item in producer_artifacts if item["samples"] == [sample]
                     and item["vcf"] == witnesses[sample]["vcf"] and item["index"] == witnesses[sample]["index"]
                     and item.get("bnd_bundle_manifest") == witnesses[sample].get("bnd_bundle_manifest")]
            if not pairs:
                raise InputValidationError("Consumer source pair is not the selected generated producer output.")
            for pair in pairs:
                if pair not in selected_producers:
                    selected_producers.append(pair)
        elif origin == "EXTERNAL":
            if any(info is not None and frozen_inputs.get(info["path"]) != info for info in witnesses[sample].values()):
                raise InputValidationError("Consumer external input differs from its frozen launch identity.")
        else:
            raise InputValidationError("Consumer input has no proven materialized source origin.")
    bundle_bindings = {}
    for sample, roles in selected.items():
        if 'bnd_bundle_manifest' in roles:
            if group != 'annotation' or unit_id != 'sv':
                raise InputValidationError('BND derived roles are explicit SV annotation sources only.')
            from hifivar.junction_consumers import validate_bnd_consumer_bundle
            from hifivar.reference import ReferenceGenome
            ref = identity['config']['reference']
            genome = ReferenceGenome.from_fasta(Path(ref['fasta']), build=ref['build'])
            bundle_bindings[sample] = validate_bnd_consumer_bundle(roles['bnd_bundle_manifest'], roles['vcf'], genome, sample)
    return {"schema_version": 1, "mode": "MATERIALIZED_PLAN" if generated else "EXTERNAL_READY",
            "execution_stage": _stage(context_path, identity, track, live=live),
            "frozen_control_sha256": frozen_controls, "effective_config": effective_identity,
            "materialization": materialization_identity, "input_manifest": artifact_identity(expected_manifest),
            "selection": {"track": track, "samples": list(samples), "unit_id": unit_id},
            "source_inputs": witnesses, "producer_receipts": selected_producers,
            "junction_bundle_bindings": bundle_bindings,
            "source_validation_scope": "ACTUAL_SELECTED_SOURCE_ROLES_ONLY_OTHER_OUTPUT_BYTES_NOT_REHASHED"}


def _commands(records, execution_kind: str, *, reading: bool = False) -> list[dict]:
    if execution_kind not in {"NATIVE", "IN_PROCESS"}:
        raise InputValidationError("Unknown receipt execution kind.")
    if execution_kind == "NATIVE" and not records:
        raise OutputValidationError("Native completion requires observed CommandRunner results.")
    result = []
    for observed in records:
        if not reading and not isinstance(observed, ObservedCommand):
            raise OutputValidationError("Native results must originate from CommandRunner observation.")
        record = observed if reading else observed.to_dict()
        if not isinstance(record, dict) or record.get("executed") is not True or type(record.get("returncode")) is not int or record["returncode"] != 0:
            raise OutputValidationError("Planned, failed or unobserved commands cannot certify completion.")
        duration = record.get("duration_seconds")
        if isinstance(duration, bool) or not isinstance(duration, (int, float)) or not math.isfinite(duration) or duration < 0:
            raise OutputValidationError("Invalid observed command duration.")
        result.append({key: record.get(key) for key in (
            "args", "returncode", "duration_seconds", "cwd", "executed", "stdout_path", "stderr_path"
        )})
    if execution_kind == "NATIVE" and not any(
        isinstance(item.get("args"), list) and len(item["args"]) > 1
        and not any(arg in {"--version", "-V", "--help", "-h"} for arg in item["args"][1:])
        for item in result
    ):
        raise OutputValidationError("Version/help probes alone cannot certify native analysis.")
    return result


@dataclass
class StepReceipt:
    """A single step's explicit, validated output declaration."""

    context_path: Path | None
    track: str
    samples: tuple[str, ...]
    commands: CommandObservation = field(default_factory=CommandObservation)
    completed: bool = False
    expected_outputs: dict[str, Path] = field(default_factory=dict)
    unit_id: str | None = None
    source_inputs: Mapping[str, Mapping[str, Path | None]] | None = None
    input_manifest: Path | None = None
    consumer_binding: dict | None = None

    def complete(
        self, *, outputs: Mapping[str, Path], tool_versions: Mapping[str, str],
        validation: Mapping[str, object], metadata: Mapping[str, object] | None = None,
        execution_kind: str = "NATIVE",
    ) -> None:
        """Seal only after native execution and the tool's validators returned."""
        if self.context_path is None:
            return
        if self.completed:
            raise OutputValidationError("A step receipt cannot be completed twice.")
        if not isinstance(self.commands, CommandObservation) or not is_current_observation(self.commands):
            raise OutputValidationError("Command results must originate from this active step observation.")
        context, identity, root = _context(self.context_path, live=True)
        if context["dry_run"]:
            raise OutputValidationError("A dry-run cannot emit a completion receipt.")
        _allowed(identity, self.track, self.samples)
        _stage(self.context_path, identity, self.track, live=True)
        current_binding = _consumer_binding(
            self.context_path, context, identity, root, self.track, self.samples, self.unit_id,
            inputs=self.source_inputs, input_manifest=self.input_manifest, live=True,
        )
        if current_binding != self.consumer_binding:
            raise InputValidationError("Consumer source/control binding changed during execution.")
        if self.track == "benchmark" and execution_kind != current_binding['expected_execution_kind']:
            raise InputValidationError('Benchmark receipt execution kind contradicts its enabled engines.')
        if self.track == "benchmark" and {role:str(Path(path).absolute()) for role,path in outputs.items()} != current_binding["expected_outputs"]:
            raise InputValidationError("Benchmark receipt must retain all declared native metrics and assignment outputs.")
        if self.track == "review":
            if ({role:str(Path(path).absolute()) for role,path in outputs.items()} != current_binding["expected_outputs"]
                or (metadata or {}).get("human_review") != "NOT_REVIEWED"
                or (metadata or {}).get("review_ids") != current_binding["selection"]["review_ids"]):
                raise InputValidationError("Review receipt must retain every bound output and unreviewed candidate.")
        if not outputs or not tool_versions or validation.get("validated") is not True:
            raise OutputValidationError("Completion requires explicit outputs, versions and wrapper validation.")
        if any(not isinstance(key, str) or not _TOKEN.fullmatch(key) for key in outputs):
            raise InputValidationError("Unsafe receipt output role.")
        normalized_outputs = {role: Path(path).absolute() for role, path in outputs.items()}
        if not self.expected_outputs or normalized_outputs != self.expected_outputs:
            raise OutputValidationError("Completion outputs differ from the fresh pre-execution declaration.")
        if not required_output_roles(self.track).issubset(outputs):
            raise OutputValidationError("Receipt is missing required output roles for this module.")
        if execution_kind == "IN_PROCESS" and self.track not in {"annotation.overlap", "cohort.sv", "cohort.tr", "benchmark"}:
            raise OutputValidationError("A native caller cannot be certified as an in-process operation.")
        if any(not isinstance(key, str) or not key.strip() or not isinstance(value, str) or not value.strip() for key, value in tool_versions.items()):
            raise InputValidationError("Receipt tool versions must be nonempty recorded text.")
        command_records = _commands(self.commands, execution_kind)
        if self.track.startswith("assembly_sv.") and identity.get("config", {}).get("assembly_sv", {}).get("input_manifest"):
            if len(self.samples) != 1:
                raise InputValidationError("Ready assembly-SV receipts require exactly one sample.")
            from hifivar.assembly_source_roles import _verify_binding
            _verify_binding((metadata or {}).get("assembly_role_binding"), self.samples[0], self.context_path)
        artifacts = {role: artifact_identity(_path(Path(path).absolute(), root=root)) for role, path in outputs.items()}
        if len({item["path"] for item in artifacts.values()}) != len(artifacts):
            raise OutputValidationError("Multiple output roles alias the same artifact.")
        inodes = [(Path(item["path"]).stat().st_dev, Path(item["path"]).stat().st_ino) for item in artifacts.values()]
        if len(set(inodes)) != len(inodes):
            raise OutputValidationError("Multiple output roles hardlink the same artifact.")
        record = {
            "schema_version": 1, "kind": "hifivar_step_execution_receipt", "status": "COMPLETE",
            "attempt_id": context["attempt_id"], "producer_attempt_id": context["attempt_id"],
            "execution_kind": execution_kind, "track": self.track, "samples": list(self.samples),
            "unit_id": self.unit_id,
            "context_sha256": artifact_identity(self.context_path)["sha256"],
            "manifest_sha256": context["manifest_sha256"], "runtime": context["runtime"],
            "consumer_binding": current_binding,
            "commands": command_records, "outputs": artifacts,
            "tool_versions": dict(tool_versions), "validation": dict(validation),
            "version_interpretation": "RECORDED_WITH_EXPLICIT_SOURCES_NOT_AUTOMATICALLY_PROBED",
            "metadata": redact_sensitive_data(standardize_data(metadata or {})),
            "created": utc_now_iso8601(), "scientific_accuracy": "NOT_ASSESSED",
            "provenance_class": "THIS_ATTEMPT_EXECUTED_NOT_REUSED",
        }
        target = self.context_path.parent / "receipts" / receipt_name(self.track, self.samples, self.unit_id)
        _write_new(target, redact_sensitive_data(record))
        self.completed = True


@contextmanager
def attempt_step(track: str, samples: Sequence[str], *, outputs: Mapping[str, Path] | None = None,
                 unit_id: str | None = None, inputs: Mapping[str, Mapping[str, Path | None]] | None = None,
                 input_manifest: Path | None = None) -> Iterator[StepReceipt]:
    """Observe an owned bridge; independent Snakemake runs remain unlinked."""
    samples = _samples(samples)
    receipt_name(track, samples, unit_id)
    path = os.environ.get(ATTEMPT_CONTEXT_ENV)
    if not path:
        yield StepReceipt(None, track, samples, unit_id=unit_id)
        return
    context_path = Path(path)
    context, identity, root = _context(context_path, live=True)
    _allowed(identity, track, samples)
    _stage(context_path, identity, track, live=True)
    binding = _consumer_binding(context_path, context, identity, root, track, samples, unit_id,
                                inputs=inputs, input_manifest=input_manifest, live=True)
    if not outputs or any(not isinstance(role, str) or not _TOKEN.fullmatch(role) for role in outputs):
        raise InputValidationError("An owned step requires explicit expected output roles before execution.")
    expected = {role: _path(Path(value).absolute(), root=root) for role, value in outputs.items()}
    if any(not any(path.is_relative_to(root / area) for area in ("results", "work")) for path in expected.values()):
        raise InputValidationError("Step outputs must be under owned results/work, not controls or logs.")
    if len(set(expected.values())) != len(expected) or any(path.exists() for path in expected.values()):
        raise OutputValidationError("Owned-step outputs must be fresh and distinct; reuse needs a separate receipt.")
    with observe_commands() as commands:
        if track.startswith("annotation.") and unit_id not in {"small", "sv", "tr"}:
            raise InputValidationError("Annotation receipts require an explicit variant-category unit.")
        receipt = StepReceipt(context_path, track, samples, commands, expected_outputs=expected, unit_id=unit_id,
                              source_inputs=inputs, input_manifest=input_manifest, consumer_binding=binding)
        yield receipt
        if not receipt.completed:
            raise OutputValidationError("Owned workflow bridge did not seal its execution receipt.")


def read_step_receipt(root: Path, attempt_id: str, track: str, samples: Sequence[str], *, verify_outputs: bool = True, unit_id: str | None = None) -> dict | None:
    """Read one explicitly selected unit; no latest-attempt or output globbing."""
    root = _path(root)
    if not isinstance(attempt_id, str) or not _TOKEN.fullmatch(attempt_id):
        raise InputValidationError("Invalid receipt attempt id.")
    samples = _samples(samples)
    context_path = _path(root / "attempts" / attempt_id / "execution_context.json", root=root)
    if not context_path.exists():
        return None
    context, identity, _ = _context(context_path, live=False)
    if track == "alignment":
        if unit_id is not None or len(samples) != 1:
            raise InputValidationError("Alignment preparation requires exactly one sample and no category unit.")
        if context.get("alignment_handoff") is None:
            return None
        from hifivar.alignment_preparation import read_preparation
        return read_preparation(root, context_path.parent, samples[0], verify_outputs=verify_outputs)
    target = context_path.parent / "receipts" / receipt_name(track, samples, unit_id)
    if not target.exists():
        return None
    record = _read(target)
    _allowed(identity, track, samples)
    expected = {
        "schema_version": 1, "kind": "hifivar_step_execution_receipt", "status": "COMPLETE",
        "attempt_id": attempt_id, "producer_attempt_id": attempt_id, "track": track,
        "samples": list(samples), "runtime": context["runtime"],
        "unit_id": unit_id,
        "manifest_sha256": context["manifest_sha256"],
        "context_sha256": artifact_identity(context_path)["sha256"],
        "scientific_accuracy": "NOT_ASSESSED", "provenance_class": "THIS_ATTEMPT_EXECUTED_NOT_REUSED",
    }
    if context["dry_run"] or any(record.get(key) != value for key, value in expected.items()):
        raise InputValidationError("Receipt is not linked to the selected native attempt.")
    binding = record.get("consumer_binding")
    if _consumer(track) or track == "sv.harmonization":
        if not isinstance(binding, dict) or not isinstance(binding.get("source_inputs"), dict):
            raise InputValidationError("Consumer receipt lacks an independently verified source/control binding.")
        source_inputs = {}
        for sample, roles in binding["source_inputs"].items():
            if not isinstance(roles, dict):
                raise InputValidationError("Consumer receipt input witness is malformed.")
            source_inputs[sample] = {role: Path(_recorded_artifact(info)["path"]) if info is not None else None
                                     for role, info in roles.items()}
        manifest_input = Path(_recorded_artifact(binding.get("input_manifest"))["path"])
        current = _consumer_binding(context_path, context, identity, root, track, samples, unit_id,
                                    inputs=source_inputs, input_manifest=manifest_input, live=False)
        if binding != current:
            raise InputValidationError("Consumer receipt source/control binding differs from the selected stage.")
        if track == 'benchmark' and record.get('execution_kind') != current['expected_execution_kind']:
            raise InputValidationError('Benchmark receipt cannot relabel native/in-process comparison.')
        if track == "benchmark" and {role:item.get("path") for role,item in record.get("outputs",{}).items()} != current["expected_outputs"]:
            raise InputValidationError("Benchmark receipt omitted a declared native metric/assignment artifact.")
        if track == "review":
            if ({role:item.get("path") for role,item in record.get("outputs",{}).items()} != current["expected_outputs"]
                or record.get("metadata",{}).get("human_review") != "NOT_REVIEWED"
                or record.get("metadata",{}).get("review_ids") != current["selection"]["review_ids"]):
                raise InputValidationError("Review receipt omitted a selected output or changed its human status.")
    elif binding is not None:
        raise InputValidationError("Producer receipt unexpectedly contains a consumer binding.")
    _commands(record.get("commands", ()), record.get("execution_kind", ""), reading=True)
    if track.startswith("assembly_sv.") and identity.get("config", {}).get("assembly_sv", {}).get("input_manifest"):
        if len(samples) != 1:
            raise InputValidationError("Ready assembly-SV receipts require exactly one sample.")
        from hifivar.assembly_source_roles import _verify_binding
        _verify_binding(record.get("metadata", {}).get("assembly_role_binding"), samples[0], context_path)
    if record.get("execution_kind") == "IN_PROCESS" and track not in {"annotation.overlap", "cohort.sv", "cohort.tr", "benchmark"}:
        raise InputValidationError("A native caller cannot have an in-process receipt.")
    if not isinstance(record.get("validation"), dict) or record["validation"].get("validated") is not True:
        raise InputValidationError("Receipt lacks explicit wrapper validation.")
    versions = record.get("tool_versions")
    if not isinstance(versions, dict) or not versions or any(
        not isinstance(key, str) or not key.strip() or not isinstance(value, str) or not value.strip()
        for key, value in versions.items()
    ):
        raise InputValidationError("Receipt lacks explicit tool versions.")
    if not isinstance(record.get("outputs"), dict) or not record["outputs"]:
        raise InputValidationError("Receipt lacks outputs.")
    if not required_output_roles(track).issubset(record["outputs"]):
        raise InputValidationError("Receipt is missing required output roles for this module.")
    for role, item in record["outputs"].items():
        if not isinstance(role, str) or not _TOKEN.fullmatch(role) or not isinstance(item, dict) or not isinstance(item.get("sha256"), str) or not _HASH.fullmatch(item["sha256"]):
            raise InputValidationError("Receipt output identity is malformed.")
        path = _path(item.get("path", ""), root=root)
        if not any(path.is_relative_to(root / area) for area in ("results", "work")):
            raise InputValidationError("Receipt output is outside owned results/work.")
        if verify_outputs and artifact_identity(path) != item:
            raise OutputValidationError("Selected-attempt output bytes differ from their execution receipt.")
    paths = [Path(item["path"]) for item in record["outputs"].values()]
    if len(set(paths)) != len(paths):
        raise InputValidationError("Receipt output roles alias a path.")
    if verify_outputs:
        inodes = [(path.stat().st_dev, path.stat().st_ino) for path in paths]
        if len(set(inodes)) != len(inodes):
            raise InputValidationError("Receipt output roles hardlink the same artifact.")
    return record


__all__ = ["ATTEMPT_CONTEXT_ENV", "artifact_identity", "attempt_step", "create_attempt_context", "read_step_receipt", "receipt_name", "required_output_roles"]
