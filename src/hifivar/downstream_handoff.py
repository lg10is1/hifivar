"""Explicit post-calling handoff from attempt receipts to downstream manifests.

Plans name a producer, never a guessed future filename. Only a validated
receipt from the selected attempt may supply a generated path and version.
Small variants, SV and TR remain independent artifacts throughout.
"""

from __future__ import annotations

import copy
import csv
import gzip
import io
from dataclasses import asdict
from pathlib import Path

import yaml

from hifivar.attempt_receipts import artifact_identity, read_step_receipt, receipt_name
from hifivar.downstream_plan import DownstreamPlanRow, read_downstream_plan
from hifivar.exceptions import ConfigurationError, InputValidationError, OutputValidationError
from hifivar.serialization import standardize_data, write_json_atomic, write_yaml_atomic


def _manifest(fields: tuple[str, ...], rows: list[dict]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def _plan_document(plan) -> dict:
    return {
        "schema_version": 1,
        "annotation": [standardize_data(asdict(row)) for row in plan.annotation],
        "cohort": [standardize_data(asdict(row)) for row in plan.cohort],
    }


def prepare_downstream_intent(config: dict, context, base: Path, root: Path) -> tuple[dict[str, str], list[Path], dict]:
    """Validate/freeze explicit selectors before creating any run directory."""
    from hifivar.launcher import _absolute, _audit_ready_gvcf, _audit_ready_region_bed, _ready_vcf, _vcf_header_markers
    supplied = config.get("workflow", {}).get("downstream_source_plan")
    if not supplied:
        raise ConfigurationError("No downstream_source_plan was provided.")
    if any(config.get(group, {}).get("input_manifest") for group in ("annotation", "cohort")):
        raise ConfigurationError("Use a downstream source plan OR legacy ready manifests, not both.")
    source = Path(_absolute(supplied, base))
    for key in ("catalog",):
        if config.get("tr", {}).get(key):
            config["tr"][key] = _absolute(config["tr"][key], base)
    plan = read_downstream_plan(source, sample_ids=context.sample_ids, reference_build=context.reference.build, config=config)
    if config.get("annotation", {}).get("enabled") and not plan.annotation:
        raise InputValidationError("Enabled annotation requires explicit source-plan rows.")
    inputs: list[Path] = [source]
    if config.get("tr", {}).get("catalog"):
        inputs.append(Path(config["tr"]["catalog"]))
    audits = []
    annotation = config.get("annotation", {})
    for key in ("annovar_database_root", "vep_cache_directory", "gene_bed", "exon_bed", "regulatory_bed", "repeat_bed", "segdup_bed"):
        if annotation.get(key):
            annotation[key] = _absolute(annotation[key], base)
            if key.endswith("_bed"):
                inputs.append(Path(annotation[key]))
                if annotation.get("overlap_enabled"):
                    _audit_ready_region_bed(Path(annotation[key]), context)
    for group, rows in (("annotation", plan.annotation), ("cohort", plan.cohort)):
        for row in rows:
            if group == "annotation" and annotation.get("overlap_enabled") and row.track in {"sv", "tr"} and not row.variant_ids:
                raise InputValidationError("Region-overlap source plans require explicit candidate variant_ids.")
            if row.path is None:
                continue
            if not str(row.path).lower().endswith(".vcf.gz") or row.index_path != Path(str(row.path) + ".tbi"):
                raise InputValidationError("Downstream sources require BGZF VCF and its conventional TBI.")
            header = _ready_vcf(row.path, row.sample_id, context)
            if group == "annotation":
                from hifivar.annotation import AnnotationInput, VariantCategory
                item = AnnotationInput(row.sample_id, VariantCategory(row.track), row.path,
                                       row.source_tool, context.reference, row.variant_ids,
                                       bnd_bundle_manifest=row.bnd_bundle_manifest)
                if item.bnd_bundle_manifest is not None:
                    binding = dict(item.source_identity)["bnd_bundle_binding"]
                    inputs.append(Path(binding["manifest"]["path"]))
                    inputs.extend(Path(info["path"]) for info in binding["artifacts"].values())
            if group == "cohort" and row.track == "small_variants":
                _audit_ready_gvcf(row.path, context)
            if group == "cohort" and row.track in {"sv", "tr"}:
                required = {"sv": {"SVTYPE"}, "tr": {"TRID", "MOTIFS", "STRUC"}}[row.track]
                if not required.issubset(_vcf_header_markers(row.path)):
                    raise InputValidationError("External cohort source lacks category-specific header declarations.")
            inputs.extend((row.path, row.index_path))
            audits.append({"sample_id": row.sample_id, "track": row.track, "header_audit": header,
                           "version_source": "USER_DECLARED_NOT_REPROBED"})
    config["workflow"]["downstream_source_plan"] = str(root / "downstream_plan.yaml")
    controls = {"downstream_plan.yaml": yaml.safe_dump(_plan_document(plan), allow_unicode=True, sort_keys=False)}
    return controls, inputs, {
        "scope": "EXPLICIT_GENERATED_OR_EXTERNAL_SOURCE_INTENT_NOT_CALLED_OUTPUTS",
        "original_plan_sha256": plan.sha256, "external_inputs": audits,
        "generated_rows": sum(row.origin == "GENERATED" for row in (*plan.annotation, *plan.cohort)),
        "accuracy": "NOT_ASSESSED", "functional_backend": "USER_DEFERRED",
    }


def producer_config(config: dict) -> dict:
    """Disable downstream consumers only; calling algorithms stay unchanged."""
    result = copy.deepcopy(config)
    for group in ("annotation", "cohort"):
        result[group]["enabled"] = False
    result["workflow"]["downstream_source_plan"] = None
    return result


def consumer_config(config: dict) -> dict:
    """Disable every producer, including independently switched harmonization.

    Disabling only sv.enabled is insufficient: standalone harmonization has its
    own enabled flag and representation guards. The ready bundle remains an
    explicit annotation input, not a request to produce it again. Original
    producer config/policies stay frozen in the launch identity and receipts.
    """
    result = copy.deepcopy(config)
    for group in ("alignment", "small", "sv", "tr", "phasing", "assembly", "assembly_sv"):
        result[group]["enabled"] = False
    # Ready FASTAs belong to the frozen producer inputs, not a consumer-stage
    # request to perform assembly-SV again. Original launch identity is retained.
    result["assembly_sv"]["input_manifest"] = None
    # Review is a separately bound final consumer stage, never reactivated by
    # annotation/cohort materialization or an unresolved source plan.
    result["review"].update(enabled=False, source_plan=None, selection_file=None)
    result["benchmark"].update(enabled=False, source_plan=None)
    for track in ("small_variants", "sv", "assembly_sv", "tr"):
        result["benchmark"][track]["enabled"] = False
    integration = result["sv"].setdefault("harmonization", {})
    integration["enabled"] = False
    integration["bnd_representation"] = {"enabled": False}
    integration["input_identity"] = {"enabled": False}
    result["workflow"]["downstream_source_plan"] = None
    return result


def _generated_source(root: Path, attempt: Path, row: DownstreamPlanRow, reference) -> tuple[Path, Path, str, str, dict, Path | None]:
    selectors = {
        "small.deepvariant": ("small", "deepvariant"),
        "tr.trgt": ("tr", "trgt"), "phasing.hiphase": ("phasing", "hiphase"),
        **{"sv." + name: ("sv." + name, "jasmine" if name == "harmonization" else name)
           for name in ("sawfish", "sniffles2", "pbsv", "cutesv", "harmonization")},
        **{"assembly_sv." + name: ("assembly_sv." + name, name) for name in ("pav", "svim_asm")},
    }
    track, tool = selectors[row.producer]
    receipt = read_step_receipt(root, attempt.name, track, (row.sample_id,))
    if receipt is None:
        raise InputValidationError(f"Generated downstream source has no selected-attempt execution receipt: {track}/{row.sample_id}.")
    role = row.output_role
    index_role = {"gvcf": "gvcf_index", "bnd_vcf": "bnd_index"}.get(role, "index")
    artifacts = receipt["outputs"]
    if role not in artifacts or index_role not in artifacts:
        raise OutputValidationError("Generated source receipt lacks its explicitly selected VCF/index roles.")
    version_keys = [key for key in receipt["tool_versions"] if key.casefold().replace("-", "_") == tool]
    if len(version_keys) != 1:
        raise InputValidationError("Generated source has no unique recorded producer version.")
    if receipt["metadata"].get("reference_build") != row.reference_build:
        raise InputValidationError("Generated source receipt reference build differs from the source plan.")
    if row.track == "tr" and receipt["metadata"].get("catalog_id") != row.catalog_id:
        raise InputValidationError("Generated TR catalog checksum differs from the source plan.")
    bundle_manifest = None
    if role == "bnd_vcf":
        from hifivar.downstream_bnd import generated_bnd_roles
        bound = generated_bnd_roles(receipt, reference, row.sample_id)
        bundle_manifest = Path(bound["bnd_bundle_manifest"]["path"])
    receipt_path = attempt / "receipts" / receipt_name(track, (row.sample_id,))
    linkage = {
        **artifact_identity(receipt_path), "producer_attempt_id": receipt["producer_attempt_id"],
        "track": track, "samples": [row.sample_id], "output_role": role, "index_role": index_role,
        "bnd_bundle_manifest_role": "bnd_manifest" if bundle_manifest else None,
        "source_tool": tool, "source_version": receipt["tool_versions"][version_keys[0]],
        "version_source": receipt["metadata"].get("tool_version_sources", {}).get(version_keys[0], "WRAPPER_RECORDED_NOT_INDEPENDENTLY_REPROBED"),
    }
    return (Path(artifacts[role]["path"]), Path(artifacts[index_role]["path"]), tool,
            linkage["source_version"], linkage, bundle_manifest)


def _has_records(path: Path) -> bool:
    """Distinguish a valid header-only source from an absent/failed artifact."""
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        while True:
            line = handle.readline(1024 * 1024 + 1)
            if not line:
                return False
            if len(line) > 1024 * 1024:
                raise InputValidationError("VCF record exceeds the bounded 1 MiB handoff line limit.")
            if line and not line.startswith("#") and line.strip():
                return True


def materialize_downstream(config: dict, context, root: Path, attempt: Path, *, calling_config: dict | None = None) -> tuple[Path, dict]:
    """Create ready manifests only after every generated producer is proven."""
    from hifivar.launcher import _prepare_ready_downstream
    plan_path = root / "downstream_plan.yaml"
    plan = read_downstream_plan(plan_path, sample_ids=context.sample_ids, reference_build=context.reference.build, config=config)
    destination = attempt / "downstream"
    if destination.exists():
        raise OutputValidationError("Downstream materialization already exists; it is not silently replaced.")
    annotation_rows, cohort_rows, links = [], [], []
    for group, rows in (("annotation", plan.annotation), ("cohort", plan.cohort)):
        for row in rows:
            path, index, tool, version = row.path, row.index_path, row.source_tool, row.source_version
            state = row.state
            bundle_manifest = row.bnd_bundle_manifest
            if row.origin == "GENERATED":
                path, index, tool, version, link, bundle_manifest = _generated_source(root, attempt, row, context.reference)
                links.append(link)
                if group == "cohort":
                    state = "CALLED" if _has_records(path) else "NO_CALLS"
            if group == "annotation":
                annotation_rows.append({"sample": row.sample_id, "variant_category": row.track,
                                        "source_vcf": str(path), "source_tool": tool,
                                        "source_variant_ids": ";".join(row.variant_ids)})
                if bundle_manifest is not None:
                    annotation_rows[-1]["bnd_bundle_manifest"] = str(bundle_manifest)
            else:
                cohort_rows.append({"sample": row.sample_id, "track": row.track, "state": state,
                                    "source_path": str(path) if path else "", "index_path": str(index) if index else "",
                                    "source_tool": tool or "", "source_version": version or "",
                                    "reference_build": row.reference_build, "catalog_id": row.catalog_id or ""})
    # No writes happen until all selected producer receipts have passed.
    destination.mkdir(parents=True, exist_ok=False)
    effective = copy.deepcopy(calling_config or config)
    for group in ("annotation", "cohort"):
        effective[group] = copy.deepcopy(config[group])
    annotation_fields = ("sample", "variant_category", "source_vcf", "source_tool", "source_variant_ids")
    if any("bnd_bundle_manifest" in row for row in annotation_rows):
        annotation_fields += ("bnd_bundle_manifest",)
    for group, fields, rows in (
        ("annotation", annotation_fields, annotation_rows),
        ("cohort", ("sample", "track", "state", "source_path", "index_path", "source_tool", "source_version", "reference_build", "catalog_id"), cohort_rows),
    ):
        if effective[group].get("enabled"):
            request_path = destination / (group + ".request.tsv")
            with request_path.open("x", encoding="utf-8", newline="") as handle:
                handle.write(_manifest(fields, rows))
            effective[group]["input_manifest"] = str(request_path)
    controls, _, audit = _prepare_ready_downstream(effective, context, root, root)
    generated_cohort = {(row.sample_id, row.track): row for row in plan.cohort if row.origin == "GENERATED"}
    for item in audit.get("cohort", []):
        if (item["sample_id"], item["track"]) in generated_cohort:
            item["ready_validator_version_label"] = item["version_source"]
            item["version_source"] = "SELECTED_ATTEMPT_PRODUCER_RECEIPT"
            item["state_source"] = "VALIDATED_GENERATED_VCF_RECORD_PRESENCE_NOT_BIOLOGICAL_ABSENCE"
            item["source_origin"] = "GENERATED"
        else:
            item["source_origin"] = "EXTERNAL"
    generated_annotation = {(row.sample_id, row.track) for row in plan.annotation if row.origin == "GENERATED"}
    for item in audit.get("annotation", []):
        item["source_origin"] = "GENERATED" if (item["sample_id"], item["category"]) in generated_annotation else "EXTERNAL"
    owned_controls = {}
    for name, text in controls.items():
        target = destination / name
        with target.open("x", encoding="utf-8", newline="") as handle:
            handle.write(text)
        owned_controls[name] = artifact_identity(target)
    for group in ("annotation", "cohort"):
        if effective[group].get("enabled"):
            effective[group]["input_manifest"] = str(destination / (group + ".inputs.tsv"))
    effective = consumer_config(effective)
    target_config = destination / "effective.yaml"
    write_yaml_atomic(effective, target_config)
    owned_controls["effective.yaml"] = artifact_identity(target_config)
    record = {
        "schema_version": 1, "kind": "hifivar_downstream_materialization", "attempt_id": attempt.name,
        "context_sha256": artifact_identity(attempt / "execution_context.json")["sha256"],
        "plan_sha256": artifact_identity(plan_path)["sha256"], "controls": owned_controls,
        "source_receipts": links, "external_source_inputs": "FROZEN_LAUNCH_MANIFEST_INPUT_IDENTITIES",
        "source_validation": audit, "accuracy": "NOT_ASSESSED",
        "ordinary_vcf_substituted_for_gvcf": False, "implicit_caller_choice": False,
        "cross_sample_sv_clustering": False,
    }
    write_json_atomic(record, destination / "materialization.json")
    return target_config, record


__all__ = ["prepare_downstream_intent", "producer_config", "consumer_config", "materialize_downstream"]
