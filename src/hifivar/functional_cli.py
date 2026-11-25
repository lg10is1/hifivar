"""Explicit user-key AlphaGenome plan execution after local REF validation."""
from __future__ import annotations
from dataclasses import asdict
import hashlib
from pathlib import Path
import yaml
from hifivar import __version__
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.functional import read_functional_selection
from hifivar.serialization import write_json_atomic, utc_now_iso8601, standardize_data

class _UniqueLoader(yaml.SafeLoader):
    """Reject duplicate or non-string keys before authorization is interpreted."""

def _unique_mapping(loader, node):
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=True)
        if not isinstance(key, str) or key in result:
            raise InputValidationError("AlphaGenome YAML contains duplicate or non-string mapping keys.")
        result[key] = loader.construct_object(value_node, deep=True)
    return result

_UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _unique_mapping)

def _mapping(value, allowed, name):
    if not isinstance(value, dict) or set(value) - set(allowed):
        raise InputValidationError(f"{name} requires a mapping with only documented fields.")
    return value

def _path(value, base, name):
    if not isinstance(value, str) or not value.strip():
        raise InputValidationError(f"{name} requires a path.")
    p = Path(value).expanduser()
    if not p.is_absolute():
        p = base / p
    p = p.resolve()
    if not p.is_file():
        raise InputValidationError(f"{name} must be an existing regular file.")
    return p

def _sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for data in iter(lambda: f.read(1024 * 1024), b""):
            h.update(data)
    return h.hexdigest()

def _stat(path):
    s = path.stat()
    return (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)

def _same_inputs(before):
    try:
        return all(_stat(p) == state for p, state in before.items())
    except OSError:
        return False

def _validate_ref_alleles(fasta, fai, variants, *, sequence_length=None):
    entries = {}
    try:
        for line in fai.read_text(encoding="ascii").splitlines():
            fields = line.split("\t")
            if len(fields) < 5 or fields[0] in entries:
                raise ValueError
            length, offset, bases, width = map(int, fields[1:5])
            if length <= 0 or offset < 0 or bases <= 0 or width < bases:
                raise ValueError
            entries[fields[0]] = (length, offset, bases, width)
        verified = []
        with fasta.open("rb") as handle:
            for variant in variants:
                if variant.contig not in entries:
                    raise InputValidationError("Selected contig is absent from the local FASTA index; contigs are never renamed.")
                length, offset, bases, width = entries[variant.contig]
                pos = variant.position - 1
                remaining = len(variant.reference_bases)
                if pos < 0 or pos + remaining > length:
                    raise InputValidationError("Selected REF extends beyond its indexed contig.")
                context = None
                if sequence_length is not None:
                    # Official SDK v0.10.0 positive/unstranded Interval.center/resize.
                    center = pos + (remaining + 1) // 2
                    context = (center - (sequence_length + 1) // 2, center + sequence_length // 2)
                    if context[0] < 0 or context[1] > length:
                        raise InputValidationError("The full model context exceeds the indexed contig; padded inference is not enabled.")
                pieces = []
                while remaining:
                    within = pos % bases
                    take = min(remaining, bases - within)
                    handle.seek(offset + (pos // bases) * width + within)
                    part = handle.read(take).decode("ascii").upper()
                    if len(part) != take:
                        raise InputValidationError("Unable to read the complete indexed REF allele.")
                    pieces.append(part)
                    pos += take
                    remaining -= take
                if "".join(pieces) != variant.reference_bases.upper():
                    raise InputValidationError("Selected REF does not match the local reference; no normalization or build conversion was performed.")
                verified.append({"sample": variant.sample_id, "source_variant_id": variant.source_variant_id,
                                 "contig": variant.contig, "position": variant.position,
                                 "ref_length": len(variant.reference_bases), "status": "MATCH",
                                 "model_context_0based_half_open": list(context) if context else None})
        return verified
    except InputValidationError:
        raise
    except (OSError, UnicodeError, ValueError):
        raise InputValidationError("Unable to validate selected REF alleles against the existing FASTA/FAI.") from None

def run_alphagenome_plan(plan_path, output_dir, *, dry_run=False, allow_network=False, backend=None):
    """Execute an explicit selection only; no credential in plan, argv or report."""
    from hifivar.alphagenome_backend import AlphaGenomeSettings, AlphaGenomeRequest, AlphaGenomeBackend
    if not isinstance(dry_run, bool) or not isinstance(allow_network, bool):
        raise InputValidationError("Execution flags must be booleans.")
    plan = Path(plan_path).expanduser().resolve()
    if not plan.is_file():
        raise InputValidationError("AlphaGenome plan must be an existing file.")
    try:
        plan_bytes = plan.read_bytes()
        plan_digest = hashlib.sha256(plan_bytes).hexdigest()
        raw = yaml.load(plan_bytes.decode("utf-8"), Loader=_UniqueLoader)
    except (OSError, UnicodeError, yaml.YAMLError):
        raise InputValidationError("Unable to read the AlphaGenome YAML plan.") from None
    raw = _mapping(raw, {"schema_version", "selection_file", "reference", "alphagenome", "data_sharing_consent"}, "AlphaGenome plan")
    if type(raw.get("schema_version")) is not int or raw["schema_version"] != 1:
        raise InputValidationError("AlphaGenome schema_version must be integer 1.")
    consent = raw.get("data_sharing_consent", False)
    if not isinstance(consent, bool):
        raise InputValidationError("data_sharing_consent must be boolean.")
    ref = _mapping(raw.get("reference"), {"fasta", "build"}, "reference")
    if ref.get("build") != "GRCh38":
        raise InputValidationError("This AlphaGenome route requires explicit human GRCh38; no liftover is performed.")
    settings_raw = _mapping(raw.get("alphagenome"), {
        "model_version", "requested_modalities", "ontology_terms", "sequence_length", "max_variants",
        "max_rows_per_variant", "max_allele_length", "credential_env", "connect_timeout_seconds", "rpc_timeout_seconds"
    }, "alphagenome")
    if not settings_raw.get("model_version") or not settings_raw.get("requested_modalities"):
        raise InputValidationError("Explicit AlphaGenome model_version and requested_modalities are required.")
    for key in ("requested_modalities", "ontology_terms"):
        if key in settings_raw and not isinstance(settings_raw[key], list):
            raise InputValidationError(f"alphagenome.{key} must be a list.")
    kw = dict(settings_raw)
    for key in ("requested_modalities", "ontology_terms"):
        if key in kw:
            kw[key] = tuple(kw[key])
    settings = AlphaGenomeSettings(reference_build=ref["build"], **kw)
    fasta = _path(ref.get("fasta"), plan.parent, "reference.fasta")
    fai = Path(str(fasta) + ".fai")
    if not fai.is_file():
        raise InputValidationError("An existing reference.fasta.fai is required; no index is generated.")
    selection = _path(raw.get("selection_file"), plan.parent, "selection_file")
    before = {p: _stat(p) for p in (plan, selection, fasta, fai)}
    selection_digest, fai_digest = _sha(selection), _sha(fai)
    variants = read_functional_selection(selection)
    request = AlphaGenomeRequest(variants, settings)
    before.update({v.source_annotation: _stat(v.source_annotation) for v in variants})
    protected = tuple(before)
    def unchanged():
        return (_same_inputs(before) and _sha(plan) == plan_digest
                and _sha(selection) == selection_digest and _sha(fai) == fai_digest)
    verified = _validate_ref_alleles(fasta, fai, variants, sequence_length=settings.sequence_length)
    root = Path(output_dir).expanduser().resolve()
    if root.exists() or any(root == p or root in p.parents for p in protected):
        raise OutputValidationError("AlphaGenome output directory must be new and must not contain protected inputs.")
    installation = Path(__file__).resolve().parent
    if root == installation or installation in root.parents:
        raise OutputValidationError("AlphaGenome outputs cannot be written into the package installation.")
    report = {
        "schema_version": 1, "hifivar_version": __version__, "created_at": utc_now_iso8601(),
        "status": "PLANNED" if dry_run else "RUNNING",
        "plan_sha256": plan_digest, "selection_sha256": selection_digest,
        "reference": {"build": "GRCh38", "fasta": str(fasta), "fai_sha256": fai_digest,
                      "local_selected_ref_check": verified, "whole_fasta_sha256": "NOT_COMPUTED"},
        "settings": standardize_data(asdict(settings), context="AlphaGenome settings"), "selected_variant_count": len(variants),
        "data_sharing_consent": consent, "network_execution": "NOT_RUN_DRY_RUN" if dry_run else "NOT_STARTED",
        "credential_value_recorded": False,
        "scientific_policy": {"functional_impact_is_call_confidence": False,
                              "functional_impact_is_pathogenicity": False,
                              "raw_calls_modified": False},
    }
    if not unchanged():
        raise InputValidationError("A protected input changed during planning.")
    if dry_run:
        return report
    if not allow_network or not consent:
        raise InputValidationError("Real AlphaGenome requests require --allow-network and data_sharing_consent: true. Selected variant coordinates and alleles are sent to the official service.")
    root.mkdir(parents=True, exist_ok=False)
    destination = root / "alphagenome.result.json"
    write_json_atomic(report, destination, overwrite=False, artifact_name="AlphaGenome execution")
    try:
        report["network_execution"] = "BACKEND_DISPATCH_ATTEMPTED"
        result = (backend or AlphaGenomeBackend()).score(request, allow_network=True)
        if not unchanged():
            raise InputValidationError("A protected input changed during AlphaGenome execution.")
        payload = result.to_dict()
        if not isinstance(payload, dict):
            raise InputValidationError("AlphaGenome backend returned an invalid result.")
        report.update(status="COMPLETE", network_execution="BACKEND_RETURNED_SUCCESS", result=payload,
                      source_immutability="UNCHANGED_STAT_AND_SELECTED_REF_VERIFIED_NOT_FULL_REFERENCE_DIGEST")
    except Exception as error:
        # Never persist provider/SDK exception text or cause chains.
        report.update(status="FAILED", error_type=type(error).__name__,
                      error="AlphaGenome request did not complete; no successful prediction is asserted.",
                      source_immutability="UNCHANGED_STAT_AND_SMALL_INPUT_DIGESTS" if unchanged() else "FAILED")
        write_json_atomic(report, destination, overwrite=True, artifact_name="AlphaGenome execution")
        raise InputValidationError("AlphaGenome execution failed; a sanitized failure report was retained.") from None
    report["completed_at"] = utc_now_iso8601()
    write_json_atomic(report, destination, overwrite=True, artifact_name="AlphaGenome execution")
    return report
