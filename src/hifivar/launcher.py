"""Bounded Linux launcher with optional raw-input alignment preparation.

It delegates alignment and calling, never tunes callers, recovers stale locks,
or certifies biological accuracy. It owns an explicitly selected run directory.
"""
from __future__ import annotations

import csv
import copy
import gzip
import hashlib
import io
import json
import os
import platform
import re
import shlex
import sys
from pathlib import Path
from uuid import uuid4

from hifivar import __version__
from hifivar.command import CommandRunner
from hifivar.attempt_receipts import ATTEMPT_CONTEXT_ENV, create_attempt_context
from hifivar.config import HiFiVarConfig, validate_config
from hifivar.context import AnalysisContext
from hifivar.exceptions import ConfigurationError, InputValidationError, OutputValidationError, WorkflowError
from hifivar.launch_inputs import IndexedBamPreflight, unique_bam_index
from hifivar.logging_utils import get_logger
from hifivar.package_resources import installed_workflow_root
from hifivar.sample import InputType
from hifivar.sample_sheet import SampleSheet
from hifivar.serialization import utc_now_iso8601, write_json_atomic, write_yaml_atomic

_LOGGER = get_logger(__name__)
_UNSUPPORTED = ("low_depth",)
_GVCF_RECORD_PREFIX_LIMIT = 4096
_GVCF_LINE_CHARACTER_LIMIT = 1024 * 1024


def file_identity(path: Path) -> dict[str, object]:
    """Stream a complete checksum and reject a file changing during the read."""
    path = path.resolve(strict=True)
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise InputValidationError(f"Input changed while hashing: {path}")
    return {"path": str(path), "bytes": after.st_size, "sha256": digest.hexdigest()}


def runtime_identity(workflow_root: Path) -> dict[str, object]:
    """Fingerprint actual Python/config/workflow bytes, not a version label alone."""
    paths = {}
    for prefix, root in (("package", Path(__file__).parent), ("workflow", workflow_root)):
        for path in sorted(root.rglob("*")):
            if path.is_file() and (path.suffix in {".py", ".yaml", ".smk", ".json"} or path.name == "Snakefile"):
                paths[f"{prefix}/{path.relative_to(root).as_posix()}"] = file_identity(path)["sha256"]
    digest = hashlib.sha256(json.dumps(paths, sort_keys=True).encode("utf-8")).hexdigest()
    return {"hifivar_version": __version__, "payload_sha256": digest,
            "python": str(Path(sys.executable).resolve()), "python_version": platform.python_version()}


def _absolute(value: str, base: Path) -> str:
    path = Path(value).expanduser()
    return str((path if path.is_absolute() else base / path).resolve())


def _normalize_executables(section: dict, base: Path) -> None:
    for key, value in section.items():
        if isinstance(value, dict):
            _normalize_executables(value, base)
        elif key.endswith("executable") and isinstance(value, str) and ("/" in value or "\\" in value):
            section[key] = _absolute(value, base)


def _enabled_tracks(config: dict) -> list[tuple[str, dict]]:
    from hifivar.joint_phasing import validate_phasing_producers
    validate_phasing_producers(config)
    for name in _UNSUPPORTED:
        if config.get(name, {}).get("enabled", False):
            raise ConfigurationError(f"hifivar run does not yet route {name}; use its documented workflow separately.")
    tracks = [(name, config[name]) for name in ("small", "tr", "phasing", "review") if config.get(name, {}).get("enabled", False)]
    if config.get("review", {}).get("enabled") and (not config["review"].get("source_plan") or config["review"].get("selection_file")):
        raise ConfigurationError("hifivar run review requires explicit review.source_plan; selection_file is standalone-only.")
    if config.get("assembly", {}).get("enabled", False):
        tracks.append(("assembly", config["assembly"]))
    annotation = config.get("annotation", {})
    if annotation.get("enabled"):
        if annotation.get("functional_enabled"):
            raise ConfigurationError("AlphaGenome is available through 'hifivar functional alphagenome' with an explicit plan and user credential; annotation.functional_enabled is not an automatic main-run route.")
        tracks.append(("annotation", annotation))
    cohort = config.get("cohort", {})
    if cohort.get("enabled"):
        if cohort.get("overwrite"):
            raise ConfigurationError("Launcher requires cohort.overwrite=false.")
        enabled_cohort = []
        for name in ("small_variants", "sv", "tr"):
            subsection = cohort.get(name, {})
            if not subsection.get("enabled"):
                continue
            values = dict(subsection)
            if name == "small_variants":
                if not isinstance(values.get("memory_gb"), int) or isinstance(values.get("memory_gb"), bool) or values["memory_gb"] <= 0:
                    raise ConfigurationError("cohort.small_variants.memory_gb must be an explicit positive integer.")
                values["memory_mb"] = values["memory_gb"] * 1024
            elif values["threads"] != 1:
                raise ConfigurationError(f"cohort.{name}.threads must be 1 for the native table builder.")
            enabled_cohort.append((f"cohort.{name}", values))
        if not enabled_cohort:
            raise ConfigurationError("cohort.enabled requires at least one enabled cohort track.")
        tracks.extend(enabled_cohort)
    if config.get("benchmark", {}).get("enabled"):
        if not config["benchmark"].get("source_plan") or config["benchmark"].get("overwrite"):
            raise ConfigurationError("hifivar run benchmark requires source_plan and overwrite=false; legacy query paths are standalone-only.")
        tracks.extend(("benchmark."+name, config["benchmark"][name]) for name in ("small_variants","sv","assembly_sv","tr") if config["benchmark"][name].get("enabled"))
    assembly_sv = config.get("assembly_sv", {})
    if assembly_sv.get("enabled"):
        if not config.get("assembly", {}).get("enabled") and not assembly_sv.get("input_manifest"):
            raise ConfigurationError("assembly_sv requires assembly.enabled and independent HiFi FASTQ input.")
        if assembly_sv.get("overwrite"):
            raise ConfigurationError("Launcher requires assembly_sv.overwrite=false.")
        assembly_callers = [(f"assembly_sv.{name}", assembly_sv[name]) for name in ("pav", "svim_asm") if assembly_sv.get(name, {}).get("enabled", False)]
        if not assembly_callers:
            raise ConfigurationError("assembly_sv.enabled requires at least one enabled assembly caller.")
        tracks.extend(assembly_callers)
    if config.get("phasing", {}).get("enabled") and not config.get("small", {}).get("enabled"):
        raise ConfigurationError("phasing requires small.enabled for this launcher.")
    sv = config.get("sv", {})
    if sv.get("enabled"):
        callers = [(f"sv.{name}", sv[name]) for name in ("sawfish", "sniffles2", "pbsv", "cutesv") if sv.get(name, {}).get("enabled", False)]
        if not callers:
            raise ConfigurationError("sv.enabled requires at least one enabled caller.")
        tracks.extend(callers)
    harmonization = sv.get("harmonization", {})
    if harmonization.get("enabled"):
        from hifivar.harmonization_inputs import generated_callers
        if not generated_callers(config) and not harmonization.get("input_vcfs"):
            raise ConfigurationError("Harmonization requires generated or external evidence.")
        if harmonization["distance_type"] not in {"linear", "nonlinear"}:
            raise ConfigurationError("sv.harmonization.distance_type must be linear or nonlinear.")
        minimum = harmonization.get("minimum_caller_support", 1)
        if type(minimum) is not int or not 1 <= minimum <= 6:
            raise ConfigurationError("sv.harmonization.minimum_caller_support must be an integer between 1 and 6.")
        tracks.append(("sv.harmonization", harmonization))
    if not tracks:
        raise ConfigurationError("Enable at least one calling track; an empty phase0 run is not batch analysis.")
    if sv.get("enabled") and sv.get("overwrite"):
        raise ConfigurationError("Launcher requires sv.overwrite=false.")
    for name, values in tracks:
        if values.get("overwrite"):
            raise ConfigurationError(f"Launcher requires {name}.overwrite=false.")
        if values["threads"] > config["runtime"]["threads"]:
            raise ConfigurationError(f"{name}.threads exceeds runtime.threads; no silent thread scaling.")
        if values["memory_mb"] > config["workflow"]["memory_mb"]:
            raise ConfigurationError(f"{name}.memory_mb exceeds workflow.memory_mb.")
    return tracks


def build_assembly_context(config, primary_context: AnalysisContext | None = None) -> AnalysisContext:
    """Bind original/explicit FASTQ to assembly independently of calling BAMs.

    No BAM/uBAM extraction or derived-read guessing is performed. An independent
    sample sheet is a user-declared dataset, not a verified biological identity.
    Its sample set and declared metadata must match the primary context exactly.
    """
    primary = primary_context or AnalysisContext.from_config(config)
    section = config.get("assembly", {})
    supplied = section.get("sample_sheet")
    if supplied:
        source = config.sources.get("user") if isinstance(config, HiFiVarConfig) else None
        value = Path(supplied).expanduser()
        if not value.is_absolute() and source is None:
            raise ConfigurationError("Relative assembly.sample_sheet requires a loaded user config source.")
        sheet = SampleSheet.from_tsv(value if value.is_absolute() else Path(source).resolve().parent / value)
        context = AnalysisContext(primary.reference, sheet.records, primary.config, sheet.source_path)
    else:
        context = primary
    if any(record.sample.input.input_type is not InputType.FASTQ for record in context.samples):
        raise InputValidationError(
            "Assembly requires primary HiFi FASTQ or an explicit assembly.sample_sheet "
            "containing FASTQ; implicit BAM/uBAM/CRAM extraction is disabled."
        )
    original = {record.sample.sample_id: record for record in primary.samples}
    observed = {record.sample.sample_id: record for record in context.samples}
    if original.keys() != observed.keys():
        raise InputValidationError("ASSEMBLY_SAMPLE_SET_MISMATCH: primary and assembly sample sets differ.")
    for sample_id, record in observed.items():
        expected = original[sample_id]
        if any(getattr(record, key) != getattr(expected, key) for key in ("sex", "father", "mother", "phenotype", "group")):
            raise InputValidationError(f"ASSEMBLY_SAMPLE_METADATA_MISMATCH: '{sample_id}' metadata differ.")
    seen: list[Path] = []
    for record in context.samples:
        for path in record.sample.input.files:
            resolved = path.resolve()
            if any(resolved.samefile(previous) for previous in seen):
                raise InputValidationError("Assembly FASTQ is assigned more than once, including hardlink/symlink aliases.")
            seen.append(resolved)
    return context


def _calling_enabled(config: dict) -> bool:
    return any(config.get(name, {}).get("enabled", False) for name in ("small", "sv", "tr", "phasing", "review"))


def _ready_vcf(path: Path, sample_id: str, context: AnalysisContext) -> dict:
    """Check a ready external VCF header/index, never infer future success."""
    from hifivar.annotation import audit_annotation_vcf_identity
    from hifivar.small import _validate_bgzf, _validate_tabix_index
    audit = audit_annotation_vcf_identity(path, sample_id, context.reference)
    if str(path).lower().endswith(".vcf.gz"):
        # Existing downstream DAGs declare conventional TBI dependencies; CSI
        # alone would pass some generic readers but would not satisfy the DAG.
        _validate_bgzf(path, "Downstream VCF")
        _validate_tabix_index(Path(f"{path}.tbi"))
    return audit


def _vcf_header_markers(path: Path) -> set[str]:
    """Read only header declarations needed to distinguish gVCF/TR/SV inputs."""
    markers = set()
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    try:
        with opener(path, "rt", encoding="utf-8", newline="") as handle:
            for raw in handle:
                if raw.startswith("#CHROM\t"):
                    break
                for marker in ("NON_REF", "SVTYPE", "TRID", "MOTIFS", "STRUC"):
                    if raw.startswith((f"##ALT=<ID={marker},", f"##ALT=<ID={marker}>", f"##INFO=<ID={marker},")):
                        markers.add(marker)
    except (OSError, EOFError, UnicodeError) as error:
        raise InputValidationError(f"Unable to inspect ready VCF declarations '{path}': {error}") from error
    return markers


def _audit_ready_gvcf(path: Path, context: AnalysisContext) -> dict[str, object]:
    """Check typed metadata and a bounded prefix for a real reference block.

    The preceding ``_ready_vcf`` guard still owns exact sample/build/contig
    header identity and index presence. DeepVariant uses the reserved ``<*>``
    allele without an ALT declaration; older gVCFs use declared ``<NON_REF>``.
    Neither an ALT header alone nor an ordinary variant is reference-block
    evidence. This read-only check is not full-record or REF-sequence validation.
    Lack of a witness within the prefix is explicitly *unverified*, not proof
    that the entire file is invalid. GLnexus still validates its complete input.
    """
    lengths = {contig.name: contig.length for contig in context.reference.contigs}
    definitions: dict[tuple[str, str], tuple[str | None, str | None]] = {}
    alt_ids: set[str] = set()
    required = {("INFO", "END"): ("1", "Integer"),
                ("FORMAT", "GT"): ("1", "String"),
                ("FORMAT", "PL"): ("G", "Integer")}
    examined = blocks = 0
    first_block = None
    exhausted = False
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    try:
        with opener(path, "rt", encoding="utf-8", newline="") as handle:
            while True:
                raw = handle.readline(_GVCF_LINE_CHARACTER_LIMIT + 1)
                if len(raw) > _GVCF_LINE_CHARACTER_LIMIT:
                    raise InputValidationError("GVCF_STRUCTURE_UNVERIFIED: header/record exceeds the bounded line limit.")
                if not raw:
                    raise InputValidationError("GVCF_STRUCTURE_UNVERIFIED: missing genotyped VCF header.")
                line = raw.rstrip("\r\n")
                if line.startswith(("##INFO=<", "##FORMAT=<", "##ALT=<")):
                    kind = line[2:line.index("=")]
                    if not line.endswith(">"):
                        raise InputValidationError("Malformed gVCF metadata declaration.")
                    lexer = shlex.shlex(line[line.index("<") + 1:-1], posix=True)
                    lexer.whitespace, lexer.whitespace_split, lexer.commenters = ",", True, ""
                    attributes = {}
                    for token in lexer:
                        key, separator, value = token.partition("=")
                        if not separator or not key or not value or key in attributes:
                            raise InputValidationError("Malformed or duplicate gVCF metadata attributes.")
                        attributes[key] = value
                    field_id = attributes.get("ID", "")
                    if kind == "ALT":
                        alt_ids.add(field_id)
                    elif (kind, field_id) in required:
                        key = kind, field_id
                        if key in definitions:
                            raise InputValidationError("Duplicate gVCF END/GT/PL declaration.")
                        definitions[key] = attributes.get("Number"), attributes.get("Type")
                if line.startswith("#CHROM\t"):
                    break
            missing = [f"{kind}/{field_id} Number={number},Type={value_type}"
                       for (kind, field_id), (number, value_type) in required.items()
                       if definitions.get((kind, field_id)) != (number, value_type)]
            if missing:
                raise InputValidationError(f"Cohort small_variants input lacks required declarations: {missing}.")
            while examined < _GVCF_RECORD_PREFIX_LIMIT:
                raw = handle.readline(_GVCF_LINE_CHARACTER_LIMIT + 1)
                if not raw:
                    exhausted = True
                    break
                if len(raw) > _GVCF_LINE_CHARACTER_LIMIT:
                    raise InputValidationError("GVCF_STRUCTURE_UNVERIFIED: record exceeds the bounded line limit.")
                if not raw.strip():
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                examined += 1
                if len(fields) != 10:
                    raise InputValidationError("Malformed single-sample gVCF record in inspected prefix.")
                if fields[0] not in lengths:
                    raise InputValidationError("REFERENCE_CONTIG_MISMATCH in inspected gVCF record prefix.")
                position = fields[1]
                if not position.isascii() or not position.isdecimal() or not 1 <= int(position) <= lengths[fields[0]]:
                    raise InputValidationError("Invalid gVCF record position in inspected prefix.")
                keys, values = fields[8].split(":"), fields[9].split(":")
                if len(keys) != len(set(keys)) or len(values) != len(keys) or not {"GT", "PL"}.issubset(keys):
                    raise InputValidationError("Inspected gVCF records require unambiguous GT and PL FORMAT values.")
                alt = fields[4]
                if alt not in {"<*>", "<NON_REF>"}:
                    continue
                if alt == "<NON_REF>" and "NON_REF" not in alt_ids:
                    raise InputValidationError("gVCF NON_REF reference block lacks its ALT declaration.")
                ends = [item.partition("=")[2] for item in fields[7].split(";") if item.startswith("END=")]
                if (len(fields[3]) != 1 or len(ends) != 1 or not ends[0].isascii()
                        or not ends[0].isdecimal() or not int(position) <= int(ends[0]) <= lengths[fields[0]]):
                    raise InputValidationError("Malformed gVCF reference-block END/span in inspected prefix.")
                genotype = re.split(r"[/|]", values[keys.index("GT")])
                if not genotype or any(allele not in {"0", "."} for allele in genotype):
                    raise InputValidationError("gVCF reference-block GT cannot call the symbolic non-reference placeholder.")
                likelihoods = values[keys.index("PL")].split(",")
                if (not likelihoods or any(not item.isascii() or not item.isdecimal() for item in likelihoods)):
                    if likelihoods != ["."] or "." not in genotype:
                        raise InputValidationError("Malformed gVCF reference-block PL values in inspected prefix.")
                blocks += 1
                if first_block is None:
                    first_block = {"prefix_record_number": examined, "contig": fields[0],
                                   "position": int(position), "end": int(ends[0]), "alt": alt}
    except (OSError, EOFError, UnicodeError, ValueError) as error:
        raise InputValidationError(f"Unable to inspect ready gVCF structure '{path}': {error}") from error
    if not blocks:
        raise InputValidationError(
            f"GVCF_STRUCTURE_UNVERIFIED: no NON_REF or <*> reference block with END/GT/PL "
            f"was observed in the inspected prefix ({examined}/{_GVCF_RECORD_PREFIX_LIMIT} records); "
            "ordinary variant-only VCFs are not accepted as ready gVCFs."
        )
    return {"scope": "TYPED_HEADER_AND_BOUNDED_RECORD_PREFIX", "records_examined": examined,
            "record_prefix_limit": _GVCF_RECORD_PREFIX_LIMIT,
            "line_character_limit": _GVCF_LINE_CHARACTER_LIMIT,
            "reached_eof": exhausted, "reference_blocks_in_prefix": blocks,
            "first_reference_block": first_block, "full_record_validation": "NOT_RUN",
            "reference_allele_sequence_check": "NOT_RUN", "input_modified": False}


def _manifest_text(fields: tuple[str, ...], rows: list[dict]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def _audit_ready_region_bed(path: Path, context: AnalysisContext) -> dict:
    """Check ready overlap BED bounds without assuming sorted input or rewriting it."""
    lengths = {contig.name: contig.length for contig in context.reference.contigs}
    count = 0
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    try:
        with opener(path, "rt", encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, 1):
                if not raw.strip() or raw.startswith(("#", "track", "browser")):
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) < 3 or fields[0] not in lengths:
                    raise InputValidationError(f"REFERENCE_CONTIG_MISMATCH or malformed annotation BED row: {path}:{line_number}.")
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError as error:
                    raise InputValidationError(f"Non-integer annotation BED coordinates: {path}:{line_number}.") from error
                if not 0 <= start < end <= lengths[fields[0]]:
                    raise InputValidationError(f"Annotation BED interval outside reference bounds: {path}:{line_number}.")
                count += 1
    except (OSError, EOFError, UnicodeError) as error:
        raise InputValidationError(f"Unable to inspect annotation BED '{path}': {error}") from error
    if not count:
        raise InputValidationError(f"Annotation BED has no intervals: {path}.")
    return {"path": str(path), "record_count": count, "coordinates": "ZERO_BASED_HALF_OPEN", "sorted_required": False}


def _prepare_ready_downstream(config: dict, context: AnalysisContext, base: Path,
                              root: Path) -> tuple[dict[str, str], list[Path], dict]:
    """Freeze explicitly selected ready inputs; return owned manifest contents.

    ``config`` is an owned effective copy and receives absolute deployment paths
    and owned manifest paths. No file is created here. There is deliberately no
    generated-source or caller auto-selection mode.
    """
    controls: dict[str, str] = {}
    inputs: list[Path] = []
    audit = {"scope": "READY_EXTERNAL_INPUTS_ONLY", "annotation": [], "cohort": [],
             "backend_failure_isolation": "NOT_IMPLEMENTED_BY_NATIVE_DAG", "functional_backend": "SEPARATE_EXPLICIT_API_CLI_NOT_THIS_DAG",
             "accuracy": "NOT_ASSESSED", "reference_sequence_check": "NOT_RUN_HEADER_COMPATIBILITY_ONLY"}
    annotation = config.get("annotation", {})
    if annotation.get("enabled"):
        from hifivar.annotation import read_annotation_inputs, read_annotation_input_rows
        from hifivar.annovar import AnnovarRequest
        from hifivar.vep import VepRequest
        source = Path(_absolute(annotation["input_manifest"], base))
        rows = read_annotation_input_rows(source)
        if not rows:
            raise InputValidationError("Annotation manifest must select at least one ready input.")
        if any(row["sample"].strip() not in context.sample_ids for row in rows):
            raise InputValidationError("Annotation manifest references an unknown sample.")
        selected = read_annotation_inputs(source, reference=context.reference)
        region_audits = []
        for key in ("annovar_database_root", "vep_cache_directory", "gene_bed", "exon_bed", "regulatory_bed", "repeat_bed", "segdup_bed"):
            if annotation.get(key):
                annotation[key] = _absolute(annotation[key], base)
                dependency = Path(annotation[key])
                if key.endswith("_bed"):
                    inputs.append(dependency)
                    if annotation.get("overlap_enabled"):
                        region_audits.append(_audit_ready_region_bed(dependency, context))
        audit["annotation_region_databases"] = region_audits
        normalized = []
        has_bnd_bundle_column = any('bnd_bundle_manifest' in row for row in rows)
        for item, row in zip(selected, rows, strict=True):
            if item.sample_id not in context.sample_ids:
                raise InputValidationError(f"Annotation manifest references unknown sample '{item.sample_id}'.")
            check = _ready_vcf(item.source_vcf.resolve(), item.sample_id, context)
            category = item.variant_category.value
            prefix = root / "results" / "annotation" / item.sample_id / category
            if annotation.get("annovar_enabled"):
                AnnovarRequest(item, Path(annotation["annovar_database_root"]), annotation["annovar_database_version"],
                               tuple(annotation["annovar_protocols"]), tuple(annotation["annovar_operations"]),
                               prefix / f"{item.sample_id}.{category}.annovar", annotation["annovar_version"])
            if annotation.get("vep_enabled"):
                VepRequest(item, Path(annotation["vep_cache_directory"]), annotation["vep_cache_version"],
                           annotation["vep_species"], annotation["vep_assembly"], prefix / f"{item.sample_id}.{category}.vep.tsv", annotation["threads"])
            tool_enabled = annotation.get("annovar_enabled") or annotation.get("vep_enabled")
            if annotation.get("overlap_enabled") and category in {"sv", "tr"}:
                if not item.source_variant_ids:
                    raise InputValidationError("Region-overlap launch requires explicit source_variant_ids for every SV/TR input.")
            elif not tool_enabled:
                raise InputValidationError("The selected annotation category has no enabled applicable backend.")
            normalized.append({"sample": item.sample_id, "variant_category": category,
                               "source_vcf": str(item.source_vcf.resolve()), "source_tool": item.source_tool,
                               "source_variant_ids": row["source_variant_ids"].strip()})
            if has_bnd_bundle_column:
                normalized[-1]['bnd_bundle_manifest'] = str(item.bnd_bundle_manifest.resolve()) if item.bnd_bundle_manifest else ''
            if item.bnd_bundle_manifest:
                from hifivar.junction_consumers import validate_bnd_consumer_bundle
                bundle = validate_bnd_consumer_bundle(item.bnd_bundle_manifest, item.source_vcf, context.reference, item.sample_id)
                inputs.append(Path(bundle['manifest']['path']))
                inputs.extend(Path(info['path']) for info in bundle['artifacts'].values())
            inputs.append(item.source_vcf.resolve())
            if str(item.source_vcf).lower().endswith(".vcf.gz"):
                inputs.append(Path(f"{item.source_vcf.resolve()}.tbi"))
            audit["annotation"].append({"sample_id": item.sample_id, "category": category, "source_tool": item.source_tool,
                                        "header_audit": check, "tool_annotation_scope": "WHOLE_SOURCE_VCF_NOT_FILTERED_BY_ID_LIST",
                                        "database_identity": "CONFIGURED_PATH_AND_VERSION_ONLY_NO_DIRECTORY_CHECKSUM"})
        fields = ("sample", "variant_category", "source_vcf", "source_tool", "source_variant_ids")
        if has_bnd_bundle_column:
            fields = (*fields, 'bnd_bundle_manifest')
        controls["annotation.inputs.tsv"] = _manifest_text(fields, normalized)
        inputs.append(source)
        annotation["input_manifest"] = str(root / "annotation.inputs.tsv")
    cohort_section = config.get("cohort", {})
    if cohort_section.get("enabled"):
        from hifivar.cohort import CohortDefinition, CohortSampleInput, CohortTrack, SampleCallState, validate_track_inputs
        cohort = CohortDefinition(cohort_section["cohort_id"], context.sample_ids, context.reference)
        source = Path(_absolute(cohort_section["input_manifest"], base))
        fields = ("sample", "track", "state", "source_path", "index_path", "source_tool", "source_version", "reference_build", "catalog_id")
        try:
            with source.open("r", encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if set(reader.fieldnames or ()) != set(fields) or len(reader.fieldnames or ()) != len(fields):
                    raise InputValidationError("Ready cohort manifest requires exactly the documented nine columns.")
                raw_rows = list(reader)
        except (OSError, UnicodeError, csv.Error) as error:
            raise InputValidationError(f"Unable to read ready cohort manifest '{source}': {error}") from error
        grouped: dict[str, dict[str, dict]] = {track.value: {} for track in CohortTrack}
        for row in raw_rows:
            if None in row or any(not isinstance(value, str) for value in row.values()):
                raise InputValidationError("Malformed cohort manifest row.")
            row = {key: value.strip() for key, value in row.items()}
            sample, track = row["sample"], row["track"]
            if track not in grouped or sample not in context.sample_ids:
                raise InputValidationError("Cohort manifest contains an unknown sample or track.")
            if sample in grouped[track]:
                raise InputValidationError("Cohort manifest contains duplicate sample/track rows.")
            grouped[track][sample] = row
        normalized = []
        for track in CohortTrack:
            if not cohort_section[track.value].get("enabled"):
                continue
            if set(grouped[track.value]) != set(context.sample_ids):
                raise InputValidationError(f"Cohort {track.value} requires exactly every declared sample; no samples are dropped.")
            items = []
            for sample in context.sample_ids:
                row = dict(grouped[track.value][sample])
                try:
                    state = SampleCallState(row["state"].upper())
                except ValueError as error:
                    raise InputValidationError("Unknown cohort sample state.") from error
                path = Path(_absolute(row["source_path"], source.parent)) if row["source_path"] else None
                index = Path(_absolute(row["index_path"], source.parent)) if row["index_path"] else None
                item = CohortSampleInput(sample, state, path, index, row["source_tool"] or None,
                                         row["source_version"] or None, row["reference_build"] or None, row["catalog_id"] or None)
                if track is CohortTrack.SMALL_VARIANTS and not item.callable:
                    raise InputValidationError("GLnexus requires a ready gVCF for every cohort sample; failed/missing samples are not dropped.")
                if item.callable:
                    if item.reference_build != context.reference.build or not item.source_tool or not item.source_version:
                        raise InputValidationError("Ready callable cohort inputs require explicit matching build, source tool and source version.")
                    assert path is not None
                    if not str(path).lower().endswith(".vcf.gz") or index != Path(f"{path}.tbi"):
                        raise InputValidationError("Ready cohort inputs require BGZF .vcf.gz and its conventional .tbi.")
                    check = _ready_vcf(path, sample, context)
                    structure = {}
                    if track is CohortTrack.SMALL_VARIANTS:
                        structure["gvcf_structure_audit"] = _audit_ready_gvcf(path, context)
                    else:
                        markers = _vcf_header_markers(path)
                        required = {CohortTrack.SV: {"SVTYPE"}, CohortTrack.TR: {"TRID", "MOTIFS", "STRUC"}}[track]
                        if not required.issubset(markers):
                            raise InputValidationError(f"Cohort {track.value} input lacks required declarations: {sorted(required - markers)}.")
                    inputs.extend((path, index))
                    audit["cohort"].append({"sample_id": sample, "track": track.value, "state": state.value,
                                            "header_audit": check, "source_version": item.source_version,
                                            "version_source": "USER_DECLARED_NOT_REPROBED", "state_source": "USER_DECLARED_NOT_RECORD_COUNT_INFERRED", **structure})
                elif path is not None or index is not None:
                    raise InputValidationError("Noncallable cohort samples must not attach unvalidated/future source files.")
                row.update(state=state.value, source_path=str(path) if path else "", index_path=str(index) if index else "")
                normalized.append(row)
                items.append(item)
            validate_track_inputs(cohort, tuple(items))
            if not any(item.callable for item in items):
                raise InputValidationError(f"Cohort {track.value} has no ready callable inputs.")
            if track is CohortTrack.TR:
                catalogs = {item.catalog_id for item in items if item.callable}
                if None in catalogs or len(catalogs) != 1:
                    raise InputValidationError("Ready TR cohort inputs require one identical explicit catalog_id.")
        controls["cohort.inputs.tsv"] = _manifest_text(fields, normalized)
        inputs.append(source)
        cohort_section["input_manifest"] = str(root / "cohort.inputs.tsv")
        audit["cohort_sv_cross_sample_clustering"] = False
    return controls, inputs, audit


def prepare_launch(config: HiFiVarConfig, run_dir: Path) -> tuple[dict, AnalysisContext, str, list[Path]]:
    """Resolve user-relative inputs and validate scope without creating any outputs."""
    validate_config(config, require_complete=True)
    effective = config.to_dict()
    _enabled_tracks(effective)
    raw_enabled = effective["alignment"].get("enabled", False)
    if raw_enabled:
        from hifivar.launch_alignment import validate_alignment_budget
        validate_alignment_budget(effective)
    context = AnalysisContext.from_config(config)
    from hifivar.calling_region import validate_calling_region
    validate_calling_region(effective["small"].get("region"), context.reference)
    if raw_enabled and effective["alignment"]["tool"] == "minimap2" and any(c.length >= 2**29 for c in context.reference.contigs):
        raise ConfigurationError("minimap2 launcher requires a BAI-compatible reference.")
    source = config.sources.get("user")
    base = Path(source).resolve().parent if source else Path.cwd()
    assembly_context = build_assembly_context(config, context) if effective["assembly"].get("enabled") else None
    for section, key, destination in (
        ("paths", "workdir", "work"), ("paths", "outdir", "results"),
        ("runtime", "tmpdir", "tmp"), ("logging", "file", "logs/hifivar.log"),
    ):
        value = effective[section].get(key)
        expected = run_dir / destination
        if value is not None and Path(_absolute(value, run_dir)) != expected:
            raise ConfigurationError(f"{section}.{key} must be null or {destination!r} for the owned run directory.")
        effective[section][key] = str(expected)
    effective["reference"]["fasta"] = str(context.reference.fasta.resolve())
    effective["samples"]["sheet"] = str(run_dir / "samples.tsv")
    inputs = [context.reference.fasta.resolve(), context.reference.fai.resolve(), context.source_sample_sheet.resolve()]
    if effective["assembly_sv"].get("input_manifest"):
        from hifivar.ready_assemblies import load_ready_assembly_manifest
        ready_path = Path(_absolute(effective["assembly_sv"]["input_manifest"], base))
        ready = load_ready_assembly_manifest(ready_path, expected_samples=context.sample_ids)
        effective["assembly_sv"]["input_manifest"] = str(ready.path)
        inputs.extend((ready.path, *(item.path for item in ready.assemblies)))
    if assembly_context is not None:
        if assembly_context.source_sample_sheet is None:
            raise InputValidationError("Assembly requires a traceable FASTQ sample sheet.")
        assembly_sheet = assembly_context.source_sample_sheet.resolve()
        effective["assembly"]["sample_sheet"] = str(assembly_sheet)
        inputs.append(assembly_sheet)
        inputs.extend(path.resolve() for record in assembly_context.samples for path in record.sample.input.files)
    if source:
        inputs.append(Path(source).resolve())
    for name in ("workflow", "alignment", "small", "sv", "tr", "phasing", "assembly", "assembly_sv", "annotation", "cohort", "review", "benchmark"):
        _normalize_executables(effective[name], base)
    if effective['sv']['harmonization'].get('enabled'):
        from hifivar.harmonization_inputs import harmonization_sources
        from hifivar.attempt_receipts import artifact_identity
        external = effective['sv']['harmonization'].get('input_vcfs', {})
        if set(external) - set(context.sample_ids):
            raise ConfigurationError('External harmonization contains undeclared sample IDs.')
        for sample, callers in external.items():
            if not isinstance(callers, dict):
                raise ConfigurationError('External harmonization sample must map callers to VCF paths.')
            for caller, value in callers.items():
                if not isinstance(value, str) or not value.strip():
                    raise ConfigurationError('External harmonization VCF path must be nonempty text.')
                raw_path = Path(value).expanduser()
                absolute = raw_path if raw_path.is_absolute() else base / raw_path
                # Check aliases before resolve() can hide a redirected input.
                artifact_identity(absolute.absolute())
                callers[caller] = str(absolute.resolve(strict=True))
        for sample in context.sample_ids:
            selected_sources = harmonization_sources(effective, sample, run_dir)
            if effective['sv']['harmonization'].get('minimum_caller_support', 1) > len(selected_sources):
                raise ConfigurationError(f"minimum_caller_support exceeds the actual input callers for sample '{sample}'.")
            for selected in selected_sources:
                if selected['origin'] == 'EXTERNAL':
                    _ready_vcf(selected['vcf'], sample, context)
                    artifact_identity(selected['index'])
                    inputs.extend((selected['vcf'], selected['index']))
    if effective["workflow"].get("downstream_source_plan"):
        from hifivar.downstream_handoff import prepare_downstream_intent
        _, downstream_inputs, _ = prepare_downstream_intent(effective, context, base, run_dir)
    else:
        _, downstream_inputs, _ = _prepare_ready_downstream(effective, context, base, run_dir)
    inputs.extend(downstream_inputs)
    if effective.get("benchmark", {}).get("enabled"):
        from hifivar.benchmark_handoff import prepare_intent
        _, benchmark_inputs, _ = prepare_intent(effective, context, base, run_dir)
        inputs.extend(benchmark_inputs)
    if effective.get("review", {}).get("enabled"):
        from hifivar.review_handoff import prepare_review_intent
        _, review_inputs, _ = prepare_review_intent(effective, context, base, run_dir)
        inputs.extend(review_inputs)
    if effective["assembly_sv"].get("enabled") and effective["assembly_sv"]["pav"].get("enabled"):
        pav = effective["assembly_sv"]["pav"]
        pav["snakefile"] = _absolute(pav["snakefile"], base)
        inputs.append(Path(pav["snakefile"]))
    if effective["tr"]["enabled"]:
        effective["tr"]["catalog"] = _absolute(effective["tr"]["catalog"], base)
        inputs.append(Path(effective["tr"]["catalog"]))
    if effective["small"]["enabled"] and effective["small"]["execution_mode"] == "apptainer":
        effective["small"]["deepvariant_image"] = _absolute(effective["small"]["deepvariant_image"], base)
        inputs.append(Path(effective["small"]["deepvariant_image"]))
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(["sample_id", "input", "input_type", "sex", "father", "mother", "phenotype", "group"])
    seen = []
    for record in context.samples:
        sample = record.sample
        if _calling_enabled(effective) and not raw_enabled and (sample.input.input_type is not InputType.BAM or len(sample.input.files) != 1):
            raise InputValidationError("hifivar run currently accepts one indexed aligned BAM per sample; FASTQ/uBAM/CRAM routing is deferred.")
        if sample.input.input_type is InputType.CRAM:
            raise InputValidationError("CRAM routing is deferred; supply FASTQ, uBAM or an indexed aligned BAM.")
        paths = [p.resolve() for p in sample.input.files]
        for path in paths:
            if any(path.samefile(previous) for previous in seen):
                raise InputValidationError("The same BAM/FASTQ is assigned more than once, including a symlink/hardlink alias.")
            seen.append(path)
            inputs.append(path)
        if sample.input.input_type is InputType.BAM and (_calling_enabled(effective) or raw_enabled):
            bam = paths[0]
            indexes = {Path(str(bam) + s) for s in (".bai", ".csi")} | {bam.with_suffix(s) for s in (".bai", ".csi")}
            if not raw_enabled or any(p.exists() for p in indexes):
                inputs.append(unique_bam_index(bam))
        writer.writerow([sample.sample_id, ";".join(map(str, paths)), sample.input.input_type.value, record.sex.value if record.sex else "",
                         record.father or "", record.mother or "", record.phenotype or "", record.group or ""])
    for path in inputs:
        if not path.is_file() or path.stat().st_size == 0:
            raise InputValidationError(f"Required input missing or empty: {path}")
        if path.resolve().is_relative_to(run_dir):
            raise InputValidationError("Run directory must not contain source inputs.")
    validate_config(effective, require_complete=True)
    return effective, context, output.getvalue(), list(dict.fromkeys(inputs))


def build_workflow_command(workflow_root: Path, run_dir: Path, config: dict,
                           *, dry_run: bool, rerun_incomplete: bool,
                           config_path: Path | None = None) -> list[str]:
    """Use the current interpreter's installed workflow, with an explicit all target."""
    command = [sys.executable, "-I", "-B", "-m", "snakemake", "--snakefile", str(workflow_root / "Snakefile"),
               "--configfile", str(config_path or run_dir / "effective.yaml"), "--cores", str(config["runtime"]["threads"]),
               "--resources", f"mem_mb={config['workflow']['memory_mb']}",
               f"deepvariant_slots={config['small']['max_concurrent_samples']}", "--scheduler", "greedy",
               "--printshellcmds", "--shared-fs-usage", "input-output", "persistence", "software-deployment",
               "software-deployment-cache", "sources", "storage-local-copies"]
    if dry_run:
        command.append("--dry-run")
    if rerun_incomplete:
        command.append("--rerun-incomplete")
    return [*command, "--", "all"]


def _reject_links(root: Path) -> None:
    """Refuse redirected output trees on resume (never follow symlinks/junctions)."""
    for current, dirs, files in os.walk(root, followlinks=False):
        for name in (*dirs, *files):
            path = Path(current) / name
            if path.is_symlink() or (hasattr(path, "is_junction") and path.is_junction()):
                raise OutputValidationError(f"Run tree contains a redirected path: {path}")


class WorkflowLauncher:
    """Launch once or explicitly resume; preserve every attempt and fail closed."""

    def __init__(self, runner: CommandRunner | None = None) -> None:
        self.runner = runner or CommandRunner()

    def run(self, config: HiFiVarConfig, *, run_dir: Path, resume: bool = False,
            rerun_incomplete: bool = False, dry_run: bool = False) -> str:
        """Run a Linux DAG under an exclusive directory lock and immutable identity."""
        if platform.system() != "Linux":
            raise WorkflowError("hifivar run requires Linux; Windows supports config/unit checks only.")
        if rerun_incomplete and not resume:
            raise ConfigurationError("--rerun-incomplete requires --resume.")
        supplied = run_dir.expanduser().absolute()
        root = supplied.resolve()
        if supplied != root or root == Path(root.anchor) or root == Path.home().resolve():
            raise OutputValidationError("Use a real, non-root, non-symlink run directory.")
        if root.exists() != resume:
            raise OutputValidationError("A fresh run needs a new directory; --resume needs an existing launcher directory.")
        effective, context, sheet, inputs = prepare_launch(config, root)
        source = config.sources.get("user")
        base = Path(source).resolve().parent if source else Path.cwd()
        downstream_config = config.to_dict()
        if downstream_config["workflow"].get("downstream_source_plan"):
            from hifivar.downstream_handoff import prepare_downstream_intent
            downstream_controls, _, downstream_audit = prepare_downstream_intent(downstream_config, context, base, root)
        else:
            downstream_controls, _, downstream_audit = _prepare_ready_downstream(downstream_config, context, base, root)
        if downstream_config.get("review", {}).get("enabled"):
            from hifivar.review_handoff import prepare_review_intent
            review_controls, _, review_audit = prepare_review_intent(downstream_config, context, base, root)
            downstream_controls.update(review_controls)
            downstream_audit["review"] = review_audit
        if downstream_config.get("benchmark", {}).get("enabled"):
            from hifivar.benchmark_handoff import prepare_intent
            benchmark_controls, _, benchmark_audit = prepare_intent(downstream_config, context, base, root)
            downstream_controls.update(benchmark_controls)
            downstream_audit["benchmark"] = benchmark_audit
        wfroot = installed_workflow_root().resolve()
        if not (wfroot / "Snakefile").is_file():
            raise WorkflowError("Installed Snakefile is missing.")
        if resume:
            _reject_links(root)
            if not (root / "launch_manifest.json").is_file():
                raise OutputValidationError("Not an owned launcher directory (missing launch_manifest.json).")
        else:
            root.mkdir(parents=True, exist_ok=False)
        lock = root / ".hifivar-launch.lock"
        try:
            with lock.open("x", encoding="utf-8") as stream:
                json.dump({"pid": os.getpid(), "host": platform.node(), "created": utc_now_iso8601()}, stream)
        except FileExistsError as error:
            raise OutputValidationError("Run is locked; do not clear a stale lock until scheduler/children are independently confirmed stopped.") from error
        interrupted = False
        try:
            return self._execute(root, config=effective, context=context, sheet=sheet, inputs=inputs,
                                 workflow_root=wfroot, resume=resume, dry_run=dry_run, rerun_incomplete=rerun_incomplete,
                                 downstream_controls=downstream_controls, downstream_audit=downstream_audit)
        except (OSError, ValueError) as error:
            raise WorkflowError(f"Launcher filesystem/metadata error: {error}") from error
        except BaseException as error:
            # Preserve the lock on Ctrl-C/SystemExit: descendants may still exist.
            interrupted = not isinstance(error, Exception)
            raise
        finally:
            if not interrupted:
                lock.unlink(missing_ok=True)

    def _execute(self, root: Path, *, config: dict, context: AnalysisContext, sheet: str,
                 inputs: list[Path], workflow_root: Path, resume: bool,
                 dry_run: bool, rerun_incomplete: bool,
                 downstream_controls: dict[str, str] | None = None, downstream_audit: dict | None = None) -> str:
        identity = {"schema": 2, "kind": "hifivar_two_stage_launcher", "run_dir": str(root),
                    "config": config, "sample_sheet": sheet, "runtime": runtime_identity(workflow_root),
                    "inputs": [file_identity(path) for path in inputs]}
        manifest = root / "launch_manifest.json"
        if resume:
            previous = json.loads(manifest.read_text(encoding="utf-8"))
            if previous["identity"] != identity:
                raise InputValidationError("RESUME_IDENTITY_MISMATCH: config, input bytes, interpreter or runtime source changed; use a new run directory.")
            for path, checksum in previous["control_sha256"].items():
                if file_identity(root / path)["sha256"] != checksum:
                    raise InputValidationError("RESUME_CONTROL_FILE_CHANGED: effective.yaml/samples.tsv changed.")
            for name, text in (downstream_controls or {}).items():
                if (root / name).read_text(encoding="utf-8") != text:
                    raise InputValidationError("RESUME_DOWNSTREAM_MANIFEST_CHANGED: normalized source selections changed.")
        else:
            write_yaml_atomic(config, root / "effective.yaml")
            with (root / "samples.tsv").open("x", encoding="utf-8", newline="") as stream:
                stream.write(sheet)
            for name, text in (downstream_controls or {}).items():
                with (root / name).open("x", encoding="utf-8", newline="") as stream:
                    stream.write(text)
            write_json_atomic({"identity": identity, "created": utc_now_iso8601(),
                               "control_sha256": {name: file_identity(root / name)["sha256"] for name in ("effective.yaml", "samples.tsv", *(downstream_controls or {}))}}, manifest)
        for name in ("logs", "tmp", "work", "results", "cache", "attempts", "config"):
            (root / name).mkdir(exist_ok=True)
        attempt = root / "attempts" / ("attempt_" + uuid4().hex)
        attempt.mkdir()
        status = {"started": utc_now_iso8601(), "dry_run": dry_run, "resume": resume,
                  "rerun_incomplete": rerun_incomplete, "status": "PREFLIGHT_RUNNING"}
        write_json_atomic(status, attempt / "status.json")
        try:
            if downstream_controls:
                write_json_atomic(downstream_audit or {}, attempt / "downstream_input_audit.json")
            calling_config = None
            alignment_inputs = []
            if config["assembly"].get("enabled"):
                assembly_context = build_assembly_context(config, context)
                write_json_atomic({
                    "kind": "explicit_assembly_fastq_context", "implicit_bam_extraction": False,
                    "primary_sample_order": list(context.sample_ids),
                    "assembly_sample_order": list(assembly_context.sample_ids),
                    "sample_sheet": str(assembly_context.source_sample_sheet),
                    "samples": [record.to_dict() for record in assembly_context.samples],
                    "biological_identity": "USER_DECLARED_METADATA_MATCH_NOT_INDEPENDENTLY_VERIFIED",
                    "native_pav_container_identity": "NOT_ROUTED_BY_NATIVE_ADAPTER",
                }, attempt / "assembly_input_context.json")
            if config["alignment"].get("enabled", False):
                from hifivar.launch_alignment import prepare_aligned_samples
                status["status"] = "ALIGNMENT_PREPARATION"
                write_json_atomic(status, attempt / "status.json", overwrite=True)
                calling_config = prepare_aligned_samples(context, config, root, attempt, self.runner,
                                                         dry_run=dry_run, rerun_incomplete=rerun_incomplete)
                if dry_run:
                    if identity["inputs"] != [file_identity(path) for path in inputs]:
                        raise InputValidationError("Input bytes changed during alignment planning.")
                    status.update(status="RAW_INPUT_PLAN_ONLY", input_immutability="PASS", calling_dag="NOT_RUN")
                    return str(status["status"])
                from hifivar.alignment_preparation import verify_handoff
                alignment_binding = file_identity(attempt / "alignment_handoff.json")
                alignment_inputs = verify_handoff(root, attempt, alignment_binding, verify_bams=True)["selected_artifacts"]
                if identity["inputs"] != [file_identity(path) for path in inputs]:
                    raise InputValidationError("Input bytes changed during alignment preparation.")
            elif _calling_enabled(config):
                preflight = IndexedBamPreflight(self.runner, config["workflow"]["samtools_executable"])
                checks = []
                for record in context.samples:
                    bam = record.sample.input.files[0].resolve()
                    checks.append(preflight.run(bam, unique_bam_index(bam), record.sample.sample_id, context.reference, attempt))
                write_json_atomic(checks, attempt / "input_preflight.json")
            environment = {key: value for key, value in os.environ.items()
                           if key in {"PATH", "HOME", "USER", "LOGNAME", "LANG", "LC_ALL", "TZ", "CONDA_PREFIX", "CONDA_EXE"}}
            environment.update(TMPDIR=str(root / "tmp"), TMP=str(root / "tmp"), TEMP=str(root / "tmp"),
                               XDG_CACHE_HOME=str(root / "cache"), XDG_CONFIG_HOME=str(root / "config"),
                               PYTHONDONTWRITEBYTECODE="1")
            version = self.runner.run([sys.executable, "-I", "-B", "-m", "snakemake", "--version"], timeout=60,
                                      env=environment, inherit_env=False)
            status["snakemake_version"] = (version.stdout or "").strip()
            if not status["snakemake_version"]:
                raise WorkflowError("Unable to resolve Snakemake version in the active interpreter.")
            command = build_workflow_command(workflow_root, root, config, dry_run=dry_run, rerun_incomplete=rerun_incomplete,
                                             config_path=calling_config)
            generated_plan = bool(config["workflow"].get("downstream_source_plan"))
            calling_base = config
            if calling_config is not None:
                import yaml
                calling_base = yaml.safe_load(calling_config.read_text(encoding='utf-8'))
            if config.get('review',{}).get('enabled') or config.get('benchmark',{}).get('enabled') or generated_plan:
                from hifivar.review_handoff import without_review
                from hifivar.downstream_handoff import producer_config
                producer_path = attempt / "producer.effective.yaml"
                base_config = calling_base
                stage_config = producer_config(base_config) if generated_plan else copy.deepcopy(base_config)
                stage_config = without_review(stage_config)
                from hifivar.benchmark_handoff import without_benchmark
                stage_config = without_benchmark(stage_config)
                write_yaml_atomic(stage_config, producer_path)
                producer_config_identity = file_identity(producer_path)
                command = build_workflow_command(workflow_root, root, config, dry_run=dry_run,
                                                 rerun_incomplete=rerun_incomplete, config_path=producer_path)
            status.update(status="PRODUCER_DAG_RUNNING" if generated_plan else "DAG_RUNNING", command=command)
            write_json_atomic(status, attempt / "status.json", overwrite=True)
            if config["assembly_sv"].get("input_manifest"):
                from hifivar.assembly_reuse import prepare_ready_handoff
                prepare_ready_handoff(root, attempt)
            execution_context = create_attempt_context(root, attempt, dry_run=dry_run)
            environment[ATTEMPT_CONTEXT_ENV] = str(execution_context)
            status.update(receipt_schema_version=1,
                          execution_context_sha256=file_identity(execution_context)["sha256"])
            write_json_atomic(status, attempt / "status.json", overwrite=True)
            _LOGGER.info("Launching %d samples; progress log: %s", context.n_samples, attempt / "snakemake.stderr.log")
            result = self.runner.run(command, cwd=root, check=False,
                                     env=environment, inherit_env=False,
                                     stdout_path=attempt / "snakemake.stdout.log", stderr_path=attempt / "snakemake.stderr.log")
            status.update(returncode=result.returncode, duration_seconds=result.duration_seconds)
            if result.returncode != 0:
                raise WorkflowError(f"Snakemake failed (exit {result.returncode}); inspect {attempt}. No automatic retry performed.")
            if (config.get("review", {}).get("enabled") or config.get("benchmark", {}).get("enabled")) and file_identity(producer_path) != producer_config_identity:
                raise InputValidationError("Review-disabled producer configuration changed during execution.")
            if generated_plan:
                if file_identity(producer_path) != producer_config_identity:
                    raise InputValidationError("Producer-stage config changed during execution.")
                if dry_run:
                    status.update(status="GENERATED_SOURCE_PLAN_ONLY", input_immutability="PASS",
                                  downstream_dag="NOT_RUN_WAITING_FOR_GENERATED_RECEIPTS")
                    if identity["inputs"] != [file_identity(path) for path in inputs]:
                        raise InputValidationError("Input bytes changed during generated-source planning.")
                    return str(status["status"])
                from hifivar.downstream_handoff import materialize_downstream
                # Do not begin a consumer stage with changed upstream inputs,
                # even when a filename/header still looks compatible.
                if identity["inputs"] != [file_identity(path) for path in inputs]:
                    raise InputValidationError("Source inputs changed before downstream execution.")
                if any(file_identity(Path(item["path"])) != item for item in alignment_inputs):
                    raise InputValidationError("Derived alignment changed before downstream execution.")
                if calling_config is not None:
                    verify_handoff(root, attempt, alignment_binding)
                expected_manifest_sha = json.loads(execution_context.read_text(encoding="utf-8"))["manifest_sha256"]
                if file_identity(manifest)["sha256"] != expected_manifest_sha:
                    raise InputValidationError("Launcher manifest changed before downstream execution.")
                for name, digest in json.loads(manifest.read_text(encoding="utf-8"))["control_sha256"].items():
                    if file_identity(root / name)["sha256"] != digest:
                        raise InputValidationError("Frozen launcher controls changed before downstream execution.")
                status["status"] = "DOWNSTREAM_PREPARATION"
                write_json_atomic(status, attempt / "status.json", overwrite=True)
                downstream_path, handoff = materialize_downstream(config, context, root, attempt, calling_config=base_config)
                downstream_command = build_workflow_command(workflow_root, root, config, dry_run=False,
                                                            rerun_incomplete=False, config_path=downstream_path)
                status.update(status="DOWNSTREAM_DAG_RUNNING", producer_returncode=result.returncode,
                              downstream_command=downstream_command,
                              downstream_materialization_sha256=file_identity(attempt / "downstream" / "materialization.json")["sha256"])
                write_json_atomic(status, attempt / "status.json", overwrite=True)
                downstream_result = self.runner.run(
                    downstream_command, cwd=root, check=False, env=environment, inherit_env=False,
                    stdout_path=attempt / "downstream.snakemake.stdout.log",
                    stderr_path=attempt / "downstream.snakemake.stderr.log",
                )
                status.update(returncode=downstream_result.returncode,
                              downstream_duration_seconds=downstream_result.duration_seconds)
                if downstream_result.returncode != 0:
                    raise WorkflowError(f"Downstream Snakemake failed (exit {downstream_result.returncode}); no automatic retry performed.")
                if file_identity(attempt / "downstream" / "materialization.json")["sha256"] != status["downstream_materialization_sha256"]:
                    raise InputValidationError("Downstream materialization record changed during execution.")
                for item in handoff["controls"].values():
                    if file_identity(Path(item["path"])) != item:
                        raise InputValidationError("Materialized downstream controls changed during execution.")
                for link in handoff["source_receipts"]:
                    recorded = {key: link[key] for key in ("path", "bytes", "sha256")}
                    if file_identity(Path(link["path"])) != recorded:
                        raise InputValidationError("A producer receipt changed during downstream execution.")
                    from hifivar.attempt_receipts import read_step_receipt
                    read_step_receipt(root, attempt.name, link["track"], tuple(link["samples"]))
            if config.get('benchmark',{}).get('enabled'):
                from hifivar.attempt_receipts import _context as verify_benchmark_context, read_step_receipt
                verify_benchmark_context(execution_context,live=False)
                if identity['inputs'] != [file_identity(path) for path in inputs]:
                    raise InputValidationError('Benchmark inputs changed before consumer preparation.')
                if dry_run:
                    status.update(status='BENCHMARK_SOURCE_PLAN_ONLY',benchmark_dag='NOT_RUN_WAITING_FOR_QUERY_RECEIPTS',input_immutability='PASS')
                    return str(status['status'])
                from hifivar.benchmark_handoff import materialize,validate_materialization
                status['status']='BENCHMARK_PREPARATION'
                write_json_atomic(status,attempt/'status.json',overwrite=True)
                benchmark_path,benchmark_handoff=materialize(config,context,root,attempt,calling_config=calling_base)
                benchmark_command=build_workflow_command(workflow_root,root,config,dry_run=False,rerun_incomplete=False,config_path=benchmark_path)
                status.update(status='BENCHMARK_DAG_RUNNING',benchmark_command=benchmark_command,
                    benchmark_materialization_sha256=file_identity(attempt/'benchmark/materialization.json')['sha256'])
                write_json_atomic(status,attempt/'status.json',overwrite=True)
                benchmark_result=self.runner.run(benchmark_command,cwd=root,check=False,env=environment,inherit_env=False,
                    stdout_path=attempt/'benchmark.snakemake.stdout.log',stderr_path=attempt/'benchmark.snakemake.stderr.log')
                status.update(returncode=benchmark_result.returncode,benchmark_returncode=benchmark_result.returncode,benchmark_duration_seconds=benchmark_result.duration_seconds)
                if benchmark_result.returncode!=0:
                    raise WorkflowError(f'Benchmark Snakemake failed (exit {benchmark_result.returncode}); no automatic retry performed.')
                validate_materialization(root,attempt)
                if read_step_receipt(root,attempt.name,'benchmark',(config['benchmark']['sample_id'],)) is None:
                    raise OutputValidationError('Benchmark scheduler success has no selected-attempt native benchmark receipt.')
            if config.get('review',{}).get('enabled'):
                from hifivar.attempt_receipts import _context as verify_review_context
                verify_review_context(execution_context, live=False)
                if dry_run:
                    if identity['inputs'] != [file_identity(path) for path in inputs]:
                        raise InputValidationError('Input bytes changed during review planning.')
                    status.update(status='REVIEW_SOURCE_PLAN_ONLY',review_dag='NOT_RUN_WAITING_FOR_SOURCE_RECEIPTS',input_immutability='PASS')
                    return str(status['status'])
                from hifivar.review_handoff import materialize_review, validate_review_materialization
                if identity['inputs'] != [file_identity(path) for path in inputs]:
                    raise InputValidationError('Input bytes changed before review preparation.')
                if calling_config is not None:
                    verify_handoff(root,attempt,alignment_binding)
                status['status'] = 'REVIEW_PREPARATION'
                write_json_atomic(status,attempt/'status.json',overwrite=True)
                review_path, review_handoff = materialize_review(config,context,root,attempt,calling_config=calling_base)
                review_command = build_workflow_command(workflow_root,root,config,dry_run=False,
                                                       rerun_incomplete=False,config_path=review_path)
                status.update(status='REVIEW_DAG_RUNNING',review_command=review_command,
                              review_materialization_sha256=file_identity(attempt/'review/materialization.json')['sha256'])
                write_json_atomic(status,attempt/'status.json',overwrite=True)
                review_result = self.runner.run(review_command,cwd=root,check=False,env=environment,inherit_env=False,
                    stdout_path=attempt/'review.snakemake.stdout.log',stderr_path=attempt/'review.snakemake.stderr.log')
                status.update(returncode=review_result.returncode,review_returncode=review_result.returncode,
                              review_duration_seconds=review_result.duration_seconds)
                if review_result.returncode != 0:
                    raise WorkflowError(f'Review Snakemake failed (exit {review_result.returncode}); no automatic retry performed.')
                _,_,_,review_samples = validate_review_materialization(root,attempt)
                from hifivar.attempt_receipts import read_step_receipt
                if read_step_receipt(root,attempt.name,'review',review_samples) is None:
                    raise OutputValidationError('Review scheduler succeeded without a selected-attempt native review receipt.')
            if identity["inputs"] != [file_identity(path) for path in inputs]:
                raise InputValidationError("Input bytes changed during workflow execution.")
            if any(file_identity(Path(item["path"])) != item for item in alignment_inputs):
                raise InputValidationError("Derived alignment bytes changed during calling.")
            if calling_config is not None:
                verify_handoff(root, attempt, alignment_binding)
            frozen_controls = json.loads(manifest.read_text(encoding="utf-8"))["control_sha256"]
            for name, digest in frozen_controls.items():
                if file_identity(root / name)["sha256"] != digest:
                    raise InputValidationError("Frozen launcher controls changed during workflow execution.")
            status.update(status="DRY_RUN_PASS" if dry_run else "WORKFLOW_SUCCEEDED", input_immutability="PASS")
        except BaseException as error:
            status.update(status="FAILED" if isinstance(error, Exception) else "INTERRUPTED_LOCK_RETAINED",
                          error_type=type(error).__name__, error=str(error))
            raise
        finally:
            status["finished"] = utc_now_iso8601()
            write_json_atomic(status, attempt / "status.json", overwrite=True)
        return str(status["status"])


__all__ = ["WorkflowLauncher", "build_assembly_context", "build_workflow_command", "prepare_launch"]
