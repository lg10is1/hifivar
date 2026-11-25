"""Read explicit downstream source intent without running or certifying callers.

Generated rows select a producer and artifact role, not a future file that is
already called. External rows establish only existing file paths and declared
provenance; launcher/backend validation still has to inspect their contents.
"""

from __future__ import annotations

import hashlib
import os
import re
import unicodedata
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml
from yaml.events import AliasEvent

from hifivar.cohort import SampleCallState
from hifivar.exceptions import InputValidationError, ReferenceError
from hifivar.reference import normalize_reference_build
from hifivar.sample import validate_sample_id


_MAX_BYTES = 1024 * 1024
_MAX_LINES = 4096
_MAX_ROWS = 4096
_TOP_KEYS = frozenset({"schema_version", "annotation", "cohort"})
_ROW_KEYS = frozenset({
    "sample_id", "track", "origin", "producer", "output_role", "path",
    "index_path", "source_tool", "source_version", "reference_build",
    "catalog_id", "state", "variant_ids", "bnd_bundle_manifest",
})
_CALLABLE_STATES = frozenset({SampleCallState.CALLED.value, SampleCallState.NO_CALLS.value})
_CATALOG_SHA = re.compile(r"(?:sha256:)?([0-9a-fA-F]{64})\Z")
_PRODUCERS = {
    "small.deepvariant": ("small", frozenset({"vcf", "gvcf"})),
    "sv.sawfish": ("sv", frozenset({"vcf"})),
    "sv.sniffles2": ("sv", frozenset({"vcf"})),
    "sv.pbsv": ("sv", frozenset({"vcf"})),
    "sv.cutesv": ("sv", frozenset({"vcf"})),
    "sv.harmonization": ("sv", frozenset({"vcf", "bnd_vcf"})),
    "tr.trgt": ("tr", frozenset({"vcf"})),
    "phasing.hiphase": ("small", frozenset({"vcf"})),
    "assembly_sv.pav": ("sv", frozenset({"vcf"})),
    "assembly_sv.svim_asm": ("sv", frozenset({"vcf"})),
}


@dataclass(frozen=True, slots=True)
class DownstreamPlanRow:
    """One immutable source request; not an execution or readiness receipt."""

    sample_id: str
    track: str
    origin: str
    producer: str | None
    output_role: str | None
    path: Path | None
    index_path: Path | None
    source_tool: str | None
    source_version: str | None
    reference_build: str
    catalog_id: str | None
    state: str = "CALLED"
    variant_ids: tuple[str, ...] = ()
    bnd_bundle_manifest: Path | None = None


@dataclass(frozen=True, slots=True)
class DownstreamSourcePlan:
    """Frozen bounded plan and SHA256 of its exact original UTF-8 bytes."""

    schema_version: int
    path: Path
    sha256: str
    annotation: tuple[DownstreamPlanRow, ...]
    cohort: tuple[DownstreamPlanRow, ...]


class _PlanLoader(yaml.SafeLoader):
    """Reject duplicate keys, merge keys and aliases instead of guessing intent."""

    def compose_node(self, parent: Any, index: Any) -> Any:
        if self.check_event(AliasEvent):
            raise InputValidationError("Downstream plan YAML aliases are not supported.")
        return super().compose_node(parent, index)

    def construct_mapping(self, node: Any, deep: bool = False) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key_node, value_node in node.value:
            if key_node.tag == "tag:yaml.org,2002:merge":
                raise InputValidationError("Downstream plan YAML merge keys are not supported.")
            key = self.construct_object(key_node, deep=deep)
            if not isinstance(key, str):
                raise InputValidationError("Downstream plan mapping keys must be strings.")
            _text(key, "mapping key")
            if key in result:
                raise InputValidationError("Downstream plan contains a duplicate mapping key.")
            result[key] = self.construct_object(value_node, deep=deep)
        return result


def read_downstream_plan(
    path: Path,
    *,
    sample_ids: Sequence[str],
    reference_build: str,
    config: Mapping[str, Any],
) -> DownstreamSourcePlan:
    """Read strict schema 1 source selectors without creating any files.

    External paths are resolved beside the plan and must already be nonempty
    regular files. Their VCF/index contents are deliberately not certified here.
    A generated row retains only an explicitly enabled producer and output role;
    its final source, version and CALLED/NO_CALLS outcome require a later receipt.

    Configured relative TR catalog paths follow a retained user-config source
    when available, otherwise the plan directory. Catalog IDs are content SHA256
    values, never filenames or arbitrary equivalent names.
    """
    samples = _sample_ids(sample_ids)
    build = _build(reference_build, "reference_build")
    if not isinstance(config, Mapping):
        raise InputValidationError("Downstream plan config must be a mapping.")
    source, raw = _plan_bytes(path)
    try:
        text = raw.decode("utf-8-sig")
        if len(text.splitlines()) > _MAX_LINES:
            raise InputValidationError("Downstream plan exceeds the 4096-line limit.")
        document = yaml.load(text, Loader=_PlanLoader)
    except (UnicodeError, yaml.YAMLError, RecursionError) as error:
        raise InputValidationError("Downstream plan must be bounded valid UTF-8 YAML.") from error
    document = _mapping(document, "Downstream plan")
    _unknown_keys(document, _TOP_KEYS, "Downstream plan")
    if type(document.get("schema_version")) is not int or document["schema_version"] != 1:
        raise InputValidationError("Downstream plan schema_version must be integer 1.")
    groups: dict[str, list[Any]] = {}
    for name in ("annotation", "cohort"):
        value = document.get(name, [])
        if not isinstance(value, list):
            raise InputValidationError(f"Downstream plan {name} must be a list.")
        groups[name] = value
    if sum(map(len, groups.values())) > _MAX_ROWS:
        raise InputValidationError("Downstream plan exceeds the 4096-row limit.")
    has_tr = any(isinstance(row, Mapping) and row.get("track") == "tr"
                 for rows in groups.values() for row in rows)
    catalog = _configured_catalog(config, source.parent) if has_tr else None
    if has_tr:
        tr = config.get("tr", {})
        if not isinstance(tr, Mapping):
            raise InputValidationError("config.tr must be a mapping.")
        if tr.get("catalog_reference_build") is not None and _build(
            tr["catalog_reference_build"], "TR catalog reference build"
        ) != build:
            raise InputValidationError("Configured TR catalog reference build mismatch.")
    parsed: dict[str, tuple[DownstreamPlanRow, ...]] = {}
    for group, values in groups.items():
        rows = tuple(
            _row(value, group, samples, build, config, source.parent, catalog)
            for value in values
        )
        keys = [(row.sample_id, row.track) for row in rows]
        if len(keys) != len(set(keys)):
            raise InputValidationError(f"Downstream plan {group} contains duplicate sample/track rows.")
        parsed[group] = rows
    _cohort_coverage(parsed["cohort"], samples, config)
    tr_ids = {row.catalog_id for rows in parsed.values() for row in rows if row.track == "tr"}
    if None in tr_ids or len(tr_ids) > 1:
        raise InputValidationError("All TR rows require the same explicit catalog SHA256 identity.")
    return DownstreamSourcePlan(
        schema_version=1,
        path=source,
        sha256=hashlib.sha256(raw).hexdigest(),
        annotation=parsed["annotation"],
        cohort=parsed["cohort"],
    )


def _row(
    value: Any,
    group: str,
    samples: tuple[str, ...],
    build: str,
    config: Mapping[str, Any],
    base: Path,
    catalog: str | None,
) -> DownstreamPlanRow:
    row = _mapping(value, f"Downstream {group} row")
    _unknown_keys(row, _ROW_KEYS, f"Downstream {group} row")
    sample = _text(row.get("sample_id"), "sample_id")
    validate_sample_id(sample)
    if sample not in samples:
        raise InputValidationError("Downstream row identifies an unknown sample.")
    track = _text(row.get("track"), "track")
    allowed = {"small", "sv", "tr"} if group in {"annotation", "review", "benchmark"} else {"small_variants", "sv", "tr"}
    if track not in allowed:
        raise InputValidationError(f"Unsupported downstream {group} track.")
    _destination_enabled(config, group, track)
    origin = _text(row.get("origin"), "origin")
    if origin not in {"GENERATED", "EXTERNAL"}:
        raise InputValidationError("Downstream origin must be GENERATED or EXTERNAL.")
    state = _text(row.get("state", "CALLED"), "state")
    try:
        SampleCallState(state)
    except ValueError as error:
        raise InputValidationError("Unknown downstream cohort sample state.") from error
    if group in {"annotation", "review", "benchmark"} and state != "CALLED":
        raise InputValidationError(f"{group.title()} source rows require state CALLED.")
    row_build = _build(row.get("reference_build", build), "row reference_build")
    if row_build != build:
        raise InputValidationError("Downstream row reference build mismatch.")
    ids = _variant_ids(row.get("variant_ids", []))
    producer = role = tool = version = None
    vcf = index = bundle = None
    if origin == "GENERATED":
        if state != "CALLED":
            raise InputValidationError("Generated source requests require state CALLED; execution is not yet certified.")
        if any(row.get(key) is not None for key in ("path", "index_path", "source_tool", "source_version", "bnd_bundle_manifest")):
            raise InputValidationError("Generated rows cannot declare future paths or tool/version receipts.")
        producer = _text(row.get("producer"), "producer")
        role = _text(row.get("output_role"), "output_role")
        _generated_selector(producer, role, track, group, config)
    elif state in _CALLABLE_STATES:
        if any(row.get(key) is not None for key in ("producer", "output_role")):
            raise InputValidationError("External rows must not declare a generated producer/output role.")
        if "reference_build" not in row:
            raise InputValidationError("Callable external rows require an explicit reference_build.")
        vcf = _existing_file(row.get("path"), base, "external source path")
        index = _existing_file(row.get("index_path"), base, "external index path")
        try:
            if vcf == index or os.path.samefile(vcf, index):
                raise InputValidationError("External source and index must be distinct files.")
        except OSError as error:
            raise InputValidationError("External source/index changed during path validation.") from error
        tool = _text(row.get("source_tool"), "source_tool")
        version = _text(row.get("source_version"), "source_version")
        if row.get("bnd_bundle_manifest") is not None:
            if group not in {"annotation", "review"} or track != "sv":
                raise InputValidationError("A BND bundle source is only an explicit SV annotation or review role.")
            bundle = _existing_file(row["bnd_bundle_manifest"], base, "BND bundle manifest")
    else:
        if any(row.get(key) is not None for key in (
            "producer", "output_role", "path", "index_path", "source_tool", "source_version", "bnd_bundle_manifest",
        )):
            raise InputValidationError("Noncallable cohort rows must not carry source paths or execution metadata.")
        if ids:
            raise InputValidationError("Noncallable cohort rows must not carry variant IDs.")
    identity = _catalog_id(row.get("catalog_id"))
    if track != "tr" and identity is not None:
        raise InputValidationError("catalog_id is only valid for TR rows.")
    if track == "tr":
        if origin == "GENERATED":
            if catalog is None:
                raise InputValidationError("Generated TR requires readable config.tr.catalog.")
            if identity is not None and identity != catalog:
                raise InputValidationError("TR catalog identity differs from config.tr.catalog SHA256.")
            identity = catalog
        elif state in _CALLABLE_STATES and identity is None:
            raise InputValidationError("Callable external TR rows require an explicit catalog_id SHA256.")
        elif identity is None:
            identity = catalog
        if catalog is not None and identity != catalog:
            raise InputValidationError("TR catalog identity differs from config.tr.catalog SHA256.")
    return DownstreamPlanRow(
        sample, track, origin, producer, role, vcf, index, tool, version,
        row_build, identity, state, ids, bundle,
    )


def _generated_selector(
    producer: str, role: str, track: str, group: str, config: Mapping[str, Any],
) -> None:
    if producer not in _PRODUCERS:
        raise InputValidationError("Unknown downstream generated producer.")
    category, roles = _PRODUCERS[producer]
    expected = "small" if track == "small_variants" else track
    if category != expected or role not in roles:
        raise InputValidationError("Generated producer/output role is incompatible with the requested category.")
    if (group == "cohort" and track == "small_variants") != (role == "gvcf"):
        raise InputValidationError("gVCF is required only for cohort small_variants, never substituted by an ordinary VCF.")
    section, backend = producer.split(".")
    if producer == "sv.harmonization" and role == "bnd_vcf":
        from hifivar.breakends import BndPolicy
        if group not in {"annotation", "review"} or track != "sv":
            raise InputValidationError("Generated BND derived roles require the SV annotation or review destination, not cohort.")
        if BndPolicy.from_config(config.get("sv", {}).get("harmonization", {}).get("bnd_representation", {})) is None:
            raise InputValidationError("Generated bnd_vcf requires explicitly enabled fixed BND representation policy.")
    if producer != "sv.harmonization" and not _enabled(config, section):
        raise InputValidationError("Requested generated producer is disabled in config.")
    if section in {"sv", "assembly_sv"} and not _enabled(config, section, backend):
        raise InputValidationError("Requested generated caller/backend is disabled in config.")
    if section == "phasing" and not _enabled(config, "small"):
        raise InputValidationError("Generated HiPhase source requires enabled small-variant calling.")
    if section == "assembly_sv":
        # GENERATED describes the new caller output, not necessarily a new
        # hifiasm execution. Ready FASTAs have their own launcher validation;
        # source intent must still select exactly one assembly input mode.
        ready_manifest = config["assembly_sv"].get("input_manifest")
        new_assembly = _enabled(config, "assembly")
        if ready_manifest is not None:
            _text(ready_manifest, "assembly_sv.input_manifest")
            if new_assembly:
                raise InputValidationError("Generated assembly-SV source cannot combine ready inputs and enabled assembly.")
        elif not new_assembly:
            raise InputValidationError("Generated assembly-SV source requires enabled assembly or an explicit ready input_manifest.")


def _destination_enabled(config: Mapping[str, Any], group: str, track: str) -> None:
    if not _enabled(config, group):
        raise InputValidationError(f"Downstream {group} destination is disabled in config.")
    if group == "cohort" and not _enabled(config, "cohort", track):
        raise InputValidationError("Downstream cohort track is disabled in config.")


def _cohort_coverage(
    rows: tuple[DownstreamPlanRow, ...], samples: tuple[str, ...], config: Mapping[str, Any],
) -> None:
    if not _enabled(config, "cohort"):
        return
    for track in ("small_variants", "sv", "tr"):
        if _enabled(config, "cohort", track):
            observed = {row.sample_id for row in rows if row.track == track}
            if observed != set(samples):
                raise InputValidationError("Every enabled cohort track must explicitly cover every declared sample.")


def _enabled(config: Mapping[str, Any], section: str, child: str | None = None) -> bool:
    value = config.get(section, {})
    if not isinstance(value, Mapping):
        raise InputValidationError(f"config.{section} must be a mapping.")
    if child is not None:
        value = value.get(child, {})
        if not isinstance(value, Mapping):
            raise InputValidationError(f"config.{section}.{child} must be a mapping.")
    enabled = value.get("enabled", False)
    if not isinstance(enabled, bool):
        raise InputValidationError("Configuration enabled flags must be booleans.")
    return enabled


def _configured_catalog(config: Mapping[str, Any], base: Path) -> str | None:
    tr = config.get("tr", {})
    if not isinstance(tr, Mapping):
        raise InputValidationError("config.tr must be a mapping.")
    value = tr.get("catalog")
    if value is None:
        return None
    sources = getattr(config, "sources", None)
    if isinstance(sources, Mapping) and sources.get("user") is not None:
        base = Path(sources["user"]).resolve().parent
    path = _existing_file(value, base, "config.tr.catalog")
    digest = hashlib.sha256()
    try:
        before = path.stat()
        with path.open("rb") as source:
            for chunk in iter(lambda: source.read(1024 * 1024), b""):
                digest.update(chunk)
        if _stat_identity(before) != _stat_identity(path.stat()):
            raise InputValidationError("TR catalog changed while its identity was read.")
    except OSError as error:
        raise InputValidationError("Unable to read configured TR catalog.") from error
    return "sha256:" + digest.hexdigest()


def _catalog_id(value: Any) -> str | None:
    if value is None:
        return None
    match = _CATALOG_SHA.fullmatch(_text(value, "catalog_id"))
    if match is None:
        raise InputValidationError("catalog_id must be a content SHA256, not a catalog name/path.")
    return "sha256:" + match.group(1).lower()


def _variant_ids(value: Any) -> tuple[str, ...]:
    if not isinstance(value, list):
        raise InputValidationError("variant_ids must be a list of strings.")
    result = tuple(_text(item, "variant ID") for item in value)
    if any(";" in item for item in result) or len(result) != len(set(result)):
        raise InputValidationError("Variant IDs must be unique and must not contain a manifest delimiter.")
    return result


def _sample_ids(value: Sequence[str]) -> tuple[str, ...]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise InputValidationError("sample_ids must be a nonempty ordered sequence.")
    samples = tuple(value)
    for sample in samples:
        validate_sample_id(sample)
    if len(samples) != len(set(samples)):
        raise InputValidationError("sample_ids must not contain duplicates.")
    return samples


def _build(value: Any, label: str) -> str:
    text = _text(value, label)
    try:
        result = normalize_reference_build(text)
    except ReferenceError as error:
        raise InputValidationError(f"Invalid {label}.") from error
    assert result is not None
    return result


def _text(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value or value != value.strip():
        raise InputValidationError(f"{label} must be a nonempty string without surrounding whitespace.")
    if any(unicodedata.category(char).startswith("C") or char in "\u2028\u2029" for char in value):
        raise InputValidationError(f"{label} must not contain control or invisible format characters.")
    return value


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or any(not isinstance(key, str) for key in value):
        raise InputValidationError(f"{label} must be a string-keyed mapping.")
    return value


def _unknown_keys(value: Mapping[str, Any], allowed: frozenset[str], label: str) -> None:
    if set(value).difference(allowed):
        raise InputValidationError(f"{label} contains unknown keys.")


def _existing_file(value: Any, base: Path, label: str) -> Path:
    text = _text(value, label)
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = base / path
    try:
        path = path.resolve(strict=True)
        if not path.is_file() or path.stat().st_size == 0:
            raise InputValidationError(f"{label} must identify a nonempty existing regular file.")
        return path
    except (OSError, RuntimeError) as error:
        raise InputValidationError(f"{label} must identify a readable existing file.") from error


def _plan_bytes(value: Path) -> tuple[Path, bytes]:
    if not isinstance(value, (str, Path)):
        raise InputValidationError("Downstream plan path must be a string or Path.")
    path = _existing_file(str(value), Path.cwd(), "Downstream plan path")
    try:
        before = path.stat()
        if before.st_size > _MAX_BYTES:
            raise InputValidationError("Downstream plan exceeds the 1 MiB limit.")
        with path.open("rb") as source:
            raw = source.read(_MAX_BYTES + 1)
        if len(raw) > _MAX_BYTES or _stat_identity(before) != _stat_identity(path.stat()):
            raise InputValidationError("Downstream plan is oversized or changed while being read.")
        return path, raw
    except OSError as error:
        raise InputValidationError("Unable to read downstream plan.") from error


def _stat_identity(stat: os.stat_result) -> tuple[int, int, int, int]:
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns


__all__ = ["DownstreamPlanRow", "DownstreamSourcePlan", "read_downstream_plan"]
