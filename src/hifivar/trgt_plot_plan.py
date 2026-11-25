"""Strict external-input YAML plans for the independent TRGT review consumer.

The seven supplied SHA256 values bind file bytes, not a native producer receipt.
This entry never starts calling, changes a workflow configuration, or completes
a human review. Plotting and its output safety belong to ``trgt_plot``.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math
import os
from pathlib import Path
import re
import stat
from typing import TYPE_CHECKING, Any

import yaml

from hifivar.exceptions import InputValidationError
from hifivar.reference import ReferenceGenome, normalize_reference_build
from hifivar.sample import validate_sample_id
from hifivar.tr import TandemRepeatCatalog

if TYPE_CHECKING:
    from hifivar.trgt_plot import TrgtPlotRequest


_MAX_BYTES = 256 * 1024
_MAX_LINES = 4096
_MAX_DEPTH = 24
_ROLES = frozenset({"vcf", "tbi", "spanning_bam", "bai", "catalog", "fasta", "fai"})
_TOP_KEYS = frozenset({
    "schema_version", "sample_id", "selector", "reference", "catalog", "source",
    "expected_sha256", "runtime", "plot", "execution", "output",
})
_LOCUS = re.compile(r"[A-Za-z0-9_.:-]+\Z")


class _PlotPlanLoader(yaml.SafeLoader):
    """Reject duplicate/non-string keys, aliases, merges and excessive nesting."""

    def compose_node(self, parent: Any, index: Any) -> Any:
        if self.check_event(yaml.AliasEvent):
            raise InputValidationError("TR plot plan YAML aliases are not supported.")
        depth = getattr(self, "_plot_depth", 0)
        if depth >= _MAX_DEPTH:
            raise InputValidationError("TR plot plan YAML is too deeply nested.")
        self._plot_depth = depth + 1
        try:
            return super().compose_node(parent, index)
        finally:
            self._plot_depth = depth

    def construct_mapping(self, node: Any, deep: bool = False) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key_node, value_node in node.value:
            if key_node.tag == "tag:yaml.org,2002:merge":
                raise InputValidationError("TR plot plan YAML merge keys are not supported.")
            key = self.construct_object(key_node, deep=deep)
            if not isinstance(key, str) or not key or any(ord(c) < 32 for c in key):
                raise InputValidationError("TR plot plan mapping keys must be nonempty strings.")
            if key in result:
                raise InputValidationError("TR plot plan contains a duplicate mapping key.")
            result[key] = self.construct_object(value_node, deep=deep)
        return result


@dataclass(frozen=True, slots=True)
class TrgtPlotPlan:
    """A parsed, exact-byte plan and its normalized wrapper request."""

    path: Path
    sha256: str
    executable: str
    request: TrgtPlotRequest
    source_mode: str = "EXTERNAL_FILES"


def _text(value: object, label: str) -> str:
    if (not isinstance(value, str) or not value or value.strip() != value
            or any(ord(character) < 32 for character in value)):
        raise InputValidationError(f"TR plot plan {label} requires a nonempty single-line string.")
    return value


def _mapping(value: object, keys: frozenset[str], label: str) -> dict[str, Any]:
    if not isinstance(value, dict) or any(not isinstance(key, str) for key in value):
        raise InputValidationError(f"TR plot plan {label} requires a mapping with string keys.")
    if set(value) != keys:
        raise InputValidationError(f"TR plot plan {label} has unknown or missing fields.")
    return value


def _path(value: object, base: Path, label: str) -> Path:
    text = _text(value, label)
    try:
        path = Path(text).expanduser()
        return (path if path.is_absolute() else base / path).absolute()
    except (OSError, ValueError, RuntimeError) as error:
        raise InputValidationError(f"TR plot plan {label} is not a usable local path.") from error


def _read_bytes(path: Path) -> bytes:
    """Read bounded regular plan bytes and reject observable concurrent changes."""
    try:
        if not path.is_file():
            raise InputValidationError("TR plot plan must be an existing regular file.")
        with path.open("rb") as handle:
            before = os.fstat(handle.fileno())
            if not stat.S_ISREG(before.st_mode):
                raise InputValidationError("TR plot plan must be a regular file.")
            raw = handle.read(_MAX_BYTES + 1)
            after = os.fstat(handle.fileno())
        current = path.stat()
    except OSError as error:
        raise InputValidationError("Unable to read the TR plot plan.") from error
    signature = lambda item: (item.st_dev, item.st_ino, item.st_size, item.st_mtime_ns)
    if signature(before) != signature(after) or signature(after) != signature(current):
        raise InputValidationError("TR plot plan changed while it was read.")
    if not raw or len(raw) > _MAX_BYTES:
        raise InputValidationError("TR plot plan must be nonempty and at most 256 KiB.")
    return raw


def read_trgt_plot_plan(path: str | Path) -> TrgtPlotPlan:
    """Resolve one strict external-input plan, without invoking a tool or writing."""
    from hifivar.trgt_plot import TrgtPlotRequest

    if not isinstance(path, (str, Path)) or not str(path).strip():
        raise InputValidationError("TR plot plan path must be explicit and nonempty.")
    try:
        source = Path(_text(str(path), "path")).expanduser().absolute()
    except (OSError, ValueError, RuntimeError) as error:
        raise InputValidationError("TR plot plan path is not usable.") from error
    raw = _read_bytes(source)
    digest = hashlib.sha256(raw).hexdigest()
    try:
        text = raw.decode("utf-8-sig")
        if len(text.splitlines()) > _MAX_LINES:
            raise InputValidationError("TR plot plan exceeds the 4096-line limit.")
        document = yaml.load(text, Loader=_PlotPlanLoader)
    except (UnicodeError, yaml.YAMLError, RecursionError, ValueError, OverflowError) as error:
        raise InputValidationError("TR plot plan requires bounded valid UTF-8 YAML.") from error
    document = _mapping(document, _TOP_KEYS, "root")
    if type(document["schema_version"]) is not int or document["schema_version"] != 1:
        raise InputValidationError("TR plot plan schema_version must be integer 1.")
    sample_id = _text(document["sample_id"], "sample_id")
    validate_sample_id(sample_id)
    selector = _text(document["selector"], "selector")
    if (not selector.startswith("TRID:") or not _LOCUS.fullmatch(selector[5:])
            or selector[5:] == "." or selector[5:].startswith("TRID:")):
        raise InputValidationError("TR plot selector requires one explicit TRID:<raw catalog ID>.")
    reference = _mapping(document["reference"], frozenset({"fasta", "fai", "build"}), "reference")
    catalog = _mapping(document["catalog"], frozenset({"path", "build"}), "catalog")
    inputs = _mapping(document["source"], frozenset({"mode", "vcf", "tbi", "spanning_bam", "bai"}), "source")
    if inputs["mode"] != "EXTERNAL_FILES":
        raise InputValidationError("TR plot plans currently support EXTERNAL_FILES only; producer receipts are not verified.")
    hashes = _mapping(document["expected_sha256"], _ROLES, "expected_sha256")
    if any(not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None
           for value in hashes.values()):
        raise InputValidationError("TR plot plan requires exactly seven lowercase SHA256 strings.")
    runtime = _mapping(document["runtime"], frozenset({"executable", "supported_version"}), "runtime")
    if runtime["supported_version"] != "5.1.0":
        raise InputValidationError("TR plot plan supports TRGT version 5.1.0 only.")
    executable = _text(runtime["executable"], "runtime.executable")
    if "/" in executable or "\\" in executable:
        executable = str(_path(executable, source.parent, "runtime.executable"))
    elif any(character.isspace() for character in executable) or executable.startswith("-"):
        raise InputValidationError("TR plot executable must be a single executable name or path, not shell arguments.")
    settings = _mapping(document["plot"], frozenset({"type", "show", "format", "max_allele_reads"}), "plot")
    if (settings["type"], settings["show"], settings["format"]) != ("allele", "motifs", "svg"):
        raise InputValidationError("TR plot plan supports allele/motifs/SVG only.")
    reads = settings["max_allele_reads"]
    if type(reads) is not int or reads <= 0:
        raise InputValidationError("TR plot max_allele_reads must be a positive integer.")
    execution = _mapping(document["execution"], frozenset({"timeout_seconds"}), "execution")
    timeout = execution["timeout_seconds"]
    if (isinstance(timeout, bool) or not isinstance(timeout, (int, float))
            or not 0 < timeout <= 1200 or not math.isfinite(timeout)):
        raise InputValidationError("TR plot timeout_seconds must be finite, positive and no more than 1200.")
    output = _mapping(document["output"], frozenset({"directory"}), "output")
    base = source.parent
    fasta = _path(reference["fasta"], base, "reference.fasta")
    fai = _path(reference["fai"], base, "reference.fai")
    if fai != Path(str(fasta) + ".fai"):
        raise InputValidationError("TR plot reference.fai must be the conventional FASTA.fai.")
    reference_build = _text(reference["build"], "reference.build")
    catalog_build = _text(catalog["build"], "catalog.build")
    if normalize_reference_build(reference_build) != normalize_reference_build(catalog_build):
        raise InputValidationError("TR plot catalog and reference builds must match explicitly.")
    genome = ReferenceGenome.from_fasta(fasta, build=reference_build)
    request = TrgtPlotRequest(
        sample_id=sample_id, repeat_id=selector[5:], reference=genome,
        catalog=TandemRepeatCatalog(_path(catalog["path"], base, "catalog.path"), genome.build),
        vcf=_path(inputs["vcf"], base, "source.vcf"), tbi=_path(inputs["tbi"], base, "source.tbi"),
        spanning_bam=_path(inputs["spanning_bam"], base, "source.spanning_bam"),
        bai=_path(inputs["bai"], base, "source.bai"), expected_sha256=dict(hashes),
        output_dir=_path(output["directory"], base, "output.directory"),
        max_allele_reads=reads, timeout_seconds=timeout, plan_path=source, plan_sha256=digest,
    )
    if hashlib.sha256(_read_bytes(source)).hexdigest() != digest:
        raise InputValidationError("TR plot plan changed while its request was constructed.")
    return TrgtPlotPlan(source, digest, executable, request)


def run_trgt_plot_plan(path: str | Path, *, dry_run: bool = False) -> dict:
    """Run only the explicit plotting consumer, never calling or human review."""
    from hifivar.trgt_plot import TrgtPlotWrapper

    if not isinstance(dry_run, bool):
        raise InputValidationError("TR plot dry_run must be boolean.")
    plan = read_trgt_plot_plan(path)
    if hashlib.sha256(_read_bytes(plan.path)).hexdigest() != plan.sha256:
        raise InputValidationError("TR plot plan changed before execution.")
    return TrgtPlotWrapper(command_prefix=(plan.executable,)).run(plan.request, dry_run=dry_run)


__all__ = ["TrgtPlotPlan", "read_trgt_plot_plan", "run_trgt_plot_plan"]
