"""Fresh-directory, single-locus native TRGT 5.1.0 SVG rendering.

This consumer does not call/genotype repeats, establish caller lineage, assess
accuracy or perform a human review. A completion receipt is published last;
failed directories and partial files remain available and are never recycled.
"""
from __future__ import annotations

from dataclasses import dataclass
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import stat
from types import MappingProxyType
from typing import Mapping, Sequence
import xml.etree.ElementTree as ET

from hifivar import __version__
from hifivar.command import CommandRunner
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar.logging_utils import get_logger
from hifivar.reference import ReferenceGenome
from hifivar.sample import validate_sample_id
from hifivar.tr import TandemRepeatCatalog, _validate_bgzf, _validate_tabix_index, _validate_trgt_vcf_header


_ROLES = frozenset({"vcf", "tbi", "spanning_bam", "bai", "catalog", "fasta", "fai"})
_SHA = re.compile(r"[0-9a-f]{64}")
_LOCUS = re.compile(r"[A-Za-z0-9_.:-]+")
_VERSION = re.compile(r"^[ \t]*trgt\b[ \t]*(?:version[ \t]+)?v?([^\r\n]*)\r?$", re.I | re.M)
_NAMES = ("plot.svg", "completion.json", "command.stderr.log", "command.stdout.log",
          "effective_input.json", "version.stderr.log", "version.stdout.log", ".completion.pending.json", "failure.json")
_PathInfo = tuple[Path, tuple[int, int] | None]
_LOGGER = get_logger(__name__)


def _path(value: str | Path, label: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise InputValidationError(f"TR plot {label} requires a nonempty path.")
    return Path(value).expanduser().absolute()


def _fingerprint(path: Path, *, allow_empty: bool = False) -> dict[str, str | int]:
    try:
        before = path.stat()
        if not stat.S_ISREG(before.st_mode) or (not allow_empty and before.st_size == 0):
            raise InputValidationError(f"TR plot input/tool is not a nonempty regular file: '{path}'.")
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
                digest.update(block)
        after = path.stat()
        if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns):
            raise InputValidationError("TR plot protected file changed during hashing.")
        return {"path": str(path), "bytes": after.st_size, "sha256": digest.hexdigest()}
    except OSError as error:
        raise InputValidationError(f"Cannot fingerprint TR plot input/tool '{path}': {error}") from error


@dataclass(frozen=True)
class TrgtPlotRequest:
    """Explicit selected external files; seven hashes are not producer proof."""

    sample_id: str
    repeat_id: str
    reference: ReferenceGenome
    catalog: TandemRepeatCatalog
    vcf: Path
    tbi: Path
    spanning_bam: Path
    bai: Path
    expected_sha256: Mapping[str, str]
    output_dir: Path
    max_allele_reads: int = 20
    timeout_seconds: float = 120
    plan_path: Path | None = None
    plan_sha256: str | None = None

    def __post_init__(self) -> None:
        validate_sample_id(self.sample_id)
        if not isinstance(self.repeat_id, str) or not _LOCUS.fullmatch(self.repeat_id) or self.repeat_id.startswith("TRID:"):
            raise InputValidationError("TR plot repeat_id must be one raw catalog ID, without the TRID: selector prefix.")
        if not isinstance(self.reference, ReferenceGenome) or not isinstance(self.catalog, TandemRepeatCatalog):
            raise InputValidationError("TR plot requires explicit reference and catalog objects.")
        if not self.reference.build or self.catalog.reference_build != self.reference.build:
            raise InputValidationError("TR plot reference/catalog require the same declared build.")
        if self.reference.fai.expanduser().absolute() != Path(str(self.reference.fasta.expanduser().absolute()) + ".fai"):
            raise InputValidationError("TR plot FAI must be the conventional FASTA.fai used by native TRGT.")
        if isinstance(self.max_allele_reads, bool) or not isinstance(self.max_allele_reads, int) or self.max_allele_reads <= 0:
            raise InputValidationError("TR plot max_allele_reads must be a positive integer.")
        if (isinstance(self.timeout_seconds, bool) or not isinstance(self.timeout_seconds, (int, float))
                or not 0 < self.timeout_seconds <= 1200 or not math.isfinite(self.timeout_seconds)):
            raise InputValidationError("TR plot timeout_seconds must be finite, positive and no more than 1200.")
        for name in ("vcf", "tbi", "spanning_bam", "bai", "output_dir"):
            object.__setattr__(self, name, _path(getattr(self, name), name))
        if (not self.vcf.name.endswith(".vcf.gz") or self.tbi != Path(str(self.vcf) + ".tbi")
                or not self.spanning_bam.name.endswith(".bam") or self.bai != Path(str(self.spanning_bam) + ".bai")):
            raise InputValidationError("TR plot requires explicit conventional VCF/TBI and spanning BAM/BAI paths.")
        if not isinstance(self.expected_sha256, Mapping):
            raise InputValidationError("TR plot expected_sha256 must be a seven-role mapping.")
        values = dict(self.expected_sha256)
        if set(values) != _ROLES or any(not isinstance(v, str) or not _SHA.fullmatch(v) for v in values.values()):
            raise InputValidationError("TR plot requires exactly seven explicit lowercase SHA256 bindings.")
        object.__setattr__(self, "expected_sha256", MappingProxyType(values))
        if (self.plan_path is None) != (self.plan_sha256 is None):
            raise InputValidationError("TR plot plan_path and parsed raw-byte plan_sha256 must be supplied together.")
        if self.plan_path is not None:
            object.__setattr__(self, "plan_path", _path(self.plan_path, "plan"))
            if not isinstance(self.plan_sha256, str) or not _SHA.fullmatch(self.plan_sha256):
                raise InputValidationError("TR plot plan_sha256 must be one lowercase SHA256.")

    @property
    def inputs(self) -> dict[str, Path]:
        """Seven immutable file roles, deliberately separate from plan identity."""
        return {"vcf": self.vcf, "tbi": self.tbi, "spanning_bam": self.spanning_bam, "bai": self.bai,
                "catalog": self.catalog.path.expanduser().absolute(),
                "fasta": self.reference.fasta.expanduser().absolute(), "fai": self.reference.fai.expanduser().absolute()}


def _inspect(path: Path) -> _PathInfo:
    try:
        resolved = path.resolve()
        try:
            value = path.stat()
        except FileNotFoundError:
            if path.is_symlink():
                raise OutputValidationError("TR_PLOT_PATH_PROTECTION: dangling symlink.")
            return resolved, None
        return resolved, (value.st_dev, value.st_ino)
    except (OSError, RuntimeError, ValueError) as error:
        raise OutputValidationError(f"TR_PLOT_PATH_PROTECTION: cannot inspect path: {error}") from error


def _overlap(a: _PathInfo, b: _PathInfo) -> bool:
    return a[0] == b[0] or a[0] in b[0].parents or b[0] in a[0].parents or (a[1] is not None and a[1] == b[1])


def _protect(request: TrgtPlotRequest, tools: Mapping[str, Path], *, owner: _PathInfo | None = None) -> None:
    directory = request.output_dir
    for path in (directory, *directory.parents):
        if path.is_symlink() or getattr(path, "is_junction", lambda: False)():
            raise OutputValidationError("TR_PLOT_PATH_PROTECTION: output ancestry must not contain symlinks or junctions.")
    if not directory.parent.is_dir():
        raise OutputValidationError("TR plot output parent must already exist; only the requested fresh directory is created.")
    if owner is None:
        if directory.exists():
            raise OutputValidationError("TR_PLOT_OUTPUT_EXISTS: require a nonexistent output directory; no overwrite or resume.")
    elif not directory.is_dir() or _inspect(directory) != owner:
        raise OutputValidationError("TR_PLOT_PATH_PROTECTION: the owned output directory changed.")
    outputs = [_inspect(directory), *(_inspect(directory / name) for name in _NAMES)]
    resources = [*request.inputs.values(), *tools.values()]
    if request.plan_path is not None:
        resources.append(request.plan_path)
    for resource in resources:
        if any(_overlap(_inspect(resource), output) for output in outputs):
            raise OutputValidationError("TR_PLOT_PATH_PROTECTION: output directory/file overlaps an immutable input, plan or tool.")


def _reference_slice(reference: ReferenceGenome, contig: str, start: int, end: int) -> str:
    """Read [start,end], 1-based inclusive: VCF anchor plus BED repeat interval."""
    entries = []
    with reference.fai.open("r", encoding="utf-8") as stream:
        for line in stream:
            fields = line.rstrip("\r\n").split("\t")
            if fields[0] == contig:
                if len(fields) != 5:
                    raise InputValidationError("TR plot selected FAI entry requires exactly five fields.")
                entries.append(tuple(int(x) for x in fields[1:]))
    if len(entries) != 1:
        raise InputValidationError("TR plot selected FASTA contig has no unique FAI entry.")
    length, offset, line_bases, line_bytes = entries[0]
    known = {item.name: item.length for item in reference.contigs}
    if (length != known.get(contig) or not 1 <= start <= end <= length or offset < 0
            or line_bases <= 0 or line_bytes < line_bases or line_bytes - line_bases not in {0, 1, 2}):
        raise InputValidationError("TR plot selected FAI/anchor interval is inconsistent.")
    pieces = []
    current = start - 1
    with reference.fasta.open("rb") as stream:
        while current < end:
            row, column = divmod(current, line_bases)
            count = min(end - current, line_bases - column)
            stream.seek(offset + row * line_bytes + column)
            part = stream.read(count)
            if len(part) != count or not re.fullmatch(b"[ACGTNacgtn]+", part):
                raise InputValidationError("TR plot FASTA span is truncated or has unsupported bases/index geometry.")
            pieces.append(part)
            current += count
    return b"".join(pieces).decode("ascii")


def _selection(request: TrgtPlotRequest) -> dict[str, str | int]:
    request.catalog.validate(request.reference)
    selected = []
    with request.catalog.path.open("r", encoding="utf-8") as stream:
        for line in stream:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            pairs = [x.split("=", 1) for x in fields[3].split(";") if "=" in x]
            tags = dict(pairs)
            if tags["ID"] == request.repeat_id:
                if len(pairs) != len(tags):
                    raise InputValidationError("TR plot selected catalog contains duplicate annotation keys.")
                selected.append((fields[0], int(fields[1]), int(fields[2]), tags["MOTIFS"], tags["STRUC"]))
    if len(selected) != 1:
        raise InputValidationError("TR plot requires exactly one matching catalog locus.")
    contig, bed_start, end, motifs, structure = selected[0]
    if bed_start < 1:
        raise InputValidationError("TR plot requires a real preceding VCF anchor base.")
    reference_bases = _reference_slice(request.reference, contig, bed_start, end)
    _validate_bgzf(request.vcf, "TR plot VCF")
    _validate_tabix_index(request.tbi)
    _validate_trgt_vcf_header(request.vcf, request.sample_id, set(request.reference.contig_names))
    matches = []
    header = False
    with gzip.open(request.vcf, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#CHROM\t"):
                if header or line.rstrip("\r\n").split("\t")[9:] != [request.sample_id]:
                    raise InputValidationError("TR plot requires exactly its one declared VCF sample.")
                header = True
            elif not line.startswith("#") and line.strip():
                if not header:
                    raise InputValidationError("TR plot VCF records precede the sample header.")
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 10:
                    raise InputValidationError("TR plot single-sample VCF row width mismatch.")
                pairs = [x.split("=", 1) for x in fields[7].split(";") if "=" in x]
                info = dict(pairs)
                if info.get("TRID") != request.repeat_id:
                    continue
                if (len(pairs) != len(info)
                        or (fields[0], int(fields[1]), int(info.get("END", -1)), info.get("MOTIFS")) != (contig, bed_start, end, motifs)):
                    raise InputValidationError("TR plot selected VCF/catalog coordinates or motifs differ.")
                if not info.get("STRUC") or len(fields[3]) != end - bed_start + 1:
                    raise InputValidationError("TR plot REF must contain the repeat plus one preceding anchor.")
                if info["STRUC"] != structure:
                    raise InputValidationError("TR plot selected VCF STRUC differs from catalog STRUC.")
                if fields[3].upper() != reference_bases.upper():
                    raise InputValidationError("TR plot selected VCF REF differs from actual FASTA, including its anchor.")
                keys, values = fields[8].split(":"), fields[9].split(":")
                if len(set(keys)) != len(keys) or "GT" not in keys or len(values) != len(keys):
                    raise InputValidationError("TR plot selected repeat needs an unambiguous GT field.")
                gt = values[keys.index("GT")]
                if not re.fullmatch(r"\d+(?:[/|]\d+)?", gt):
                    raise InputValidationError("TR plot needs a fully called haploid/diploid GT.")
                alleles = [fields[3], *fields[4].split(",")]
                indexes = [int(x) for x in re.split(r"[/|]", gt)]
                if any(i >= len(alleles) or not re.fullmatch("[ACGTNacgtn]+", alleles[i]) for i in indexes):
                    raise InputValidationError("TR plot selected GT needs explicit available sequence alleles.")
                matches.append({"repeat_id": request.repeat_id, "contig": contig, "bed_start": bed_start,
                                "bed_end": end, "GT": gt, "MOTIFS": motifs, "STRUC": structure,
                                "reference_validation": "PASS_SELECTED_ANCHOR_PLUS_REPEAT_FASTA_SEQUENCE",
                                "selected_reference_sha256": hashlib.sha256(reference_bases.upper().encode("ascii")).hexdigest()})
    if not header or len(matches) != 1:
        raise InputValidationError("TR plot VCF must contain exactly one selected TRID record.")
    return matches[0]


def _validate_svg(path: Path) -> None:
    """Conservative static SVG only; this is not rendered-image or manual QA."""
    if not path.is_file() or path.is_symlink() or not 0 < path.stat().st_size <= 16 * 1024 * 1024:
        raise OutputValidationError("TR plot SVG must be a nonempty regular file no larger than 16 MiB.")
    raw = path.read_bytes()
    if b"<!DOCTYPE" in raw.upper() or b"<!ENTITY" in raw.upper():
        raise OutputValidationError("TR plot SVG entity/doctype declarations are prohibited.")
    without_declaration = re.sub(br"^\s*<\?xml\s+[^?]*\?>", b"", raw, count=1)
    if b"<?" in without_declaration:
        raise OutputValidationError("TR plot SVG processing instructions are prohibited.")
    try:
        root = ET.fromstring(raw)
    except ET.ParseError as error:
        raise OutputValidationError("TR plot output is not parseable SVG.") from error
    namespace = "{http://www.w3.org/2000/svg}"
    tags = {"svg", "g", "rect", "line", "polygon", "polyline", "path", "circle", "ellipse", "text", "tspan", "title", "desc"}
    attributes = {"width", "height", "viewBox", "preserveAspectRatio", "x", "y", "x1", "y1", "x2", "y2",
                  "points", "d", "cx", "cy", "r", "rx", "ry", "dx", "dy", "transform", "fill", "stroke",
                  "stroke-width", "stroke-linecap", "stroke-linejoin", "stroke-dasharray", "stroke-miterlimit",
                  "fill-rule", "clip-rule", "opacity", "fill-opacity", "stroke-opacity", "font-family", "font-weight",
                  "font-size", "font-style", "text-anchor", "dominant-baseline", "alignment-baseline", "letter-spacing", "id"}
    if root.tag != namespace + "svg":
        raise OutputValidationError("TR plot root is not namespace-qualified SVG.")
    for node in root.iter():
        if node.tag not in {namespace + name for name in tags} or set(node.attrib).difference(attributes):
            raise OutputValidationError("TR plot SVG contains unsupported/active elements or attributes.")
        for key, value in node.attrib.items():
            if re.search(r"url\s*\(|@import|javascript:|data:", value, re.I):
                raise OutputValidationError("TR plot SVG resource/CSS references are prohibited.")
            if key in {"fill", "stroke"} and not re.fullmatch(r"#[0-9A-Fa-f]{3,8}|[A-Za-z]+|rgba?\([0-9.,%\s]+\)", value):
                raise OutputValidationError("TR plot SVG paint value is not a literal color.")
    if not any(node.tag in {namespace + tag for tag in tags - {"svg", "g", "title", "desc"}} for node in root.iter()):
        raise OutputValidationError("TR plot SVG contains no graphical content.")


def _json_bytes(data: Mapping[str, object]) -> bytes:
    return (json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")


def _write_new(path: Path, data: Mapping[str, object]) -> _PathInfo:
    encoded = _json_bytes(data)
    try:
        with path.open("xb") as out:
            out.write(encoded)
            out.flush()
            os.fsync(out.fileno())
            created = os.fstat(out.fileno())
    except OSError as error:
        raise OutputValidationError(f"TR plot refuses an existing/unwritable evidence destination '{path}': {error}") from error
    return path.resolve(), (created.st_dev, created.st_ino)


def _reserve(path: Path) -> _PathInfo:
    try:
        with path.open("xb"):
            pass
    except OSError as error:
        raise OutputValidationError(f"TR plot refuses an existing/late output or log '{path}': {error}") from error
    return _inspect(path)


def _owned(path: Path, expected: _PathInfo) -> None:
    if path.is_symlink() or not path.is_file() or _inspect(path) != expected:
        raise OutputValidationError("TR_PLOT_PATH_PROTECTION: owned output/log identity changed.")


class TrgtPlotWrapper:
    """Render one selected repeat into a fresh directory; never overwrite/retry."""

    def __init__(self, command_prefix: Sequence[str] = ("trgt",), *, runner: CommandRunner | None = None) -> None:
        if (not isinstance(command_prefix, Sequence) or isinstance(command_prefix, (str, bytes)) or not command_prefix
                or any(not isinstance(x, str) or not x.strip() or "\x00" in x for x in command_prefix)):
            raise InputValidationError("TR plot command_prefix requires a nonempty argument sequence, not a shell command.")
        self.prefix = tuple(command_prefix)
        self.runner = runner or CommandRunner()

    def _tools(self) -> dict[str, Path]:
        files = {}
        for index, item in enumerate(self.prefix):
            found = shutil.which(item) if index == 0 else None
            path = Path(found or item).expanduser().absolute()
            if path.is_file() or path.is_symlink() or Path(item).expanduser().is_absolute():
                files[f"prefix_file_{index}"] = path
        return files

    def command(self, request: TrgtPlotRequest) -> list[str]:
        """Official 5.1.0 plot interface; no threads or arbitrary native flags."""
        return [*self.prefix, "plot", "--genome", str(request.inputs["fasta"]), "--repeats", str(request.inputs["catalog"]),
                "--vcf", str(request.vcf), "--spanning-reads", str(request.spanning_bam), "--repeat-id", request.repeat_id,
                "--image", str(request.output_dir / "plot.svg"), "--plot-type", "allele", "--show", "motifs",
                "--max-allele-reads", str(request.max_allele_reads)]

    def _raw_snapshots(self, request: TrgtPlotRequest, tools: Mapping[str, Path]) -> dict[str, object]:
        return {"inputs": {role: _fingerprint(path) for role, path in request.inputs.items()},
                "tools": {name: _fingerprint(path) for name, path in tools.items() if path.exists() or path.is_symlink()},
                "plan": _fingerprint(request.plan_path) if request.plan_path is not None else None}

    def _snapshots(self, request: TrgtPlotRequest, tools: Mapping[str, Path]) -> dict[str, object]:
        snapshot = self._raw_snapshots(request, tools)
        if any(snapshot["inputs"][role]["sha256"] != expected for role, expected in request.expected_sha256.items()):
            raise InputValidationError("TR plot selected seven input SHA256 identities differ from the explicit bindings.")
        plan = snapshot["plan"]
        if plan is not None and plan["sha256"] != request.plan_sha256:
            raise InputValidationError("TR plot plan bytes changed after parsing; refusing a mismatched effective request.")
        return snapshot

    def run(self, request: TrgtPlotRequest, *, dry_run: bool = False) -> dict[str, object]:
        """Execute via CommandRunner; no completion exists on a precommit failure.

        The exclusively claimed directory is not an OS sandbox against hostile
        concurrent mutation. There is no paired-file overwrite transaction,
        recursive cleanup, automatic retry or SIGKILL recovery guarantee. A
        fully committed receipt can survive an interruption before API return;
        that is not a normal successful return and needs independent readback.
        """
        if not isinstance(request, TrgtPlotRequest) or not isinstance(dry_run, bool):
            raise InputValidationError("TR plot requires TrgtPlotRequest and a boolean dry_run.")
        tools = self._tools()
        _protect(request, tools)
        before = self._snapshots(request, tools)
        try:
            selection = _selection(request)
            _validate_bgzf(request.spanning_bam, "TR plot spanning BAM")
            with gzip.open(request.spanning_bam, "rb") as bam:
                if bam.read(4) != b"BAM\x01":
                    raise InputValidationError("TR plot spanning BAM has invalid BAM magic.")
            with request.bai.open("rb") as index:
                if index.read(4) != b"BAI\x01":
                    raise InputValidationError("TR plot spanning index has invalid BAI magic.")
        except (OSError, EOFError, UnicodeError, KeyError, ValueError) as error:
            raise InputValidationError("Cannot validate the selected TR plot catalog/VCF/reference/BAM record.") from error
        argv = self.command(request)
        effective = {"schema_version": 1, "track": "tr.plot", "status": "PLANNED_NOT_COMPLETE",
                     "sample_id": request.sample_id, "repeat_id": request.repeat_id, "selection": selection,
                     "command": argv, "input_identity": before, "reference_build": request.reference.build,
                     "plot_type": "allele", "show": "motifs", "max_allele_reads": request.max_allele_reads,
                     "timeout_seconds": request.timeout_seconds,
                     "source_origin": "EXTERNAL_FILES_NATIVE_PRODUCER_RECEIPT_NOT_VERIFIED",
                     "manual_review": "NOT_REVIEWED", "accuracy": "NOT_ESTABLISHED",
                     "allele_panel_order": "NATIVE_DISPLAY_ORDER_NOT_PARENTAL_OR_PHASE_ACCURACY"}
        if dry_run:
            self.runner.run(argv, dry_run=True, timeout=request.timeout_seconds)
            return {**effective, "status": "PLANNED", "executed": False, "completion_path": str(request.output_dir / "completion.json")}
        self.runner.require_executable(self.prefix[0])
        _protect(request, tools)
        if self._tools() != tools or self._snapshots(request, tools) != before:
            raise InputValidationError("TR plot protected files/tool resolution changed before claiming output.")
        try:
            request.output_dir.mkdir(exist_ok=False)
        except OSError as error:
            raise OutputValidationError(f"TR plot requires a fresh exclusively created output directory: {error}") from error
        owner = _inspect(request.output_dir)
        try:
            return self._execute(request, tools, before, selection, effective, argv, owner)
        except BaseException as error:
            # Failure diagnostics are not cleanup/recovery. The original error,
            # including KeyboardInterrupt/SystemExit, is always re-raised.
            failure = {"schema_version": 1, "track": "tr.plot", "status": "FAILED_OR_INTERRUPTED_NO_NORMAL_RETURN",
                       "failure_class": type(error).__name__, "completion_present": (request.output_dir / "completion.json").exists()}
            try:
                after = self._raw_snapshots(request, self._tools())
                failure.update(protected_after=after, protected_identity="PASS" if after == before else "FAIL")
            except Exception as diagnostic:
                failure.update(protected_identity="NOT_VERIFIABLE", posthash_failure_class=type(diagnostic).__name__)
            try:
                _protect(request, tools, owner=owner)
                _write_new(request.output_dir / "failure.json", failure)
            except Exception as diagnostic:
                _LOGGER.error("TR plot failure evidence could not be appended (%s); original %s remains authoritative.",
                              type(diagnostic).__name__, type(error).__name__)
            raise

    def _execute(
        self, request: TrgtPlotRequest, tools: Mapping[str, Path], before: dict[str, object],
        selection: dict[str, str | int], effective: dict[str, object], argv: list[str], owner: _PathInfo,
    ) -> dict[str, object]:
        _write_new(request.output_dir / "effective_input.json", effective)
        version_out, version_err = (request.output_dir / name for name in ("version.stdout.log", "version.stderr.log"))
        owned_version = {path: _reserve(path) for path in (version_out, version_err)}
        _protect(request, tools, owner=owner)
        for path, expected in owned_version.items():
            _owned(path, expected)
        version = self.runner.run([*self.prefix, "--version"], timeout=request.timeout_seconds,
                                  stdout_path=version_out, stderr_path=version_err)
        if not version.executed or version.returncode != 0:
            raise OutputValidationError("TR plot version probe did not execute successfully.")
        for path, expected in owned_version.items():
            _owned(path, expected)
        if any(path.stat().st_size > 65536 for path in owned_version):
            raise ToolVersionError("TR plot version output exceeds the bounded 64 KiB parser.")
        version_text = "\n".join(path.read_text(encoding="utf-8", errors="replace") for path in (version_out, version_err))
        reported_versions = {value.strip() for value in _VERSION.findall(version_text)}
        if len(reported_versions) != 1 or not reported_versions <= {"5.1.0", "5.1.0-ec66463"}:
            raise ToolVersionError("TR plot supports only the exact TRGT version reports 5.1.0 or 5.1.0-ec66463.")
        reported_version = next(iter(reported_versions))
        if self._snapshots(request, tools) != before or self._tools() != tools:
            raise InputValidationError("TR plot protected input/tool identity changed during the version probe.")
        _protect(request, tools, owner=owner)
        svg, stdout, stderr = (request.output_dir / name for name in ("plot.svg", "command.stdout.log", "command.stderr.log"))
        owned_command = {path: _reserve(path) for path in (svg, stdout, stderr)}
        _protect(request, tools, owner=owner)
        for path, expected in owned_command.items():
            _owned(path, expected)
        result = self.runner.run(argv, timeout=request.timeout_seconds, stdout_path=stdout, stderr_path=stderr)
        if not result.executed or result.returncode != 0:
            raise OutputValidationError("TR plot command did not execute successfully.")
        _protect(request, tools, owner=owner)
        for path, expected in {**owned_version, **owned_command}.items():
            _owned(path, expected)
        _validate_svg(svg)
        after = self._snapshots(request, tools)
        if after != before or self._tools() != tools:
            raise InputValidationError("TR plot protected input/tool/plan identity changed during execution.")
        outputs = {name: _fingerprint(request.output_dir / name) for name in ("plot.svg", "effective_input.json")}
        # Empty tool logs are valid; record their raw byte identities separately.
        logs = {}
        for name in ("command.stderr.log", "command.stdout.log", "version.stderr.log", "version.stdout.log"):
            path = request.output_dir / name
            logs[name] = _fingerprint(path, allow_empty=True)
        receipt = {"schema_version": 1, "track": "tr.plot", "status": "COMPLETE_TECHNICAL_SVG",
                   "hifivar_version": __version__, "sample_id": request.sample_id, "selection": selection,
                   "tool_version": "5.1.0", "tool_version_reported": reported_version,
                   "command": argv, "version_command": list(version.args),
                   "returncode": result.returncode, "version_returncode": version.returncode,
                   "duration_seconds": result.duration_seconds, "inputs": before["inputs"], "tool_files": before["tools"],
                   "plan": before["plan"], "source_immutability": "PASS_EXPLICIT_SEVEN_FILES_PLAN_AND_TOOL_SHA256",
                   "source_origin": effective["source_origin"], "spanning_bam_scope": "PINNED_FILE_BYTES_NATIVE_PARSER_CHECKS_NOT_PRODUCER_LINEAGE_PROOF",
                   "manual_review": "NOT_REVIEWED", "accuracy": "NOT_ESTABLISHED", "allele_panel_order": effective["allele_panel_order"],
                   "visual_validation": "STRUCTURAL_SVG_ONLY_NOT_PIXEL_OR_MANUAL_QA", "plot_type": "allele", "show": "motifs",
                   "max_allele_reads": request.max_allele_reads, "outputs": outputs, "logs": logs,
                   "completion_path": str(request.output_dir / "completion.json"), "fresh_directory_only": True}
        pending = request.output_dir / ".completion.pending.json"
        pending_bytes = _json_bytes(receipt)
        pending_sha256 = hashlib.sha256(pending_bytes).hexdigest()
        pending_owner = _write_new(pending, receipt)
        _protect(request, tools, owner=owner)
        if self._snapshots(request, tools) != before or self._tools() != tools:
            raise InputValidationError("TR plot protected identities changed before final receipt publication.")
        for name, expected in {**outputs, **logs}.items():
            path = request.output_dir / name
            if path.is_symlink() or _fingerprint(path, allow_empty=True) != expected:
                raise OutputValidationError("TR plot output/log bytes changed before final receipt publication.")
        _owned(pending, pending_owner)
        sealed = _fingerprint(pending)
        if sealed["bytes"] != len(pending_bytes) or sealed["sha256"] != pending_sha256:
            raise OutputValidationError("TR plot pending completion bytes changed before final receipt publication.")
        # One no-overwrite commit point followed only by readback, never cleanup
        # or rollback. Pending evidence remains. A later interrupt/failure is
        # not a normal successful return even if a committed file is present.
        try:
            os.link(pending, request.output_dir / "completion.json")
        except OSError as error:
            raise OutputValidationError(f"TR plot completion destination already exists or cannot be published: {error}") from error
        for path in (pending, request.output_dir / "completion.json"):
            _owned(path, (path.resolve(), pending_owner[1]))
            published = _fingerprint(path)
            if published["bytes"] != len(pending_bytes) or published["sha256"] != pending_sha256:
                raise OutputValidationError("TR plot final completion readback differs from the sealed canonical receipt.")
        for name, expected in {**outputs, **logs}.items():
            path = request.output_dir / name
            if path.is_symlink() or _fingerprint(path, allow_empty=True) != expected:
                raise OutputValidationError("TR plot final receipt output/log linkage changed during publication.")
        if self._snapshots(request, tools) != before or self._tools() != tools:
            raise InputValidationError("TR plot protected identities changed during final completion readback.")
        return receipt


__all__ = ["TrgtPlotRequest", "TrgtPlotWrapper"]
