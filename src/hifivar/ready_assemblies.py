"""Explicit, validation-only reuse of existing haplotype FASTAs.

Declared sample/role/source metadata is not biological phasing truth and does
not establish that hifiasm ran in the current attempt. No command is executed,
no input is converted, and no native producer receipt is issued here.
"""
from __future__ import annotations

import csv
import hashlib
import io
import os
import re
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Sequence

from hifivar.assembly import AssemblyRole
from hifivar.exceptions import InputValidationError
from hifivar.sample import validate_sample_id

_DNA = b"ACGTRYSWKMBDHVNacgtryswkmbdhvn"
_ROLES = (AssemblyRole.HAPLOTYPE1, AssemblyRole.HAPLOTYPE2)
_SHA = re.compile(r"[0-9a-fA-F]{64}\Z")
_REQUIRED = {"sample_id", "role", "fasta", "sha256"}
_OPTIONAL = {"source_tool", "source_version"}
_MANIFEST_LIMIT = 8 * 1024 * 1024


def _identity(stat: os.stat_result) -> tuple[int, int, int, int, int]:
    # On Windows Python 3.12 stat()/fstat() expose different legacy ctime
    # meanings (creation versus metadata-change time). Content hashing plus
    # device/inode/size/mtime remain the cross-platform identity checks.
    ctime = stat.st_ctime_ns if os.name != "nt" else 0
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, ctime


def _path(path: Path) -> Path:
    try:
        result = Path(path).expanduser().resolve(strict=True)
        if not result.is_file():
            raise InputValidationError(f"Ready assembly input is not a regular file: '{path}'.")
        return result
    except (OSError, RuntimeError, ValueError, TypeError) as error:
        raise InputValidationError(f"Cannot read ready assembly input: '{path}'.") from error


class _FastaScan:
    """Bounded fragments, including chromosome-length unwrapped sequences."""

    def __init__(self) -> None:
        self.line_start = True
        self.header = False
        self.header_bytes = bytearray()
        self.pending_cr = False
        self.ids: set[bytes] = set()
        self.current_bases = 0
        self.bases = 0

    def feed(self, data: bytes) -> None:
        if not data:
            return
        if self.pending_cr:
            raise InputValidationError("FASTA contains CR outside a CRLF line ending.")
        if data.endswith(b"\r"):
            self.pending_cr = True
            data = data[:-1]
        if b"\r" in data:
            raise InputValidationError("FASTA contains an internal CR character.")
        if not data:
            return
        if self.line_start:
            self.header = data.startswith(b">")
            self.line_start = False
        if self.header:
            self.header_bytes.extend(data)
            if len(self.header_bytes) > 65536:
                raise InputValidationError("FASTA header exceeds the 64 KiB limit.")
        else:
            if not self.ids or data.translate(None, _DNA):
                raise InputValidationError("FASTA sequence precedes a header or contains non-IUPAC DNA.")
            self.current_bases += len(data)
            self.bases += len(data)

    def finish_line(self) -> None:
        if self.line_start:
            raise InputValidationError("Ready assembly FASTA contains an empty line.")
        if self.header:
            header = bytes(self.header_bytes)
            if any(c < 32 and c != 9 or c > 126 for c in header):
                raise InputValidationError("FASTA headers must contain printable ASCII identifiers.")
            fields = header[1:].split()
            if not fields or header[1:2].isspace() or fields[0] in self.ids:
                raise InputValidationError("FASTA contains an empty or duplicate sequence identifier.")
            if self.ids and self.current_bases == 0:
                raise InputValidationError("FASTA record has no sequence.")
            self.ids.add(fields[0])
            self.current_bases = 0
            self.header_bytes.clear()
        self.line_start = True
        self.pending_cr = False

    def finish(self) -> None:
        if self.pending_cr:
            raise InputValidationError("FASTA ends with CR rather than a CRLF line ending.")
        if not self.line_start:
            self.finish_line()
        if not self.ids or not self.current_bases:
            raise InputValidationError("Ready assembly FASTA contains no sequence or an empty final record.")


def _scan_fasta(path: Path, *, chunk_size: int = 1024 * 1024) -> tuple[str, int, int, int]:
    """Validate every byte and fingerprint the same file descriptor."""
    if type(chunk_size) is not int or not 1 <= chunk_size <= 8 * 1024 * 1024:
        raise InputValidationError("FASTA validation chunks must be bounded positive integers.")
    scanner = _FastaScan()
    digest = hashlib.sha256()
    try:
        before = path.stat()
        with path.open("rb") as stream:
            opened = os.fstat(stream.fileno())
            if _identity(before) != _identity(opened):
                raise InputValidationError("Ready assembly changed while opening.")
            while block := stream.read(chunk_size):
                digest.update(block)
                pieces = block.split(b"\n")
                for part in pieces[:-1]:
                    scanner.feed(part)
                    scanner.finish_line()
                scanner.feed(pieces[-1])
            scanner.finish()
            if _identity(opened) != _identity(os.fstat(stream.fileno())):
                raise InputValidationError("Ready assembly changed during validation.")
        if _identity(opened) != _identity(path.stat()):
            raise InputValidationError("Ready assembly path changed during validation.")
    except OSError as error:
        raise InputValidationError(f"Cannot validate ready assembly: '{path}'.") from error
    return digest.hexdigest(), opened.st_size, len(scanner.ids), scanner.bases


@dataclass(frozen=True, slots=True)
class ReadyAssemblyArtifact:
    """Content-validated external FASTA with explicitly user-declared roles."""

    sample_id: str
    role: AssemblyRole
    path: Path
    sha256: str
    file_size: int
    sequence_count: int
    base_count: int
    source_tool: str | None = None
    source_version: str | None = None
    manifest_path: Path | None = None
    manifest_sha256: str | None = None

    def __post_init__(self) -> None:
        validate_sample_id(self.sample_id)
        if not isinstance(self.role, AssemblyRole) or self.role not in _ROLES:
            raise InputValidationError("Ready assembly requires an explicit haplotype1/haplotype2 role.")
        if not isinstance(self.sha256, str) or _SHA.fullmatch(self.sha256) is None:
            raise InputValidationError("Ready assembly requires a full expected SHA256.")
        for value in (self.file_size, self.sequence_count, self.base_count):
            if type(value) is not int or value <= 0:
                raise InputValidationError("Ready assembly counts must be positive integers.")
        for value in (self.source_tool, self.source_version):
            if value is not None and (not isinstance(value, str) or not value.strip()
                                      or any(ord(c) < 32 for c in value)):
                raise InputValidationError("Declared source metadata must be nonempty printable text.")
        if (self.manifest_path is None) != (self.manifest_sha256 is None):
            raise InputValidationError("Ready assembly manifest path and SHA256 must be supplied together.")
        if self.manifest_sha256 is not None and (not isinstance(self.manifest_sha256, str)
                                                or _SHA.fullmatch(self.manifest_sha256) is None):
            raise InputValidationError("Ready assembly manifest SHA256 is invalid.")
        object.__setattr__(self, "path", Path(self.path))
        object.__setattr__(self, "sha256", self.sha256.lower())
        if self.manifest_path is not None:
            object.__setattr__(self, "manifest_path", Path(self.manifest_path))
            object.__setattr__(self, "manifest_sha256", self.manifest_sha256.lower())

    def verify(self) -> None:
        """Revalidate content and manifest before reuse or accepting new outputs."""
        actual = _scan_fasta(_path(self.path))
        if actual != (self.sha256, self.file_size, self.sequence_count, self.base_count):
            raise InputValidationError("Ready assembly content differs from the validated declaration.")
        if self.manifest_path is not None:
            data, digest = _read_manifest_bytes(self.manifest_path)
            if digest != self.manifest_sha256:
                raise InputValidationError("Ready assembly manifest changed after validation.")
            if not _matches_manifest_row(data, self):
                raise InputValidationError("Ready assembly sample/role/source differs from its manifest row.")

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id, "role": self.role.value, "path": str(self.path),
            "sha256": self.sha256, "file_size": self.file_size,
            "sequence_count": self.sequence_count, "base_count": self.base_count,
            "input_mode": "EXTERNAL_READY", "execution_kind": "VALIDATION_ONLY",
            "role_source": "USER_DECLARED_NOT_PHASE_TRUTH",
            "source_tool": self.source_tool, "source_version": self.source_version,
            "source_version_source": "USER_DECLARED" if self.source_version else "UNKNOWN",
            "source_gfa": None, "hifiasm_version": None, "reference_independent": None,
            "manifest_path": str(self.manifest_path) if self.manifest_path else None,
            "manifest_sha256": self.manifest_sha256,
        }


def validate_ready_assembly(sample_id: str, role: AssemblyRole, fasta: Path,
                            expected_sha256: str, *, source_tool: str | None = None,
                            source_version: str | None = None) -> ReadyAssemblyArtifact:
    """Read-only full validation; only plain FASTA is currently supported."""
    validate_sample_id(sample_id)
    if not isinstance(role, AssemblyRole) or role not in _ROLES:
        raise InputValidationError("Ready assembly requires an explicit haplotype1/haplotype2 role.")
    if not isinstance(expected_sha256, str) or _SHA.fullmatch(expected_sha256) is None:
        raise InputValidationError("Ready assembly expected SHA256 must contain 64 hexadecimal characters.")
    path = _path(fasta)
    if path.suffix.lower() not in {".fa", ".fasta", ".fna"}:
        raise InputValidationError("Ready assembly requires uncompressed .fa/.fasta/.fna; no implicit conversion.")
    digest, size, records, bases = _scan_fasta(path)
    if digest != expected_sha256.lower():
        raise InputValidationError("Ready assembly SHA256 does not match its explicit declaration.")
    return ReadyAssemblyArtifact(sample_id, role, path, digest, size, records, bases,
                                 source_tool, source_version)


def _read_manifest_bytes(path: Path) -> tuple[bytes, str]:
    try:
        with _path(path).open("rb") as stream:
            before = os.fstat(stream.fileno())
            data = stream.read(_MANIFEST_LIMIT + 1)
            if len(data) > _MANIFEST_LIMIT:
                raise InputValidationError("Ready assembly manifest exceeds the 8 MiB limit.")
            if _identity(before) != _identity(os.fstat(stream.fileno())) or _identity(before) != _identity(Path(path).stat()):
                raise InputValidationError("Ready assembly manifest changed during validation.")
    except OSError as error:
        raise InputValidationError("Cannot read ready assembly manifest.") from error
    return data, hashlib.sha256(data).hexdigest()


def _matches_manifest_row(data: bytes, item: ReadyAssemblyArtifact) -> bool:
    """Bind roles/source fields as well as bytes; do not infer from filenames."""
    try:
        reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig"), newline=""), delimiter="\t", strict=True)
        matches = 0
        for row in reader:
            if row.get("sample_id") != item.sample_id or row.get("role") != item.role.value:
                continue
            raw = Path(row["fasta"])
            path = raw if raw.is_absolute() else Path(item.manifest_path).parent / raw
            if (_path(path) == _path(item.path) and row.get("sha256", "").lower() == item.sha256
                    and (row.get("source_tool") or None) == item.source_tool
                    and (row.get("source_version") or None) == item.source_version):
                matches += 1
        return matches == 1
    except (UnicodeError, csv.Error, KeyError, TypeError):
        return False


@dataclass(frozen=True, slots=True)
class ReadyAssemblyManifest:
    path: Path
    sha256: str
    declared_sample_order: tuple[str, ...]
    assemblies: tuple[ReadyAssemblyArtifact, ...]

    def for_sample(self, sample_id: str) -> tuple[ReadyAssemblyArtifact, ...]:
        """Return canonical h1/h2 argument order, without rewriting inputs."""
        items = tuple(item for item in self.assemblies if item.sample_id == sample_id)
        if not items:
            raise InputValidationError("Sample is not present in the ready assembly manifest.")
        return tuple(sorted(items, key=lambda item: _ROLES.index(item.role)))

    def to_dict(self) -> dict[str, object]:
        return {"schema_version": 1, "path": str(self.path), "sha256": self.sha256,
                "execution_kind": "VALIDATION_ONLY", "declared_sample_order": list(self.declared_sample_order),
                "declared_rows": [item.to_dict() for item in self.assemblies],
                "caller_role_order": [role.value for role in _ROLES], "new_hifiasm_execution": False}


def load_ready_assembly_manifest(path: Path, *, expected_samples: Sequence[str]) -> ReadyAssemblyManifest:
    """Strict TSV; relative paths are anchored at the manifest, never the CWD."""
    if not isinstance(expected_samples, Sequence) or isinstance(expected_samples, (str, bytes)) or not expected_samples:
        raise InputValidationError("Expected ready assembly samples must be a nonempty sequence.")
    samples = tuple(expected_samples)
    for value in samples:
        validate_sample_id(value)
    if len(set(samples)) != len(samples):
        raise InputValidationError("Expected ready assembly samples must be unique.")
    manifest_path = _path(path)
    data, digest = _read_manifest_bytes(manifest_path)
    try:
        reader = csv.reader(io.StringIO(data.decode("utf-8-sig"), newline=""), delimiter="\t", strict=True)
        columns = next(reader)
        if len(set(columns)) != len(columns) or not _REQUIRED <= set(columns) or set(columns) - _REQUIRED - _OPTIONAL:
            raise InputValidationError("Ready assembly manifest has missing, duplicate or unknown columns.")
        entries: list[ReadyAssemblyArtifact] = []
        seen_roles: set[tuple[str, AssemblyRole]] = set()
        seen_files: set[tuple[int, int]] = set()
        order: list[str] = []
        for fields in reader:
            if len(fields) != len(columns) or len(entries) >= 10000:
                raise InputValidationError("Ready assembly manifest row is malformed or exceeds the row limit.")
            row = dict(zip(columns, fields))
            sample = row["sample_id"]
            validate_sample_id(sample)
            try:
                role = AssemblyRole(row["role"])
            except ValueError as error:
                raise InputValidationError("Unknown ready assembly role.") from error
            if (sample, role) in seen_roles:
                raise InputValidationError("Ready assembly manifest contains duplicate sample roles.")
            if sample not in samples:
                raise InputValidationError("Ready assembly manifest has an unexpected sample.")
            raw_path = Path(row["fasta"])
            if not row["fasta"] or row["fasta"] != row["fasta"].strip():
                raise InputValidationError("Ready assembly FASTA path is empty or has surrounding whitespace.")
            resolved = _path(raw_path if raw_path.is_absolute() else manifest_path.parent / raw_path)
            stat = resolved.stat()
            physical = (stat.st_dev, stat.st_ino)
            if resolved == manifest_path or physical in seen_files:
                raise InputValidationError("Ready assembly roles/samples must not alias one physical file.")
            item = validate_ready_assembly(sample, role, resolved, row["sha256"],
                                           source_tool=row.get("source_tool") or None,
                                           source_version=row.get("source_version") or None)
            entries.append(replace(item, manifest_path=manifest_path, manifest_sha256=digest))
            seen_roles.add((sample, role))
            seen_files.add(physical)
            if sample not in order:
                order.append(sample)
        if set(order) != set(samples):
            raise InputValidationError("Ready assembly manifest is missing expected samples.")
        if _read_manifest_bytes(manifest_path)[1] != digest:
            raise InputValidationError("Ready assembly manifest changed while its FASTAs were validated.")
        return ReadyAssemblyManifest(manifest_path, digest, tuple(order), tuple(entries))
    except (UnicodeError, csv.Error, StopIteration) as error:
        raise InputValidationError("Ready assembly manifest must be a nonempty UTF-8 TSV.") from error


def protect_ready_assembly_inputs(assemblies: Sequence[object], *, outputs: Sequence[Path] = (),
                                  work_directory: Path | None = None, verify: bool = True,
                                  require_fresh_work_directory: bool = False) -> None:
    """Reject write-path aliases before any native/version command is issued."""
    inputs: list[Path] = []
    physical: set[tuple[int, int]] = set()
    for item in assemblies:
        if not isinstance(item, ReadyAssemblyArtifact):
            continue
        if verify:
            item.verify()
        path = _path(item.path)
        stat = path.stat()
        identity = (stat.st_dev, stat.st_ino)
        if identity in physical:
            raise InputValidationError("Ready assembly roles alias one physical file.")
        physical.add(identity)
        inputs.append(path)
        if item.manifest_path is not None:
            inputs.append(_path(item.manifest_path))
    if inputs and require_fresh_work_directory and work_directory is not None:
        directory = Path(work_directory)
        if directory.exists() and (not directory.is_dir() or any(directory.iterdir())):
            raise InputValidationError("Ready assembly reuse requires a new or empty caller work directory.")
    for source in inputs:
        if work_directory is not None and source.is_relative_to(Path(work_directory).resolve()):
            raise InputValidationError("Ready assembly input cannot reside within the caller work directory.")
        for output in outputs:
            target = Path(output).resolve()
            if target == source or target.exists() and os.path.samefile(source, target):
                raise InputValidationError("Caller output/log path aliases a ready assembly input.")


__all__ = ["ReadyAssemblyArtifact", "ReadyAssemblyManifest", "load_ready_assembly_manifest",
           "validate_ready_assembly", "protect_ready_assembly_inputs"]
