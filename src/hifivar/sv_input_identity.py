"""Explicit, reversible SV source-ID views; never normalize or deduplicate calls.

Canonical sources remain authoritative. SQLite bounds memory, JSONL records every
source ordinal, and a success manifest is published only after round-trip and
native BGZF/index validation. Ambiguous mate references are rejected, not guessed.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import gzip
import hashlib
import io
from importlib.resources import files
from itertools import zip_longest
import json
import os
from pathlib import Path
import re
import sqlite3
from typing import TYPE_CHECKING

from hifivar.command import CommandRunner
from hifivar.exceptions import ConfigurationError, InputValidationError, OutputValidationError
from hifivar.reference import ReferenceGenome
from hifivar.sample import validate_sample_id
from hifivar.logging_utils import get_logger

if TYPE_CHECKING:
    from hifivar.harmonization import SVEvidenceSourceArtifact

POLICY = json.loads(files("hifivar").joinpath("resources/policies/sv_source_identity_v1.json").read_text(encoding="utf-8"))
POLICY_SHA256 = hashlib.sha256(json.dumps(POLICY, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")).hexdigest()
SOURCE_IDENTITY_MIN_HTSLIB = (1, 23, 1)
_NAMESPACE = POLICY["analysis_namespace"]
_LOGGER = get_logger(__name__)


@dataclass(frozen=True, slots=True)
class SourceIdentityPolicy:
    policy_id: str
    policy_version: int
    policy_sha256: str
    mode: str
    timeout_seconds: int = 300

    def __post_init__(self) -> None:
        if (self.policy_id != POLICY["policy_id"] or type(self.policy_version) is not int
                or self.policy_version != POLICY["policy_version"] or self.policy_sha256 != POLICY_SHA256
                or self.mode != POLICY["mode"]):
            raise ConfigurationError("SV source identity requires the exact supported policy id/version/hash/mode.")
        if type(self.timeout_seconds) is not int or self.timeout_seconds <= 0:
            raise ConfigurationError("SV source identity timeout_seconds must be positive.")

    @classmethod
    def from_config(cls, section: object) -> SourceIdentityPolicy | None:
        if not isinstance(section, Mapping):
            raise ConfigurationError("sv.harmonization.input_identity must be a mapping.")
        expected = {"enabled", "policy_id", "policy_version", "policy_sha256", "mode", "timeout_seconds"}
        if set(section) - expected:
            raise ConfigurationError("Unknown SV source identity configuration key.")
        enabled, timeout = section.get("enabled", False), section.get("timeout_seconds", 300)
        if type(enabled) is not bool or type(timeout) is not int or timeout <= 0:
            raise ConfigurationError("SV source identity enabled must be boolean and timeout_seconds positive.")
        if not enabled:
            if any(section.get(k) is not None for k in expected - {"enabled", "timeout_seconds"}):
                raise ConfigurationError("Disabled SV source identity must not declare an active policy.")
            return None
        try:
            return cls(section["policy_id"], section["policy_version"], section["policy_sha256"], section["mode"], timeout)
        except KeyError as error:
            raise ConfigurationError("Enabled SV source identity requires explicit policy id/version/hash/mode.") from error

    def to_dict(self) -> dict[str, object]:
        return {"enabled": True, "policy_id": self.policy_id, "policy_version": self.policy_version,
                "policy_sha256": self.policy_sha256, "mode": self.mode, "timeout_seconds": self.timeout_seconds}


def file_identity(path: Path) -> dict[str, object]:
    """Hash a stable regular file, without following a last-component symlink."""
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise InputValidationError(f"SV source identity requires a regular non-symlink file: '{path}'.")
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as reader:
        for block in iter(lambda: reader.read(1024 * 1024), b""):
            digest.update(block)
    after = path.stat()
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns):
        raise InputValidationError("SV source identity input changed while hashing.")
    return {"path": str(path.resolve()), "sha256": digest.hexdigest(), "bytes": after.st_size,
            "device": after.st_dev, "inode": after.st_ino, "mtime_ns": after.st_mtime_ns}


def _json(path: Path, value: object) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as writer:
        json.dump(value, writer, ensure_ascii=True, sort_keys=True, indent=2)
        writer.write("\n")


def _reader(path: Path):
    return gzip.open(path, "rb") if str(path).endswith(".gz") else path.open("rb")


def _split(raw: bytes) -> tuple[list[bytes], bytes]:
    ending = b"\r\n" if raw.endswith(b"\r\n") else b"\n" if raw.endswith(b"\n") else b""
    parts = (raw[:-len(ending)] if ending else raw).split(b"\t")
    if len(parts) != 10 or not parts[2] or any(x in parts[2] for x in (b",", b";", b" ", b"\r", b"\n")):
        raise InputValidationError("SV ID view requires one sample and an unambiguous single ID token (or '.').")
    return parts, ending


def _records(path: Path):
    with _reader(path) as reader:
        for line in reader:
            if not line.startswith(b"#"):
                yield line


def _inventory(db: sqlite3.Connection, sources: Sequence[SVEvidenceSourceArtifact],
               sample_id: str, identities: Sequence[dict[str, object]],
               reference_lengths: Mapping[str, int]) -> None:
    db.executescript("""
        CREATE TABLE records (source INTEGER, ordinal INTEGER, original_id TEXT, analysis_id TEXT,
          native_sha TEXT, view_sha TEXT, emissions INTEGER DEFAULT 0, PRIMARY KEY(source,ordinal));
        CREATE INDEX original_ids ON records(original_id);
        CREATE UNIQUE INDEX analysis_ids ON records(analysis_id);
        CREATE TABLE links (source INTEGER, target TEXT);
    """)
    for i, source in enumerate(sources):
        chrom_header = False
        with _reader(Path(source.vcf_path)) as reader:
            ordinal = 0
            for raw in reader:
                if raw.startswith(b"#"):
                    if ordinal or raw.startswith(b"##hifivarSourceIdentity="):
                        raise InputValidationError("Source view already marked, or header appears after records.")
                    if raw.startswith(b"#CHROM\t"):
                        if chrom_header or raw.rstrip(b"\r\n").split(b"\t")[9:] != [sample_id.encode("utf-8")]:
                            raise InputValidationError("Source VCF sample header is not the declared single sample.")
                        chrom_header = True
                    if raw.startswith(b"##contig=<"):
                        declaration = raw.decode("utf-8").strip()
                        name = re.search(r"(?:<|,)ID=([^,>]+)", declaration)
                        length = re.search(r"(?:<|,)length=([0-9]+)(?:,|>)", declaration)
                        if (name is None or name[1] not in reference_lengths
                                or length is not None and int(length[1]) != reference_lengths[name[1]]):
                            raise InputValidationError("REFERENCE_CONTIG_MISMATCH: source contig name/length differs.")
                    continue
                if not chrom_header:
                    raise InputValidationError("Source VCF is missing its sample header.")
                ordinal += 1
                parts, _ = _split(raw)
                contig = parts[0].decode("utf-8")
                try:
                    position = int(parts[1])
                except ValueError as error:
                    raise InputValidationError("Source VCF POS is not an integer.") from error
                if contig not in reference_lengths or not 0 <= position <= reference_lengths[contig] + 1:
                    raise InputValidationError("REFERENCE_CONTIG_MISMATCH: source position exceeds exact 0..N+1 bounds.")
                identifier = parts[2].decode("utf-8")
                if identifier.startswith(_NAMESPACE):
                    raise InputValidationError("Source ID collides with the reserved analysis namespace.")
                db.execute("INSERT INTO records(source,ordinal,original_id,native_sha) VALUES(?,?,?,?)",
                           (i, ordinal, identifier, hashlib.sha256(raw).hexdigest()))
                for token in parts[7].split(b";"):
                    key, sep, value = token.partition(b"=")
                    if sep and key in {b"MATEID", b"PARID"}:
                        for target in value.decode("utf-8").split(","):
                            if target and target != ".":
                                db.execute("INSERT INTO links VALUES(?,?)", (i, target))
        if not chrom_header:
            raise InputValidationError("Source VCF lacks a #CHROM sample header.")
    db.execute("CREATE TABLE collisions AS SELECT original_id FROM records WHERE original_id!='.' GROUP BY original_id HAVING count(*)>1")
    db.execute("CREATE UNIQUE INDEX collision_ids ON collisions(original_id)")
    if db.execute("SELECT 1 FROM links JOIN collisions ON target=original_id LIMIT 1").fetchone():
        raise InputValidationError("Ambiguous MATEID/PARID references a duplicated original ID; graph mapping is required.")
    for i, source in enumerate(sources):
        source_sha = str(identities[i]["sha256"])
        db.execute("UPDATE records SET analysis_id=CASE WHEN original_id='.' OR original_id IN (SELECT original_id FROM collisions) "
                   "THEN ? || ordinal ELSE original_id END WHERE source=?", (f"{_NAMESPACE}{source.caller}_{source_sha}_", i))
    db.commit()


def _derive(db: sqlite3.Connection, source: Path, destination: Path, ledger: Path,
            source_index: int, caller: str, source_sha: str, policy: SourceIdentityPolicy) -> dict[str, int]:
    changed = count = 0
    marker = (f"##hifivarSourceIdentity=<ID={policy.policy_id},Version=1,PolicySHA256={policy.policy_sha256},Caller={caller},SourceSHA256={source_sha}>\n").encode("ascii")
    with _reader(source) as reader, destination.open("xb") as writer, ledger.open("x", encoding="utf-8", newline="\n") as book:
        for raw in reader:
            if raw.startswith(b"#"):
                if raw.startswith(b"#CHROM\t"):
                    writer.write(marker)
                writer.write(raw)
                continue
            count += 1
            parts, ending = _split(raw)
            entry = db.execute("SELECT original_id,analysis_id,native_sha FROM records WHERE source=? AND ordinal=?", (source_index, count)).fetchone()
            if entry is None or entry[0] != parts[2].decode("utf-8") or entry[2] != hashlib.sha256(raw).hexdigest():
                raise OutputValidationError("Source record changed between identity passes.")
            parts[2] = entry[1].encode("utf-8")
            represented = b"\t".join(parts) + ending
            digest = hashlib.sha256(represented).hexdigest()
            db.execute("UPDATE records SET view_sha=? WHERE source=? AND ordinal=?", (digest, source_index, count))
            row = {"source_index": source_index, "caller": caller, "source_vcf_sha256": source_sha,
                   "source_record_ordinal": count, "original_record_id": entry[0], "analysis_record_id": entry[1],
                   "source_record_key": f"caller:{caller}:sha256:{source_sha}:record:{count}",
                   "native_record_sha256": entry[2], "view_record_sha256": digest,
                   "ID_changed": entry[0] != entry[1]}
            book.write(json.dumps(row, sort_keys=True, ensure_ascii=True) + "\n")
            writer.write(represented)
            changed += row["ID_changed"]
    if count != db.execute("SELECT count(*) FROM records WHERE source=?", (source_index,)).fetchone()[0]:
        raise OutputValidationError("Source record count changed between identity passes.")
    db.commit()
    return {"records": count, "changed_ID_records": changed}


def verify_source_view(source: Path, view: Path, ledger: Path, *, caller: str | None = None) -> dict[str, object]:
    """Independently reconstruct every canonical row, including original EOL."""
    source_sha = file_identity(source)["sha256"]
    def headers(path: Path):
        with _reader(path) as reader:
            for raw in reader:
                if not raw.startswith(b"#"):
                    break
                yield raw
    native_headers, view_headers = list(headers(source)), list(headers(view))
    markers = [r for r in view_headers if r.startswith(b"##hifivarSourceIdentity=")]
    if len(markers) != 1 or [r for r in view_headers if not r.startswith(b"##hifivarSourceIdentity=")] != native_headers:
        raise OutputValidationError("Source identity view modified non-policy header bytes.")
    match = re.fullmatch(
        rb"##hifivarSourceIdentity=<ID=" + re.escape(POLICY["policy_id"].encode())
        + rb",Version=1,PolicySHA256=" + POLICY_SHA256.encode()
        + rb",Caller=([A-Za-z0-9_]+),SourceSHA256=" + str(source_sha).encode() + rb">\n", markers[0])
    if match is None or caller is not None and match[1].decode() != caller:
        raise OutputValidationError("Source identity policy/source/caller marker differs.")
    caller = match[1].decode()
    count = changed = 0
    with ledger.open(encoding="utf-8") as book:
        for count, (native, represented, line) in enumerate(zip_longest(_records(source), _records(view), book), 1):
            if native is None or represented is None or line is None:
                raise OutputValidationError("Source/view/ledger cardinality mismatch.")
            row = json.loads(line)
            left, eol = _split(native)
            right, right_eol = _split(represented)
            if (row["source_record_ordinal"] != count or row["source_vcf_sha256"] != source_sha
                    or row["caller"] != caller
                    or row["original_record_id"] != left[2].decode("utf-8")
                    or row["analysis_record_id"] != right[2].decode("utf-8")
                    or row["native_record_sha256"] != hashlib.sha256(native).hexdigest()
                    or row["view_record_sha256"] != hashlib.sha256(represented).hexdigest()
                    or type(row["ID_changed"]) is not bool
                    or row["ID_changed"] != (left[2] != right[2]) or eol != right_eol
                    or row["source_record_key"] != f"caller:{caller}:sha256:{source_sha}:record:{count}"):
                raise OutputValidationError("Source identity ledger disagrees with actual record bytes/order.")
            right[2] = left[2]
            if b"\t".join(right) + right_eol != native:
                raise OutputValidationError("A non-ID record field was modified.")
            changed += row["ID_changed"]
    return {"status": "PASS_EXACT_SOURCE_RECORD_ROUND_TRIP", "records": count, "changed_ID_records": changed}


@dataclass(frozen=True, slots=True)
class SourceIdentityResult:
    directory: Path
    manifest_path: Path
    database_path: Path
    views: tuple[Path, ...]
    ledgers: tuple[Path, ...]
    source_order: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return {"manifest": str(self.manifest_path), "directory": str(self.directory),
                "database": str(self.database_path), "views": [str(x) for x in self.views],
                "ledgers": [str(x) for x in self.ledgers], "source_order": list(self.source_order),
                "canonical_sources_unchanged": True, "standalone_view_is_lossless": False}

    def validate(self) -> dict[str, object]:
        payload = json.loads(self.manifest_path.read_text(encoding="utf-8"))
        if payload.get("status") != "completed" or payload.get("source_order") != list(self.source_order):
            raise OutputValidationError("SV source identity success manifest or order differs.")
        if (payload.get("views") != [str(p.relative_to(self.directory)) for p in self.views]
                or payload.get("ledgers") != [str(p.relative_to(self.directory)) for p in self.ledgers]):
            raise OutputValidationError("SV source identity representation order differs.")
        if ([s.get("caller") for s in payload["canonical_sources"]] != list(self.source_order)
                or len(self.views) != len(self.source_order) or len(self.ledgers) != len(self.source_order)):
            raise OutputValidationError("SV source identity canonical source cardinality/order differs.")
        if SourceIdentityPolicy.from_config(payload["policy"]) is None:
            raise OutputValidationError("Source identity manifest must bind an active policy.")
        expected = {self.database_path, *self.views, *(Path(str(p) + ".tbi") for p in self.views), *self.ledgers}
        if set(payload["artifact_identities"]) != {str(p.relative_to(self.directory)) for p in expected}:
            raise OutputValidationError("Source identity manifest does not bind all required artifacts.")
        for relative, identity in payload["artifact_identities"].items():
            path = (self.directory / relative).resolve()
            if not path.is_relative_to(self.directory.resolve()) or file_identity(path) != identity:
                raise OutputValidationError("SV source identity artifact changed or escaped its bundle.")
        for source in payload["canonical_sources"]:
            for role in ("vcf", "index"):
                if file_identity(Path(source[role]["path"])) != source[role]:
                    raise OutputValidationError("Canonical SV source/index changed after identity preparation.")
        for identity in payload["reference_inputs"]:
            if file_identity(Path(identity["path"])) != identity:
                raise OutputValidationError("Source identity reference changed.")
        return payload


class SourceIdentityWrapper:
    """Prepare explicit BGZF views using the unified command runner."""
    def __init__(self, *, bgzip_executable: str = "bgzip", tabix_executable: str = "tabix",
                 runner: CommandRunner | None = None) -> None:
        self.bgzip, self.tabix = bgzip_executable, tabix_executable
        self.runner = runner or CommandRunner()

    def run(self, sources: Sequence[SVEvidenceSourceArtifact], *, sample_id: str,
            reference: ReferenceGenome, output_directory: Path, policy: SourceIdentityPolicy,
            dry_run: bool = False) -> SourceIdentityResult | dict[str, object]:
        validate_sample_id(sample_id)
        if not isinstance(policy, SourceIdentityPolicy) or not isinstance(reference, ReferenceGenome):
            raise InputValidationError("Source identity requires explicit policy and reference objects.")
        directory = Path(output_directory)
        if not directory.is_absolute() or directory.exists() or directory.is_symlink():
            raise OutputValidationError("SV source identity requires a fresh absolute directory; overwrite is unsupported.")
        if not sources or len({s.caller for s in sources}) != len(sources):
            raise InputValidationError("Source identity requires unique ordered callers.")
        for s in sources:
            if (s.sample_id != sample_id or re.fullmatch(r"[A-Za-z0-9_]+", s.caller) is None
                    or s.vcf_path is None or s.index_path is None
                    or s.status.value not in {"completed", "no_calls"}):
                raise InputValidationError("Invalid source identity sample/caller/VCF/index.")
            if Path(s.index_path).resolve() != Path(str(s.vcf_path) + ".tbi").resolve():
                raise InputValidationError("SV source identity requires its declared adjacent .tbi index.")
            for p in (s.vcf_path, s.index_path, reference.fasta, reference.fai):
                if Path(p).resolve().is_relative_to(directory.resolve()):
                    raise InputValidationError("Source identity output must not contain canonical inputs.")
        if dry_run:
            return {"status": "planned", "policy": policy.to_dict(), "source_order": [s.caller for s in sources],
                    "directory": str(directory), "executed": False, "will_create_outputs": False}
        from hifivar.sv import validate_sv_vcf
        before = []
        reference_inputs = [file_identity(reference.fasta), file_identity(reference.fai)]
        seen = set()
        for s in sources:
            validate_sv_vcf(s.vcf_path, s.index_path, sample_id=sample_id, reference=reference)
            vcf, index = file_identity(s.vcf_path), file_identity(s.index_path)
            key = (vcf["device"], vcf["inode"])
            if key in seen:
                raise InputValidationError("Different callers must not alias the same physical VCF.")
            seen.add(key)
            before.append({"caller": s.caller, "vcf": vcf, "index": index})
        directory.mkdir(parents=True, exist_ok=False, mode=0o700)
        work, logs = directory / "work", directory / "logs"
        work.mkdir(mode=0o700)
        logs.mkdir(mode=0o700)
        db_path, db = directory / "records.sqlite3", None
        commands: list[dict[str, object]] = []
        try:
            env = {k: os.environ[k] for k in ("PATH", "SYSTEMROOT", "WINDIR") if k in os.environ}
            for key, name in (("HOME", "home"), ("TMPDIR", "tmp"), ("XDG_CONFIG_HOME", "config"), ("XDG_CACHE_HOME", "cache")):
                p = work / name
                p.mkdir(mode=0o700)
                env[key] = str(p)
            env.update(TMP=env["TMPDIR"], TEMP=env["TMPDIR"], LC_ALL="C", TZ="UTC")
            def command(args: list[str], step: str, stdout_path: Path | None = None):
                record = {"step": step, "args": args, "status": "started", "stderr_path": str(logs / f"{step}.stderr.log")}
                commands.append(record)
                try:
                    result = self.runner.run(args, cwd=work, env=env, inherit_env=False,
                                             timeout=policy.timeout_seconds, stdout_path=stdout_path,
                                             stderr_path=logs / f"{step}.stderr.log")
                except Exception as error:
                    record.update(status="failed", error_type=type(error).__name__)
                    raise
                commands[-1] = {"step": step, "status": "completed", **result.to_dict()}
                if not result.executed or result.returncode != 0:
                    raise OutputValidationError("Source identity command did not execute successfully.")
                return result
            tools, versions = {}, {}
            for name, exe in (("bgzip", self.bgzip), ("tabix", self.tabix)):
                resolved = self.runner.require_executable(exe, env=env, inherit_env=False, cwd=work)
                tools[name] = file_identity(resolved)
                result = command([str(resolved), "--version"], f"{name}_version")
                text = Path(result.stdout_path).read_text(encoding="utf-8") if result.stdout_path else result.stdout or ""
                found = re.search(r"(?<!\d)(\d+)\.(\d+)(?:\.(\d+))?", text)
                if found is None or tuple(int(x or 0) for x in found.groups()) < SOURCE_IDENTITY_MIN_HTSLIB:
                    raise OutputValidationError("SV source identity requires bgzip/tabix >=1.23.1 for mixed VCF4.4 inputs.")
                versions[name] = found.group(0)
            db = sqlite3.connect(db_path)
            for i, s in enumerate(sources):
                command([str(tools["tabix"]["path"]), "-l", str(Path(s.vcf_path).absolute())], f"canonical_tabix_list_{i}")
            _inventory(db, sources, sample_id, [s["vcf"] for s in before],
                       {c.name: c.length for c in reference.contigs})
            views, ledgers, summaries = [], [], []
            for i, s in enumerate(sources):
                stem = f"{i:02d}.{s.caller}.analysis"
                plain, ledger, view = work / f"{stem}.vcf", directory / f"{stem}.records.jsonl", directory / f"{stem}.vcf.gz"
                summary = _derive(db, Path(s.vcf_path), plain, ledger, i, s.caller, str(before[i]["vcf"]["sha256"]), policy)
                command([str(tools["bgzip"]["path"]), "-c", str(plain)], f"bgzip_{i}", view)
                command([str(tools["tabix"]["path"]), "-p", "vcf", str(view)], f"tabix_{i}")
                validate_sv_vcf(view, Path(str(view) + ".tbi"), sample_id=sample_id, reference=reference)
                verify_source_view(Path(s.vcf_path), view, ledger, caller=s.caller)
                command([str(tools["tabix"]["path"]), "-l", str(view)], f"tabix_list_{i}")
                first = next(_records(view), None)
                if first is not None:
                    parts, _ = _split(first)
                    contig, pos = parts[0].decode("utf-8"), int(parts[1])
                    lengths = {c.name: c.length for c in reference.contigs}
                    if 0 < pos <= lengths[contig] and ":" not in contig:
                        fetched = command([str(tools["tabix"]["path"]), str(view), f"{contig}:{pos}-{pos}"], f"tabix_query_{i}")
                        output = Path(fetched.stdout_path).open("rb") if fetched.stdout_path else io.BytesIO((fetched.stdout or "").encode("utf-8"))
                        with output as reader:
                            if not any(first.rstrip(b"\r\n") == r.rstrip(b"\r\n") for r in reader):
                                raise OutputValidationError("Source view index query did not recover the selected record.")
                summaries.append(summary)
                views.append(view)
                ledgers.append(ledger)
            db.close()
            db = None
            for source in before:
                for role in ("vcf", "index"):
                    if file_identity(Path(source[role]["path"])) != source[role]:
                        raise OutputValidationError("Canonical source/index changed during identity preparation.")
            for identity in tools.values():
                if file_identity(Path(identity["path"])) != identity:
                    raise OutputValidationError("Source identity tool changed during execution.")
            for identity in reference_inputs:
                if file_identity(Path(identity["path"])) != identity:
                    raise OutputValidationError("Source identity reference changed during execution.")
            artifacts = (db_path, *views, *(Path(str(p) + ".tbi") for p in views), *ledgers)
            payload = {"schema_version": 1, "status": "completed", "policy": policy.to_dict(),
                       "source_order": [s.caller for s in sources], "canonical_sources": before,
                       "reference_inputs": reference_inputs,
                       "views": [str(p.relative_to(directory)) for p in views],
                       "ledgers": [str(p.relative_to(directory)) for p in ledgers],
                       "source_summaries": summaries, "tools": tools, "tool_versions": versions, "commands": commands,
                       "canonical_sources_unchanged": True, "standalone_view_is_lossless": False,
                       "artifact_identities": {str(p.relative_to(directory)): file_identity(p) for p in artifacts}}
            temporary = directory / "manifest.partial.json"
            _json(temporary, payload)
            temporary.rename(directory / "manifest.json")
            _LOGGER.info("SV source identity completed for sample %s (%d ordered sources); canonical calls unchanged.", sample_id, len(sources))
            return SourceIdentityResult(directory, directory / "manifest.json", db_path, tuple(views), tuple(ledgers), tuple(s.caller for s in sources))
        except Exception as error:
            if db is not None:
                db.close()
            _json(directory / "FAILED.json", {"status": "failed", "error_type": type(error).__name__,
                  "success_manifest_published": False, "commands": commands, "canonical_sources": before})
            _LOGGER.error("SV source identity failed for sample %s (%s); owned evidence preserved.", sample_id, type(error).__name__)
            raise


__all__ = ["SourceIdentityPolicy", "SourceIdentityResult", "SourceIdentityWrapper", "POLICY", "POLICY_SHA256", "verify_source_view"]
