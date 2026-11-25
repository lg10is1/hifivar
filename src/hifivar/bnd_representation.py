"""Opt-in, no-overwrite native-plus-derived BND representation bundles.

The native Jasmine artifact retains its original output role. This wrapper
does not select derived calls for any benchmark, annotation or review consumer.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import time

from hifivar.bnd_records import derive_record_stream, verify_record_round_trip
from hifivar.breakends import BndPolicy
from hifivar.command import CommandRunner
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar.jasmine import JasmineResult, JasmineResultStatus, VCF44_MIN_TABIX
from hifivar.logging_utils import get_logger


_LOGGER = get_logger(__name__)


def _identity(path: Path) -> dict[str, object]:
    if not path.is_file() or path.is_symlink():
        raise InputValidationError(f"BND input must be an existing regular, non-symlink file: '{path}'.")
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as reader:
        for chunk in iter(lambda: reader.read(1024 * 1024), b""):
            digest.update(chunk)
    after = path.stat()
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns):
        raise InputValidationError("BND input changed while its identity was being recorded.")
    return {"path": str(path.resolve()), "sha256": digest.hexdigest(), "bytes": after.st_size,
            "device": after.st_dev, "inode": after.st_ino, "mtime_ns": after.st_mtime_ns}


def _write_json(path: Path, payload: object) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as writer:
        json.dump(payload, writer, indent=2, sort_keys=True)
        writer.write("\n")


def _copy_exact(source: Path, target: Path, identity: dict[str, object]) -> None:
    with source.open("rb") as reader, target.open("xb") as writer:
        shutil.copyfileobj(reader, writer, length=1024 * 1024)
    copied = _identity(target)
    if (copied["sha256"], copied["bytes"]) != (identity["sha256"], identity["bytes"]):
        raise OutputValidationError("Native binary copy differs from its bound input identity.")


@dataclass(frozen=True, slots=True)
class BndRepresentationRequest:
    native_result: JasmineResult
    output_directory: Path
    policy: BndPolicy

    def __post_init__(self) -> None:
        result = self.native_result
        if not isinstance(result, JasmineResult) or result.status != JasmineResultStatus.COMPLETED or result.artifact is None:
            raise InputValidationError("BND derivation requires a completed native JasmineResult, not a planned output path.")
        if not isinstance(self.policy, BndPolicy):
            raise InputValidationError("BND derivation requires an explicit frozen BndPolicy.")
        artifact, request = result.artifact, result.request
        if (artifact.jasmine_version != "1.1.5" or not result.tool_versions or result.tool_versions.get("jasmine") != "1.1.5"
                or artifact.sample_id != request.sample_id or artifact.vcf_path.resolve() != request.output_vcf.resolve()
                or artifact.index_path.resolve() != request.output_index.resolve()
                or artifact.evidence_table.resolve() != request.evidence_table.resolve()
                or artifact.reference_fasta.resolve() != request.reference.fasta.resolve()
                or artifact.reference_build != request.reference.build or artifact.sources != request.sources
                or not result.commands or artifact.commands != tuple(item.args for item in result.commands)):
            raise InputValidationError("Native Jasmine result, commands, source order or artifact identity disagree.")
        directory = Path(self.output_directory)
        if not directory.is_absolute():
            raise InputValidationError("BND output directory must be an explicit absolute path.")
        object.__setattr__(self, "output_directory", directory)
        if directory.exists() or directory.is_symlink():
            raise OutputValidationError("BND representation requires a fresh directory; overwrite is unsupported.")

    @property
    def native_raw(self) -> Path:
        matches = [Path(path) for path in self.native_result.artifact.intermediate_files
                   if Path(path).name == f"{self.native_result.request.sample_id}.jasmine.raw.vcf"]
        if len(matches) != 1:
            raise InputValidationError("BND bundle requires the unique native Jasmine raw VCF provenance artifact.")
        return matches[0]

    def inputs(self) -> tuple[Path, ...]:
        result = self.native_result
        request, artifact = result.request, result.artifact
        paths = [Path(artifact.vcf_path), Path(artifact.index_path), self.native_raw,
                 Path(artifact.evidence_table), request.reference.fasta, request.reference.fai]
        for source in request.runnable_sources:
            if source.vcf_path is None or source.index_path is None:
                raise InputValidationError("Completed source lacks its VCF or index identity.")
            paths.extend([Path(source.vcf_path), Path(source.index_path)])
        if result.source_identity is not None:
            result.source_identity.validate()
            identity = result.source_identity
            paths.extend([identity.manifest_path, identity.database_path, *identity.views,
                          *(Path(str(p) + ".tbi") for p in identity.views), *identity.ledgers])
        output = self.output_directory.resolve()
        for path in paths:
            if path.resolve().is_relative_to(output):
                raise InputValidationError("BND bundle must not contain or replace an upstream input.")
        return tuple(dict.fromkeys(paths))


@dataclass(frozen=True, slots=True)
class BndRepresentationResult:
    output_directory: Path
    manifest_path: Path
    derived_vcf: Path
    derived_index: Path
    ledger_path: Path
    junction_path: Path
    summary: dict[str, object]

    def to_dict(self) -> dict[str, object]:
        return {"status": "completed", "native_output_role_unchanged": True,
                "consumer_opt_in_required": True, "standalone_derived_VCF_is_lossless": False,
                "output_directory": str(self.output_directory), "manifest": str(self.manifest_path),
                "derived_vcf": str(self.derived_vcf), "derived_index": str(self.derived_index),
                "record_ledger": str(self.ledger_path), "junction_table": str(self.junction_path),
                "summary": self.summary}


class BndRepresentationWrapper:
    """Stream a narrowly validated derived copy, then verify compression/index.

    Failure preserves the owned directory and records FAILED.json, with no
    manifest.json success sentinel. Existing/native files are never replaced.
    """
    def __init__(self, *, bgzip_executable: str = "bgzip", tabix_executable: str = "tabix",
                 runner: CommandRunner | None = None) -> None:
        self.bgzip, self.tabix = bgzip_executable, tabix_executable
        self.runner = runner or CommandRunner()

    def plan(self, request: BndRepresentationRequest) -> dict[str, object]:
        directory, sample = request.output_directory, request.native_result.request.sample_id
        return {"status": "planned", "policy": request.policy.to_dict(), "native_result_status": "completed",
                "native_vcf": str(request.native_result.artifact.vcf_path),
                "output_directory": str(directory), "native_output_role_unchanged": True,
                "commands": [[self.bgzip, "-c", str(directory / "work" / "derived.vcf")],
                             [self.tabix, "-p", "vcf", str(directory / f"{sample}.bndrepr-v1.sv.vcf.gz")]],
                "will_create_outputs": False}

    def run(self, request: BndRepresentationRequest, *, dry_run: bool = False) -> BndRepresentationResult | dict[str, object]:
        if dry_run:
            return self.plan(request)
        started = time.perf_counter()
        directory, native = request.output_directory, request.native_result
        inputs = request.inputs()
        before = {str(path.resolve()): _identity(path) for path in inputs}
        # Exclusive mkdir also rejects a directory created since request construction.
        directory.mkdir(parents=True, exist_ok=False, mode=0o700)
        work, logs = directory / "work", directory / "logs"
        work.mkdir(mode=0o700)
        logs.mkdir(mode=0o700)
        commands: list[dict[str, object]] = []
        try:
            env = {key: os.environ[key] for key in ("PATH", "SYSTEMROOT", "WINDIR") if key in os.environ}
            for key, name in (("HOME", "home"), ("TMPDIR", "tmp"), ("XDG_CONFIG_HOME", "config"), ("XDG_CACHE_HOME", "cache")):
                owned = work / name
                owned.mkdir(mode=0o700)
                env[key] = str(owned)
            env.update(TEMP=env["TMPDIR"], TMP=env["TMPDIR"], LC_ALL="C", TZ="UTC")

            def command(args: list[str], step: str, *, stdout_path: Path | None = None):
                record: dict[str, object] = {"step": step, "args": args, "status": "started",
                    "stderr_path": str(logs / f"{step}.stderr.log")}
                commands.append(record)
                try:
                    result = self.runner.run(args, cwd=work, env=env, inherit_env=False,
                        timeout=request.policy.timeout_seconds, stdout_path=stdout_path,
                        stderr_path=logs / f"{step}.stderr.log")
                except Exception as error:
                    record.update(status="failed", error_type=type(error).__name__)
                    raise
                record = result.to_dict()
                record["step"] = step
                record["status"] = "completed"
                commands[-1] = record
                if not result.executed or result.returncode != 0:
                    raise OutputValidationError("BND command did not complete successfully.")
                return result

            versions, tool_identities = {}, {}
            for name, executable in (("bgzip", self.bgzip), ("tabix", self.tabix)):
                resolved = self.runner.require_executable(executable, env=env, inherit_env=False, cwd=work)
                tool_identities[name] = _identity(resolved)
                result = command([str(resolved), "--version"], f"{name}_version")
                output = Path(result.stdout_path).read_text(encoding="utf-8") if result.stdout_path else result.stdout or ""
                match = re.search(r"(?<![0-9])(\d+)\.(\d+)(?:\.(\d+))?", output)
                if match is None:
                    raise ToolVersionError(f"Unable to parse {name} version.")
                version = tuple(int(value or 0) for value in match.groups())
                if version < VCF44_MIN_TABIX:
                    raise ToolVersionError("Explicit BND VCFv4.4 representation requires bgzip/tabix >=1.23.1.")
                versions[name] = ".".join(str(value) for value in version)
            native_dir = directory / "native"
            native_dir.mkdir(mode=0o700)
            copies = {}
            for role, path in (("vcf", native.artifact.vcf_path), ("index", native.artifact.index_path),
                               ("raw", request.native_raw), ("evidence", native.artifact.evidence_table)):
                target = native_dir / Path(path).name
                _copy_exact(Path(path), target, before[str(Path(path).resolve())])
                copies[role] = str(target.relative_to(directory))
            _write_json(directory / "native_result.json", native.to_dict())
            plain, ledger, junctions = work / "derived.vcf", directory / "records.jsonl", directory / "junctions.tsv"
            summary = derive_record_stream(native.artifact.vcf_path, plain, ledger, junctions,
                reference_lengths={contig.name: contig.length for contig in native.request.reference.contigs},
                sample_id=native.request.sample_id, native_sha256=before[str(native.artifact.vcf_path.resolve())]["sha256"],
                policy=request.policy, work_directory=work)
            compressed = directory / f"{native.request.sample_id}.bndrepr-v1.sv.vcf.gz"
            command([str(tool_identities["bgzip"]["path"]), "-c", str(plain)], "bgzip", stdout_path=compressed)
            command([str(tool_identities["tabix"]["path"]), "-p", "vcf", str(compressed)], "tabix_index")
            index = Path(f"{compressed}.tbi")
            if not index.is_file() or index.stat().st_size == 0:
                raise OutputValidationError("BND output index is missing or empty.")
            with compressed.open("rb") as reader:
                prefix = reader.read(18)
                if prefix[:4] != b"\x1f\x8b\x08\x04" or prefix[12:16] != b"BC\x02\x00":
                    raise OutputValidationError("Derived VCF is not standard BGZF output.")
            summary["compressed_round_trip"] = verify_record_round_trip(native.artifact.vcf_path, compressed, ledger)
            command([str(tool_identities["tabix"]["path"]), "-l", str(compressed)], "tabix_list")
            # Full sequential equality is distinct from the bounded local index probes.
            local_queries = []
            with ledger.open(encoding="utf-8") as reader:
                for line in reader:
                    row = json.loads(line)
                    junction = row["junction"]
                    if not junction or len(local_queries) >= 4:
                        continue
                    position = junction["position"]
                    # Virtual boundaries are retained, not converted into ordinary intervals.
                    if junction["local_virtual"] or ":" in junction["contig"]:
                        local_queries.append({"ordinal": row["record_ordinal"], "status": "NOT_TESTED_NONORDINARY_QUERY", "junction": junction})
                        continue
                    region = f"{junction['contig']}:{position}-{position}"
                    result = command([str(tool_identities["tabix"]["path"]), str(compressed), region], f"tabix_query_{row['record_ordinal']}")
                    expected_hash = row["derived_record_sha256"]
                    with Path(result.stdout_path).open("rb") as records:
                        ending = b"\r\n" if row["eol"] == "CRLF" else b"\n"
                        found = any(hashlib.sha256(record.rstrip(b"\r\n") + ending).hexdigest() == expected_hash for record in records)
                    if not found:
                        raise OutputValidationError("Local index query did not return the expected BND record bytes.")
                    local_queries.append({"ordinal": row["record_ordinal"], "region": region, "status": "PASS"})
            summary["local_index_queries"] = local_queries
            after = {str(path.resolve()): _identity(path) for path in inputs}
            if after != before:
                raise OutputValidationError("An upstream BND/native/source/reference input changed.")
            for name, identity in tool_identities.items():
                if _identity(Path(identity["path"])) != identity:
                    raise OutputValidationError(f"{name} identity changed during BND derivation.")
            payload = {"schema_version": 1, "status": "completed", "policy": request.policy.to_dict(),
                "native_output_role_unchanged": True, "consumer_opt_in_required": True,
                "standalone_derived_VCF_is_lossless": False, "native_binary_copies": copies,
                "native_inputs_before": before, "native_inputs_after": after,
                "source_order": [source.caller for source in native.request.runnable_sources],
                "tools": tool_identities, "tool_versions": versions, "commands": commands,
                "summary": summary, "runtime_seconds": time.perf_counter() - started,
                "artifact_identities": {str(path.relative_to(directory)): _identity(path) for path in
                    (compressed, index, ledger, junctions, directory / "native_result.json", *(native_dir / name for name in [Path(p).name for p in copies.values()]))}}
            manifest = directory / "manifest.json"
            temporary = directory / "manifest.partial.json"
            _write_json(temporary, payload)
            temporary.rename(manifest)
            _LOGGER.info("BND representation completed for sample %s; native artifacts unchanged.", native.request.sample_id)
            return BndRepresentationResult(directory, manifest, compressed, index, ledger, junctions, summary)
        except Exception as error:
            _write_json(directory / "FAILED.json", {"status": "failed", "error_type": type(error).__name__,
                "success_manifest_published": False, "commands_completed": commands,
                "native_inputs_before": before, "no_existing_file_replaced_or_removed": True})
            raise


__all__ = ["BndRepresentationRequest", "BndRepresentationResult", "BndRepresentationWrapper"]
