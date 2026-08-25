"""Jasmine primary SV clustering wrapper without truth inference."""

from __future__ import annotations

import gzip
import heapq
import shutil
import os
import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hifivar.command import CommandRunner, format_command
from hifivar.exceptions import OutputValidationError, ToolVersionError
from hifivar.harmonization import HarmonizedSvArtifact, SVHarmonizationRequest, validate_harmonized_artifact, write_evidence_table
from hifivar.validation import validate_output_file


_VERSION = re.compile(r"(\d+(?:\.\d+)+)")
_CONTIG_HEADER = re.compile(r"^##contig=<ID=([^,>]+)")
_FORMAT_HEADER = re.compile(r"^##FORMAT=<ID=([^,>]+)")
_SORT_CHUNK_RECORDS = 100_000


class JasmineResultStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class JasmineCommandPlan:
    step: str
    args: tuple[str, ...]
    stdout_path: Path | None = None

    def to_dict(self) -> dict[str, object]:
        return {"step": self.step, "args": list(self.args), "display_command": format_command(self.args), "stdout_path": str(self.stdout_path) if self.stdout_path else None, "shell": False}


@dataclass(frozen=True, slots=True)
class JasmineResult:
    request: SVHarmonizationRequest
    status: JasmineResultStatus
    commands: tuple[JasmineCommandPlan, ...]
    tool_versions: dict[str, str] | None = None
    runtime_seconds: float = 0.0
    artifact: HarmonizedSvArtifact | None = None

    def to_dict(self) -> dict[str, object]:
        return {"request": self.request.to_dict(), "status": self.status.value,
                "commands": [item.to_dict() for item in self.commands],
                "tool_versions": self.tool_versions, "runtime_seconds": self.runtime_seconds,
                "artifact": self.artifact.to_dict() if self.artifact else None}


class JasmineWrapper:
    def __init__(self, *, executable="jasmine", bgzip_executable="bgzip",
                 tabix_executable="tabix", runner: CommandRunner | None = None) -> None:
        self.executable, self.bgzip, self.tabix = executable, bgzip_executable, tabix_executable
        self.runner = runner or CommandRunner()

    def plan_commands(
        self,
        request: SVHarmonizationRequest,
        *,
        jasmine_prefix: tuple[str, ...] | None = None,
    ) -> tuple[JasmineCommandPlan, ...]:
        raw = request.work_directory / f"{request.sample_id}.jasmine.raw.vcf"
        sorted_vcf = request.work_directory / f"{request.sample_id}.jasmine.sorted.vcf"
        listing = request.work_directory / "jasmine.inputs.txt"
        prefix = jasmine_prefix or (self.executable,)
        args = [
            *prefix,
            f"file_list={listing.absolute()}",
            f"out_file={raw.absolute()}",
            f"genome_file={request.reference.fasta.absolute()}",
            f"max_dist={request.max_dist}",
        ]
        if request.distance_type == "nonlinear":
            args.append("--nonlinear_dist")
        return (
            JasmineCommandPlan("jasmine", tuple(args)),
            JasmineCommandPlan(
                "bgzip",
                (self.bgzip, "-c", str(sorted_vcf.absolute())),
                request.output_vcf,
            ),
            JasmineCommandPlan(
                "tabix",
                (self.tabix, "-p", "vcf", str(request.output_vcf.absolute())),
            ),
        )

    def detect_versions(
        self,
        *,
        jasmine_prefix: tuple[str, ...] | None = None,
    ) -> dict[str, str]:
        prefix = jasmine_prefix or self._resolve_jasmine_prefix()
        commands = {
            "jasmine": (*prefix, "--version"),
            "bgzip": (self.bgzip, "--version"),
            "tabix": (self.tabix, "--version"),
        }
        versions: dict[str, str] = {}
        for name, command in commands.items():
            self.runner.require_executable(command[0])
            result = self.runner.run(command)
            output = "\n".join(
                item for item in (result.stdout, result.stderr) if item
            )
            match = _VERSION.search(output)
            if match is None:
                raise ToolVersionError(f"Unable to parse {name} version.")
            versions[name] = match.group(1)
        return versions

    def _resolve_jasmine_prefix(self) -> tuple[str, ...]:
        resolved = self.runner.require_executable(self.executable)
        if not resolved.is_file():
            return (self.executable,)
        try:
            with resolved.open("rb") as handle:
                prefix = handle.read(512)
        except OSError as error:
            raise ToolVersionError(
                f"Unable to inspect Jasmine launcher '{resolved}': {error}"
            ) from error
        try:
            prefix.decode("utf-8")
            is_text_script = b"\x00" not in prefix
        except UnicodeDecodeError:
            is_text_script = False
        if is_text_script and not prefix.startswith(b"#!"):
            self.runner.require_executable("bash")
            return ("bash", str(resolved.absolute()))
        return (self.executable,)

    def run(
        self,
        request: SVHarmonizationRequest,
        *,
        dry_run: bool = False,
        stderr_path: Path | None = None,
    ) -> JasmineResult:
        if dry_run:
            commands = self.plan_commands(request)
            for item in commands:
                self.runner.run(
                    item.args,
                    dry_run=True,
                    stdout_path=item.stdout_path,
                )
            return JasmineResult(
                request,
                JasmineResultStatus.PLANNED,
                commands,
            )

        jasmine_prefix = self._resolve_jasmine_prefix()
        commands = self.plan_commands(
            request,
            jasmine_prefix=jasmine_prefix,
        )
        versions = self.detect_versions(jasmine_prefix=jasmine_prefix)
        unzipped_inputs = self._prepare(request)
        raw = request.work_directory / f"{request.sample_id}.jasmine.raw.vcf"
        sorted_vcf = request.work_directory / f"{request.sample_id}.jasmine.sorted.vcf"

        runtime = 0.0
        jasmine_result = self.runner.run(
            commands[0].args,
            stderr_path=_step_log(stderr_path, commands[0].step),
        )
        runtime += jasmine_result.duration_seconds
        validate_output_file(raw)
        _sort_jasmine_vcf(
            raw,
            sorted_vcf,
            unzipped_inputs,
            overwrite=request.overwrite,
        )
        for item in commands[1:]:
            result = self.runner.run(
                item.args,
                stdout_path=item.stdout_path,
                stderr_path=_step_log(stderr_path, item.step),
            )
            runtime += result.duration_seconds

        write_evidence_table(request, request.output_vcf)
        artifact = HarmonizedSvArtifact(
            sample_id=request.sample_id,
            reference_fasta=request.reference.fasta,
            reference_build=request.reference.build,
            vcf_path=request.output_vcf,
            index_path=request.output_index,
            evidence_table=request.evidence_table,
            sources=request.sources,
            jasmine_version=versions["jasmine"],
            commands=tuple(item.args for item in commands),
            intermediate_files=(
                request.work_directory / "jasmine.inputs.txt",
                *unzipped_inputs,
                raw,
                sorted_vcf,
            ),
        )
        validated = validate_harmonized_artifact(
            artifact,
            request.reference,
        )
        return JasmineResult(
            request,
            JasmineResultStatus.COMPLETED,
            commands,
            versions,
            runtime,
            validated,
        )

    @staticmethod
    def _prepare(
        request: SVHarmonizationRequest,
    ) -> tuple[Path, ...]:
        request.work_directory.mkdir(parents=True, exist_ok=True)
        request.output_vcf.parent.mkdir(parents=True, exist_ok=True)
        listing = request.work_directory / "jasmine.inputs.txt"
        listing_temporary = listing.with_name(
            f".{listing.name}.hifivar.tmp"
        )
        raw = request.work_directory / f"{request.sample_id}.jasmine.raw.vcf"
        sorted_vcf = (
            request.work_directory
            / f"{request.sample_id}.jasmine.sorted.vcf"
        )
        for owned in (listing, listing_temporary, raw, sorted_vcf):
            _prepare_owned_file(owned, overwrite=request.overwrite)

        unzipped: list[Path] = []
        for index, source in enumerate(request.runnable_sources):
            if source.vcf_path is None or source.index_path is None:
                raise OutputValidationError(
                    f"Runnable caller '{source.caller}' lacks VCF/index."
                )
            validate_output_file(source.vcf_path)
            validate_output_file(source.index_path)
            destination = (
                request.work_directory
                / f"{index:02d}.{source.caller}.unzipped.vcf"
            )
            _decompress_vcf(
                source.vcf_path,
                destination,
                overwrite=request.overwrite,
            )
            unzipped.append(destination)

        listing_temporary.write_text(
            "\n".join(str(path.absolute()) for path in unzipped) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        os.replace(listing_temporary, listing)
        return tuple(unzipped)


def _prepare_owned_file(path: Path, *, overwrite: bool) -> None:
    if path.exists() and path.is_dir():
        raise OutputValidationError(
            f"Jasmine owned file is a directory: '{path}'."
        )
    if path.exists() and not overwrite:
        raise OutputValidationError(
            f"Jasmine owned file already exists: '{path}'."
        )
    if path.exists():
        path.unlink()


def _decompress_vcf(
    source: Path,
    destination: Path,
    *,
    overwrite: bool,
) -> None:
    _prepare_owned_file(destination, overwrite=overwrite)
    temporary = destination.with_name(
        f".{destination.name}.hifivar.tmp"
    )
    _prepare_owned_file(temporary, overwrite=overwrite)
    try:
        with gzip.open(source, "rb") as reader:
            with temporary.open("xb") as writer:
                shutil.copyfileobj(reader, writer, length=1024 * 1024)
        validate_output_file(temporary)
        os.replace(temporary, destination)
    except (OSError, EOFError) as error:
        if temporary.exists():
            temporary.unlink()
        raise OutputValidationError(
            f"Unable to stream-decompress Jasmine input '{source}': {error}"
        ) from error


def _source_format_headers(
    sources: tuple[Path, ...],
) -> dict[str, str]:
    definitions: dict[str, str] = {}
    for source in sources:
        try:
            with source.open("r", encoding="utf-8", newline="") as handle:
                for raw_line in handle:
                    line = raw_line.rstrip("\r\n")
                    match = _FORMAT_HEADER.match(line)
                    if match and match.group(1) not in definitions:
                        definitions[match.group(1)] = line
                    if line.startswith("#CHROM\t"):
                        break
        except (OSError, UnicodeError) as error:
            raise OutputValidationError(
                f"Unable to read decompressed Jasmine input '{source}': {error}"
            ) from error
    return definitions


def _sort_jasmine_vcf(
    raw_vcf: Path,
    destination: Path,
    source_vcfs: tuple[Path, ...],
    *,
    overwrite: bool,
) -> None:
    """Sort Jasmine records without loading the complete VCF into memory.

    Jasmine 1.1.5 may emit chromosome blocks out of reference order and may
    omit FORMAT declarations copied from its inputs.  Both conditions make a
    biologically valid merge impossible to index.  Preserve the raw output,
    restore only definitions present in the source headers, and external-sort
    records using the contig order declared by Jasmine's output header.
    """

    _prepare_owned_file(destination, overwrite=overwrite)
    temporary = destination.with_name(f".{destination.name}.hifivar.tmp")
    _prepare_owned_file(temporary, overwrite=overwrite)
    source_formats = _source_format_headers(source_vcfs)
    headers: list[str] = []
    chrom_header: str | None = None
    contig_order: dict[str, int] = {}
    declared_formats: set[str] = set()
    used_formats: set[str] = set()
    chunks: list[Path] = []
    records: list[str] = []

    try:
        with raw_vcf.open("r", encoding="utf-8", newline="") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\r\n")
                if line.startswith("##"):
                    headers.append(line)
                    contig_match = _CONTIG_HEADER.match(line)
                    if contig_match and contig_match.group(1) not in contig_order:
                        contig_order[contig_match.group(1)] = len(contig_order)
                    format_match = _FORMAT_HEADER.match(line)
                    if format_match:
                        declared_formats.add(format_match.group(1))
                    continue
                if line.startswith("#CHROM\t"):
                    chrom_header = line
                    continue
                if line.startswith("#"):
                    headers.append(line)
                    continue
                if not line:
                    continue
                if chrom_header is None:
                    raise OutputValidationError(
                        "Jasmine VCF contains records before the #CHROM header."
                    )
                fields = line.split("\t")
                if len(fields) < 8:
                    raise OutputValidationError("Malformed Jasmine VCF record.")
                _record_sort_key(line, contig_order)
                if len(fields) > 8 and fields[8] not in ("", "."):
                    used_formats.update(fields[8].split(":"))
                records.append(line + "\n")
                if len(records) >= _SORT_CHUNK_RECORDS:
                    chunk = destination.with_name(
                        f".{destination.name}.sort.{len(chunks):06d}.tmp"
                    )
                    _write_sort_chunk(
                        records,
                        chunk,
                        contig_order,
                        overwrite=overwrite,
                    )
                    chunks.append(chunk)
                    records = []

        if chrom_header is None:
            raise OutputValidationError("Jasmine VCF lacks a #CHROM header.")
        if not any(line.startswith("##fileformat=") for line in headers):
            raise OutputValidationError("Jasmine VCF lacks ##fileformat.")
        missing_formats = used_formats - declared_formats
        unknown_formats = missing_formats - source_formats.keys()
        if unknown_formats:
            names = ", ".join(sorted(unknown_formats))
            raise OutputValidationError(
                "Jasmine used FORMAT identifiers without source definitions: "
                f"{names}."
            )

        with temporary.open("x", encoding="utf-8", newline="\n") as output:
            for header in headers:
                output.write(header + "\n")
            for name, definition in source_formats.items():
                if name in missing_formats:
                    output.write(definition + "\n")
            output.write(chrom_header + "\n")
            if chunks:
                if records:
                    chunk = destination.with_name(
                        f".{destination.name}.sort.{len(chunks):06d}.tmp"
                    )
                    _write_sort_chunk(
                        records,
                        chunk,
                        contig_order,
                        overwrite=overwrite,
                    )
                    chunks.append(chunk)
                iterators = tuple(_iter_chunk(path) for path in chunks)
                for record in heapq.merge(
                    *iterators,
                    key=lambda item: _record_sort_key(item, contig_order),
                ):
                    output.write(record)
            else:
                records.sort(key=lambda item: _record_sort_key(item, contig_order))
                output.writelines(records)
        validate_output_file(temporary)
        os.replace(temporary, destination)
    except (OSError, UnicodeError) as error:
        raise OutputValidationError(
            f"Unable to sort Jasmine VCF '{raw_vcf}': {error}"
        ) from error
    finally:
        if temporary.exists():
            temporary.unlink()
        for chunk in chunks:
            if chunk.exists():
                chunk.unlink()


def _record_sort_key(
    line: str,
    contig_order: dict[str, int],
) -> tuple[int, int, str]:
    fields = line.split("\t", 3)
    if len(fields) < 3:
        raise OutputValidationError("Malformed Jasmine VCF record.")
    contig = fields[0]
    if contig not in contig_order:
        raise OutputValidationError(
            f"Jasmine record contig lacks a header declaration: '{contig}'."
        )
    try:
        position = int(fields[1])
    except ValueError as error:
        raise OutputValidationError(
            f"Jasmine VCF position is not an integer: '{fields[1]}'."
        ) from error
    return contig_order[contig], position, line


def _write_sort_chunk(
    records: list[str],
    path: Path,
    contig_order: dict[str, int],
    *,
    overwrite: bool,
) -> None:
    _prepare_owned_file(path, overwrite=overwrite)
    records.sort(key=lambda line: _record_sort_key(line, contig_order))
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        handle.writelines(records)


def _iter_chunk(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        yield from handle


def _step_log(path: Path | None, step: str) -> Path | None:
    return None if path is None else path.with_name(f"{path.stem}.{step}{path.suffix}")


__all__ = ["JasmineCommandPlan", "JasmineResult", "JasmineResultStatus", "JasmineWrapper"]
