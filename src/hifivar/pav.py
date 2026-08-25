"""PAV workflow execution adapter for independent assembly-derived SV evidence."""

from __future__ import annotations

import json
import os
import re
import shutil
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hifivar.assembly import AssemblyRole
from hifivar.assembly_sv import AssemblySvArtifact, AssemblySvCaller, AssemblySvRequest, create_assembly_sv_artifact
from hifivar.command import CommandRunner, format_command
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar.validation import validate_output_file
from hifivar.sv import validate_sv_vcf


_VERSION = re.compile(r"\d+(?:\.\d+)+")
_PENDING_VERSION = "VERSION_PENDING_LINUX_VERIFICATION"


class PavResultStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class PavCommandPlan:
    args: tuple[str, ...]
    display_command: str
    backend: str = "native_snakemake_adapter"

    def to_dict(self) -> dict[str, object]:
        return {"args": list(self.args), "display_command": self.display_command, "backend": self.backend, "shell": False}


@dataclass(frozen=True, slots=True)
class PavResult:
    request: AssemblySvRequest
    status: PavResultStatus
    command: PavCommandPlan
    pav_version: str
    adapter_version: str | None = None
    runtime_seconds: float = 0.0
    artifact: AssemblySvArtifact | None = None
    pav_version_source: str = "config"

    def to_dict(self) -> dict[str, object]:
        return {
            "request": self.request.to_dict(), "status": self.status.value,
            "command": self.command.to_dict(), "pav_version": self.pav_version,
            "adapter_version": self.adapter_version, "runtime_seconds": self.runtime_seconds,
            "artifact": self.artifact.to_dict() if self.artifact else None,
            "pav_version_source": self.pav_version_source,
        }


class PavWrapper:
    """Run a configured PAV Snakefile as a transparent workflow adapter."""

    def __init__(self, *, snakefile: str | Path, executable: str = "snakemake",
                 pav_version: str = _PENDING_VERSION,
                 runner: CommandRunner | None = None) -> None:
        self.snakefile = Path(snakefile)
        self.executable = executable
        self.pav_version = pav_version
        self.runner = runner or CommandRunner()
    def plan_command(self, request: AssemblySvRequest) -> PavCommandPlan:
        if request.caller is not AssemblySvCaller.PAV:
            raise InputValidationError("PAV wrapper requires caller=pav.")
        self._validate_inputs(request)
        args = (
            self.executable, "--snakefile", str(self.snakefile.absolute()),
            "--directory", str(request.work_directory.absolute()),
            "--cores", str(request.resources.threads),
        )
        return PavCommandPlan(args, format_command(args))

    def detect_adapter_version(self) -> str:
        self.runner.require_executable(self.executable)
        result = self.runner.run([self.executable, "--version"])
        output = "\n".join(item for item in (result.stdout, result.stderr) if item)
        match = _VERSION.search(output)
        if match is None:
            raise ToolVersionError("Unable to parse PAV Snakemake adapter version.")
        return match.group(0)

    def run(self, request: AssemblySvRequest, *, dry_run: bool = False,
            stderr_path: Path | None = None) -> PavResult:
        command = self.plan_command(request)
        if dry_run:
            self.runner.run(command.args, dry_run=True)
            version_source = "unresolved" if self.pav_version == _PENDING_VERSION else "config"
            return PavResult(
                request,
                PavResultStatus.PLANNED,
                command,
                self.pav_version,
                pav_version_source=version_source,
            )
        pav_version = self._require_pav_version()
        adapter_version = self.detect_adapter_version()
        self._prepare(request)
        result = self.runner.run(command.args, cwd=request.work_directory, stderr_path=stderr_path)
        raw_vcf = request.work_directory / f"{request.sample_id}.vcf.gz"
        raw_index = Path(f"{raw_vcf}.tbi")
        validate_sv_vcf(
            raw_vcf,
            raw_index,
            sample_id=request.sample_id,
            reference=request.reference,
        )
        _copy_atomic(raw_vcf, request.output_vcf, overwrite=request.overwrite)
        _copy_atomic(raw_index, request.output_index, overwrite=request.overwrite)
        artifact = create_assembly_sv_artifact(
            request, raw_vcf=raw_vcf,
            intermediate_files=(
                request.work_directory / "config.json",
                request.work_directory / "assemblies.tsv",
                raw_vcf,
                raw_index,
            ),
            caller_version=pav_version, backend=command.backend,
            commands=(command.args,),
        )
        return PavResult(
            request,
            PavResultStatus.COMPLETED,
            command,
            pav_version,
            adapter_version,
            result.duration_seconds,
            artifact,
            "config",
        )

    def _require_pav_version(self) -> str:
        version = self.pav_version.strip()
        if version == _PENDING_VERSION or _VERSION.fullmatch(version) is None:
            raise ToolVersionError(
                "PAV version must be configured as an explicit numeric release "
                "(for example, '2.4.6') before real execution."
            )
        return version

    def _validate_inputs(self, request: AssemblySvRequest) -> None:
        validate_output_file(request.reference.fasta)
        validate_output_file(request.reference.fai)
        if not self.snakefile.exists() or not self.snakefile.is_file():
            raise InputValidationError(f"PAV Snakefile is missing: '{self.snakefile}'.")
        for assembly in request.assemblies:
            validate_output_file(assembly.path)

    def _prepare(self, request: AssemblySvRequest) -> None:
        if request.output_vcf.exists() and not request.overwrite:
            raise OutputValidationError(f"PAV output already exists: '{request.output_vcf}'.")
        request.work_directory.mkdir(parents=True, exist_ok=True)
        request.output_vcf.parent.mkdir(parents=True, exist_ok=True)
        config = request.work_directory / "config.json"
        table = request.work_directory / "assemblies.tsv"
        _write_atomic(
            config,
            json.dumps({"reference": str(request.reference.fasta.absolute())}, indent=2) + "\n",
            overwrite=request.overwrite,
        )
        columns = ["NAME", *[f"HAP_{'h1' if item.role is AssemblyRole.HAPLOTYPE1 else 'h2'}" for item in request.assemblies]]
        values = [request.sample_id, *[str(item.path.absolute()) for item in request.assemblies]]
        _write_atomic(
            table,
            "\t".join(columns) + "\n" + "\t".join(values) + "\n",
            overwrite=request.overwrite,
        )

def _write_atomic(path: Path, text: str, *, overwrite: bool) -> None:
    temporary = path.with_name(f".{path.name}.hifivar.tmp")
    for owned in (path, temporary):
        if owned.exists() and owned.is_dir():
            raise OutputValidationError(f"PAV owned config is a directory: '{owned}'.")
        if owned.exists() and not overwrite:
            raise OutputValidationError(f"PAV owned config already exists: '{owned}'.")
    if temporary.exists():
        temporary.unlink()
    try:
        temporary.write_text(text, encoding="utf-8", newline="\n")
        os.replace(temporary, path)
    except OSError as error:
        if temporary.exists():
            temporary.unlink()
        raise OutputValidationError(
            f"Unable to write PAV owned config '{path}': {error}"
        ) from error


def _copy_atomic(source: Path, destination: Path, *, overwrite: bool) -> None:
    validate_output_file(source)
    temporary = destination.with_name(f".{destination.name}.hifivar.tmp")
    for owned in (destination, temporary):
        if owned.exists() and owned.is_dir():
            raise OutputValidationError(
                f"PAV owned output is a directory: '{owned}'."
            )
        if owned.exists() and not overwrite:
            raise OutputValidationError(
                f"PAV output already exists: '{owned}'."
            )
        if owned.exists():
            owned.unlink()
    try:
        with source.open("rb") as reader:
            with temporary.open("xb") as writer:
                shutil.copyfileobj(reader, writer, length=1024 * 1024)
        validate_output_file(temporary)
        os.replace(temporary, destination)
    except OSError as error:
        if temporary.exists():
            temporary.unlink()
        raise OutputValidationError(
            f"Unable to capture PAV output '{source}' as '{destination}': {error}"
        ) from error


__all__ = ["PavCommandPlan", "PavResult", "PavResultStatus", "PavWrapper"]
