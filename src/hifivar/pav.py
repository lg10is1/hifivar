"""PAV workflow execution adapter for independent assembly-derived SV evidence."""

from __future__ import annotations

import gzip
import json
import os
import re
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
_PAV_SOURCE = re.compile(r"^##source=PAV (\d+(?:\.\d+)+)$")
_PENDING_VERSION = "VERSION_PENDING_LINUX_VERIFICATION"
_VALIDATED_PAV_VERSION = (2, 4, 6)
_PAV_SV_MIN_LENGTH = 50
_PAV_SV_SELECTION_POLICY = (
    "PAV_2.4.6_VARTYPE:SVTYPE=INV_OR_ABS(SVLEN)>=50"
)


class PavResultStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class PavCommandPlan:
    args: tuple[str, ...]
    display_command: str
    backend: str = "native_snakemake_adapter"
    step: str = "pav"
    stdout_path: Path | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "step": self.step,
            "args": list(self.args),
            "display_command": self.display_command,
            "backend": self.backend,
            "stdout_path": str(self.stdout_path) if self.stdout_path else None,
            "shell": False,
        }


@dataclass(frozen=True, slots=True)
class PavSvSelection:
    source_vcf: Path
    plain_vcf: Path
    policy: str
    total_records: int | None = None
    selected_records: int | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "source_vcf": str(self.source_vcf),
            "plain_vcf": str(self.plain_vcf),
            "policy": self.policy,
            "total_records": self.total_records,
            "selected_records": self.selected_records,
        }


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
    finalization_commands: tuple[PavCommandPlan, ...] = ()
    finalizer_versions: dict[str, str] | None = None
    selection: PavSvSelection | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "request": self.request.to_dict(), "status": self.status.value,
            "command": self.command.to_dict(), "pav_version": self.pav_version,
            "adapter_version": self.adapter_version, "runtime_seconds": self.runtime_seconds,
            "artifact": self.artifact.to_dict() if self.artifact else None,
            "pav_version_source": self.pav_version_source,
            "finalization_commands": [
                item.to_dict() for item in self.finalization_commands
            ],
            "finalizer_versions": self.finalizer_versions,
            "selection": self.selection.to_dict() if self.selection else None,
        }


class PavWrapper:
    """Run a configured PAV Snakefile as a transparent workflow adapter."""

    def __init__(self, *, snakefile: str | Path, executable: str = "snakemake",
                 bgzip_executable: str = "bgzip", tabix_executable: str = "tabix",
                 pav_version: str = _PENDING_VERSION,
                 runner: CommandRunner | None = None) -> None:
        self.snakefile = Path(snakefile)
        self.executable = executable
        self.bgzip = bgzip_executable
        self.tabix = tabix_executable
        self.pav_version = pav_version
        self.runner = runner or CommandRunner()

    def plan_command(self, request: AssemblySvRequest) -> PavCommandPlan:
        return self.plan_commands(request)[0]

    def plan_commands(
        self,
        request: AssemblySvRequest,
    ) -> tuple[PavCommandPlan, PavCommandPlan, PavCommandPlan]:
        if request.caller is not AssemblySvCaller.PAV:
            raise InputValidationError("PAV wrapper requires caller=pav.")
        self._validate_inputs(request)
        pav_args = (
            self.executable, "--snakefile", str(self.snakefile.absolute()),
            "--directory", str(request.work_directory.absolute()),
            "--cores", str(request.resources.threads),
        )
        plain_vcf = _selection_vcf(request)
        bgzip_args = (self.bgzip, "-c", str(plain_vcf.absolute()))
        tabix_args = (
            self.tabix,
            "-p",
            "vcf",
            str(request.output_vcf.absolute()),
        )
        return (
            PavCommandPlan(pav_args, format_command(pav_args)),
            PavCommandPlan(
                bgzip_args,
                format_command(bgzip_args),
                "native_snakemake_adapter",
                "bgzip_sv_only",
                request.output_vcf,
            ),
            PavCommandPlan(
                tabix_args,
                format_command(tabix_args),
                "native_snakemake_adapter",
                "tabix_sv_only",
            ),
        )

    def detect_adapter_version(self) -> str:
        self.runner.require_executable(self.executable)
        result = self.runner.run([self.executable, "--version"])
        output = "\n".join(item for item in (result.stdout, result.stderr) if item)
        match = _VERSION.search(output)
        if match is None:
            raise ToolVersionError("Unable to parse PAV Snakemake adapter version.")
        return match.group(0)

    def detect_finalizer_versions(self) -> dict[str, str]:
        versions: dict[str, str] = {}
        for name, executable in (("bgzip", self.bgzip), ("tabix", self.tabix)):
            self.runner.require_executable(executable)
            result = self.runner.run([executable, "--version"])
            output = "\n".join(
                item for item in (result.stdout, result.stderr) if item
            )
            match = _VERSION.search(output)
            if match is None:
                raise ToolVersionError(f"Unable to parse {name} version.")
            versions[name] = match.group(0)
        return versions

    def run(self, request: AssemblySvRequest, *, dry_run: bool = False,
            stderr_path: Path | None = None) -> PavResult:
        command, *finalization_commands = self.plan_commands(request)
        finalizers = tuple(finalization_commands)
        planned_selection = PavSvSelection(
            request.work_directory / f"{request.sample_id}.vcf.gz",
            _selection_vcf(request),
            _PAV_SV_SELECTION_POLICY,
        )
        if dry_run:
            self.runner.run(command.args, dry_run=True)
            for item in finalizers:
                self.runner.run(
                    item.args,
                    dry_run=True,
                    stdout_path=item.stdout_path,
                )
            version_source = "unresolved" if self.pav_version == _PENDING_VERSION else "config"
            return PavResult(
                request,
                PavResultStatus.PLANNED,
                command,
                self.pav_version,
                pav_version_source=version_source,
                finalization_commands=finalizers,
                selection=planned_selection,
            )
        pav_version = self._require_pav_version()
        adapter_version = self.detect_adapter_version()
        finalizer_versions = self.detect_finalizer_versions()
        self._prepare(request)
        result = self.runner.run(command.args, cwd=request.work_directory, stderr_path=stderr_path)
        runtime_seconds = result.duration_seconds
        raw_vcf = request.work_directory / f"{request.sample_id}.vcf.gz"
        raw_index = Path(f"{raw_vcf}.tbi")
        validate_sv_vcf(
            raw_vcf,
            raw_index,
            sample_id=request.sample_id,
            reference=request.reference,
        )
        _validate_pav_source_version(raw_vcf, pav_version)
        selection = _write_pav_sv_only_vcf(
            raw_vcf,
            _selection_vcf(request),
            overwrite=request.overwrite,
        )
        bgzip_result = self.runner.run(
            finalizers[0].args,
            stdout_path=request.output_vcf,
            stderr_path=_step_log(stderr_path, finalizers[0].step),
        )
        runtime_seconds += bgzip_result.duration_seconds
        validate_output_file(request.output_vcf)
        tabix_result = self.runner.run(
            finalizers[1].args,
            stderr_path=_step_log(stderr_path, finalizers[1].step),
        )
        runtime_seconds += tabix_result.duration_seconds
        artifact = create_assembly_sv_artifact(
            request, raw_vcf=raw_vcf,
            intermediate_files=(
                request.work_directory / "config.json",
                request.work_directory / "assemblies.tsv",
                raw_vcf,
                raw_index,
                selection.plain_vcf,
            ),
            caller_version=pav_version, backend=command.backend,
            commands=(command.args, *(item.args for item in finalizers)),
        )
        return PavResult(
            request,
            PavResultStatus.COMPLETED,
            command,
            pav_version,
            adapter_version,
            runtime_seconds,
            artifact,
            "config",
            finalizers,
            finalizer_versions,
            selection,
        )

    def _require_pav_version(self) -> str:
        version = self.pav_version.strip()
        if version == _PENDING_VERSION or _VERSION.fullmatch(version) is None:
            raise ToolVersionError(
                "PAV version must be configured as an explicit numeric release "
                "(for example, '2.4.6') before real execution."
            )
        parsed = tuple(int(item) for item in version.split("."))
        if parsed[:3] != _VALIDATED_PAV_VERSION:
            raise ToolVersionError(
                "PAV SV-only selection is validated only for PAV 2.4.6; "
                f"configured version is '{version}'."
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
        for output in (request.output_vcf, request.output_index):
            if output.exists() and not request.overwrite:
                raise OutputValidationError(f"PAV output already exists: '{output}'.")
            if output.exists() and output.is_dir():
                raise OutputValidationError(f"PAV output is a directory: '{output}'.")
            if output.exists():
                output.unlink()
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


def _selection_vcf(request: AssemblySvRequest) -> Path:
    return request.work_directory / f"{request.sample_id}.pav.sv_only.vcf"


def _validate_pav_source_version(raw_vcf: Path, configured_version: str) -> None:
    observed: str | None = None
    try:
        with gzip.open(raw_vcf, "rt", encoding="utf-8", newline="") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\r\n")
                match = _PAV_SOURCE.match(line)
                if match:
                    observed = match.group(1)
                    break
                if line.startswith("#CHROM\t") or not line.startswith("#"):
                    break
    except (OSError, EOFError, UnicodeError) as error:
        raise OutputValidationError(
            f"Unable to inspect PAV source version in '{raw_vcf}': {error}"
        ) from error
    if observed is None:
        raise OutputValidationError("PAV VCF header lacks an exact ##source=PAV version.")
    configured = tuple(int(item) for item in configured_version.split("."))
    actual = tuple(int(item) for item in observed.split("."))
    if configured[:3] != actual[:3]:
        raise OutputValidationError(
            f"PAV VCF version '{observed}' differs from configured version "
            f"'{configured_version}'."
        )


def _write_pav_sv_only_vcf(
    raw_vcf: Path,
    destination: Path,
    *,
    overwrite: bool,
) -> PavSvSelection:
    """Stream PAV's mixed VCF into its own 2.4.6 VARTYPE=SV subset."""

    temporary = destination.with_name(f".{destination.name}.hifivar.tmp")
    for owned in (destination, temporary):
        if owned.exists() and owned.is_dir():
            raise OutputValidationError(
                f"PAV SV-only intermediate is a directory: '{owned}'."
            )
        if owned.exists() and not overwrite:
            raise OutputValidationError(
                f"PAV SV-only intermediate already exists: '{owned}'."
            )
        if owned.exists():
            owned.unlink()

    total = 0
    selected = 0
    has_chrom_header = False
    try:
        with gzip.open(raw_vcf, "rt", encoding="utf-8", newline="") as source:
            with temporary.open("x", encoding="utf-8", newline="\n") as output:
                for raw_line in source:
                    line = raw_line.rstrip("\r\n")
                    if line.startswith("##"):
                        output.write(line + "\n")
                        continue
                    if line.startswith("#CHROM\t"):
                        output.write(
                            "##hifivar_pav_sv_selection="
                            f"{_PAV_SV_SELECTION_POLICY}\n"
                        )
                        output.write(line + "\n")
                        has_chrom_header = True
                        continue
                    if not line:
                        continue
                    if not has_chrom_header:
                        raise OutputValidationError(
                            "PAV VCF contains records before the #CHROM header."
                        )
                    total += 1
                    if _is_pav_structural_variant(line):
                        output.write(line + "\n")
                        selected += 1
        if not has_chrom_header:
            raise OutputValidationError("PAV VCF lacks a #CHROM header.")
        validate_output_file(temporary)
        os.replace(temporary, destination)
    except OutputValidationError:
        if temporary.exists():
            temporary.unlink()
        raise
    except (OSError, EOFError, UnicodeError) as error:
        if temporary.exists():
            temporary.unlink()
        raise OutputValidationError(
            f"Unable to derive PAV SV-only VCF from '{raw_vcf}': {error}"
        ) from error

    return PavSvSelection(
        raw_vcf,
        destination,
        _PAV_SV_SELECTION_POLICY,
        total,
        selected,
    )


def _is_pav_structural_variant(record: str) -> bool:
    fields = record.split("\t")
    if len(fields) < 8:
        raise OutputValidationError("Malformed PAV VCF record.")
    info: dict[str, str | None] = {}
    for item in fields[7].split(";"):
        if not item:
            continue
        key, separator, value = item.partition("=")
        info[key] = value if separator else None
    svtype = info.get("SVTYPE")
    if svtype == "SNV":
        return False
    if svtype == "INV":
        return True
    if svtype not in {"INS", "DEL"}:
        raise OutputValidationError(
            f"Unsupported PAV SVTYPE in mixed VCF: '{svtype}'."
        )
    svlen = info.get("SVLEN")
    if svlen is None or "," in svlen:
        raise OutputValidationError(
            f"PAV {svtype} record '{fields[2]}' lacks one scalar SVLEN."
        )
    try:
        length = abs(int(svlen))
    except ValueError as error:
        raise OutputValidationError(
            f"PAV {svtype} record '{fields[2]}' has invalid SVLEN '{svlen}'."
        ) from error
    return length >= _PAV_SV_MIN_LENGTH


def _step_log(path: Path | None, step: str) -> Path | None:
    if path is None:
        return None
    return path.with_name(f"{path.stem}.{step}{path.suffix}")


__all__ = [
    "PavCommandPlan",
    "PavResult",
    "PavResultStatus",
    "PavSvSelection",
    "PavWrapper",
]
