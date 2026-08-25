"""Common structural-variant artifacts and lightweight VCF finalization."""

from __future__ import annotations

import gzip
import re
from collections.abc import Collection, Sequence
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hifivar.command import CommandRunner, format_command
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar.reference import ReferenceGenome
from hifivar.validation import validate_output_file


_CONTIG_PATTERN = re.compile(r"^##contig=<ID=([^,>]+)")
_VERSION_PATTERN = re.compile(r"(\d+(?:\.\d+)+)")


class SvCaller(str, Enum):
    SAWFISH = "sawfish"
    SNIFFLES2 = "sniffles2"
    PBSV = "pbsv"
    CUTESV = "cutesv"


@dataclass(frozen=True, slots=True)
class StructuralVariantArtifact:
    """One independently validated caller VCF; not a harmonized callset."""

    caller: SvCaller
    sample_id: str
    reference_fasta: Path
    reference_build: str | None
    vcf_path: Path
    index_path: Path
    caller_version: str | None
    commands: tuple[tuple[str, ...], ...]
    reference_compatibility: str = "header_contigs_subset_of_declared_reference"
    harmonized: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.caller, SvCaller):
            raise InputValidationError("SV artifact caller must be SvCaller.")
        if not isinstance(self.sample_id, str) or not self.sample_id:
            raise InputValidationError("SV artifact sample_id must be non-empty.")
        expected = f".{self.caller.value}.sv.vcf.gz"
        if not str(self.vcf_path).lower().endswith(expected):
            raise InputValidationError(f"{self.caller.value} artifact must end with '{expected}'.")
        if self.index_path != Path(f"{self.vcf_path}.tbi"):
            raise InputValidationError("SV artifact index must be '<vcf>.tbi'.")
        if self.harmonized:
            raise InputValidationError("Raw caller artifacts cannot be marked harmonized in Phase 4.")

    def to_dict(self) -> dict[str, object]:
        return {
            "caller": self.caller.value,
            "sample_id": self.sample_id,
            "reference_fasta": str(self.reference_fasta),
            "reference_build": self.reference_build,
            "vcf_path": str(self.vcf_path),
            "index_path": str(self.index_path),
            "caller_version": self.caller_version,
            "commands": [list(command) for command in self.commands],
            "reference_compatibility": self.reference_compatibility,
            "harmonized": self.harmonized,
        }


def validate_sv_vcf(
    vcf_path: str | Path,
    index_path: str | Path,
    *,
    sample_id: str,
    reference: ReferenceGenome,
) -> None:
    """Validate one indexed SV VCF without interpreting or normalizing records."""
    vcf = Path(vcf_path)
    index = Path(index_path)
    validate_output_file(vcf)
    validate_output_file(index)
    _validate_bgzf(vcf, "SV VCF")
    observed_sample, observed_contigs, has_svtype = _read_vcf_header(vcf)
    if observed_sample != sample_id:
        raise OutputValidationError(
            f"SV VCF sample '{observed_sample}' does not match '{sample_id}'."
        )
    unknown = sorted(observed_contigs.difference({item.name for item in reference.contigs}))
    if unknown:
        raise OutputValidationError(f"SV VCF contains contigs absent from the reference: {unknown!r}.")
    if not has_svtype:
        raise OutputValidationError("SV VCF header lacks an INFO/SVTYPE definition.")
    _validate_tabix_index(index)


def validate_structural_variant_artifact(artifact: StructuralVariantArtifact, reference: ReferenceGenome) -> StructuralVariantArtifact:
    """Stream only VCF headers and validate BGZF, sample, contigs, and TBI."""
    if not isinstance(artifact, StructuralVariantArtifact):
        raise InputValidationError("SV validation requires StructuralVariantArtifact.")
    if not isinstance(reference, ReferenceGenome):
        raise InputValidationError("SV validation requires ReferenceGenome.")
    if artifact.reference_fasta.absolute() != reference.fasta.absolute() or artifact.reference_build != reference.build:
        raise OutputValidationError("SV artifact reference metadata differs from the declared reference.")
    validate_sv_vcf(
        artifact.vcf_path,
        artifact.index_path,
        sample_id=artifact.sample_id,
        reference=reference,
    )
    return artifact


def create_structural_variant_artifact(
    *,
    caller: SvCaller,
    sample_id: str,
    reference: ReferenceGenome,
    vcf_path: str | Path,
    caller_version: str | None,
    commands: Sequence[Sequence[str]],
) -> StructuralVariantArtifact:
    vcf = Path(vcf_path).expanduser()
    artifact = StructuralVariantArtifact(
        caller=caller,
        sample_id=sample_id,
        reference_fasta=reference.fasta,
        reference_build=reference.build,
        vcf_path=vcf,
        index_path=Path(f"{vcf}.tbi"),
        caller_version=caller_version,
        commands=tuple(tuple(str(arg) for arg in command) for command in commands),
    )
    return validate_structural_variant_artifact(artifact, reference)


@dataclass(frozen=True, slots=True)
class VcfFinalizeRequest:
    caller: SvCaller
    sample_id: str
    reference: ReferenceGenome
    source_vcf: Path
    output_vcf: Path
    caller_version: str | None = None
    caller_commands: tuple[tuple[str, ...], ...] = ()
    overwrite: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.caller, SvCaller) or not isinstance(self.reference, ReferenceGenome):
            raise InputValidationError("VCF finalization requires SvCaller and ReferenceGenome.")
        if not isinstance(self.sample_id, str) or not self.sample_id:
            raise InputValidationError("VCF finalization sample_id must be non-empty.")
        source = _coerce_path(self.source_vcf, "source VCF")
        output = _coerce_path(self.output_vcf, "final SV VCF")
        if source.suffix.lower() != ".vcf":
            raise InputValidationError("VCF finalization source must be a plain '.vcf' file.")
        expected = f".{self.caller.value}.sv.vcf.gz"
        if not str(output).lower().endswith(expected):
            raise InputValidationError(f"Final SV VCF must end with '{expected}'.")
        if source.absolute() == output.absolute():
            raise InputValidationError("Source and final SV VCF paths must differ.")
        if not isinstance(self.overwrite, bool):
            raise InputValidationError("VCF finalization overwrite must be boolean.")
        object.__setattr__(self, "source_vcf", source)
        object.__setattr__(self, "output_vcf", output)

    @property
    def output_index(self) -> Path:
        return Path(f"{self.output_vcf}.tbi")


@dataclass(frozen=True, slots=True)
class VcfFinalizeCommandPlan:
    args: tuple[str, ...]
    display_command: str
    step: str
    shell: bool = False

    def to_dict(self) -> dict[str, object]:
        return {"step": self.step, "args": list(self.args), "display_command": self.display_command, "shell": self.shell}


class VcfFinalizeStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class VcfFinalizeResult:
    request: VcfFinalizeRequest
    status: VcfFinalizeStatus
    commands: tuple[VcfFinalizeCommandPlan, VcfFinalizeCommandPlan]
    tool_versions: dict[str, str] | None = None
    artifact: StructuralVariantArtifact | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "status": self.status.value,
            "commands": [command.to_dict() for command in self.commands],
            "tool_versions": self.tool_versions,
            "artifact": self.artifact.to_dict() if self.artifact else None,
        }


class BgzipTabixWrapper:
    """Finalize a native plain caller VCF without altering its VCF records."""

    def __init__(self, *, runner: CommandRunner | None = None, bgzip_executable: str = "bgzip", tabix_executable: str = "tabix") -> None:
        self.runner = runner or CommandRunner()
        self.bgzip_executable = bgzip_executable
        self.tabix_executable = tabix_executable

    def build_commands(self, request: VcfFinalizeRequest) -> tuple[list[str], list[str]]:
        if not isinstance(request, VcfFinalizeRequest):
            raise InputValidationError("BGZF finalizer requires VcfFinalizeRequest.")
        return (
            [self.bgzip_executable, "-c", str(request.source_vcf.absolute())],
            [self.tabix_executable, "-p", "vcf", str(request.output_vcf.absolute())],
        )

    def plan_commands(self, request: VcfFinalizeRequest, *, redact_values: Collection[str] | None = None) -> tuple[VcfFinalizeCommandPlan, VcfFinalizeCommandPlan]:
        bgzip, tabix = self.build_commands(request)
        return (
            VcfFinalizeCommandPlan(tuple(bgzip), format_command(bgzip, redact_values=redact_values), "bgzip"),
            VcfFinalizeCommandPlan(tuple(tabix), format_command(tabix, redact_values=redact_values), "tabix"),
        )

    def detect_versions(self) -> dict[str, str]:
        versions: dict[str, str] = {}
        for name, executable in (("bgzip", self.bgzip_executable), ("tabix", self.tabix_executable)):
            self.runner.require_executable(executable)
            result = self.runner.run([executable, "--version"])
            output = "\n".join(part for part in (result.stdout, result.stderr) if part)
            match = _VERSION_PATTERN.search(output)
            if match is None:
                raise ToolVersionError(f"Unable to parse {name} version from {output.strip()!r}.")
            versions[name] = match.group(1)
        return versions

    def run(self, request: VcfFinalizeRequest, *, dry_run: bool = False, timeout: float | None = None, redact_values: Collection[str] | None = None, stderr_path: str | Path | None = None) -> VcfFinalizeResult:
        if not dry_run and (not request.source_vcf.exists() or not request.source_vcf.is_file()):
            raise InputValidationError(f"Native caller VCF is missing: '{request.source_vcf}'.")
        commands = self.plan_commands(request, redact_values=redact_values)
        if dry_run:
            self.runner.run(commands[0].args, dry_run=True, timeout=timeout, redact_values=redact_values, stdout_path=request.output_vcf, stderr_path=_step_log(stderr_path, "bgzip"))
            self.runner.run(commands[1].args, dry_run=True, timeout=timeout, redact_values=redact_values, stderr_path=_step_log(stderr_path, "tabix"))
            return VcfFinalizeResult(request, VcfFinalizeStatus.PLANNED, commands)
        versions = self.detect_versions()
        self._prepare_outputs(request)
        self.runner.run(commands[0].args, timeout=timeout, redact_values=redact_values, stdout_path=request.output_vcf, stderr_path=_step_log(stderr_path, "bgzip"))
        validate_output_file(request.output_vcf)
        self.runner.run(commands[1].args, timeout=timeout, redact_values=redact_values, stderr_path=_step_log(stderr_path, "tabix"))
        artifact = create_structural_variant_artifact(
            caller=request.caller,
            sample_id=request.sample_id,
            reference=request.reference,
            vcf_path=request.output_vcf,
            caller_version=request.caller_version,
            commands=(*request.caller_commands, *(command.args for command in commands)),
        )
        return VcfFinalizeResult(request, VcfFinalizeStatus.COMPLETED, commands, versions, artifact)

    @staticmethod
    def _prepare_outputs(request: VcfFinalizeRequest) -> None:
        for output in (request.output_vcf, request.output_index):
            if output.exists() and output.is_dir():
                raise OutputValidationError(f"Final SV output is a directory: '{output}'.")
            if output.exists() and not request.overwrite:
                raise OutputValidationError(f"Final SV output already exists: '{output}'.")
            if output.exists():
                output.unlink()
        request.output_vcf.parent.mkdir(parents=True, exist_ok=True)


def _read_vcf_header(path: Path) -> tuple[str, set[str], bool]:
    contigs: set[str] = set()
    has_fileformat = False
    has_svtype = False
    try:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\r\n")
                if line.startswith("##fileformat=VCFv"):
                    has_fileformat = True
                elif line.startswith("##INFO=<ID=SVTYPE,"):
                    has_svtype = True
                match = _CONTIG_PATTERN.match(line)
                if match:
                    contigs.add(match.group(1))
                if line.startswith("#CHROM\t"):
                    columns = line.split("\t")
                    if len(columns) != 10:
                        raise OutputValidationError("SV VCF must contain exactly one sample column.")
                    if not has_fileformat:
                        raise OutputValidationError("SV VCF header lacks ##fileformat.")
                    return columns[9], contigs, has_svtype
                if line and not line.startswith("#"):
                    break
    except OutputValidationError:
        raise
    except (OSError, EOFError, UnicodeError) as error:
        raise OutputValidationError(f"Unable to stream SV VCF header '{path}': {error}") from error
    raise OutputValidationError("SV VCF header lacks a valid #CHROM line.")


def _validate_bgzf(path: Path, label: str) -> None:
    try:
        with path.open("rb") as handle:
            fixed = handle.read(12)
            if len(fixed) != 12 or fixed[:4] != b"\x1f\x8b\x08\x04":
                raise OutputValidationError(f"{label} is not BGZF-compressed: '{path}'.")
            extra = handle.read(int.from_bytes(fixed[10:12], "little"))
    except OutputValidationError:
        raise
    except OSError as error:
        raise OutputValidationError(f"Unable to read {label} '{path}': {error}") from error
    offset = 0
    while offset + 4 <= len(extra):
        subfield_id = extra[offset : offset + 2]
        length = int.from_bytes(extra[offset + 2 : offset + 4], "little")
        offset += 4
        if subfield_id == b"BC" and length == 2:
            return
        offset += length
    raise OutputValidationError(f"{label} is not BGZF-compressed: '{path}'.")


def _validate_tabix_index(path: Path) -> None:
    _validate_bgzf(path, "tabix index")
    try:
        with gzip.open(path, "rb") as handle:
            magic = handle.read(4)
    except (OSError, EOFError) as error:
        raise OutputValidationError(f"Unable to read tabix index '{path}': {error}") from error
    if magic != b"TBI\x01":
        raise OutputValidationError(f"Invalid tabix index header: '{path}'.")


def _coerce_path(value: str | Path, label: str) -> Path:
    if not isinstance(value, (str, Path)) or (isinstance(value, str) and not value.strip()):
        raise InputValidationError(f"{label} must be a non-empty string or Path.")
    return Path(value).expanduser()


def _step_log(path: str | Path | None, step: str) -> Path | None:
    if path is None:
        return None
    log = Path(path)
    return log.with_name(f"{log.stem}.{step}{log.suffix}")


__all__ = [
    "BgzipTabixWrapper",
    "StructuralVariantArtifact",
    "SvCaller",
    "VcfFinalizeCommandPlan",
    "VcfFinalizeRequest",
    "VcfFinalizeResult",
    "VcfFinalizeStatus",
    "create_structural_variant_artifact",
    "validate_sv_vcf",
    "validate_structural_variant_artifact",
]
