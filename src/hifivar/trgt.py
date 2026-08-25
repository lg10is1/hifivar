"""TRGT single-sample tandem-repeat wrapper and output finalization."""

from __future__ import annotations

import re
from collections.abc import Collection
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hifivar import validation
from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import AlignmentArtifact, validate_alignment_artifact
from hifivar.command import CommandRunner, format_command
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar.logging_utils import get_logger
from hifivar.tr import TandemRepeatArtifact, TandemRepeatCatalog, validate_tandem_repeat_outputs


_LOGGER = get_logger(__name__)
_TRGT_VERSION = re.compile(r"(?:^|\s)trgt(?:\s+version)?\s+v?(\d+(?:\.\d+)+)", re.I)
_BCFTOOLS_VERSION = re.compile(r"bcftools\s+(\d+(?:\.\d+)+)", re.I)
_SAMTOOLS_VERSION = re.compile(r"samtools\s+(\d+(?:\.\d+)+)", re.I)
_LOCUS_PROCESSING_ERROR = re.compile(
    r"\[ERROR\]\s*-\s*Locus processing:",
    re.I,
)


class TrgtPreset(str, Enum):
    WGS = "wgs"
    TARGETED = "targeted"


class TrgtResultStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class TrgtResources:
    threads: int = 8
    memory_mb: int = 16000
    runtime_minutes: int = 720

    def __post_init__(self) -> None:
        for label, value in (
            ("threads", self.threads),
            ("memory_mb", self.memory_mb),
            ("runtime_minutes", self.runtime_minutes),
        ):
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise InputValidationError(f"TRGT {label} must be a positive integer.")

    def to_dict(self) -> dict[str, int]:
        return {
            "threads": self.threads,
            "memory_mb": self.memory_mb,
            "runtime_minutes": self.runtime_minutes,
        }


@dataclass(frozen=True, slots=True)
class TrgtRequest:
    artifact: AlignmentArtifact
    catalog: TandemRepeatCatalog
    raw_output_prefix: Path
    final_vcf: Path
    final_spanning_bam: Path
    karyotype: str
    resources: TrgtResources = TrgtResources()
    preset: TrgtPreset = TrgtPreset.WGS
    overwrite: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.artifact, AlignmentArtifact):
            raise InputValidationError("TRGT request requires AlignmentArtifact.")
        if self.artifact.output_format is not AlignmentOutputFormat.BAM:
            raise InputValidationError("TRGT genotype requires indexed aligned BAM; CRAM is not accepted.")
        if not isinstance(self.catalog, TandemRepeatCatalog):
            raise InputValidationError("TRGT request requires TandemRepeatCatalog.")
        if not isinstance(self.resources, TrgtResources):
            raise InputValidationError("TRGT request resources must be TrgtResources.")
        if not isinstance(self.preset, TrgtPreset):
            raise InputValidationError("TRGT preset must be wgs or targeted.")
        if self.karyotype not in {"XX", "XY"}:
            raise InputValidationError("TRGT karyotype must be XX or XY.")
        for field_name in ("raw_output_prefix", "final_vcf", "final_spanning_bam"):
            path = Path(getattr(self, field_name)).expanduser()
            object.__setattr__(self, field_name, path)
        if self.final_vcf.name != f"{self.sample_id}.tr.vcf.gz":
            raise InputValidationError("TRGT final VCF must follow {sample}.tr.vcf.gz.")
        if self.final_spanning_bam.name != f"{self.sample_id}.tr.spanning.bam":
            raise InputValidationError("TRGT spanning BAM must follow {sample}.tr.spanning.bam.")

    @property
    def sample_id(self) -> str:
        return self.artifact.sample_id

    @property
    def raw_vcf(self) -> Path:
        return Path(f"{self.raw_output_prefix}.vcf.gz")

    @property
    def raw_spanning_bam(self) -> Path:
        return Path(f"{self.raw_output_prefix}.spanning.bam")

    @property
    def final_vcf_index(self) -> Path:
        return Path(f"{self.final_vcf}.tbi")

    @property
    def final_spanning_bam_index(self) -> Path:
        return Path(f"{self.final_spanning_bam}.bai")

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "alignment": self.artifact.to_dict(),
            "catalog": self.catalog.to_dict(),
            "raw_output_prefix": str(self.raw_output_prefix),
            "final_vcf": str(self.final_vcf),
            "final_spanning_bam": str(self.final_spanning_bam),
            "karyotype": self.karyotype,
            "preset": self.preset.value,
            "resources": self.resources.to_dict(),
            "overwrite": self.overwrite,
        }


@dataclass(frozen=True, slots=True)
class TrgtCommandPlan:
    args: tuple[str, ...]
    display_command: str

    def to_dict(self) -> dict[str, object]:
        return {"args": list(self.args), "display_command": self.display_command}


@dataclass(frozen=True, slots=True)
class TrgtResult:
    request: TrgtRequest
    status: TrgtResultStatus
    commands: tuple[TrgtCommandPlan, ...]
    tool_versions: dict[str, str] | None = None
    runtime_seconds: float | None = None
    artifact: TandemRepeatArtifact | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "request": self.request.to_dict(),
            "status": self.status.value,
            "commands": [command.to_dict() for command in self.commands],
            "tool_versions": self.tool_versions,
            "runtime_seconds": self.runtime_seconds,
            "artifact": self.artifact.to_dict() if self.artifact else None,
        }


class TrgtWrapper:
    """Run TRGT then sort/index both official unsorted outputs."""

    def __init__(
        self,
        *,
        executable: str = "trgt",
        bcftools_executable: str = "bcftools",
        samtools_executable: str = "samtools",
        runner: CommandRunner | None = None,
    ) -> None:
        self.executable = executable
        self.bcftools_executable = bcftools_executable
        self.samtools_executable = samtools_executable
        self.runner = runner or CommandRunner()

    def plan_commands(
        self,
        request: TrgtRequest,
        *,
        redact_values: Collection[str] | None = None,
    ) -> tuple[TrgtCommandPlan, ...]:
        self._validate_contract(request)
        threads = str(request.resources.threads)
        commands = (
            (
                self.executable,
                "genotype",
                "--genome",
                str(request.artifact.reference.fasta.absolute()),
                "--reads",
                str(request.artifact.path.absolute()),
                "--repeats",
                str(request.catalog.path.absolute()),
                "--output-prefix",
                str(request.raw_output_prefix.absolute()),
                "--sample-name",
                request.sample_id,
                "--karyotype",
                request.karyotype,
                "--threads",
                threads,
                "--preset",
                request.preset.value,
            ),
            (
                self.bcftools_executable,
                "sort",
                "-Oz",
                "-o",
                str(request.final_vcf.absolute()),
                str(request.raw_vcf.absolute()),
            ),
            (
                self.bcftools_executable,
                "index",
                "--tbi",
                "--threads",
                threads,
                str(request.final_vcf.absolute()),
            ),
            (
                self.samtools_executable,
                "sort",
                "-@",
                threads,
                "-o",
                str(request.final_spanning_bam.absolute()),
                str(request.raw_spanning_bam.absolute()),
            ),
            (
                self.samtools_executable,
                "index",
                "-@",
                threads,
                "-o",
                str(request.final_spanning_bam_index.absolute()),
                str(request.final_spanning_bam.absolute()),
            ),
        )
        return tuple(
            TrgtCommandPlan(command, format_command(command, redact_values=redact_values))
            for command in commands
        )

    def detect_versions(self) -> dict[str, str]:
        versions: dict[str, str] = {}
        for name, executable, pattern in (
            ("trgt", self.executable, _TRGT_VERSION),
            ("bcftools", self.bcftools_executable, _BCFTOOLS_VERSION),
            ("samtools", self.samtools_executable, _SAMTOOLS_VERSION),
        ):
            self.runner.require_executable(executable)
            result = self.runner.run([executable, "--version"])
            output = "\n".join(part for part in (result.stdout, result.stderr) if part)
            match = pattern.search(output)
            if match is None:
                raise ToolVersionError(f"Unable to parse {name} version from {output.strip()!r}.")
            versions[name] = match.group(1)
        return versions

    def run(
        self,
        request: TrgtRequest,
        *,
        dry_run: bool = False,
        timeout: float | None = None,
        redact_values: Collection[str] | None = None,
        stderr_path: str | Path | None = None,
    ) -> TrgtResult:
        self._validate_inputs(request)
        commands = self.plan_commands(request, redact_values=redact_values)
        if dry_run:
            for command in commands:
                self.runner.run(
                    command.args,
                    dry_run=True,
                    timeout=timeout,
                    redact_values=redact_values,
                    stderr_path=stderr_path,
                )
            return TrgtResult(request, TrgtResultStatus.PLANNED, commands)
        versions = self.detect_versions()
        self._prepare_outputs(request)
        duration = 0.0
        for index, command in enumerate(commands):
            command_stderr_path = _command_log_path(stderr_path, index)
            result = self.runner.run(
                command.args,
                timeout=timeout,
                redact_values=redact_values,
                stderr_path=command_stderr_path,
            )
            duration += result.duration_seconds
            if index == 0:
                locus_error_count = _count_locus_processing_errors(
                    result.stderr,
                    command_stderr_path,
                )
                if locus_error_count:
                    log_hint = (
                        f" Inspect '{command_stderr_path}'."
                        if command_stderr_path is not None
                        else ""
                    )
                    raise OutputValidationError(
                        f"TRGT reported {locus_error_count} locus-processing "
                        f"error(s) for sample '{request.sample_id}' despite exit "
                        f"code 0; refusing incomplete tandem-repeat outputs."
                        f"{log_hint}"
                    )
                validation.validate_output_file(request.raw_vcf)
                validation.validate_output_file(request.raw_spanning_bam)
        artifact = validate_tandem_repeat_outputs(
            sample_id=request.sample_id,
            reference=request.artifact.reference,
            catalog=request.catalog,
            vcf_path=request.final_vcf,
            vcf_index_path=request.final_vcf_index,
            spanning_bam_path=request.final_spanning_bam,
            spanning_bam_index_path=request.final_spanning_bam_index,
            trgt_version=versions["trgt"],
            commands=tuple(command.args for command in commands),
        )
        _LOGGER.info(
            "TRGT completed sample=%s version=%s runtime=%.3fs",
            request.sample_id,
            versions["trgt"],
            duration,
        )
        return TrgtResult(
            request,
            TrgtResultStatus.COMPLETED,
            commands,
            versions,
            duration,
            artifact,
        )

    @staticmethod
    def _validate_contract(request: TrgtRequest) -> None:
        if not isinstance(request, TrgtRequest):
            raise InputValidationError("TRGT wrapper requires TrgtRequest.")

    def _validate_inputs(self, request: TrgtRequest) -> None:
        self._validate_contract(request)
        validation.validate_fasta(request.artifact.reference.fasta, require_fai=True)
        validate_alignment_artifact(request.artifact, require_index=True)
        request.catalog.validate(request.artifact.reference)

    @staticmethod
    def _prepare_outputs(request: TrgtRequest) -> None:
        outputs = (
            request.raw_vcf,
            request.raw_spanning_bam,
            request.final_vcf,
            request.final_vcf_index,
            request.final_spanning_bam,
            request.final_spanning_bam_index,
        )
        for output in outputs:
            if output.exists() and output.is_dir():
                raise OutputValidationError(f"TRGT output is a directory: '{output}'.")
            if output.exists() and not request.overwrite:
                raise OutputValidationError(f"TRGT output already exists: '{output}'.")
        for output in outputs:
            if output.exists():
                output.unlink()
        request.raw_output_prefix.parent.mkdir(parents=True, exist_ok=True)
        request.final_vcf.parent.mkdir(parents=True, exist_ok=True)
        request.final_spanning_bam.parent.mkdir(parents=True, exist_ok=True)


def _command_log_path(base: str | Path | None, command_index: int) -> Path | None:
    if base is None:
        return None
    path = Path(base)
    labels = ("trgt", "bcftools-sort", "bcftools-index", "samtools-sort", "samtools-index")
    return path.with_name(f"{path.stem}.{labels[command_index]}{path.suffix or '.log'}")


def _count_locus_processing_errors(
    captured_stderr: str | None,
    stderr_path: Path | None,
) -> int:
    """Count TRGT's zero-exit locus failures without scanning output VCF records."""
    if captured_stderr is not None:
        return sum(
            1
            for line in captured_stderr.splitlines()
            if _LOCUS_PROCESSING_ERROR.search(line)
        )
    if stderr_path is None or not stderr_path.is_file():
        return 0
    try:
        with stderr_path.open("r", encoding="utf-8", errors="replace") as handle:
            return sum(
                1 for line in handle if _LOCUS_PROCESSING_ERROR.search(line)
            )
    except OSError as error:
        raise OutputValidationError(
            f"Unable to inspect TRGT stderr log '{stderr_path}': {error}"
        ) from error


__all__ = [
    "TrgtCommandPlan",
    "TrgtPreset",
    "TrgtRequest",
    "TrgtResources",
    "TrgtResult",
    "TrgtResultStatus",
    "TrgtWrapper",
]
