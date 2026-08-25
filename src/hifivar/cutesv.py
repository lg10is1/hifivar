"""cuteSV single-sample structural-variant wrapper."""

from __future__ import annotations

import re
import shutil
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


_LOGGER = get_logger(__name__)
_VERSION_PATTERN = re.compile(r"cuteSV(?:\s+version)?\s+v?(\d+(?:\.\d+)+)", re.I)


@dataclass(frozen=True, slots=True)
class CuteSvResources:
    threads: int = 8
    memory_mb: int = 16_000
    runtime_minutes: int = 720

    def __post_init__(self) -> None:
        for name in ("threads", "memory_mb", "runtime_minutes"):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
                raise InputValidationError(f"cuteSV {name} must be a positive integer.")

    def to_dict(self) -> dict[str, int]:
        return {name: getattr(self, name) for name in ("threads", "memory_mb", "runtime_minutes")}


@dataclass(frozen=True, slots=True)
class CuteSvRequest:
    artifact: AlignmentArtifact
    raw_vcf: Path
    work_directory: Path
    resources: CuteSvResources = CuteSvResources()
    minimum_support: int = 10
    minimum_sv_size: int = 30
    max_cluster_bias_ins: int = 1_000
    diff_ratio_merging_ins: float = 0.9
    max_cluster_bias_del: int = 1_000
    diff_ratio_merging_del: float = 0.5
    genotype: bool = True
    overwrite: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.artifact, AlignmentArtifact):
            raise InputValidationError("cuteSV artifact must be AlignmentArtifact.")
        if self.artifact.output_format is not AlignmentOutputFormat.BAM:
            raise InputValidationError("HiFiVar Phase 4 cuteSV execution requires BAM input.")
        if not isinstance(self.resources, CuteSvResources):
            raise InputValidationError("cuteSV resources must be CuteSvResources.")
        for name in ("minimum_support", "minimum_sv_size", "max_cluster_bias_ins", "max_cluster_bias_del"):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
                raise InputValidationError(f"cuteSV {name} must be a positive integer.")
        for name in ("diff_ratio_merging_ins", "diff_ratio_merging_del"):
            value = getattr(self, name)
            if not isinstance(value, (int, float)) or isinstance(value, bool) or not 0 <= float(value) <= 1:
                raise InputValidationError(f"cuteSV {name} must be between 0 and 1.")
        if not isinstance(self.genotype, bool):
            raise InputValidationError("cuteSV genotype must be boolean.")
        if not isinstance(self.overwrite, bool):
            raise InputValidationError("cuteSV overwrite must be boolean.")
        raw_vcf = _coerce_path(self.raw_vcf, "cuteSV raw VCF")
        work = _coerce_path(self.work_directory, "cuteSV work directory")
        if not str(raw_vcf).lower().endswith(".cutesv.raw.vcf"):
            raise InputValidationError("cuteSV native output must end with '.cutesv.raw.vcf'.")
        object.__setattr__(self, "raw_vcf", raw_vcf)
        object.__setattr__(self, "work_directory", work)

    @property
    def sample_id(self) -> str:
        return self.artifact.sample_id

    @classmethod
    def create(cls, artifact: AlignmentArtifact, output_directory: str | Path, work_directory: str | Path, **kwargs: object) -> "CuteSvRequest":
        output_root = _coerce_path(output_directory, "cuteSV output directory")
        work_root = _coerce_path(work_directory, "cuteSV work root")
        sample = artifact.sample_id if isinstance(artifact, AlignmentArtifact) else "sample"
        return cls(artifact, output_root / f"{sample}.cutesv.raw.vcf", work_root / sample, **kwargs)  # type: ignore[arg-type]

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "alignment": self.artifact.to_dict(),
            "raw_vcf": str(self.raw_vcf),
            "work_directory": str(self.work_directory),
            "resources": self.resources.to_dict(),
            "minimum_support": self.minimum_support,
            "minimum_sv_size": self.minimum_sv_size,
            "max_cluster_bias_ins": self.max_cluster_bias_ins,
            "diff_ratio_merging_ins": self.diff_ratio_merging_ins,
            "max_cluster_bias_del": self.max_cluster_bias_del,
            "diff_ratio_merging_del": self.diff_ratio_merging_del,
            "genotype": self.genotype,
            "overwrite": self.overwrite,
        }


@dataclass(frozen=True, slots=True)
class CuteSvCommandPlan:
    args: tuple[str, ...]
    display_command: str
    shell: bool = False

    def to_dict(self) -> dict[str, object]:
        return {"args": list(self.args), "display_command": self.display_command, "shell": self.shell}


class CuteSvResultStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class CuteSvResult:
    request: CuteSvRequest
    status: CuteSvResultStatus
    command: CuteSvCommandPlan
    tool_version: str | None = None
    duration_seconds: float = 0.0

    def to_dict(self) -> dict[str, object]:
        return {
            "request": self.request.to_dict(),
            "status": self.status.value,
            "executed": self.status is CuteSvResultStatus.COMPLETED,
            "command": self.command.to_dict(),
            "tool_version": self.tool_version,
            "duration_seconds": self.duration_seconds,
        }


class CuteSvWrapper:
    def __init__(self, *, runner: CommandRunner | None = None, executable: str = "cuteSV") -> None:
        if not isinstance(executable, str) or not executable.strip():
            raise InputValidationError("cuteSV executable must be non-empty.")
        self.runner = runner or CommandRunner()
        self.executable = executable

    def build_command(self, request: CuteSvRequest) -> list[str]:
        self._validate_contract(request)
        command = [
            self.executable,
            str(request.artifact.path.absolute()),
            str(request.artifact.reference.fasta.absolute()),
            str(request.raw_vcf.absolute()),
            str(request.work_directory.absolute()),
            "--threads", str(request.resources.threads),
            "--sample", request.sample_id,
            "--min_support", str(request.minimum_support),
            "--min_size", str(request.minimum_sv_size),
            "--max_cluster_bias_INS", str(request.max_cluster_bias_ins),
            "--diff_ratio_merging_INS", str(request.diff_ratio_merging_ins),
            "--max_cluster_bias_DEL", str(request.max_cluster_bias_del),
            "--diff_ratio_merging_DEL", str(request.diff_ratio_merging_del),
        ]
        if request.genotype:
            command.append("--genotype")
        return command

    def plan_command(self, request: CuteSvRequest, *, redact_values: Collection[str] | None = None) -> CuteSvCommandPlan:
        command = self.build_command(request)
        return CuteSvCommandPlan(tuple(command), format_command(command, redact_values=redact_values))

    def detect_version(self) -> str:
        self.runner.require_executable(self.executable)
        result = self.runner.run([self.executable, "--version"])
        output = "\n".join(part for part in (result.stdout, result.stderr) if part)
        match = _VERSION_PATTERN.search(output)
        if match is None:
            raise ToolVersionError(f"Unable to parse cuteSV version from {output.strip()!r}.")
        return match.group(1)

    def run(self, request: CuteSvRequest, *, dry_run: bool = False, timeout: float | None = None, redact_values: Collection[str] | None = None, stderr_path: str | Path | None = None) -> CuteSvResult:
        self._validate_inputs(request)
        command = self.plan_command(request, redact_values=redact_values)
        if dry_run:
            self.runner.run(command.args, dry_run=True, timeout=timeout, redact_values=redact_values, stderr_path=stderr_path)
            return CuteSvResult(request, CuteSvResultStatus.PLANNED, command)
        version = self.detect_version()
        self._prepare_outputs(request)
        result = self.runner.run(command.args, timeout=timeout, redact_values=redact_values, stderr_path=stderr_path)
        validation.validate_output_file(request.raw_vcf)
        _LOGGER.info("cuteSV completed sample=%s version=%s runtime=%.3fs", request.sample_id, version, result.duration_seconds)
        return CuteSvResult(request, CuteSvResultStatus.COMPLETED, command, version, result.duration_seconds)

    @staticmethod
    def _validate_contract(request: CuteSvRequest) -> None:
        if not isinstance(request, CuteSvRequest):
            raise InputValidationError("cuteSV wrapper requires CuteSvRequest.")

    def _validate_inputs(self, request: CuteSvRequest) -> None:
        self._validate_contract(request)
        validation.validate_fasta(request.artifact.reference.fasta, require_fai=True)
        validate_alignment_artifact(request.artifact, require_index=True)

    @staticmethod
    def _prepare_outputs(request: CuteSvRequest) -> None:
        if request.raw_vcf.exists() and request.raw_vcf.is_dir():
            raise OutputValidationError(f"cuteSV output is a directory: '{request.raw_vcf}'.")
        if request.raw_vcf.exists() and not request.overwrite:
            raise OutputValidationError(f"cuteSV output already exists: '{request.raw_vcf}'.")
        if request.work_directory.exists() and not request.overwrite:
            raise OutputValidationError(f"cuteSV work directory already exists: '{request.work_directory}'.")
        marker = _ownership_marker(request.work_directory)
        try:
            expected = _ownership_marker_content(request.work_directory)
            if marker.exists() and (
                not marker.is_file() or marker.read_text(encoding="utf-8") != expected
            ):
                raise OutputValidationError(
                    f"Refusing to replace an unrecognized cuteSV ownership marker: '{marker}'."
                )
            if request.work_directory.exists():
                if request.work_directory.is_symlink() or not request.work_directory.is_dir():
                    raise OutputValidationError(
                        f"cuteSV work path is not a removable owned directory: '{request.work_directory}'."
                    )
                if not marker.is_file() or marker.read_text(encoding="utf-8") != expected:
                    raise OutputValidationError(
                        "Refusing to clean a cuteSV work directory not marked as HiFiVar-owned: "
                        f"'{request.work_directory}'."
                    )
                shutil.rmtree(request.work_directory)
            if request.raw_vcf.exists():
                request.raw_vcf.unlink()
            request.raw_vcf.parent.mkdir(parents=True, exist_ok=True)
            request.work_directory.mkdir(parents=True, exist_ok=False)
            marker.write_text(
                _ownership_marker_content(request.work_directory), encoding="utf-8"
            )
        except OutputValidationError:
            raise
        except OSError as error:
            raise OutputValidationError(
                f"Unable to prepare cuteSV outputs for '{request.sample_id}': {error}"
            ) from error


def _coerce_path(value: str | Path, label: str) -> Path:
    if not isinstance(value, (str, Path)) or (isinstance(value, str) and not value.strip()):
        raise InputValidationError(f"{label} must be a non-empty string or Path.")
    return Path(value).expanduser()


def _ownership_marker(work_directory: Path) -> Path:
    return work_directory.parent / f".{work_directory.name}.hifivar-cutesv-owned"


def _ownership_marker_content(work_directory: Path) -> str:
    return f"hifivar-cutesv-workdir-v1\n{work_directory.absolute()}\n"


__all__ = ["CuteSvCommandPlan", "CuteSvRequest", "CuteSvResources", "CuteSvResult", "CuteSvResultStatus", "CuteSvWrapper"]
