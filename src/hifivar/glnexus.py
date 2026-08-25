"""GLnexus joint-genotyping wrapper for Phase 12 small-variant cohorts."""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hifivar.cohort import (
    CohortDefinition,
    CohortSampleInput,
    CohortTrack,
    CohortTrackResult,
    SampleCallState,
    scan_multisample_vcf,
    validate_track_inputs,
)
from hifivar.command import CommandRunner, format_command
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar.validation import validate_output_file


_VERSION = re.compile(r"(?:release\s+v?|glnexus_cli\s+v?)(\d+(?:\.\d+)+)", re.IGNORECASE)


@dataclass(frozen=True, slots=True)
class GLnexusResources:
    threads: int = 8
    memory_gb: int = 32

    def __post_init__(self) -> None:
        for name in ("threads", "memory_gb"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise InputValidationError(f"GLnexus {name} must be a positive integer.")


@dataclass(frozen=True, slots=True)
class GLnexusRequest:
    cohort: CohortDefinition
    inputs: tuple[CohortSampleInput, ...]
    work_directory: Path
    output_bcf: Path
    output_vcf: Path
    preset: str = "DeepVariantWGS"
    resources: GLnexusResources = GLnexusResources()
    overwrite: bool = False

    def __post_init__(self) -> None:
        inputs = tuple(self.inputs)
        validate_track_inputs(self.cohort, inputs)
        if any(not item.callable for item in inputs):
            states = {item.sample_id: item.state.value for item in inputs if not item.callable}
            raise InputValidationError(f"GLnexus requires a gVCF for every cohort sample; unavailable states: {states!r}.")
        if not isinstance(self.preset, str) or not self.preset.strip():
            raise InputValidationError("GLnexus preset must be non-empty.")
        object.__setattr__(self, "inputs", inputs)
        for name in ("work_directory", "output_bcf", "output_vcf"):
            object.__setattr__(self, name, Path(getattr(self, name)).expanduser())
        if self.output_bcf.name != f"{self.cohort.cohort_id}.small.bcf":
            raise InputValidationError("GLnexus BCF must follow {cohort}.small.bcf.")
        if self.output_vcf.name != f"{self.cohort.cohort_id}.small.vcf.gz":
            raise InputValidationError("GLnexus VCF must follow {cohort}.small.vcf.gz.")
        for path in self.expected_outputs:
            if path.exists() and not self.overwrite:
                raise OutputValidationError(f"GLnexus output already exists: '{path}'.")

    @property
    def bcf_index(self) -> Path:
        return Path(f"{self.output_bcf}.csi")

    @property
    def vcf_index(self) -> Path:
        return Path(f"{self.output_vcf}.tbi")

    @property
    def expected_outputs(self) -> tuple[Path, ...]:
        return self.output_bcf, self.bcf_index, self.output_vcf, self.vcf_index


@dataclass(frozen=True, slots=True)
class GLnexusCommandPlan:
    step: str
    args: tuple[str, ...]
    stdout_path: Path | None = None

    def to_dict(self) -> dict[str, object]:
        return {"step": self.step, "args": list(self.args), "display_command": format_command(self.args), "stdout_path": str(self.stdout_path) if self.stdout_path else None, "shell": False}


class GLnexusRunStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class GLnexusResult:
    request: GLnexusRequest
    status: GLnexusRunStatus
    commands: tuple[GLnexusCommandPlan, ...]
    versions: dict[str, str] | None = None
    runtime_seconds: float = 0.0
    qc: dict[str, object] | None = None

    def as_track_result(self) -> CohortTrackResult:
        state = SampleCallState.CALLED if self.status is GLnexusRunStatus.COMPLETED else SampleCallState.NOT_RUN
        return CohortTrackResult(
            CohortTrack.SMALL_VARIANTS,
            True,
            state,
            "glnexus",
            None if self.versions is None else self.versions.get("glnexus"),
            self.request.expected_outputs if self.status is GLnexusRunStatus.COMPLETED else (),
            tuple(item.args for item in self.commands),
            self.request.inputs,
            self.qc or {},
        )


class GLnexusWrapper:
    """Run GLnexus and bcftools as discrete shell-free commands."""

    def __init__(self, *, executable: str = "glnexus_cli", bcftools_executable: str = "bcftools", runner: CommandRunner | None = None) -> None:
        self.executable = executable
        self.bcftools = bcftools_executable
        self.runner = runner or CommandRunner()

    def plan_commands(self, request: GLnexusRequest) -> tuple[GLnexusCommandPlan, ...]:
        glnexus = (
            self.executable,
            "--dir", str(request.work_directory.absolute()),
            "--config", request.preset,
            "--threads", str(request.resources.threads),
            "--mem-gbytes", str(request.resources.memory_gb),
            *(str(item.source_path.absolute()) for item in request.inputs if item.source_path is not None),
        )
        force = ("--force",) if request.overwrite else ()
        return (
            GLnexusCommandPlan("glnexus", glnexus, request.output_bcf),
            GLnexusCommandPlan("bcf_index", (self.bcftools, "index", *force, "--threads", str(request.resources.threads), "--csi", str(request.output_bcf.absolute()))),
            GLnexusCommandPlan("vcf", (self.bcftools, "view", "-Oz", "-o", str(request.output_vcf.absolute()), str(request.output_bcf.absolute()))),
            GLnexusCommandPlan("vcf_index", (self.bcftools, "index", *force, "--threads", str(request.resources.threads), "--tbi", str(request.output_vcf.absolute()))),
        )

    def detect_versions(self) -> dict[str, str]:
        versions: dict[str, str] = {}
        for name, command in {
            "glnexus": (self.executable, "--help"),
            "bcftools": (self.bcftools, "--version"),
        }.items():
            self.runner.require_executable(command[0])
            result = self.runner.run(command)
            output = "\n".join(part for part in (result.stdout, result.stderr) if part)
            match = _VERSION.search(output) if name == "glnexus" else re.search(r"bcftools\s+(\d+(?:\.\d+)+)", output, re.IGNORECASE)
            if match is None:
                raise ToolVersionError(f"Unable to parse {name} version.")
            versions[name] = match.group(1)
        return versions

    def run(self, request: GLnexusRequest, *, dry_run: bool = False, log_path: Path | None = None) -> GLnexusResult:
        commands = self.plan_commands(request)
        if dry_run:
            for command in commands:
                self.runner.run(command.args, dry_run=True, stdout_path=command.stdout_path, stderr_path=_step_log(log_path, command.step))
            return GLnexusResult(request, GLnexusRunStatus.PLANNED, commands)

        _validate_gvcf_inputs(request)
        request.output_bcf.parent.mkdir(parents=True, exist_ok=True)
        request.output_vcf.parent.mkdir(parents=True, exist_ok=True)
        if request.work_directory.exists():
            raise OutputValidationError(
                f"GLnexus scratch directory already exists: '{request.work_directory}'. "
                "Use a new work directory or remove this HiFiVar-owned scratch directory after review."
            )
        else:
            request.work_directory.parent.mkdir(parents=True, exist_ok=True)
        versions = self.detect_versions()
        runtime = 0.0
        for command in commands:
            result = self.runner.run(command.args, stdout_path=command.stdout_path, stderr_path=_step_log(log_path, command.step))
            runtime += result.duration_seconds
            if command.step == "glnexus":
                validate_output_file(request.output_bcf)
        for output in request.expected_outputs:
            validate_output_file(output)
        qc = scan_multisample_vcf(request.output_vcf, request.cohort.sample_ids)
        return GLnexusResult(request, GLnexusRunStatus.COMPLETED, commands, versions, runtime, qc)


def _validate_gvcf_inputs(request: GLnexusRequest) -> None:
    reference_contigs = set(request.cohort.reference.contig_names)
    for item in request.inputs:
        assert item.source_path is not None
        try:
            validate_output_file(item.source_path)
        except OutputValidationError as error:
            raise InputValidationError(f"gVCF input is unavailable for '{item.sample_id}': {error}") from error
        if item.index_path is None:
            raise InputValidationError(f"gVCF index is missing for '{item.sample_id}'.")
        try:
            validate_output_file(item.index_path)
        except OutputValidationError as error:
            raise InputValidationError(f"gVCF index is unavailable for '{item.sample_id}': {error}") from error
        samples: tuple[str, ...] | None = None
        contigs: set[str] = set()
        try:
            with gzip.open(item.source_path, "rt", encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith("##contig=<ID="):
                        contigs.add(line.split("ID=", 1)[1].split(",", 1)[0].split(">", 1)[0])
                    elif line.startswith("#CHROM\t"):
                        samples = tuple(line.rstrip("\r\n").split("\t")[9:])
                        break
        except (OSError, EOFError, UnicodeError) as error:
            raise InputValidationError(f"Unable to read gVCF '{item.source_path}': {error}") from error
        if samples != (item.sample_id,):
            raise InputValidationError(f"gVCF sample mismatch for '{item.sample_id}': observed {samples!r}.")
        unexpected = sorted(contigs.difference(reference_contigs))
        if unexpected:
            raise InputValidationError(f"REFERENCE_CONTIG_MISMATCH in gVCF '{item.source_path}': {unexpected!r}.")


def _step_log(path: Path | None, step: str) -> Path | None:
    return None if path is None else path.with_name(f"{path.stem}.{step}{path.suffix}")


__all__ = [
    "GLnexusCommandPlan", "GLnexusRequest", "GLnexusResources", "GLnexusResult",
    "GLnexusRunStatus", "GLnexusWrapper",
]
