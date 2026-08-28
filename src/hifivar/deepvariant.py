"""DeepVariant wrapper and isolated native/container execution backends."""

from __future__ import annotations

import re
from collections.abc import Collection, Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path

from hifivar import validation
from hifivar.alignment_postprocess import validate_alignment_artifact
from hifivar.command import CommandRunner, format_command
from hifivar.exceptions import (
    ConfigurationError,
    InputValidationError,
    OutputValidationError,
    ToolVersionError,
)
from hifivar.logging_utils import get_logger
from hifivar.small import (
    DeepVariantRequest,
    SmallVariantCommandPlan,
    SmallVariantResult,
    SmallVariantResultStatus,
    validate_small_variant_outputs,
)


_LOGGER = get_logger(__name__)
_CONTAINER_EXECUTABLE = "/opt/deepvariant/bin/run_deepvariant"
_VERSION_PATTERN = re.compile(r"DeepVariant(?:\s+version)?\s+v?(\d+(?:\.\d+)+)", re.I)
_MINIMUM_FILE_DESCRIPTOR_LIMIT = 4_096
_RECOMMENDED_FILE_DESCRIPTOR_LIMIT = 65_536


class DeepVariantExecutionMode(str, Enum):
    """Supported DeepVariant launch strategies."""

    NATIVE = "native"
    DOCKER = "docker"
    APPTAINER = "apptainer"


@dataclass(frozen=True, slots=True)
class DeepVariantRuntime:
    """Keep executable/container details outside calling models and workflows."""

    mode: DeepVariantExecutionMode = DeepVariantExecutionMode.NATIVE
    executable: str = "run_deepvariant"
    image: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.mode, DeepVariantExecutionMode):
            raise ConfigurationError(
                "DeepVariant runtime mode must be native, docker, or apptainer."
            )
        if not isinstance(self.executable, str) or not self.executable.strip():
            raise ConfigurationError("DeepVariant executable must be non-empty.")
        if self.mode is DeepVariantExecutionMode.NATIVE:
            if self.image is not None:
                raise ConfigurationError("Native DeepVariant runtime cannot set an image.")
        elif not isinstance(self.image, str) or not self.image.strip():
            raise ConfigurationError(
                f"DeepVariant {self.mode.value} runtime requires an image."
            )

    @classmethod
    def from_config(cls, config: Mapping[str, object]) -> DeepVariantRuntime:
        section = config.get("small")
        if not isinstance(section, Mapping):
            return cls()
        try:
            mode = DeepVariantExecutionMode(str(section.get("execution_mode", "native")).lower())
        except ValueError as error:
            raise ConfigurationError(
                "small.execution_mode must be native, docker, or apptainer."
            ) from error
        executable = section.get("deepvariant_executable", "run_deepvariant")
        image = section.get("deepvariant_image")
        return cls(mode=mode, executable=executable, image=image)  # type: ignore[arg-type]

    @property
    def launcher(self) -> str:
        if self.mode is DeepVariantExecutionMode.NATIVE:
            return self.executable
        return self.mode.value

    def version_command(self) -> list[str]:
        if self.mode is DeepVariantExecutionMode.NATIVE:
            return [self.executable, "--version"]
        if self.mode is DeepVariantExecutionMode.DOCKER:
            return [
                "docker",
                "run",
                "--rm",
                str(self.image),
                _CONTAINER_EXECUTABLE,
                "--version",
            ]
        return [
            "apptainer",
            "exec",
            "--cleanenv",
            str(self.image),
            _CONTAINER_EXECUTABLE,
            "--version",
        ]

    def command_prefix(
        self,
        request: DeepVariantRequest,
        *,
        create_writable_mounts: bool = True,
    ) -> list[str]:
        if self.mode is DeepVariantExecutionMode.NATIVE:
            return [self.executable]
        mounts = _container_mounts(
            request,
            create_writable_sources=create_writable_mounts,
        )
        if self.mode is DeepVariantExecutionMode.DOCKER:
            prefix = ["docker", "run", "--rm"]
            for path, writable in mounts:
                mount = f"type=bind,source={path},target={path}"
                if not writable:
                    mount += ",readonly"
                prefix.extend(("--mount", mount))
            prefix.extend(("--env", f"TMPDIR={request.temporary_directory.absolute()}"))
            return prefix + [str(self.image), _CONTAINER_EXECUTABLE]
        prefix = ["apptainer", "exec", "--cleanenv"]
        for path, writable in mounts:
            binding = f"{path}:{path}" + ("" if writable else ":ro")
            prefix.extend(("--bind", binding))
        prefix.extend(("--env", f"TMPDIR={request.temporary_directory.absolute()}"))
        return prefix + [str(self.image), _CONTAINER_EXECUTABLE]

    def to_dict(self) -> dict[str, object]:
        return {
            "mode": self.mode.value,
            "launcher": self.launcher,
            "executable": self.executable,
            "image": self.image,
        }


class DeepVariantWrapper:
    """Build and execute deterministic PacBio DeepVariant commands."""

    def __init__(
        self,
        *,
        runner: CommandRunner | None = None,
        runtime: DeepVariantRuntime | None = None,
    ) -> None:
        self.runner = runner or CommandRunner()
        self.runtime = runtime or DeepVariantRuntime()

    def build_command(
        self,
        request: DeepVariantRequest,
        *,
        create_writable_mounts: bool = False,
    ) -> list[str]:
        self._validate_contract(request)
        return self.runtime.command_prefix(
            request,
            create_writable_mounts=create_writable_mounts,
        ) + [
            f"--model_type={request.model_type.value}",
            f"--ref={request.reference_fasta.absolute()}",
            f"--reads={request.alignment_path.absolute()}",
            f"--sample_name={request.sample_id}",
            f"--output_vcf={request.output_vcf.absolute()}",
            f"--output_gvcf={request.output_gvcf.absolute()}",
            f"--num_shards={request.resources.threads}",
            f"--intermediate_results_dir={request.intermediate_directory.absolute()}",
            f"--logging_dir={request.logging_directory.absolute()}",
            "--vcf_stats_report=false",
        ]

    def plan_command(
        self,
        request: DeepVariantRequest,
        *,
        redact_values: Collection[str] | None = None,
        create_writable_mounts: bool = False,
    ) -> SmallVariantCommandPlan:
        args = tuple(
            self.build_command(
                request,
                create_writable_mounts=create_writable_mounts,
            )
        )
        return SmallVariantCommandPlan(
            args=args,
            display_command=format_command(args, redact_values=redact_values),
        )

    def detect_version(self) -> str:
        self.runner.require_executable(self.runtime.launcher)
        result = self.runner.run(self.runtime.version_command())
        output = "\n".join(part for part in (result.stdout, result.stderr) if part)
        match = _VERSION_PATTERN.search(output)
        if match is None:
            raise ToolVersionError(
                f"Unable to parse DeepVariant version from {output.strip()!r}."
            )
        return match.group(1)

    def run(
        self,
        request: DeepVariantRequest,
        *,
        dry_run: bool = False,
        timeout: float | None = None,
        redact_values: Collection[str] | None = None,
        stderr_path: str | Path | None = None,
    ) -> SmallVariantResult:
        self._validate_inputs(request)
        runner_environment = (
            {"TMPDIR": str(request.temporary_directory.absolute())}
            if self.runtime.mode is DeepVariantExecutionMode.NATIVE
            else None
        )
        if dry_run:
            command = self.plan_command(request, redact_values=redact_values)
            self.runner.run(
                command.args,
                dry_run=True,
                env=runner_environment,
                timeout=timeout,
                redact_values=redact_values,
                stderr_path=stderr_path,
            )
            return SmallVariantResult(
                request=request,
                status=SmallVariantResultStatus.PLANNED,
                command=command,
            )

        if isinstance(self.runner, CommandRunner):
            _validate_file_descriptor_limit(request)
        version = self.detect_version()
        self._prepare_outputs(request)
        command = self.plan_command(
            request,
            redact_values=redact_values,
            create_writable_mounts=True,
        )
        result = self.runner.run(
            command.args,
            env=runner_environment,
            timeout=timeout,
            redact_values=redact_values,
            stderr_path=stderr_path,
        )
        try:
            artifact = validate_small_variant_outputs(request, tool_version=version)
        except OutputValidationError as error:
            quarantine = _quarantine_outputs(request, error)
            if quarantine is not None:
                raise OutputValidationError(
                    f"{error} DeepVariant outputs were moved to quarantine: "
                    f"'{quarantine}'."
                ) from error
            raise
        _LOGGER.info(
            "DeepVariant completed sample=%s version=%s runtime=%.3fs",
            request.sample_id,
            version,
            result.duration_seconds,
        )
        return SmallVariantResult(
            request=request,
            status=SmallVariantResultStatus.COMPLETED,
            command=command,
            tool_version=version,
            duration_seconds=result.duration_seconds,
            artifact=artifact,
        )

    def _validate_contract(self, request: DeepVariantRequest) -> None:
        if not isinstance(request, DeepVariantRequest):
            raise InputValidationError("DeepVariant wrapper requires DeepVariantRequest.")

    def _validate_inputs(self, request: DeepVariantRequest) -> None:
        self._validate_contract(request)
        validation.validate_fasta(request.reference_fasta, require_fai=True)
        validate_alignment_artifact(request.artifact, require_index=True)

    def _prepare_outputs(self, request: DeepVariantRequest) -> None:
        outputs = (
            request.output_vcf,
            request.output_gvcf,
            request.output_vcf_index,
            request.output_gvcf_index,
        )
        try:
            for output in outputs:
                if output.exists():
                    if output.is_dir():
                        raise OutputValidationError(
                            f"DeepVariant output is a directory: '{output}'."
                        )
                    if not request.overwrite:
                        raise OutputValidationError(
                            f"DeepVariant output already exists: '{output}'."
                        )
                    output.unlink()
            request.output_vcf.parent.mkdir(parents=True, exist_ok=True)
            request.output_gvcf.parent.mkdir(parents=True, exist_ok=True)
            request.intermediate_directory.mkdir(parents=True, exist_ok=True)
            request.logging_directory.mkdir(parents=True, exist_ok=True)
            request.temporary_directory.mkdir(parents=True, exist_ok=True)
        except OutputValidationError:
            raise
        except OSError as error:
            raise OutputValidationError(
                f"Unable to prepare DeepVariant outputs for sample "
                f"'{request.sample_id}': {error}"
            ) from error


def _container_mounts(
    request: DeepVariantRequest,
    *,
    create_writable_sources: bool,
) -> tuple[tuple[Path, bool], ...]:
    access: dict[str, tuple[Path, bool]] = {}
    readonly_paths = (
        request.reference_fasta.absolute().parent,
        request.alignment_path.absolute().parent,
    )
    writable_paths = {
        request.output_vcf.absolute().parent,
        request.output_gvcf.absolute().parent,
        request.intermediate_directory.absolute(),
        request.logging_directory.absolute(),
        request.temporary_directory.absolute(),
    }
    for path in readonly_paths:
        access[str(path)] = (path, False)
    for path in writable_paths:
        if create_writable_sources:
            try:
                path.mkdir(parents=True, exist_ok=True)
            except OSError as error:
                raise OutputValidationError(
                    f"Unable to create writable DeepVariant container mount "
                    f"source '{path}': {error}"
                ) from error
        access[str(path)] = (path, True)
    return tuple(access[key] for key in sorted(access))


def _get_file_descriptor_limits() -> tuple[int, int] | None:
    """Return POSIX RLIMIT_NOFILE or None on platforms without resource."""
    try:
        import resource
    except ImportError:  # pragma: no cover - exercised on Windows
        return None
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    return int(soft), int(hard)


def _validate_file_descriptor_limit(request: DeepVariantRequest) -> None:
    """Fail before DeepVariant when the POSIX open-file limit is unsafe."""
    limits = _get_file_descriptor_limits()
    if limits is None:
        return
    soft, hard = limits
    if soft < 0 or soft >= _MINIMUM_FILE_DESCRIPTOR_LIMIT:
        return
    raise ConfigurationError(
        f"DeepVariant sample '{request.sample_id}' cannot start safely: "
        f"open-file soft limit is {soft} (hard limit {hard}) for "
        f"{request.resources.threads} shard(s). HiFiVar requires at least "
        f"{_MINIMUM_FILE_DESCRIPTOR_LIMIT} and recommends "
        f"{_RECOMMENDED_FILE_DESCRIPTOR_LIMIT}. Raise it with "
        f"'ulimit -n <allowed-value>' before launch, or reduce small.threads "
        f"when the site hard limit cannot be raised."
    )


def _quarantine_outputs(
    request: DeepVariantRequest,
    validation_error: OutputValidationError,
) -> Path | None:
    """Move validation-failed products outside Snakemake-managed outputs."""
    outputs = (
        request.output_vcf,
        request.output_gvcf,
        request.output_vcf_index,
        request.output_gvcf_index,
    )
    existing = tuple(path for path in outputs if path.is_file())
    if not existing:
        return None
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    quarantine = (
        request.output_vcf.parent
        / "quarantine"
        / f"{request.sample_id}.{timestamp}"
    )
    try:
        quarantine.mkdir(parents=True, exist_ok=False)
        for output in existing:
            output.replace(quarantine / output.name)
        (quarantine / "VALIDATION_ERROR.txt").write_text(
            f"{validation_error}\n",
            encoding="utf-8",
        )
    except OSError as error:
        raise OutputValidationError(
            f"{validation_error} Unable to quarantine DeepVariant outputs "
            f"under '{quarantine}': {error}"
        ) from error
    _LOGGER.error(
        "DeepVariant validation failed sample=%s; outputs quarantined at %s",
        request.sample_id,
        quarantine,
    )
    return quarantine


__all__ = [
    "DeepVariantExecutionMode",
    "DeepVariantRuntime",
    "DeepVariantWrapper",
]
