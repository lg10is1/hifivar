"""SVIM-asm wrapper with explicit assembly-alignment and finalization steps."""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hifivar.assembly_sv import AssemblySvArtifact, AssemblySvCaller, AssemblySvRequest, create_assembly_sv_artifact
from hifivar.command import CommandRunner, format_command
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar.validation import validate_output_file


_VERSION = re.compile(r"(\d+(?:\.\d+)+)")


class SvimAsmResultStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class SvimAsmCommandPlan:
    step: str
    args: tuple[str, ...]
    stdout_path: Path | None = None

    def to_dict(self) -> dict[str, object]:
        return {"step": self.step, "args": list(self.args), "display_command": format_command(self.args), "stdout_path": str(self.stdout_path) if self.stdout_path else None, "shell": False}


@dataclass(frozen=True, slots=True)
class SvimAsmResult:
    request: AssemblySvRequest
    status: SvimAsmResultStatus
    commands: tuple[SvimAsmCommandPlan, ...]
    tool_versions: dict[str, str] | None = None
    runtime_seconds: float = 0.0
    artifact: AssemblySvArtifact | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "request": self.request.to_dict(), "status": self.status.value,
            "commands": [item.to_dict() for item in self.commands],
            "tool_versions": self.tool_versions, "runtime_seconds": self.runtime_seconds,
            "artifact": self.artifact.to_dict() if self.artifact else None,
        }


class SvimAsmWrapper:
    def __init__(self, *, executable: str = "svim-asm", minimap2_executable: str = "minimap2",
                 samtools_executable: str = "samtools", bgzip_executable: str = "bgzip",
                 tabix_executable: str = "tabix", runner: CommandRunner | None = None) -> None:
        self.executable = executable
        self.minimap2 = minimap2_executable
        self.samtools = samtools_executable
        self.bgzip = bgzip_executable
        self.tabix = tabix_executable
        self.runner = runner or CommandRunner()

    def plan_commands(self, request: AssemblySvRequest) -> tuple[SvimAsmCommandPlan, ...]:
        if request.caller is not AssemblySvCaller.SVIM_ASM:
            raise InputValidationError("SVIM-asm wrapper requires caller=svim_asm.")
        self._validate_inputs(request)
        native = request.work_directory / "native"
        plans: list[SvimAsmCommandPlan] = []
        bams: list[Path] = []
        for assembly in request.assemblies:
            stem = assembly.role.value
            sam = request.work_directory / f"{stem}.sam"
            bam = request.work_directory / f"{stem}.sorted.bam"
            bams.append(bam)
            plans.append(SvimAsmCommandPlan(f"minimap2_{stem}", (self.minimap2, "-a", "-x", "asm5", "--cs", "-r2k", "-t", str(request.resources.threads), str(request.reference.fasta.absolute()), str(assembly.path.absolute())), sam))
            plans.append(SvimAsmCommandPlan(f"samtools_sort_{stem}", (self.samtools, "sort", "-@", str(request.resources.threads), "-o", str(bam.absolute()), str(sam.absolute()))))
            plans.append(SvimAsmCommandPlan(f"samtools_index_{stem}", (self.samtools, "index", "-@", str(request.resources.threads), str(bam.absolute()))))
        mode = "diploid" if len(bams) == 2 else "haploid"
        svim = (self.executable, mode, "--sample", request.sample_id, str(native.absolute()), *[str(path.absolute()) for path in bams], str(request.reference.fasta.absolute()))
        plans.append(SvimAsmCommandPlan("svim_asm", svim))
        raw = native / "variants.vcf"
        plans.append(SvimAsmCommandPlan("bgzip", (self.bgzip, "-c", str(raw.absolute())), request.output_vcf))
        plans.append(SvimAsmCommandPlan("tabix", (self.tabix, "-p", "vcf", str(request.output_vcf.absolute()))))
        return tuple(plans)
    def detect_versions(self) -> dict[str, str]:
        versions: dict[str, str] = {}
        probes = {
            "svim_asm": (self.executable, "--version"),
            "minimap2": (self.minimap2, "--version"),
            "samtools": (self.samtools, "--version"),
            "bgzip": (self.bgzip, "--version"),
            "tabix": (self.tabix, "--version"),
        }
        for name, command in probes.items():
            self.runner.require_executable(command[0])
            result = self.runner.run(command)
            output = "\n".join(item for item in (result.stdout, result.stderr) if item)
            match = _VERSION.search(output)
            if match is None:
                raise ToolVersionError(f"Unable to parse {name} version.")
            versions[name] = match.group(1)
        return versions

    def run(self, request: AssemblySvRequest, *, dry_run: bool = False,
            stderr_path: Path | None = None) -> SvimAsmResult:
        commands = self.plan_commands(request)
        if dry_run:
            for item in commands:
                self.runner.run(item.args, dry_run=True, stdout_path=item.stdout_path)
            return SvimAsmResult(request, SvimAsmResultStatus.PLANNED, commands)
        versions = self.detect_versions()
        self._prepare(request)
        runtime = 0.0
        for item in commands:
            result = self.runner.run(
                item.args,
                stdout_path=item.stdout_path,
                stderr_path=_step_log(stderr_path, item.step),
            )
            runtime += result.duration_seconds
        raw = request.work_directory / "native" / "variants.vcf"
        validate_output_file(raw)
        intermediates: list[Path] = [raw]
        for assembly in request.assemblies:
            stem = assembly.role.value
            intermediates.extend((
                request.work_directory / f"{stem}.sam",
                request.work_directory / f"{stem}.sorted.bam",
                request.work_directory / f"{stem}.sorted.bam.bai",
            ))
        artifact = create_assembly_sv_artifact(
            request, raw_vcf=raw, intermediate_files=tuple(intermediates),
            caller_version=versions["svim_asm"], backend="native",
            commands=tuple(item.args for item in commands),
        )
        return SvimAsmResult(request, SvimAsmResultStatus.COMPLETED, commands, versions, runtime, artifact)

    @staticmethod
    def _validate_inputs(request: AssemblySvRequest) -> None:
        validate_output_file(request.reference.fasta)
        validate_output_file(request.reference.fai)
        for assembly in request.assemblies:
            validate_output_file(assembly.path)

    @staticmethod
    def _prepare(request: AssemblySvRequest) -> None:
        owned = [request.output_vcf, request.output_index, request.work_directory / "native" / "variants.vcf"]
        if any(path.exists() for path in owned) and not request.overwrite:
            raise OutputValidationError("SVIM-asm owned output already exists.")
        request.work_directory.mkdir(parents=True, exist_ok=True)
        (request.work_directory / "native").mkdir(parents=True, exist_ok=True)
        request.output_vcf.parent.mkdir(parents=True, exist_ok=True)


def _step_log(path: Path | None, step: str) -> Path | None:
    if path is None:
        return None
    return path.with_name(f"{path.stem}.{step}{path.suffix}")


__all__ = ["SvimAsmCommandPlan", "SvimAsmResult", "SvimAsmResultStatus", "SvimAsmWrapper"]
