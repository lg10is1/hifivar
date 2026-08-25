"""Phase 8 independent PAV and SVIM-asm orchestration."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from hifivar import __version__
from hifivar.assembly import AssemblyArtifact, AssemblyRole
from hifivar.assembly_sv import AssemblySvCaller, AssemblySvCollection, AssemblySvRequest, AssemblySvResources
from hifivar.context import AnalysisContext
from hifivar.exceptions import ConfigurationError, InputValidationError, ReferenceError
from hifivar.pav import PavResult, PavWrapper
from hifivar.serialization import standardize_data, utc_now_iso8601, write_json_atomic, write_yaml_atomic
from hifivar.svim_asm import SvimAsmResult, SvimAsmWrapper


PHASE8_REPORT_SCHEMA_VERSION: Final[str] = "1.0"


@dataclass(frozen=True, slots=True)
class Phase8Settings:
    overwrite: bool
    pav_enabled: bool
    pav_executable: str
    pav_snakefile: Path
    pav_version: str
    pav_resources: AssemblySvResources
    svim_enabled: bool
    svim_executable: str
    minimap2_executable: str
    samtools_executable: str
    bgzip_executable: str
    tabix_executable: str
    svim_resources: AssemblySvResources

    @classmethod
    def from_config(cls, config: Mapping[str, object]) -> "Phase8Settings":
        section = _mapping(config, "assembly_sv")
        pav = _mapping(section, "pav")
        svim = _mapping(section, "svim_asm")
        try:
            return cls(
                bool(section["overwrite"]), bool(pav["enabled"]), str(pav["executable"]),
                Path(str(pav["snakefile"])), str(pav["version"]),
                AssemblySvResources(int(pav["threads"]), int(pav["memory_mb"]), int(pav["runtime_minutes"])),
                bool(svim["enabled"]), str(svim["executable"]), str(svim["minimap2_executable"]),
                str(svim["samtools_executable"]), str(svim["bgzip_executable"]),
                str(svim["tabix_executable"]),
                AssemblySvResources(int(svim["threads"]), int(svim["memory_mb"]), int(svim["runtime_minutes"])),
            )
        except (KeyError, TypeError, ValueError) as error:
            raise ConfigurationError(f"Invalid Phase 8 settings: {error}") from error

    def to_dict(self) -> dict[str, object]:
        return {
            "overwrite": self.overwrite,
            "pav": {"enabled": self.pav_enabled, "executable": self.pav_executable,
                    "snakefile": str(self.pav_snakefile), "version": self.pav_version,
                    "resources": self.pav_resources.to_dict()},
            "svim_asm": {"enabled": self.svim_enabled, "executable": self.svim_executable,
                         "minimap2_executable": self.minimap2_executable,
                         "samtools_executable": self.samtools_executable,
                         "bgzip_executable": self.bgzip_executable,
                         "tabix_executable": self.tabix_executable,
                         "resources": self.svim_resources.to_dict()},
        }


@dataclass(frozen=True, slots=True)
class Phase8SampleResult:
    sample_id: str
    pav: PavResult | None
    svim_asm: SvimAsmResult | None
    collection: AssemblySvCollection

    def to_dict(self) -> dict[str, object]:
        return {"sample_id": self.sample_id, "pav": self.pav.to_dict() if self.pav else None,
                "svim_asm": self.svim_asm.to_dict() if self.svim_asm else None,
                "collection": self.collection.to_dict()}

@dataclass(frozen=True, slots=True)
class Phase8RunReport:
    context: AnalysisContext
    settings: Phase8Settings
    sample_results: tuple[Phase8SampleResult, ...]
    dry_run: bool
    created_at: str
    hifivar_version: str = __version__
    schema_version: str = PHASE8_REPORT_SCHEMA_VERSION

    def to_dict(self) -> dict[str, object]:
        payload = {
            "schema_version": self.schema_version, "hifivar_version": self.hifivar_version,
            "created_at": self.created_at, "status": "planned" if self.dry_run else "completed",
            "dry_run": self.dry_run, "context": self.context.to_dict(),
            "settings": self.settings.to_dict(),
            "sample_results": [item.to_dict() for item in self.sample_results],
        }
        value = standardize_data(payload, context="Phase 8 report value")
        if not isinstance(value, dict):
            raise InputValidationError("Phase 8 report serialization failed.")
        return value

    def write_json(self, path: Path, *, overwrite: bool = False) -> Path:
        return write_json_atomic(self.to_dict(), path, overwrite=overwrite, artifact_name="Phase 8 report")

    def write_yaml(self, path: Path, *, overwrite: bool = False) -> Path:
        return write_yaml_atomic(self.to_dict(), path, overwrite=overwrite, artifact_name="Phase 8 report")


def run_phase8(
    context: AnalysisContext, assemblies: Mapping[str, AssemblyArtifact], *,
    output_directory: str | Path, work_directory: str | Path,
    config: Mapping[str, object], dry_run: bool = False,
    pav_wrapper: PavWrapper | None = None, svim_wrapper: SvimAsmWrapper | None = None,
) -> Phase8RunReport:
    section = config.get("assembly_sv")
    if not isinstance(section, Mapping) or section.get("enabled") is not True:
        raise ConfigurationError("Phase 8 requires assembly_sv.enabled: true.")
    settings = Phase8Settings.from_config(config)
    if not settings.pav_enabled and not settings.svim_enabled:
        raise ConfigurationError("Phase 8 requires at least one enabled assembly-SV caller.")
    unknown = sorted(set(assemblies).difference(context.sample_ids))
    if unknown:
        raise InputValidationError(f"Phase 8 received assemblies for unknown samples: {unknown!r}.")
    output_root = Path(output_directory)
    work_root = Path(work_directory)
    pav_engine = pav_wrapper or PavWrapper(
        snakefile=settings.pav_snakefile, executable=settings.pav_executable,
        pav_version=settings.pav_version,
    )
    svim_engine = svim_wrapper or SvimAsmWrapper(
        executable=settings.svim_executable, minimap2_executable=settings.minimap2_executable,
        samtools_executable=settings.samtools_executable,
        bgzip_executable=settings.bgzip_executable, tabix_executable=settings.tabix_executable,
    )
    results: list[Phase8SampleResult] = []
    for record in context.samples:
        sample_id = record.sample.sample_id
        source = assemblies.get(sample_id)
        if source is None:
            raise InputValidationError(f"Phase 8 assembly is missing for sample '{sample_id}'.")
        if source.sample_id != sample_id:
            raise InputValidationError(f"Phase 8 assembly sample mismatch for '{sample_id}'.")
        haplotypes = tuple(item for item in source.assemblies if item.role in {AssemblyRole.HAPLOTYPE1, AssemblyRole.HAPLOTYPE2})
        pav_result = None
        svim_result = None
        artifacts = []
        if settings.pav_enabled:
            request = AssemblySvRequest(sample_id, AssemblySvCaller.PAV, context.reference, haplotypes,
                                        work_root / "pav" / sample_id,
                                        output_root / sample_id / f"{sample_id}.pav.assembly.sv.vcf.gz",
                                        settings.pav_resources, settings.overwrite)
            pav_result = pav_engine.run(request, dry_run=dry_run, stderr_path=work_root / "pav" / sample_id / f"{sample_id}.pav.log")
            if pav_result.artifact:
                artifacts.append(pav_result.artifact)
        if settings.svim_enabled:
            request = AssemblySvRequest(sample_id, AssemblySvCaller.SVIM_ASM, context.reference, haplotypes,
                                        work_root / "svim_asm" / sample_id,
                                        output_root / sample_id / f"{sample_id}.svim_asm.assembly.sv.vcf.gz",
                                        settings.svim_resources, settings.overwrite)
            svim_result = svim_engine.run(request, dry_run=dry_run, stderr_path=work_root / "svim_asm" / sample_id / f"{sample_id}.svim_asm.log")
            if svim_result.artifact:
                artifacts.append(svim_result.artifact)
        results.append(Phase8SampleResult(sample_id, pav_result, svim_result, AssemblySvCollection(sample_id, tuple(artifacts))))
    return Phase8RunReport(context, settings, tuple(results), dry_run, utc_now_iso8601())


def _mapping(mapping: Mapping[str, object], key: str) -> Mapping[str, object]:
    value = mapping.get(key)
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"Phase 8 requires configuration mapping {key}.")
    return value


__all__ = ["PHASE8_REPORT_SCHEMA_VERSION", "Phase8RunReport", "Phase8SampleResult", "Phase8Settings", "run_phase8"]
