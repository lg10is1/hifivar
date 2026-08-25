"""Phase 5 TRGT orchestration and provenance."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Final

from hifivar import __version__
from hifivar.alignment import AlignmentOutputFormat
from hifivar.alignment_postprocess import (
    AlignmentArtifact,
    AlignmentSortOrder,
    AlignmentSource,
    find_alignment_index,
    validate_alignment_artifact,
)
from hifivar.context import AnalysisContext
from hifivar.exceptions import ConfigurationError, InputValidationError, ReferenceError
from hifivar.sample import InputType
from hifivar.sample_sheet import Sex
from hifivar.serialization import standardize_data, utc_now_iso8601, write_json_atomic, write_yaml_atomic
from hifivar.tr import TandemRepeatCatalog
from hifivar.trgt import TrgtPreset, TrgtRequest, TrgtResources, TrgtResult, TrgtWrapper


PHASE5_REPORT_SCHEMA_VERSION: Final[str] = "1.0"


class Phase5RunStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class Phase5Settings:
    enabled: bool
    catalog: Path
    catalog_reference_build: str | None
    executable: str
    bcftools_executable: str
    samtools_executable: str
    resources: TrgtResources
    preset: TrgtPreset
    karyotype: str
    overwrite: bool

    @classmethod
    def from_config(cls, config: Mapping[str, object]) -> "Phase5Settings":
        section = config.get("tr")
        if not isinstance(section, Mapping):
            raise ConfigurationError("Phase 5 requires configuration mapping tr.")
        try:
            catalog = section["catalog"]
            if not isinstance(catalog, str) or not catalog.strip():
                raise ConfigurationError("tr.catalog must be a non-empty BED path when Phase 5 runs.")
            build = section["catalog_reference_build"]
            return cls(
                enabled=bool(section["enabled"]),
                catalog=Path(catalog).expanduser(),
                catalog_reference_build=None if build is None else str(build),
                executable=str(section["executable"]),
                bcftools_executable=str(section["bcftools_executable"]),
                samtools_executable=str(section["samtools_executable"]),
                resources=TrgtResources(
                    int(section["threads"]),
                    int(section["memory_mb"]),
                    int(section["runtime_minutes"]),
                ),
                preset=TrgtPreset(str(section["preset"]).lower()),
                karyotype=str(section["karyotype"]),
                overwrite=bool(section["overwrite"]),
            )
        except ConfigurationError:
            raise
        except (KeyError, TypeError, ValueError) as error:
            raise ConfigurationError(f"Invalid or incomplete Phase 5 tr configuration: {error}") from error

    def to_dict(self) -> dict[str, object]:
        return {
            "enabled": self.enabled,
            "catalog": str(self.catalog),
            "catalog_reference_build": self.catalog_reference_build,
            "executable": self.executable,
            "bcftools_executable": self.bcftools_executable,
            "samtools_executable": self.samtools_executable,
            "resources": self.resources.to_dict(),
            "preset": self.preset.value,
            "karyotype": self.karyotype,
            "overwrite": self.overwrite,
        }


@dataclass(frozen=True, slots=True)
class Phase5SampleResult:
    alignment: AlignmentArtifact
    karyotype: str
    trgt_result: TrgtResult

    @property
    def sample_id(self) -> str:
        return self.alignment.sample_id

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "alignment": self.alignment.to_dict(),
            "karyotype": self.karyotype,
            "trgt_result": self.trgt_result.to_dict(),
        }


@dataclass(frozen=True, slots=True)
class Phase5RunReport:
    context: AnalysisContext
    settings: Phase5Settings
    sample_results: tuple[Phase5SampleResult, ...]
    dry_run: bool
    schema_version: str = PHASE5_REPORT_SCHEMA_VERSION
    hifivar_version: str = __version__
    created_at: str = field(default_factory=utc_now_iso8601)
    status: Phase5RunStatus = field(init=False)

    def __post_init__(self) -> None:
        if tuple(result.sample_id for result in self.sample_results) != self.context.sample_ids:
            raise InputValidationError("Phase 5 result order differs from AnalysisContext sample order.")
        if not self.dry_run and any(result.trgt_result.artifact is None for result in self.sample_results):
            raise InputValidationError("Completed Phase 5 results require validated TR artifacts.")
        object.__setattr__(
            self,
            "status",
            Phase5RunStatus.PLANNED if self.dry_run else Phase5RunStatus.COMPLETED,
        )

    @property
    def tool_versions(self) -> dict[str, str | None]:
        versions: dict[str, str | None] = {"trgt": None, "bcftools": None, "samtools": None}
        for sample in self.sample_results:
            if sample.trgt_result.tool_versions:
                versions.update(sample.trgt_result.tool_versions)
        return versions

    def to_dict(self) -> dict[str, object]:
        payload = {
            "schema_version": self.schema_version,
            "hifivar_version": self.hifivar_version,
            "created_at": self.created_at,
            "status": self.status.value,
            "dry_run": self.dry_run,
            "settings": self.settings.to_dict(),
            "tool_versions": self.tool_versions,
            "analysis_context": self.context.to_dict(),
            "sample_results": [result.to_dict() for result in self.sample_results],
        }
        standardized = standardize_data(payload, context="Phase 5 report value")
        if not isinstance(standardized, dict):
            raise InputValidationError("Phase 5 report serialization failed.")
        return standardized

    def write_json(self, path: str | Path, *, overwrite: bool = False) -> Path:
        return write_json_atomic(self.to_dict(), path, overwrite=overwrite, artifact_name="Phase 5 report")

    def write_yaml(self, path: str | Path, *, overwrite: bool = False) -> Path:
        return write_yaml_atomic(self.to_dict(), path, overwrite=overwrite, artifact_name="Phase 5 report")


def run_phase5(
    context: AnalysisContext,
    output_directory: str | Path,
    *,
    alignment_artifacts: Mapping[str, AlignmentArtifact] | None = None,
    settings: Phase5Settings | None = None,
    wrapper: TrgtWrapper | None = None,
    dry_run: bool = False,
) -> Phase5RunReport:
    """Run one independent TRGT genotype workflow per sample."""
    if not isinstance(context, AnalysisContext):
        raise InputValidationError("run_phase5 requires AnalysisContext.")
    selected = settings or Phase5Settings.from_config(context.config)
    if not selected.enabled:
        raise ConfigurationError("Phase 5 requires tr.enabled: true.")
    alignments = _resolve_alignments(context, dict(alignment_artifacts or {}))
    catalog = TandemRepeatCatalog(selected.catalog, selected.catalog_reference_build)
    output_root = Path(output_directory).expanduser()
    tr_root = output_root / "tr"
    work_root = output_root / "work" / "tr"
    log_root = output_root / "logs" / "tr"
    runner = wrapper or TrgtWrapper(
        executable=selected.executable,
        bcftools_executable=selected.bcftools_executable,
        samtools_executable=selected.samtools_executable,
    )
    sample_results: list[Phase5SampleResult] = []
    records = {record.sample.sample_id: record for record in context.samples}
    for alignment in alignments:
        karyotype = resolve_trgt_karyotype(selected.karyotype, records[alignment.sample_id].sex, alignment.sample_id)
        request = TrgtRequest(
            artifact=alignment,
            catalog=catalog,
            raw_output_prefix=work_root / alignment.sample_id / f"{alignment.sample_id}.trgt",
            final_vcf=tr_root / f"{alignment.sample_id}.tr.vcf.gz",
            final_spanning_bam=tr_root / f"{alignment.sample_id}.tr.spanning.bam",
            karyotype=karyotype,
            resources=selected.resources,
            preset=selected.preset,
            overwrite=selected.overwrite,
        )
        result = runner.run(
            request,
            dry_run=dry_run,
            stderr_path=log_root / f"{alignment.sample_id}.trgt.log",
        )
        sample_results.append(Phase5SampleResult(alignment, karyotype, result))
    return Phase5RunReport(context, selected, tuple(sample_results), dry_run)


def resolve_trgt_karyotype(configured: str, sex: Sex | None, sample_id: str) -> str:
    """Resolve an explicit TRGT karyotype without inferring missing sex."""
    if configured in {"XX", "XY"}:
        return configured
    if configured != "auto":
        raise ConfigurationError("tr.karyotype must be auto, XX, or XY.")
    if sex is Sex.FEMALE:
        return "XX"
    if sex is Sex.MALE:
        return "XY"
    raise ConfigurationError(
        f"TRGT karyotype for sample '{sample_id}' cannot be inferred: declare sample sex or set tr.karyotype explicitly."
    )


def _resolve_alignments(
    context: AnalysisContext,
    supplied: Mapping[str, AlignmentArtifact],
) -> tuple[AlignmentArtifact, ...]:
    unknown = sorted(set(supplied).difference(context.sample_ids))
    if unknown:
        raise InputValidationError(f"Phase 5 received artifacts for unknown samples: {unknown!r}.")
    resolved: list[AlignmentArtifact] = []
    for record in context.samples:
        sample = record.sample
        artifact = supplied.get(sample.sample_id)
        if artifact is None:
            if sample.input.input_type is InputType.FASTQ:
                raise InputValidationError(
                    f"Phase 5 sample '{sample.sample_id}' has raw FASTQ; provide a completed Phase 2 AlignmentArtifact."
                )
            path = sample.input.files[0]
            artifact = AlignmentArtifact(
                sample.sample_id,
                path,
                AlignmentOutputFormat.BAM if sample.input.input_type is InputType.BAM else AlignmentOutputFormat.CRAM,
                context.reference,
                AlignmentSource.EXISTING,
                AlignmentSortOrder.UNKNOWN,
                find_alignment_index(path),
            )
        if artifact.sample_id != sample.sample_id:
            raise InputValidationError(f"Phase 5 alignment sample mismatch for '{sample.sample_id}'.")
        if (
            artifact.reference.fasta.absolute() != context.reference.fasta.absolute()
            or artifact.reference.build != context.reference.build
            or artifact.reference.contigs != context.reference.contigs
        ):
            raise ReferenceError(f"Phase 5 alignment reference conflicts for sample '{sample.sample_id}'.")
        if artifact.output_format is not AlignmentOutputFormat.BAM:
            raise InputValidationError("Phase 5 TRGT currently requires BAM input; CRAM is not accepted.")
        resolved.append(validate_alignment_artifact(artifact, require_index=True))
    return tuple(resolved)


__all__ = [
    "PHASE5_REPORT_SCHEMA_VERSION",
    "Phase5RunReport",
    "Phase5RunStatus",
    "Phase5SampleResult",
    "Phase5Settings",
    "resolve_trgt_karyotype",
    "run_phase5",
]
