"""Run-level reference, sample, and effective-configuration integration."""

from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

from hifivar.config import HiFiVarConfig
from hifivar.exceptions import ConfigurationError, InputValidationError, ReferenceError
from hifivar.reference import ReferenceGenome, normalize_reference_build
from hifivar.sample import InputType, Sample
from hifivar.sample_sheet import SampleRecord, SampleSheet


ConfigValue = Mapping[str, object] | HiFiVarConfig


@dataclass(frozen=True, slots=True)
class AnalysisContext:
    """Validated run-level inputs without workflow planning or tool execution."""

    reference: ReferenceGenome
    samples: tuple[SampleRecord, ...]
    config: ConfigValue
    source_sample_sheet: Path | None = None

    def __post_init__(self) -> None:
        """Detach caller data and enforce cross-object invariants."""
        if not isinstance(self.reference, ReferenceGenome):
            raise ReferenceError("AnalysisContext reference must be a ReferenceGenome.")

        records = tuple(self.samples)
        if not records:
            raise InputValidationError(
                "AnalysisContext must contain at least one sample."
            )
        if any(not isinstance(record, SampleRecord) for record in records):
            raise InputValidationError(
                "AnalysisContext samples must contain only SampleRecord objects."
            )
        object.__setattr__(self, "samples", records)

        if isinstance(self.config, HiFiVarConfig):
            detached_config: ConfigValue = HiFiVarConfig(
                self.config.to_dict(),
                self.config.sources,
            )
        elif isinstance(self.config, Mapping):
            detached_config = deepcopy(dict(self.config))
        else:
            raise ConfigurationError(
                "AnalysisContext config must be a mapping or HiFiVarConfig."
            )
        object.__setattr__(self, "config", detached_config)

        if self.source_sample_sheet is not None:
            if not isinstance(self.source_sample_sheet, (str, Path)):
                raise InputValidationError(
                    "source_sample_sheet must be a string, Path, or None."
                )
            object.__setattr__(
                self,
                "source_sample_sheet",
                Path(self.source_sample_sheet).expanduser(),
            )

        _validate_unique_records(records)
        _validate_reference_config(self.reference, detached_config)

    @classmethod
    def from_sample(
        cls,
        reference: ReferenceGenome,
        sample: Sample,
        config: ConfigValue,
    ) -> AnalysisContext:
        """Wrap one Phase 1.2 sample in a metadata-empty SampleRecord."""
        if not isinstance(sample, Sample):
            raise InputValidationError("sample must be a Sample instance.")
        return cls(
            reference=reference,
            samples=(SampleRecord(sample=sample),),
            config=config,
        )

    @classmethod
    def from_sample_sheet(
        cls,
        reference: ReferenceGenome,
        sample_sheet: SampleSheet,
        config: ConfigValue,
    ) -> AnalysisContext:
        """Create a context from an already validated SampleSheet."""
        if not isinstance(sample_sheet, SampleSheet):
            raise InputValidationError(
                "sample_sheet must be a SampleSheet instance."
            )
        return cls(
            reference=reference,
            samples=sample_sheet.records,
            config=config,
            source_sample_sheet=sample_sheet.source_path,
        )

    @classmethod
    def from_config(cls, config: ConfigValue) -> AnalysisContext:
        """Load reference and sample-sheet models from an effective config.

        Relative analysis input paths are accepted only when ``config`` retains
        a user-config source; they are then interpreted beside that YAML file.
        This factory does not add side effects to the configuration loader.
        """
        if not isinstance(config, Mapping):
            raise ConfigurationError(
                "AnalysisContext.from_config requires a configuration mapping."
            )
        reference_config = _required_section(config, "reference")
        samples_config = _required_section(config, "samples")
        fasta_value = _required_path_value(reference_config, "reference.fasta")
        sheet_value = _required_path_value(samples_config, "samples.sheet")
        fasta_path = _resolve_analysis_path(fasta_value, config, "reference.fasta")
        sheet_path = _resolve_analysis_path(sheet_value, config, "samples.sheet")

        build = reference_config.get("build")
        if build is not None and not isinstance(build, str):
            raise ConfigurationError(
                "reference.build must be a non-empty string or null."
            )
        reference = ReferenceGenome.from_fasta(fasta_path, build=build)
        sample_sheet = SampleSheet.from_tsv(sheet_path)
        return cls.from_sample_sheet(reference, sample_sheet, config)

    @property
    def sample_ids(self) -> tuple[str, ...]:
        """Return sample IDs in deterministic context order."""
        return tuple(record.sample.sample_id for record in self.samples)

    @property
    def n_samples(self) -> int:
        """Return the number of samples in the analysis."""
        return len(self.samples)

    @property
    def input_types(self) -> tuple[InputType, ...]:
        """Return input types in deterministic sample order."""
        return tuple(record.sample.input.input_type for record in self.samples)

    def validate_query_contigs(self, query_contigs: Iterable[str]) -> None:
        """Delegate exact-name contig compatibility to the reference model."""
        self.reference.validate_contigs(query_contigs)

    def reference_compatibility(self) -> list[dict[str, str]]:
        """Describe the intentionally limited Phase 1 compatibility boundary."""
        statuses: list[dict[str, str]] = []
        for record in self.samples:
            input_type = record.sample.input.input_type
            if input_type is InputType.FASTQ:
                status = "not_applicable"
                reason = "FASTQ input is unaligned; reference compatibility is not claimed."
            else:
                status = "not_checked"
                reason = (
                    "BAM/CRAM header and reference compatibility are not checked "
                    "in Phase 1.4."
                )
            statuses.append(
                {
                    "sample_id": record.sample.sample_id,
                    "status": status,
                    "reason": reason,
                }
            )
        return statuses

    def config_to_dict(self) -> dict[str, object]:
        """Return an independent standard mapping of the effective config."""
        if isinstance(self.config, HiFiVarConfig):
            return self.config.to_dict()
        return deepcopy(dict(self.config))

    def to_dict(self) -> dict[str, object]:
        """Return a JSON/YAML-friendly context summary."""
        return {
            "reference": self.reference.to_dict(include_contigs=True),
            "samples": [record.to_dict() for record in self.samples],
            "config": self.config_to_dict(),
            "source_sample_sheet": (
                str(self.source_sample_sheet)
                if self.source_sample_sheet is not None
                else None
            ),
            "reference_compatibility": self.reference_compatibility(),
        }


def _required_section(
    config: Mapping[str, object],
    section_name: str,
) -> Mapping[str, object]:
    section = config.get(section_name)
    if not isinstance(section, Mapping):
        raise ConfigurationError(
            f"Analysis context requires configuration section '{section_name}'."
        )
    return section


def _required_path_value(section: Mapping[str, object], field_name: str) -> str:
    key = field_name.rsplit(".", maxsplit=1)[1]
    value = section.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(
            f"Analysis context requires {field_name} to be a non-empty path."
        )
    return value


def _resolve_analysis_path(
    value: str,
    config: ConfigValue,
    field_name: str,
) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    if isinstance(config, HiFiVarConfig):
        user_source = config.sources.get("user")
        if user_source is not None:
            return user_source.absolute().parent / path
    raise ConfigurationError(
        f"Relative {field_name} requires a loaded user config source; "
        "otherwise provide an absolute path."
    )


def _validate_unique_records(records: tuple[SampleRecord, ...]) -> None:
    sample_ids: set[str] = set()
    input_owners: dict[str, str] = {}
    for record in records:
        sample_id = record.sample.sample_id
        if sample_id in sample_ids:
            raise InputValidationError(
                f"AnalysisContext contains duplicate sample_id '{sample_id}'."
            )
        sample_ids.add(sample_id)
        for path in record.sample.input.files:
            identity = _absolute_identity(path)
            previous_owner = input_owners.get(identity)
            if previous_owner is not None:
                raise InputValidationError(
                    f"Primary input '{path}' is reused by samples "
                    f"'{previous_owner}' and '{sample_id}' in this analysis."
                )
            input_owners[identity] = sample_id


def _validate_reference_config(
    reference: ReferenceGenome,
    config: ConfigValue,
) -> None:
    reference_config = config.get("reference")
    if reference_config is None:
        return
    if not isinstance(reference_config, Mapping):
        raise ConfigurationError("Configuration section reference must be a mapping.")

    configured_build = reference_config.get("build")
    if configured_build is not None:
        try:
            normalized_build = normalize_reference_build(configured_build)  # type: ignore[arg-type]
        except ReferenceError as error:
            raise ConfigurationError(f"Invalid reference.build: {error}") from error
        if normalized_build != reference.build:
            raise ReferenceError(
                "Configured reference.build conflicts with the ReferenceGenome: "
                f"{normalized_build!r} != {reference.build!r}."
            )

    configured_fasta = reference_config.get("fasta")
    if configured_fasta is not None:
        if not isinstance(configured_fasta, str) or not configured_fasta.strip():
            raise ConfigurationError(
                "reference.fasta must be a non-empty string path or null."
            )
        configured_path = _resolve_analysis_path(
            configured_fasta,
            config,
            "reference.fasta",
        )
        if _absolute_identity(configured_path) != _absolute_identity(reference.fasta):
            raise ReferenceError(
                "Configured reference.fasta conflicts with the ReferenceGenome: "
                f"'{configured_path}' != '{reference.fasta}'."
            )


def _absolute_identity(path: Path) -> str:
    return os.path.normcase(os.path.normpath(str(path.absolute())))


__all__ = ["AnalysisContext"]
