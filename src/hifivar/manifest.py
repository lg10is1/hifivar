"""Portable provenance manifests for validated Phase 1 analysis inputs."""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

import yaml

from hifivar import __version__, validation
from hifivar.context import AnalysisContext
from hifivar.exceptions import InputValidationError
from hifivar.serialization import (
    redact_sensitive_data,
    standardize_data,
    utc_now_iso8601,
    write_json_atomic,
    write_yaml_atomic,
)


MANIFEST_SCHEMA_VERSION: Final[str] = "1.0"
_REQUIRED_KEYS = frozenset(
    {
        "schema_version",
        "hifivar_version",
        "created_at",
        "reference",
        "samples",
        "config",
        "source_sample_sheet",
        "inputs",
    }
)


@dataclass(frozen=True, slots=True)
class RunManifest:
    """A serializable snapshot of validated inputs and effective settings."""

    schema_version: str
    hifivar_version: str
    created_at: str
    reference: dict[str, object]
    samples: list[dict[str, object]]
    config: dict[str, object]
    source_sample_sheet: str | None
    inputs: list[dict[str, object]]
    reference_compatibility: list[dict[str, object]] | None = None

    def __post_init__(self) -> None:
        """Detach nested values and reject non-portable payload types."""
        if not isinstance(self.schema_version, str) or not self.schema_version:
            raise InputValidationError("Manifest schema_version must be a string.")
        if not isinstance(self.hifivar_version, str) or not self.hifivar_version:
            raise InputValidationError("Manifest hifivar_version must be a string.")
        if not isinstance(self.created_at, str) or not self.created_at:
            raise InputValidationError("Manifest created_at must be a string.")
        if self.source_sample_sheet is not None and not isinstance(
            self.source_sample_sheet, str
        ):
            raise InputValidationError(
                "Manifest source_sample_sheet must be a string or null."
            )

        reference = _standardize(self.reference)
        samples = _standardize(self.samples)
        config = _standardize(self.config)
        inputs = _standardize(self.inputs)
        compatibility = _standardize(self.reference_compatibility)
        if not isinstance(reference, dict):
            raise InputValidationError("Manifest reference must be a mapping.")
        if not isinstance(samples, list) or not all(
            isinstance(item, dict) for item in samples
        ):
            raise InputValidationError("Manifest samples must be a list of mappings.")
        if not isinstance(config, dict):
            raise InputValidationError("Manifest config must be a mapping.")
        if not isinstance(inputs, list) or not all(
            isinstance(item, dict) for item in inputs
        ):
            raise InputValidationError("Manifest inputs must be a list of mappings.")
        if compatibility is not None and (
            not isinstance(compatibility, list)
            or not all(isinstance(item, dict) for item in compatibility)
        ):
            raise InputValidationError(
                "Manifest reference_compatibility must be a list of mappings or null."
            )
        object.__setattr__(self, "reference", reference)
        object.__setattr__(self, "samples", samples)
        object.__setattr__(self, "config", config)
        object.__setattr__(self, "inputs", inputs)
        object.__setattr__(self, "reference_compatibility", compatibility)

    @classmethod
    def from_context(
        cls,
        context: AnalysisContext,
        *,
        compute_input_checksums: bool = False,
    ) -> RunManifest:
        """Snapshot one AnalysisContext without executing tools or a workflow."""
        if not isinstance(context, AnalysisContext):
            raise InputValidationError(
                "RunManifest.from_context requires an AnalysisContext."
            )

        reference = context.reference.to_dict(include_contigs=True)
        reference["fasta"] = _absolute_path(context.reference.fasta)
        reference["fai"] = _absolute_path(context.reference.fai)

        samples: list[dict[str, object]] = []
        inputs: list[dict[str, object]] = []
        for record in context.samples:
            sample_metadata = record.to_dict()
            input_metadata = sample_metadata.get("input")
            if not isinstance(input_metadata, dict):
                raise InputValidationError(
                    f"Sample '{record.sample.sample_id}' has invalid input metadata."
                )
            input_metadata["files"] = [
                _absolute_path(path) for path in record.sample.input.files
            ]
            samples.append(sample_metadata)

            for path in record.sample.input.files:
                absolute_path = Path(_absolute_path(path))
                try:
                    size_bytes = absolute_path.stat().st_size
                except OSError as error:
                    raise InputValidationError(
                        f"Unable to inspect input file '{absolute_path}': {error}"
                    ) from error
                checksum = (
                    validation.compute_sha256(path)
                    if compute_input_checksums
                    else None
                )
                inputs.append(
                    {
                        "sample_id": record.sample.sample_id,
                        "path": str(absolute_path),
                        "input_type": record.sample.input.input_type,
                        "size_bytes": size_bytes,
                        "sha256": checksum,
                    }
                )

        source_sheet = (
            _absolute_path(context.source_sample_sheet)
            if context.source_sample_sheet is not None
            else None
        )
        return cls(
            schema_version=MANIFEST_SCHEMA_VERSION,
            hifivar_version=__version__,
            created_at=_utc_timestamp(),
            reference=reference,
            samples=samples,
            config=redact_sensitive_data(context.config_to_dict()),
            source_sample_sheet=source_sheet,
            inputs=inputs,
            reference_compatibility=context.reference_compatibility(),
        )

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> RunManifest:
        """Build a manifest data object from parsed JSON/YAML content."""
        if not isinstance(payload, Mapping):
            raise InputValidationError("Manifest root must be a mapping.")
        missing = sorted(_REQUIRED_KEYS.difference(payload))
        if missing:
            raise InputValidationError(
                f"Manifest is missing required key(s): {', '.join(missing)}."
            )
        try:
            return cls(
                schema_version=payload["schema_version"],  # type: ignore[arg-type]
                hifivar_version=payload["hifivar_version"],  # type: ignore[arg-type]
                created_at=payload["created_at"],  # type: ignore[arg-type]
                reference=payload["reference"],  # type: ignore[arg-type]
                samples=payload["samples"],  # type: ignore[arg-type]
                config=payload["config"],  # type: ignore[arg-type]
                source_sample_sheet=payload["source_sample_sheet"],  # type: ignore[arg-type]
                inputs=payload["inputs"],  # type: ignore[arg-type]
                reference_compatibility=payload.get("reference_compatibility"),  # type: ignore[arg-type]
            )
        except (KeyError, TypeError) as error:
            raise InputValidationError(f"Invalid manifest structure: {error}") from error

    @classmethod
    def from_json(cls, path: str | Path) -> RunManifest:
        """Parse one UTF-8 JSON manifest without rebuilding an analysis context."""
        manifest_path = validation.validate_file(path)
        try:
            with manifest_path.open("rt", encoding="utf-8") as handle:
                payload = json.load(handle)
        except (OSError, UnicodeError, json.JSONDecodeError) as error:
            raise InputValidationError(
                f"Unable to read JSON manifest '{manifest_path}': {error}"
            ) from error
        return cls.from_dict(payload)

    @classmethod
    def from_yaml(cls, path: str | Path) -> RunManifest:
        """Parse one UTF-8 YAML manifest without rebuilding an analysis context."""
        manifest_path = validation.validate_file(path)
        try:
            with manifest_path.open("rt", encoding="utf-8") as handle:
                payload = yaml.safe_load(handle)
        except (OSError, UnicodeError, yaml.YAMLError) as error:
            raise InputValidationError(
                f"Unable to read YAML manifest '{manifest_path}': {error}"
            ) from error
        return cls.from_dict(payload)

    def to_dict(self) -> dict[str, object]:
        """Return an independent payload containing only standard data types."""
        payload: dict[str, object] = {
            "schema_version": self.schema_version,
            "hifivar_version": self.hifivar_version,
            "created_at": self.created_at,
            "reference": self.reference,
            "samples": self.samples,
            "config": self.config,
            "source_sample_sheet": self.source_sample_sheet,
            "inputs": self.inputs,
        }
        if self.reference_compatibility is not None:
            payload["reference_compatibility"] = self.reference_compatibility
        standardized = _standardize(payload)
        if not isinstance(standardized, dict):  # pragma: no cover - invariant guard
            raise InputValidationError("Manifest serialization produced invalid data.")
        return standardized

    def write_json(self, path: str | Path, *, overwrite: bool = False) -> Path:
        """Atomically write UTF-8 JSON, refusing replacement by default."""
        return write_json_atomic(
            self.to_dict(),
            path,
            overwrite=overwrite,
            artifact_name="Manifest",
        )

    def write_yaml(self, path: str | Path, *, overwrite: bool = False) -> Path:
        """Atomically write UTF-8 YAML, refusing replacement by default."""
        return write_yaml_atomic(
            self.to_dict(),
            path,
            overwrite=overwrite,
            artifact_name="Manifest",
        )


def _utc_timestamp() -> str:
    return utc_now_iso8601()


def _absolute_path(path: str | Path) -> str:
    """Return an absolute spelling without resolving symbolic links."""
    return str(Path(path).expanduser().absolute())


def _standardize(value: object) -> Any:
    """Convert nested values to JSON/YAML primitive containers and scalars."""
    return standardize_data(value, context="Manifest value")


__all__ = ["MANIFEST_SCHEMA_VERSION", "RunManifest"]
