"""Explicit-selection functional-prioritization boundary for AlphaGenome."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Protocol

from hifivar.annotation import VariantCategory
from hifivar.exceptions import InputValidationError
from hifivar.validation import validate_file


@dataclass(frozen=True, slots=True)
class PrioritizedVariant:
    """One explicitly selected variant; never inferred from an impact score."""

    sample_id: str
    source_variant_id: str
    variant_category: VariantCategory
    contig: str
    position: int
    reference_bases: str
    alternate_bases: str
    source_annotation: Path
    selection_reason: str

    def __post_init__(self) -> None:
        for name in (
            "sample_id", "source_variant_id", "contig", "reference_bases",
            "alternate_bases", "selection_reason",
        ):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip() or any(char in value for char in "\r\n"):
                raise InputValidationError(f"Prioritized variant {name} must be non-empty and single-line.")
        if not isinstance(self.variant_category, VariantCategory):
            raise InputValidationError("Prioritized variant category is invalid.")
        if not isinstance(self.position, int) or isinstance(self.position, bool) or self.position < 1:
            raise InputValidationError("Prioritized variant position must be positive.")
        object.__setattr__(self, "source_annotation", validate_file(self.source_annotation))

    def to_dict(self) -> dict[str, object]:
        return {
            "sample": self.sample_id,
            "source_variant_id": self.source_variant_id,
            "variant_category": self.variant_category.value,
            "contig": self.contig,
            "position": self.position,
            "reference_bases": self.reference_bases,
            "alternate_bases": self.alternate_bases,
            "source_annotation": str(self.source_annotation),
            "selection_reason": self.selection_reason,
            "selection_is_explicit": True,
        }


@dataclass(frozen=True, slots=True)
class FunctionalPrioritizationRequest:
    selected_variants: tuple[PrioritizedVariant, ...]
    model_name: str
    model_version: str
    requested_modalities: tuple[str, ...]

    def __post_init__(self) -> None:
        selected = tuple(self.selected_variants)
        if not selected:
            raise InputValidationError(
                "Functional prioritization requires an explicit non-empty selection."
            )
        for name in ("model_name", "model_version"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise InputValidationError(f"Functional model {name} must be non-empty.")
        modalities = tuple(self.requested_modalities)
        if not modalities or any(not isinstance(item, str) or not item.strip() for item in modalities):
            raise InputValidationError("Functional prioritization requires explicit output modalities.")
        keys = [(item.sample_id, item.source_variant_id) for item in selected]
        if len(keys) != len(set(keys)):
            raise InputValidationError("Functional selection contains duplicate source variants.")
        object.__setattr__(self, "selected_variants", selected)
        object.__setattr__(self, "requested_modalities", modalities)

    def to_dict(self) -> dict[str, object]:
        return {
            "selected_variants": [item.to_dict() for item in self.selected_variants],
            "model_name": self.model_name,
            "model_version": self.model_version,
            "requested_modalities": list(self.requested_modalities),
            "selection_is_explicit": True,
            "whole_genome_unselected_execution": False,
        }


@dataclass(frozen=True, slots=True)
class FunctionalPrediction:
    variant: PrioritizedVariant
    modality_scores: Mapping[str, float]
    backend_record_id: str | None = None

    def __post_init__(self) -> None:
        scores = dict(self.modality_scores)
        if not scores or any(not isinstance(key, str) or not key for key in scores):
            raise InputValidationError("Functional prediction requires named modality scores.")
        if any(not isinstance(value, (int, float)) or isinstance(value, bool) for value in scores.values()):
            raise InputValidationError("Functional modality scores must be numeric.")
        object.__setattr__(self, "modality_scores", scores)

    def to_dict(self) -> dict[str, object]:
        return {
            "variant": self.variant.to_dict(),
            "modality_scores": dict(self.modality_scores),
            "backend_record_id": self.backend_record_id,
            "functional_impact_is_call_confidence": False,
            "functional_impact_is_truth": False,
        }


@dataclass(frozen=True, slots=True)
class FunctionalPrioritizationResult:
    request: FunctionalPrioritizationRequest
    predictions: tuple[FunctionalPrediction, ...]
    backend_version: str

    def __post_init__(self) -> None:
        predictions = tuple(self.predictions)
        selected = {(item.sample_id, item.source_variant_id) for item in self.request.selected_variants}
        observed = {(item.variant.sample_id, item.variant.source_variant_id) for item in predictions}
        if observed != selected or len(predictions) != len(selected):
            raise InputValidationError("Functional backend did not return exactly the explicit selection.")
        if not isinstance(self.backend_version, str) or not self.backend_version.strip():
            raise InputValidationError("Functional backend version must be non-empty.")
        object.__setattr__(self, "predictions", predictions)

    def to_dict(self) -> dict[str, object]:
        return {
            "request": self.request.to_dict(),
            "backend_version": self.backend_version,
            "predictions": [item.to_dict() for item in self.predictions],
            "scientific_policy": {
                "functional_impact_is_call_confidence": False,
                "functional_impact_is_truth": False,
            },
        }


class FunctionalBackend(Protocol):
    """Credential/backend-neutral AlphaGenome-compatible interface."""

    def predict(self, request: FunctionalPrioritizationRequest) -> FunctionalPrioritizationResult:
        ...


def run_functional_prioritization(
    request: FunctionalPrioritizationRequest,
    *,
    backend: FunctionalBackend,
) -> FunctionalPrioritizationResult:
    """Run only the already explicit selection using an injected backend."""
    result = backend.predict(request)
    if not isinstance(result, FunctionalPrioritizationResult) or result.request != request:
        raise InputValidationError("Functional backend returned an incompatible result.")
    return result


def read_functional_selection(path: str | Path) -> tuple[PrioritizedVariant, ...]:
    """Read the explicit AlphaGenome candidate list; never rank a whole VCF."""
    source = validate_file(path)
    required = (
        "sample", "source_variant_id", "variant_category", "contig", "position",
        "reference_bases", "alternate_bases", "source_annotation", "selection_reason",
    )
    selected: list[PrioritizedVariant] = []
    try:
        with source.open("rt", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = tuple(reader.fieldnames or ())
            if set(fields) != set(required):
                raise InputValidationError(
                    f"Functional selection columns must be exactly {required!r}."
                )
            for line_number, row in enumerate(reader, 2):
                annotation = Path(row["source_annotation"]).expanduser()
                if not annotation.is_absolute():
                    annotation = source.parent / annotation
                try:
                    selected.append(
                        PrioritizedVariant(
                            row["sample"], row["source_variant_id"],
                            VariantCategory(row["variant_category"].strip().lower()),
                            row["contig"], int(row["position"]), row["reference_bases"],
                            row["alternate_bases"], annotation, row["selection_reason"],
                        )
                    )
                except (KeyError, ValueError, TypeError) as error:
                    raise InputValidationError(
                        f"Invalid functional selection row {line_number}: {error}."
                    ) from error
    except (OSError, UnicodeError, csv.Error) as error:
        raise InputValidationError(f"Unable to read functional selection '{source}': {error}") from error
    return tuple(selected)


__all__ = [
    "FunctionalBackend", "FunctionalPrediction", "FunctionalPrioritizationRequest",
    "FunctionalPrioritizationResult", "PrioritizedVariant",
    "read_functional_selection", "run_functional_prioritization",
]
