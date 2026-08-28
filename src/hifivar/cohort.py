"""Phase 12 cohort contracts, explicit call states, and streaming QC."""

from __future__ import annotations

import gzip
import json
import csv
import hashlib
from collections import Counter
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import TextIO
from datetime import datetime, timezone

import yaml

from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.reference import ReferenceGenome
from hifivar.sample import validate_sample_id


class CohortTrack(str, Enum):
    SMALL_VARIANTS = "small_variants"
    SV = "sv"
    TR = "tr"


class SampleCallState(str, Enum):
    """Per-sample state; absence of a record is never inferred as hom-ref."""

    CALLED = "CALLED"
    NO_CALLS = "NO_CALLS"
    NOT_RUN = "NOT_RUN"
    FAILED = "FAILED"
    DISABLED = "DISABLED"
    MISSING_INPUT = "MISSING_INPUT"
    NOT_OBSERVED = "NOT_OBSERVED"


@dataclass(frozen=True, slots=True)
class CohortDefinition:
    cohort_id: str
    sample_ids: tuple[str, ...]
    reference: ReferenceGenome

    def __post_init__(self) -> None:
        validate_sample_id(self.cohort_id)
        samples = tuple(self.sample_ids)
        if not samples:
            raise InputValidationError("A cohort requires at least one sample.")
        for sample_id in samples:
            validate_sample_id(sample_id)
        if len(samples) != len(set(samples)):
            raise InputValidationError("A cohort cannot contain duplicate sample IDs.")
        if not isinstance(self.reference, ReferenceGenome):
            raise InputValidationError("Cohort reference must be a ReferenceGenome.")
        object.__setattr__(self, "sample_ids", samples)

    def to_dict(self) -> dict[str, object]:
        return {
            "cohort_id": self.cohort_id,
            "sample_ids": list(self.sample_ids),
            "reference": self.reference.to_dict(include_contigs=True),
        }


@dataclass(frozen=True, slots=True)
class CohortSampleInput:
    sample_id: str
    state: SampleCallState
    source_path: Path | None = None
    index_path: Path | None = None
    source_tool: str | None = None
    source_version: str | None = None
    reference_build: str | None = None
    catalog_id: str | None = None
    detail: str | None = None

    def __post_init__(self) -> None:
        validate_sample_id(self.sample_id)
        if not isinstance(self.state, SampleCallState):
            raise InputValidationError("Cohort sample state must be a SampleCallState.")
        for name in ("source_path", "index_path"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, Path(value).expanduser())
        if self.state in {SampleCallState.CALLED, SampleCallState.NO_CALLS} and self.source_path is None:
            raise InputValidationError(f"{self.state.value} sample '{self.sample_id}' requires a source artifact.")

    @property
    def callable(self) -> bool:
        return self.state in {SampleCallState.CALLED, SampleCallState.NO_CALLS}

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "state": self.state.value,
            "source_path": str(self.source_path) if self.source_path else None,
            "index_path": str(self.index_path) if self.index_path else None,
            "source_tool": self.source_tool,
            "source_version": self.source_version,
            "reference_build": self.reference_build,
            "catalog_id": self.catalog_id,
            "detail": self.detail,
        }


@dataclass(frozen=True, slots=True)
class CohortTrackResult:
    track: CohortTrack
    enabled: bool
    state: SampleCallState
    tool: str | None
    tool_version: str | None
    outputs: tuple[Path, ...] = ()
    commands: tuple[tuple[str, ...], ...] = ()
    sample_states: tuple[CohortSampleInput, ...] = ()
    metrics: dict[str, object] = field(default_factory=dict)
    message: str | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "track": self.track.value,
            "enabled": self.enabled,
            "state": self.state.value,
            "tool": self.tool,
            "tool_version": self.tool_version,
            "outputs": [str(path) for path in self.outputs],
            "commands": [list(command) for command in self.commands],
            "sample_states": [item.to_dict() for item in self.sample_states],
            "metrics": self.metrics,
            "message": self.message,
        }


@dataclass(frozen=True, slots=True)
class CohortManifest:
    cohort: CohortDefinition
    tracks: tuple[CohortTrackResult, ...]
    hifivar_version: str
    config_path: Path | None = None
    git_commit: str | None = None
    effective_config: dict[str, object] | None = None
    created_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())

    def __post_init__(self) -> None:
        tracks = tuple(self.tracks)
        identities = [item.track for item in tracks]
        if len(identities) != len(set(identities)):
            raise InputValidationError("Cohort manifest contains a duplicate track.")
        object.__setattr__(self, "tracks", tracks)
        if self.config_path is not None:
            object.__setattr__(self, "config_path", Path(self.config_path))

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": 1,
            "hifivar_version": self.hifivar_version,
            "cohort": self.cohort.to_dict(),
            "sample_order_sha256": hashlib.sha256(("\n".join(self.cohort.sample_ids) + "\n").encode("utf-8")).hexdigest(),
            "config_path": str(self.config_path) if self.config_path else None,
            "git_commit": self.git_commit,
            "effective_config": self.effective_config,
            "created_at": self.created_at,
            "tracks": [track.to_dict() for track in self.tracks],
        }

    def write(self, json_path: Path, yaml_path: Path | None = None) -> None:
        payload = self.to_dict()
        _write_text_no_overwrite(Path(json_path), json.dumps(payload, indent=2, sort_keys=True) + "\n")
        if yaml_path is not None:
            _write_text_no_overwrite(
                Path(yaml_path),
                yaml.safe_dump(payload, allow_unicode=True, sort_keys=False),
            )


def validate_track_inputs(
    cohort: CohortDefinition,
    inputs: tuple[CohortSampleInput, ...],
    *,
    require_all: bool = True,
) -> None:
    """Validate ordered sample identity and reference compatibility."""
    observed = tuple(item.sample_id for item in inputs)
    if require_all and observed != cohort.sample_ids:
        raise InputValidationError(
            f"Cohort input sample order/set mismatch: expected {cohort.sample_ids!r}, observed {observed!r}."
        )
    if len(observed) != len(set(observed)):
        raise InputValidationError("Cohort track contains duplicate sample inputs.")
    unknown = sorted(set(observed).difference(cohort.sample_ids))
    if unknown:
        raise InputValidationError(f"Cohort track contains unknown samples: {unknown!r}.")
    for item in inputs:
        if item.reference_build is not None and item.reference_build != cohort.reference.build:
            raise InputValidationError(
                f"Reference build mismatch for '{item.sample_id}': {item.reference_build!r} != {cohort.reference.build!r}."
            )


def scan_multisample_vcf(path: Path, expected_samples: tuple[str, ...]) -> dict[str, object]:
    """Stream a cohort VCF and compute bounded-memory descriptive QC."""
    variant_count = 0
    multiallelic_count = 0
    filter_counts: Counter[str] = Counter()
    non_ref = {sample: 0 for sample in expected_samples}
    missing = {sample: 0 for sample in expected_samples}
    called = {sample: 0 for sample in expected_samples}
    header_samples: tuple[str, ...] | None = None
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM\t"):
                header_samples = tuple(line.rstrip("\r\n").split("\t")[9:])
                if any(not sample for sample in header_samples):
                    raise OutputValidationError(
                        f"Cohort VCF contains an empty sample name in '{path}'."
                    )
                duplicate_samples = sorted(
                    sample for sample, count in Counter(header_samples).items() if count > 1
                )
                if duplicate_samples:
                    raise OutputValidationError(
                        f"Cohort VCF contains duplicate sample IDs {duplicate_samples!r} in '{path}'."
                    )
                expected_set = set(expected_samples)
                observed_set = set(header_samples)
                missing_samples = sorted(expected_set - observed_set)
                extra_samples = sorted(observed_set - expected_set)
                if missing_samples or extra_samples:
                    raise OutputValidationError(
                        "Cohort VCF sample set mismatch: "
                        f"missing={missing_samples!r}, extra={extra_samples!r}; "
                        f"declared order={expected_samples!r}, output order={header_samples!r}."
                    )
                continue
            if line.startswith("#") or not line.strip():
                continue
            if header_samples is None:
                raise OutputValidationError(f"VCF records precede #CHROM header in '{path}'.")
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9 + len(header_samples):
                raise OutputValidationError(f"Malformed cohort VCF record in '{path}'.")
            variant_count += 1
            multiallelic_count += int("," in fields[4])
            filter_counts[fields[6]] += 1
            format_keys = fields[8].split(":")
            try:
                gt_index = format_keys.index("GT")
            except ValueError:
                gt_index = -1
            for sample, cell in zip(header_samples, fields[9:]):
                values = cell.split(":")
                gt = values[gt_index] if gt_index >= 0 and gt_index < len(values) else "."
                alleles = gt.replace("|", "/").split("/")
                if not alleles or any(allele in {"", "."} for allele in alleles):
                    missing[sample] += 1
                else:
                    called[sample] += 1
                    if any(allele != "0" for allele in alleles):
                        non_ref[sample] += 1
    if header_samples is None:
        raise OutputValidationError(f"VCF lacks #CHROM header: '{path}'.")
    denominator = max(variant_count, 1)
    return {
        "sample_count": len(expected_samples),
        "declared_sample_order": list(expected_samples),
        "output_sample_order": list(header_samples),
        "sample_set_match": True,
        "sample_order_match": header_samples == expected_samples,
        "variant_count": variant_count,
        "multiallelic_count": multiallelic_count,
        "filter_distribution": dict(sorted(filter_counts.items())),
        "per_sample_non_ref_count": non_ref,
        "per_sample_missing_rate": {key: value / denominator for key, value in missing.items()},
        "per_sample_call_rate": {key: value / denominator for key, value in called.items()},
    }


def read_cohort_input_manifest(
    path: Path,
    cohort: CohortDefinition,
    track: CohortTrack,
) -> tuple[CohortSampleInput, ...]:
    """Read one explicit track from the long-form Phase 12 TSV manifest."""
    required = {"sample", "track", "state", "source_path", "index_path", "source_tool", "source_version", "reference_build", "catalog_id"}
    rows: list[CohortSampleInput] = []
    manifest = Path(path)
    try:
        with manifest.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = required.difference(reader.fieldnames or ())
            if missing:
                raise InputValidationError(f"Cohort input manifest is missing columns: {sorted(missing)!r}.")
            for row in reader:
                if row["track"].strip() != track.value:
                    continue
                try:
                    state = SampleCallState(row["state"].strip().upper())
                except ValueError as error:
                    raise InputValidationError(f"Unknown cohort state {row['state']!r} for sample {row['sample']!r}.") from error
                rows.append(CohortSampleInput(
                    sample_id=row["sample"].strip(),
                    state=state,
                    source_path=_optional_manifest_path(row["source_path"], manifest.parent),
                    index_path=_optional_manifest_path(row["index_path"], manifest.parent),
                    source_tool=row["source_tool"].strip() or None,
                    source_version=row["source_version"].strip() or None,
                    reference_build=row["reference_build"].strip() or None,
                    catalog_id=row["catalog_id"].strip() or None,
                ))
    except OSError as error:
        raise InputValidationError(f"Unable to read cohort input manifest '{manifest}': {error}") from error
    result = tuple(rows)
    validate_track_inputs(cohort, result)
    return result


def _open_text(path: Path) -> TextIO:
    return gzip.open(path, "rt", encoding="utf-8", newline="") if str(path).lower().endswith(".gz") else path.open("r", encoding="utf-8", newline="")


def _write_text_no_overwrite(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise OutputValidationError(f"Refusing to overwrite cohort artifact: '{path}'.")
    path.write_text(text, encoding="utf-8", newline="\n")


def _optional_manifest_path(value: str, root: Path) -> Path | None:
    text = value.strip()
    if not text:
        return None
    path = Path(text).expanduser()
    return path if path.is_absolute() else root / path


__all__ = [
    "CohortDefinition", "CohortManifest", "CohortSampleInput", "CohortTrack",
    "CohortTrackResult", "SampleCallState", "scan_multisample_vcf",
    "read_cohort_input_manifest", "validate_track_inputs",
]
