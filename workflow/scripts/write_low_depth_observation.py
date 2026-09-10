"""Combine measured coverage and hap.py metrics into one run observation."""

import json
from pathlib import Path

from hifivar.happy import parse_happy_summary
from hifivar.low_depth import LowDepthObservation, LowDepthRunStatus
from hifivar.serialization import write_json_atomic

coverage = json.loads(Path(str(snakemake.input.coverage)).read_text(encoding="utf-8"))  # type: ignore[name-defined]
metrics = parse_happy_summary(Path(str(snakemake.input.summary)))  # type: ignore[name-defined]
values = {
    f"{metric.stratum}.{metric.name}": float(metric.value)
    for metric in metrics if metric.value is not None
}
observation = LowDepthObservation(
    run_id=str(snakemake.params.run_id),  # type: ignore[name-defined]
    profile="baseline",
    target_coverage=float(snakemake.params.target_coverage),  # type: ignore[name-defined]
    replicate=int(snakemake.wildcards.replicate),  # type: ignore[name-defined]
    status=LowDepthRunStatus.COMPLETED,
    metrics=values,
    observed_coverage=float(coverage["observed_coverage"]),
)
payload = observation.to_dict()
payload["scientific_semantics"] = {
    "observed_coverage_is_measured": True,
    "accuracy_is_truth_set_specific": True,
    "minimum_depth_is_not_claimed_by_one_observation": True,
}
write_json_atomic(
    payload, Path(str(snakemake.output.json)),  # type: ignore[name-defined]
    overwrite=bool(snakemake.config["low_depth"].get("overwrite", False)),  # type: ignore[name-defined]
    artifact_name="Low-depth observation",
)
