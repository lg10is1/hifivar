from __future__ import annotations

import json
from pathlib import Path

from hifivar.low_depth import (
    AccuracyThreshold,
    LowDepthObservation,
    LowDepthRunStatus,
    ParameterProfile,
    build_low_depth_plan,
    evaluate_low_depth_study,
)


def test_phase15_plan_to_replicated_result_smoke(tmp_path: Path) -> None:
    plan = build_low_depth_plan(
        study_id="tiny-depth-study",
        sample_id="HG002",
        source_bam=tmp_path / "HG002.bam",
        baseline_coverage=20,
        target_coverages=(5, 10, 20),
        replicates=3,
        base_seed=42,
        profiles=(ParameterProfile("baseline"), ParameterProfile("candidate", {"mode": "explicit"})),
        thresholds=(AccuracyThreshold("SNP", "f1", 0.98),),
        output_root=tmp_path / "subsamples",
        minimum_successful_replicates=3,
    )
    observations = tuple(
        LowDepthObservation(
            run.run_id,
            run.profile.name,
            run.target_coverage,
            run.replicate,
            LowDepthRunStatus.COMPLETED,
            {"SNP.f1": 0.97 if run.target_coverage == 5 else 0.99},
            run.target_coverage * 0.99,
        )
        for run in plan.runs
    )
    paired = {}
    for run in plan.runs:
        paired.setdefault((run.target_coverage, run.replicate), set()).add(
            (run.seed, run.output_bam)
        )
    assert all(len(values) == 1 for values in paired.values())
    assert plan.to_dict()["scientific_semantics"]["subsamples_are_paired_across_profiles"] is True
    nested = {}
    for run in plan.runs:
        nested.setdefault(run.replicate, set()).add(run.seed)
    assert all(len(values) == 1 for values in nested.values())
    assert len({run.seed for run in plan.runs}) == 3
    assert plan.to_dict()["scientific_semantics"]["subsamples_are_nested_across_depths"] is True
    result = evaluate_low_depth_study(plan, observations)
    json_path, yaml_path = result.write(tmp_path / "result")
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert payload["minimum_supported_depth_by_profile"] == {"baseline": 10, "candidate": 10}
    assert payload["plan"]["thresholds"] == [
        {"variant_class": "SNP", "metric": "f1", "minimum": 0.98}
    ]
    assert len(payload["observations"]) == len(plan.runs)
    assert payload["assessments"][0]["observed_coverage_range"] == [4.95, 4.95]
    assert payload["scientific_semantics"]["result_is_not_a_universal_minimum_depth"] is True
    assert yaml_path.is_file()
    assert not plan.source_bam.exists()
