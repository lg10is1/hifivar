"""Native-assignment semantics with fake optional modules, no matching claims."""

from __future__ import annotations

import json

import pytest

from hifivar import sv_benchmark_stratification as stratification
from hifivar.exceptions import InputValidationError, OutputValidationError
from test_sv_benchmark_selection import make_vcf, row, runtime  # noqa: F401


EXTRA = ('##INFO=<ID=IsRefined,Number=0,Type=Flag,Description="Origin">\n'
         '##FORMAT=<ID=BD,Number=1,Type=String,Description="Decision">\n'
         '##FORMAT=<ID=BK,Number=1,Type=String,Description="Match">\n')


def assigned(pos, vid, kind, size, bd, bk="gm", flag=True):
    info = f"SVTYPE={kind};SVLEN={size}" + (";IsRefined" if flag else "")
    return row(pos, vid, f"<{kind}>", info, fmt="GT:BD:BK", value=f"0|1:{bd}:{bk}")


def native_run(tmp_path, base_body, comp_body, counts, *, sample="HG002", extra=EXTRA):
    run = tmp_path / "bench"
    run.mkdir()
    make_vcf(run, "refine.base.vcf.gz", body=base_body, sample=sample, extra=extra)
    make_vcf(run, "refine.comp.vcf.gz", body=comp_body, sample=sample, extra=extra)
    (run / "params.json").write_text(json.dumps({"bSample": sample, "cSample": sample, "dup_to_ins": False}), encoding="utf-8")
    rates = stratification._rates(counts)
    payload = {**rates["counts"], **{k: rates[k] for k in ("precision", "recall", "f1")}}
    (run / "refine.variant_summary.json").write_text(json.dumps(payload), encoding="utf-8")
    # Wrong initial-bench files are deliberately irrelevant.
    (run / "tp-base.vcf.gz").write_bytes(b"not used")
    return run


@pytest.fixture(autouse=True)
def inject_runtime(monkeypatch, runtime):
    monkeypatch.setattr(stratification, "_load_native_runtime", lambda: runtime)


def test_refined_independent_tp_counts_short_query_and_strata(tmp_path):
    base = assigned(10, "b1", "INS", 60, "TP") + assigned(20, "b2", "DEL", -500, "FN", ".", False)
    comp = (assigned(10, "c1", "INS", 40, "TP", "am") + assigned(15, "c2", "INS", 60, "TP")
            + assigned(30, "c3", "DEL", -500, "FP", "lm", False))
    run = native_run(tmp_path, base, comp, {"TP-base": 1, "TP-comp": 2, "FN": 1, "FP": 1})
    result = stratification.stratify_refined_sv_benchmark(run, tmp_path / "new", sample_id="HG002")
    assert result["global_reconciliation"] == "PASS_EXACT_NATIVE_COUNTS"
    assert result["global_metrics"]["counts"] == {"TP-base": 1, "TP-comp": 2, "FN": 1, "FP": 1, "base cnt": 2, "comp cnt": 3}
    bins = {item["stratum"]: item for item in result["strata"]}
    assert bins["SIZE:30-49bp"]["counts"]["TP-comp"] == 1
    assert bins["SIZE:30-49bp"]["precision"] is None  # empty truth side
    assert bins["SVTYPE:INS"]["precision"] == 1
    assert result["assignment_provenance_counts"]["comparison:BK=am"] == 1


@pytest.mark.parametrize("base,comp,counts,expected", [
    ("", "", {}, (None, None, None)),
    ("", assigned(10, "fp", "INS", 50, "FP", "."), {"FP": 1}, (None, None, None)),
    (assigned(10, "fn", "DEL", -50, "FN", "."), assigned(20, "fp", "INS", 50, "FP", "."), {"FN": 1, "FP": 1}, (0, 0, None)),
])
def test_native_null_rules(tmp_path, base, comp, counts, expected):
    run = native_run(tmp_path, base, comp, counts)
    result = stratification.stratify_refined_sv_benchmark(run, tmp_path / "new", sample_id="HG002")
    assert tuple(result["global_metrics"][name] for name in ("precision", "recall", "f1")) == expected


def test_reconciliation_failure_does_not_publish_strata(tmp_path):
    run = native_run(tmp_path, assigned(1, "b", "INS", 60, "TP"), "", {})
    with pytest.raises(OutputValidationError, match="RECONCILIATION_FAIL"):
        stratification.stratify_refined_sv_benchmark(run, tmp_path / "new", sample_id="HG002")
    failed = json.loads((tmp_path / "new/refined_stratification.json").read_text())
    assert failed["status"] == "FAILED" and "strata" not in failed


@pytest.mark.parametrize("side,bd,bk", [("base", "FP", "lm"), ("comparison", "FN", "lm"), ("base", "TP", "lm"), ("comparison", "FP", "gm"), ("base", "N", ".")])
def test_invalid_side_or_bk_fails_closed(tmp_path, side, bd, bk):
    body = assigned(10, "bad", "INS", 60, bd, bk)
    run = native_run(tmp_path, body if side == "base" else "", body if side == "comparison" else "", {})
    with pytest.raises(OutputValidationError, match="BD/BK"):
        stratification.stratify_refined_sv_benchmark(run, tmp_path / "new", sample_id="HG002")


def test_explicit_sample_mapping_must_match_params_and_headers(tmp_path):
    run = native_run(tmp_path, "", "", {}, sample="not_HG002")
    with pytest.raises(InputValidationError, match="bSample/cSample"):
        stratification.stratify_refined_sv_benchmark(run, tmp_path / "new", sample_id="HG002")


def test_no_initial_bench_vcf_substitution(tmp_path):
    run = native_run(tmp_path, "", "", {})
    (run / "refine.comp.vcf.gz").unlink()
    with pytest.raises(InputValidationError, match="missing"):
        stratification.stratify_refined_sv_benchmark(run, tmp_path / "new", sample_id="HG002")
    assert not (tmp_path / "new").exists()


def test_bnd_never_gets_a_biologically_false_length_bin(tmp_path):
    run = native_run(tmp_path, "", assigned(10, "bnd", "BND", 60, "FP", "."), {"FP": 1})
    result = stratification.stratify_refined_sv_benchmark(run, tmp_path / "new", sample_id="HG002")
    bins = {item["stratum"]: item for item in result["strata"]}
    assert bins["SIZE:UNCLASSIFIED_NON_INS_DEL"]["counts"]["FP"] == 1
    assert "SIZE:50-99bp" not in bins


@pytest.mark.parametrize("boundaries", [(50, 100), (30, 50, 50), (30, False, 50), (50, 30), (0, 30, 50)])
def test_size_bins_require_preserved_30_49_bucket(tmp_path, boundaries):
    with pytest.raises(InputValidationError, match="boundaries"):
        stratification.stratify_refined_sv_benchmark(tmp_path, tmp_path / "new", sample_id="HG002", size_boundaries=boundaries)


def test_dry_run_without_import_or_write(tmp_path, monkeypatch):
    run = native_run(tmp_path, "", "", {})
    monkeypatch.setattr(stratification, "_load_native_runtime", lambda: pytest.fail("dryrun imported native dependency"))
    result = stratification.stratify_refined_sv_benchmark(run, tmp_path / "new", sample_id="HG002", dry_run=True)
    assert result["status"] == "PLANNED" and not (tmp_path / "new").exists()


def test_native_assignment_headers_required_not_guessed(tmp_path):
    run = native_run(tmp_path, "", "", {}, extra="")
    with pytest.raises(OutputValidationError, match="FORMAT/BD"):
        stratification.stratify_refined_sv_benchmark(run, tmp_path / "new", sample_id="HG002")


@pytest.mark.parametrize("size,expected", [(29, "0-29bp"), (30, "30-49bp"), (49, "30-49bp"), (50, "50-99bp"), (50000, "10000-50000bp"), (50001, ">=50001bp")])
def test_size_boundaries_are_explicit_and_closed_without_refiltering(size, expected):
    assert stratification._bin(size, stratification.DEFAULT_SIZE_BOUNDARIES) == expected
