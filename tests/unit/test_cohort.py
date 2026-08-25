from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest
from importlib.resources import files

from hifivar.cohort import (
    CohortDefinition, CohortManifest, CohortSampleInput, CohortTrack,
    CohortTrackResult, SampleCallState, read_cohort_input_manifest,
    scan_multisample_vcf,
)
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.reference import Contig, ReferenceGenome


def reference(tmp_path: Path) -> ReferenceGenome:
    return ReferenceGenome(tmp_path / "参考.fa", tmp_path / "参考.fa.fai", "GRCh38", (Contig("chr1", 1000),), "a" * 64)


def cohort(tmp_path: Path) -> CohortDefinition:
    return CohortDefinition("C1", ("S2", "S1"), reference(tmp_path))


def test_cohort_identity_order_and_duplicate_safety(tmp_path: Path) -> None:
    model = cohort(tmp_path)
    assert model.sample_ids == ("S2", "S1")
    with pytest.raises(InputValidationError, match="duplicate"):
        CohortDefinition("C1", ("S1", "S1"), reference(tmp_path))
    with pytest.raises(InputValidationError, match="Invalid sample_id"):
        CohortDefinition("bad id", ("S1",), reference(tmp_path))


@pytest.mark.parametrize("state", [SampleCallState.NOT_RUN, SampleCallState.FAILED, SampleCallState.DISABLED, SampleCallState.MISSING_INPUT])
def test_non_callable_states_are_distinct_from_no_calls(state: SampleCallState) -> None:
    item = CohortSampleInput("S1", state)
    assert not item.callable
    assert item.state is not SampleCallState.NO_CALLS


def test_called_and_no_calls_require_artifact() -> None:
    for state in (SampleCallState.CALLED, SampleCallState.NO_CALLS):
        with pytest.raises(InputValidationError, match="requires a source"):
            CohortSampleInput("S1", state)


def test_manifest_reader_preserves_explicit_order_and_states(tmp_path: Path) -> None:
    path = tmp_path / "cohort.tsv"
    path.write_text(
        "sample\ttrack\tstate\tsource_path\tindex_path\tsource_tool\tsource_version\treference_build\tcatalog_id\n"
        "S2\tsv\tCALLED\t二.vcf.gz\t二.vcf.gz.tbi\tjasmine\t1.1.5\tGRCh38\t\n"
        "S1\tsv\tFAILED\t\t\tjasmine\t1.1.5\tGRCh38\t\n",
        encoding="utf-8",
    )
    rows = read_cohort_input_manifest(path, cohort(tmp_path), CohortTrack.SV)
    assert tuple(item.sample_id for item in rows) == ("S2", "S1")
    assert rows[0].source_path == tmp_path / "二.vcf.gz"
    assert rows[1].state is SampleCallState.FAILED


def test_manifest_reader_never_silently_drops_sample(tmp_path: Path) -> None:
    path = tmp_path / "cohort.tsv"
    path.write_text("sample\ttrack\tstate\tsource_path\tindex_path\tsource_tool\tsource_version\treference_build\tcatalog_id\nS2\tsv\tFAILED\t\t\t\t\tGRCh38\t\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="order/set mismatch"):
        read_cohort_input_manifest(path, cohort(tmp_path), CohortTrack.SV)


def test_streaming_qc_distinguishes_missing_homref_and_nonref(tmp_path: Path) -> None:
    path = tmp_path / "cohort.vcf.gz"
    text = (
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\tS1\n"
        "chr1\t10\tv1\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\n"
        "chr1\t20\tv2\tA\tC,G\t.\tLowQual\t.\tGT\t./.\t1/2\n"
    )
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)
    qc = scan_multisample_vcf(path, ("S2", "S1"))
    assert qc["variant_count"] == 2
    assert qc["multiallelic_count"] == 1
    assert qc["per_sample_non_ref_count"] == {"S2": 1, "S1": 1}
    assert qc["per_sample_missing_rate"]["S2"] == 0.5
    assert qc["per_sample_call_rate"]["S1"] == 1.0


def test_vcf_sample_set_must_match_exactly(tmp_path: Path) -> None:
    path = tmp_path / "bad.vcf"
    path.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n", encoding="utf-8")
    with pytest.raises(OutputValidationError, match="order/set mismatch"):
        scan_multisample_vcf(path, ("S2", "S1"))


def test_manifest_json_yaml_and_overwrite_protection(tmp_path: Path) -> None:
    track = CohortTrackResult(CohortTrack.SV, False, SampleCallState.DISABLED, None, None)
    manifest = CohortManifest(cohort(tmp_path), (track,), "0.0.1")
    json_path, yaml_path = tmp_path / "manifest.json", tmp_path / "manifest.yaml"
    manifest.write(json_path, yaml_path)
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert payload["tracks"][0]["state"] == "DISABLED"
    assert len(payload["sample_order_sha256"]) == 64
    assert payload["created_at"].endswith("+00:00")
    with pytest.raises(OutputValidationError, match="overwrite"):
        manifest.write(json_path)


def test_packaged_default_contains_phase12_config() -> None:
    text = files("hifivar").joinpath("resources/configs/default.yaml").read_text(encoding="utf-8")
    assert "cohort:" in text and "glnexus_executable: glnexus_cli" in text
