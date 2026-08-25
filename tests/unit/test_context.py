"""Tests for Phase 1.4 run-level AnalysisContext integration."""

from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path

import pytest
import yaml

from hifivar.config import HiFiVarConfig
from hifivar.context import AnalysisContext
from hifivar.exceptions import ConfigurationError, InputValidationError, ReferenceError
from hifivar.reference import ReferenceGenome
from hifivar.sample import InputDataset, InputType, Sample
from hifivar.sample_sheet import SampleRecord, SampleSheet


def write_reference(root: Path, name: str = "reference.fa") -> Path:
    root.mkdir(parents=True, exist_ok=True)
    fasta = root / name
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    return fasta


def write_fastq(root: Path, name: str = "sample.fastq") -> Path:
    root.mkdir(parents=True, exist_ok=True)
    path = root / name
    path.write_text("@read1\nACGT\n+\nIIII\n", encoding="utf-8")
    return path


def write_alignment(root: Path, name: str) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    path = root / name
    path.write_bytes(b"alignment-placeholder")
    return path


def make_reference(tmp_path: Path, *, build: str | None = "GRCh38") -> ReferenceGenome:
    return ReferenceGenome.from_fasta(write_reference(tmp_path), build=build)


def make_sample(tmp_path: Path, sample_id: str = "S1") -> Sample:
    return Sample(
        sample_id,
        InputDataset.from_files((write_fastq(tmp_path, f"{sample_id}.fastq"),)),
    )


def matching_config(reference: ReferenceGenome) -> dict[str, object]:
    return {
        "reference": {"fasta": str(reference.fasta.absolute()), "build": reference.build},
        "samples": {"sheet": None},
        "runtime": {"threads": 1},
    }


def make_context(tmp_path: Path) -> AnalysisContext:
    reference = make_reference(tmp_path)
    return AnalysisContext.from_sample(
        reference,
        make_sample(tmp_path),
        matching_config(reference),
    )


def test_direct_context_accepts_valid_record(tmp_path: Path) -> None:
    reference = make_reference(tmp_path)
    record = SampleRecord(make_sample(tmp_path))
    context = AnalysisContext(reference, (record,), matching_config(reference))
    assert context.samples == (record,)


def test_single_sample_factory_wraps_metadata_empty_record(tmp_path: Path) -> None:
    reference = make_reference(tmp_path)
    context = AnalysisContext.from_sample(
        reference, make_sample(tmp_path), matching_config(reference)
    )
    assert context.samples[0].sex is None
    assert context.source_sample_sheet is None


def test_sample_sheet_factory_preserves_order_and_source(tmp_path: Path) -> None:
    reference = make_reference(tmp_path)
    sheet_path = tmp_path / "samples.tsv"
    sheet_path.write_text("placeholder", encoding="utf-8")
    records = (
        SampleRecord(make_sample(tmp_path, "S2")),
        SampleRecord(make_sample(tmp_path, "S1")),
    )
    sheet = SampleSheet(sheet_path, records)
    context = AnalysisContext.from_sample_sheet(reference, sheet, matching_config(reference))
    assert context.sample_ids == ("S2", "S1")
    assert context.source_sample_sheet == sheet_path


def test_empty_context_is_rejected(tmp_path: Path) -> None:
    reference = make_reference(tmp_path)
    with pytest.raises(InputValidationError, match="at least one"):
        AnalysisContext(reference, (), matching_config(reference))


def test_non_record_context_member_is_rejected(tmp_path: Path) -> None:
    reference = make_reference(tmp_path)
    with pytest.raises(InputValidationError, match="SampleRecord"):
        AnalysisContext(reference, (object(),), matching_config(reference))  # type: ignore[arg-type]


def test_duplicate_sample_id_is_rejected(tmp_path: Path) -> None:
    reference = make_reference(tmp_path)
    records = (
        SampleRecord(make_sample(tmp_path / "a", "S1")),
        SampleRecord(make_sample(tmp_path / "b", "S1")),
    )
    with pytest.raises(InputValidationError, match="duplicate sample_id"):
        AnalysisContext(reference, records, matching_config(reference))


def test_primary_input_reuse_across_samples_is_rejected(tmp_path: Path) -> None:
    reference = make_reference(tmp_path)
    fastq = write_fastq(tmp_path)
    dataset = InputDataset.from_files((fastq,))
    records = (
        SampleRecord(Sample("S1", dataset)),
        SampleRecord(Sample("S2", dataset)),
    )
    with pytest.raises(InputValidationError, match="reused by samples"):
        AnalysisContext(reference, records, matching_config(reference))


def test_mixed_input_types_across_cohort_are_allowed(tmp_path: Path) -> None:
    reference = make_reference(tmp_path)
    fastq_record = SampleRecord(make_sample(tmp_path, "FASTQ1"))
    bam = write_alignment(tmp_path, "aligned.bam")
    bam_record = SampleRecord(Sample("BAM1", InputDataset.from_files((bam,))))
    context = AnalysisContext(
        reference,
        (fastq_record, bam_record),
        matching_config(reference),
    )
    assert context.input_types == (InputType.FASTQ, InputType.BAM)


def test_bam_index_is_not_required_by_generic_context(tmp_path: Path) -> None:
    reference = make_reference(tmp_path)
    bam = write_alignment(tmp_path, "sample.bam")
    context = AnalysisContext.from_sample(
        reference,
        Sample("S1", InputDataset.from_files((bam,))),
        matching_config(reference),
    )
    assert context.n_samples == 1
    assert not Path(f"{bam}.bai").exists()


def test_sample_ids_property_is_deterministic(tmp_path: Path) -> None:
    context = make_context(tmp_path)
    assert context.sample_ids == ("S1",)


def test_n_samples_property(tmp_path: Path) -> None:
    assert make_context(tmp_path).n_samples == 1


def test_input_types_property(tmp_path: Path) -> None:
    assert make_context(tmp_path).input_types == (InputType.FASTQ,)


def test_query_contig_subset_delegates_to_reference(tmp_path: Path) -> None:
    make_context(tmp_path).validate_query_contigs(("chr1",))


def test_query_contig_mismatch_is_explicit(tmp_path: Path) -> None:
    with pytest.raises(ReferenceError, match="REFERENCE_CONTIG_MISMATCH"):
        make_context(tmp_path).validate_query_contigs(("1",))


def test_fastq_reference_status_makes_no_compatibility_claim(tmp_path: Path) -> None:
    status = make_context(tmp_path).reference_compatibility()[0]
    assert status["status"] == "not_applicable"
    assert "unaligned" in status["reason"]


@pytest.mark.parametrize("suffix", ("bam", "cram"))
def test_alignment_reference_status_is_not_checked(tmp_path: Path, suffix: str) -> None:
    reference = make_reference(tmp_path)
    alignment = write_alignment(tmp_path, f"sample.{suffix}")
    context = AnalysisContext.from_sample(
        reference,
        Sample("S1", InputDataset.from_files((alignment,))),
        matching_config(reference),
    )
    status = context.reference_compatibility()[0]
    assert status["status"] == "not_checked"
    assert "header" in status["reason"]


def test_context_detaches_caller_owned_config(tmp_path: Path) -> None:
    reference = make_reference(tmp_path)
    config = matching_config(reference)
    original = deepcopy(config)
    context = AnalysisContext.from_sample(reference, make_sample(tmp_path), config)
    config["runtime"]["threads"] = 99  # type: ignore[index]
    assert context.config_to_dict() == original


def test_config_to_dict_returns_independent_copy(tmp_path: Path) -> None:
    context = make_context(tmp_path)
    exported = context.config_to_dict()
    exported["runtime"]["threads"] = 8  # type: ignore[index]
    assert context.config_to_dict()["runtime"]["threads"] == 1  # type: ignore[index]


def test_context_serialization_uses_standard_types(tmp_path: Path) -> None:
    payload = make_context(tmp_path).to_dict()
    assert json.loads(json.dumps(payload)) == payload
    assert yaml.safe_load(yaml.safe_dump(payload)) == payload


def test_reference_build_alias_matches_config(tmp_path: Path) -> None:
    reference = make_reference(tmp_path, build="GRCh38")
    config = matching_config(reference)
    config["reference"]["build"] = "hg38"  # type: ignore[index]
    AnalysisContext.from_sample(reference, make_sample(tmp_path), config)


def test_reference_build_conflict_is_rejected(tmp_path: Path) -> None:
    reference = make_reference(tmp_path, build="GRCh38")
    config = matching_config(reference)
    config["reference"]["build"] = "GRCh37"  # type: ignore[index]
    with pytest.raises(ReferenceError, match=r"reference\.build conflicts"):
        AnalysisContext.from_sample(reference, make_sample(tmp_path), config)


def test_reference_fasta_conflict_is_rejected(tmp_path: Path) -> None:
    reference = make_reference(tmp_path)
    config = matching_config(reference)
    config["reference"]["fasta"] = str((tmp_path / "other.fa").absolute())  # type: ignore[index]
    with pytest.raises(ReferenceError, match=r"reference\.fasta conflicts"):
        AnalysisContext.from_sample(reference, make_sample(tmp_path), config)


def write_sheet(root: Path, input_value: str = "sample.fastq") -> Path:
    sheet = root / "samples.tsv"
    sheet.write_text(f"sample_id\tinput\nS1\t{input_value}\n", encoding="utf-8")
    return sheet


def test_from_config_builds_reference_and_sample_sheet(tmp_path: Path) -> None:
    fasta = write_reference(tmp_path)
    write_fastq(tmp_path)
    sheet = write_sheet(tmp_path)
    config = {
        "reference": {"fasta": str(fasta.absolute()), "build": "GRCh38"},
        "samples": {"sheet": str(sheet.absolute())},
    }
    context = AnalysisContext.from_config(config)
    assert context.reference.build == "GRCh38"
    assert context.sample_ids == ("S1",)


def test_from_config_requires_reference_fasta(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match=r"reference\.fasta"):
        AnalysisContext.from_config(
            {"reference": {"fasta": None}, "samples": {"sheet": "x.tsv"}}
        )


def test_from_config_requires_samples_sheet(tmp_path: Path) -> None:
    fasta = write_reference(tmp_path)
    with pytest.raises(ConfigurationError, match=r"samples\.sheet"):
        AnalysisContext.from_config(
            {"reference": {"fasta": str(fasta.absolute())}, "samples": {"sheet": None}}
        )


def test_from_config_propagates_bad_sample_sheet(tmp_path: Path) -> None:
    fasta = write_reference(tmp_path)
    bad_sheet = tmp_path / "bad.tsv"
    bad_sheet.write_text("wrong\theader\nvalue\tvalue\n", encoding="utf-8")
    config = {
        "reference": {"fasta": str(fasta.absolute()), "build": None},
        "samples": {"sheet": str(bad_sheet.absolute())},
    }
    with pytest.raises(InputValidationError, match="missing required"):
        AnalysisContext.from_config(config)


def test_relative_analysis_paths_resolve_beside_user_config(tmp_path: Path) -> None:
    project = tmp_path / "project"
    fasta = write_reference(project, "ref.fa")
    write_fastq(project)
    write_sheet(project)
    user_config = project / "config.yaml"
    user_config.write_text("placeholder", encoding="utf-8")
    config = HiFiVarConfig(
        {
            "reference": {"fasta": "ref.fa", "build": None},
            "samples": {"sheet": "samples.tsv"},
        },
        {"user": user_config},
    )
    context = AnalysisContext.from_config(config)
    assert context.reference.fasta == fasta
    assert context.source_sample_sheet == project / "samples.tsv"


def test_relative_analysis_path_without_source_is_rejected() -> None:
    with pytest.raises(ConfigurationError, match="Relative reference.fasta"):
        AnalysisContext.from_config(
            {
                "reference": {"fasta": "ref.fa", "build": None},
                "samples": {"sheet": "samples.tsv"},
            }
        )


def test_wrong_reference_type_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(ReferenceError, match="ReferenceGenome"):
        AnalysisContext(object(), (SampleRecord(make_sample(tmp_path)),), {})  # type: ignore[arg-type]


def test_wrong_config_type_is_rejected(tmp_path: Path) -> None:
    reference = make_reference(tmp_path)
    with pytest.raises(ConfigurationError, match="mapping"):
        AnalysisContext(reference, (SampleRecord(make_sample(tmp_path)),), [])  # type: ignore[arg-type]
