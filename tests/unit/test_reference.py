"""Tests for the Phase 1.1 reference-genome data model."""

from __future__ import annotations

import gzip
import hashlib
import json
from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest
import yaml

from hifivar.exceptions import ReferenceError
from hifivar.reference import Contig, ReferenceGenome


VALID_FASTA = ">chr1\nACGT\n>chr2\nAAAAAA\n"
VALID_FAI = "chr1\t4\t6\t4\t5\nchr2\t6\t17\t6\t7\n"


def write_reference(
    root: Path,
    *,
    fasta_text: str = VALID_FASTA,
    fai_text: str = VALID_FAI,
    name: str = "reference.fa",
) -> Path:
    """Write a tiny reference and standard five-column FAI fixture."""
    root.mkdir(parents=True, exist_ok=True)
    fasta = root / name
    fasta.write_text(fasta_text, encoding="utf-8")
    Path(f"{fasta}.fai").write_text(fai_text, encoding="utf-8")
    return fasta


def test_valid_reference_loads_from_fasta_and_fai(tmp_path: Path) -> None:
    """A valid indexed FASTA should produce immutable reference metadata."""
    fasta = write_reference(tmp_path)

    reference = ReferenceGenome.from_fasta(fasta)

    assert reference.fasta == fasta
    assert reference.fai == Path(f"{fasta}.fai")
    assert reference.build is None
    assert reference.sha256 is None


def test_missing_fasta_raises_reference_error(tmp_path: Path) -> None:
    """Factory failures use the reference-specific exception boundary."""
    missing = tmp_path / "missing.fa"

    with pytest.raises(ReferenceError, match=r"Invalid reference FASTA.*missing"):
        ReferenceGenome.from_fasta(missing)


def test_missing_fai_reports_remediation_without_creating_it(
    tmp_path: Path,
) -> None:
    """FAI generation remains an explicit external preparation step."""
    fasta = tmp_path / "reference.fa"
    fasta.write_text(VALID_FASTA, encoding="utf-8")
    expected_fai = Path(f"{fasta}.fai")

    with pytest.raises(
        ReferenceError,
        match=r"FASTA index missing.*Create FASTA index before using",
    ):
        ReferenceGenome.from_fasta(fasta)

    assert not expected_fai.exists()


def test_empty_fasta_raises_reference_error(tmp_path: Path) -> None:
    """An index cannot make an empty FASTA a valid reference."""
    fasta = write_reference(tmp_path, fasta_text="")

    with pytest.raises(ReferenceError, match=r"Invalid reference FASTA.*empty"):
        ReferenceGenome.from_fasta(fasta)


def test_contigs_and_lengths_come_from_fai(tmp_path: Path) -> None:
    """The model should not rescan FASTA sequences to derive lengths."""
    reference = ReferenceGenome.from_fasta(write_reference(tmp_path))

    assert reference.contigs == (
        Contig(name="chr1", length=4),
        Contig(name="chr2", length=6),
    )
    assert reference.get_contig("chr1").length == 4


def test_contig_order_is_preserved_from_fai(tmp_path: Path) -> None:
    """FAI insertion order should remain deterministic in the model."""
    fasta = write_reference(
        tmp_path,
        fai_text=(
            "chr1\t4\t6\t4\t5\n"
            "chr2\t6\t17\t6\t7\n"
            "chrM\t2\t30\t2\t3\n"
        ),
    )

    reference = ReferenceGenome.from_fasta(fasta)

    assert reference.contig_names == ("chr1", "chr2", "chrM")


def test_duplicate_fai_contig_is_rejected(tmp_path: Path) -> None:
    """Duplicate contig names must fail before a model is returned."""
    fasta = write_reference(
        tmp_path,
        fai_text="chr1\t4\t6\t4\t5\nchr1\t6\t17\t6\t7\n",
    )

    with pytest.raises(ReferenceError, match=r"duplicate contig.*chr1"):
        ReferenceGenome.from_fasta(fasta)


@pytest.mark.parametrize("length", ("0", "-1", "abc"))
def test_invalid_fai_contig_length_is_rejected(
    tmp_path: Path,
    length: str,
) -> None:
    """Zero, negative, and non-integer FAI lengths are invalid."""
    fasta = write_reference(
        tmp_path,
        fai_text=f"chr1\t{length}\t6\t4\t5\n",
    )

    with pytest.raises(ReferenceError, match=r"coordinate|non-integer"):
        ReferenceGenome.from_fasta(fasta)


def test_explicit_build_is_saved_without_inference(tmp_path: Path) -> None:
    """Reference build comes only from explicit metadata."""
    reference = ReferenceGenome.from_fasta(
        write_reference(tmp_path),
        build="GRCh38",
    )

    assert reference.build == "GRCh38"


def test_build_none_is_allowed(tmp_path: Path) -> None:
    """An unspecified build should remain unknown rather than guessed."""
    reference = ReferenceGenome.from_fasta(write_reference(tmp_path))

    assert reference.build is None


@pytest.mark.parametrize(
    "alias,canonical",
    (("hg38", "GRCh38"), ("grch38", "GRCh38"), ("hg19", "GRCh37")),
)
def test_unambiguous_build_aliases_are_canonicalized(
    tmp_path: Path,
    alias: str,
    canonical: str,
) -> None:
    """Only an explicit alias table should canonicalize build labels."""
    reference = ReferenceGenome.from_fasta(
        write_reference(tmp_path),
        build=alias,
    )

    assert reference.build == canonical


def test_custom_build_is_preserved(tmp_path: Path) -> None:
    """Non-human and project-specific build labels remain supported."""
    reference = ReferenceGenome.from_fasta(
        write_reference(tmp_path),
        build="my_species_v1",
    )

    assert reference.build == "my_species_v1"


def test_checksum_is_not_computed_by_default(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Ordinary construction must not stream the complete reference."""
    fasta = write_reference(tmp_path)

    def fail_if_called(path: object) -> str:
        raise AssertionError(f"checksum unexpectedly called for {path}")

    monkeypatch.setattr(
        "hifivar.reference.validation.compute_sha256",
        fail_if_called,
    )

    reference = ReferenceGenome.from_fasta(fasta)

    assert reference.sha256 is None


def test_checksum_can_be_computed_explicitly(tmp_path: Path) -> None:
    """Explicit checksum mode should use the Phase 0 streaming helper."""
    fasta = write_reference(tmp_path)

    reference = ReferenceGenome.from_fasta(fasta, compute_checksum=True)

    assert reference.sha256 == hashlib.sha256(fasta.read_bytes()).hexdigest()


def test_with_checksum_returns_new_immutable_reference(tmp_path: Path) -> None:
    """Lazy provenance enrichment should not mutate an existing instance."""
    reference = ReferenceGenome.from_fasta(write_reference(tmp_path))

    checksummed = reference.with_checksum()

    assert reference.sha256 is None
    assert checksummed.sha256 == hashlib.sha256(
        reference.fasta.read_bytes()
    ).hexdigest()
    assert checksummed is not reference


def test_provided_checksum_is_recorded_without_rescanning(tmp_path: Path) -> None:
    """A future manifest checksum can be attached without reading the FASTA."""
    reference = ReferenceGenome.from_fasta(
        write_reference(tmp_path),
        sha256="A" * 64,
    )

    assert reference.sha256 == "a" * 64


def test_mismatched_provided_checksum_is_rejected(tmp_path: Path) -> None:
    """Explicit verification should detect stale manifest metadata."""
    fasta = write_reference(tmp_path)

    with pytest.raises(ReferenceError, match=r"SHA256 does not match"):
        ReferenceGenome.from_fasta(
            fasta,
            sha256="0" * 64,
            compute_checksum=True,
        )


def test_summary_serialization_is_json_and_yaml_friendly(tmp_path: Path) -> None:
    """Default provenance metadata should contain only standard data types."""
    reference = ReferenceGenome.from_fasta(
        write_reference(tmp_path),
        build="GRCh38",
    )

    metadata = reference.to_dict()

    assert metadata == {
        "build": "GRCh38",
        "fasta": str(reference.fasta),
        "fai": str(reference.fai),
        "sha256": None,
        "contig_count": 2,
    }
    assert json.loads(json.dumps(metadata)) == metadata
    assert yaml.safe_load(yaml.safe_dump(metadata)) == metadata


def test_full_serialization_includes_plain_contig_mappings(tmp_path: Path) -> None:
    """Detailed metadata remains serializable when explicitly requested."""
    reference = ReferenceGenome.from_fasta(write_reference(tmp_path))

    metadata = reference.to_dict(include_contigs=True)

    assert metadata["contigs"] == [
        {"name": "chr1", "length": 4},
        {"name": "chr2", "length": 6},
    ]
    json.dumps(metadata)
    yaml.safe_dump(metadata)


def test_reference_contig_subset_is_compatible(tmp_path: Path) -> None:
    """Future query files may contain a subset of reference contigs."""
    reference = ReferenceGenome.from_fasta(write_reference(tmp_path))

    reference.validate_contigs(["chr1"])


def test_reference_contig_mismatch_is_not_normalized(tmp_path: Path) -> None:
    """chr-prefixed and unprefixed names remain scientifically distinct."""
    reference = ReferenceGenome.from_fasta(write_reference(tmp_path))

    with pytest.raises(ReferenceError, match=r"REFERENCE_CONTIG_MISMATCH.*\b1\b"):
        reference.validate_contigs(["1"])

    assert reference.contig_names[0] == "chr1"


def test_reference_and_contig_are_immutable(tmp_path: Path) -> None:
    """Callers cannot change build or contig metadata during a run."""
    reference = ReferenceGenome.from_fasta(write_reference(tmp_path))

    with pytest.raises(FrozenInstanceError):
        reference.build = "GRCh37"  # type: ignore[misc]
    with pytest.raises(FrozenInstanceError):
        reference.contigs[0].name = "1"  # type: ignore[misc]


def test_unicode_reference_path_is_supported(tmp_path: Path) -> None:
    """Reference loading should work in Unicode directories on both platforms."""
    fasta = write_reference(tmp_path / "参考基因组")

    reference = ReferenceGenome.from_fasta(fasta)

    assert reference.fasta == fasta
    assert reference.contig_names == ("chr1", "chr2")


def test_compressed_reference_is_explicitly_unsupported(tmp_path: Path) -> None:
    """Phase 0 gzip validation does not imply workflow-reference support."""
    fasta = tmp_path / "reference.fa.gz"
    with gzip.open(fasta, "wt", encoding="utf-8", newline="") as handle:
        handle.write(VALID_FASTA)

    with pytest.raises(ReferenceError, match=r"Compressed FASTA.*unsupported"):
        ReferenceGenome.from_fasta(fasta)


@pytest.mark.parametrize("suffix", (".fa", ".fasta", ".fna"))
def test_supported_uncompressed_fasta_suffixes(
    tmp_path: Path,
    suffix: str,
) -> None:
    """All documented primary-reference suffixes should load."""
    fasta = write_reference(tmp_path, name=f"reference{suffix}")

    assert ReferenceGenome.from_fasta(fasta).fasta == fasta


def test_invalid_contig_model_is_rejected() -> None:
    """Direct metadata construction retains positive-length invariants."""
    with pytest.raises(ReferenceError, match=r"positive integer"):
        Contig(name="chr1", length=0)
