from __future__ import annotations

from dataclasses import replace
import gzip
from pathlib import Path
import struct
import zlib

import pytest

from hifivar.assembly_sv import SVEvidenceSource
from hifivar.command import CommandResult
from hifivar.exceptions import CommandExecutionError, InputValidationError, OutputValidationError, ToolVersionError
from hifivar.harmonization import EvidenceRunStatus, HarmonizedSvArtifact, SVEvidenceSourceArtifact, SVHarmonizationRequest, write_evidence_table
from hifivar.reference import ReferenceGenome
from hifivar.sv_support_selection import SupportSelectionRequest, SupportSelectionWrapper, support_selection_paths


def bgzf(path: Path, payload: bytes):
    compressor = zlib.compressobj(level=6, wbits=-15)
    data = compressor.compress(payload) + compressor.flush()
    size = 18 + len(data) + 8
    header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", size - 1)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(header + data + struct.pack("<II", zlib.crc32(payload), len(payload)))


HEADER = (
    '##fileformat=VCFv4.4\n##contig=<ID=chr1,length=1000>\n'
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
    '##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Sources">\n'
    '##INFO=<ID=SUPP,Number=1,Type=Integer,Description="Count">\n'
    '##INFO=<ID=IDLIST,Number=.,Type=String,Description="IDs">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
).encode()


def record(position, vector, *, name="same", kind="DEL", alt="<DEL>"):
    n = vector.count("1")
    return (f"chr1\t{position}\t{name}\tA\t{alt}\t17\tLowQual\t"
            f"SVTYPE={kind};END={position + 10};SUPP_VEC={vector};SUPP={n};IDLIST=x,y\tGT\t0|1\n").encode()


@pytest.fixture
def full(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + "A" * 1000 + "\n")
    Path(f"{fasta}.fai").write_text("chr1\t1000\t6\t1000\t1001\n")
    reference = ReferenceGenome.from_fasta(fasta, build="GRCh38")
    sources = []
    for caller, kind in (("sawfish", SVEvidenceSource.READ), ("sniffles2", SVEvidenceSource.READ), ("pav", SVEvidenceSource.ASSEMBLY)):
        path = tmp_path / f"{caller}.sv.vcf.gz"
        bgzf(path, HEADER)
        bgzf(Path(f"{path}.tbi"), b"TBI\x01")
        sources.append(SVEvidenceSourceArtifact("S1", kind, caller, path, Path(f"{path}.tbi"), EvidenceRunStatus.COMPLETED))
    directory = tmp_path / "native" / "S1"
    vcf = directory / "S1.harmonized.sv.vcf.gz"
    bgzf(vcf, HEADER + record(10, "100") + record(20, "110") + record(30, "111"))
    bgzf(Path(f"{vcf}.tbi"), b"TBI\x01")
    table = directory / "S1.sv.evidence.tsv"
    request = SVHarmonizationRequest("S1", reference, tuple(sources), tmp_path / "old_work", vcf, table, overwrite=True)
    write_evidence_table(request, vcf)
    artifact = HarmonizedSvArtifact("S1", reference.fasta, reference.build, vcf, Path(f"{vcf}.tbi"), table, tuple(sources), "1.1.5", (("jasmine", "max_dist=1000"),))
    return artifact, reference


def request(full, tmp_path, k=2):
    artifact, reference = full
    return SupportSelectionRequest(artifact, reference, k, support_selection_paths(tmp_path / "out", "S1", k)["directory"])


class FakeRunner:
    def __init__(self, *, fail=None, version="1.24"):
        self.calls = []
        self.fail = fail
        self.version = version

    def run(self, command, **kwargs):
        args = tuple(str(value) for value in command)
        self.calls.append((args, kwargs))
        dry = kwargs.get("dry_run", False)
        if dry:
            return CommandResult(args, None, None, None, 0.0, None, False)
        if args[-1] == "--version":
            return CommandResult(args, 0, f"{args[0]} (htslib) {self.version}", "", 0.01, None, True)
        if self.fail and args[:2] == self.fail:
            raise CommandExecutionError("synthetic native failure")
        if args[:2] == ("bgzip", "-c"):
            bgzf(Path(kwargs["stdout_path"]), Path(args[-1]).read_bytes())
        elif args[:2] == ("bgzip", "-t"):
            with gzip.open(args[-1], "rb") as handle:
                handle.read()
        elif args[:2] == ("tabix", "-p"):
            bgzf(Path(f"{args[-1]}.tbi"), b"TBI\x01")
        return CommandResult(args, 0, "chr1\n" if args[:2] == ("tabix", "-l") else "", "", 0.01, None, True)


@pytest.mark.parametrize("k", [0, -1, True, False, 1.0, "2", None])
def test_support_threshold_rejects_non_positive_integer(full, tmp_path, k):
    with pytest.raises(InputValidationError):
        request(full, tmp_path, k)


def test_threshold_counts_only_actual_runnable_sources(full, tmp_path):
    artifact, reference = full
    disabled = SVEvidenceSourceArtifact("S1", SVEvidenceSource.READ, "pbsv", None, None, EvidenceRunStatus.DISABLED)
    artifact = replace(artifact, sources=(*artifact.sources, disabled))
    with pytest.raises(InputValidationError, match="actual Jasmine"):
        SupportSelectionRequest(artifact, reference, 4, tmp_path / "out" / "S1" / "support-ge-4")
    artifact = replace(artifact, sources=tuple(replace(source, status=EvidenceRunStatus.NO_CALLS) for source in artifact.sources[:3]))
    item = SupportSelectionRequest(artifact, reference, 3, tmp_path / "out" / "S1" / "support-ge-3")
    assert len(item.runnable_sources) == 3


def test_selected_view_preserves_full_source_and_exact_retained_rows(full, tmp_path):
    item = request(full, tmp_path)
    original_vcf = item.artifact.vcf_path.read_bytes()
    original_table = item.artifact.evidence_table.read_bytes()
    runner = FakeRunner()
    result = SupportSelectionWrapper(runner=runner).run(item, stderr_path=tmp_path / "logs" / "selected.log")
    assert result.status == "completed"
    assert result.selected_records == 2 and result.source_records == 3
    assert gzip.decompress(result.vcf_path.read_bytes()) == HEADER + record(20, "110") + record(30, "111")
    original_rows = original_table.splitlines(keepends=True)
    assert result.evidence_table.read_bytes() == original_rows[0] + original_rows[2] + original_rows[3]
    assert item.artifact.vcf_path.read_bytes() == original_vcf
    assert item.artifact.evidence_table.read_bytes() == original_table
    payload = result.validate()
    assert payload["request"]["source_order"] == ["sawfish", "sniffles2", "pav"]
    assert payload["policy"]["caller_count_is_confidence"] is False
    assert payload["policy"]["canonical_output_role_unchanged"] is True
    assert payload["policy"]["bnd_mate_completeness"] == "NOT_ASSESSED"
    commands = [command for command, _ in runner.calls]
    assert commands[-1] == ("tabix", "-l", str(item.paths["directory"] / ".working" / item.paths["vcf"].name))
    for _, kwargs in runner.calls[2:]:
        assert kwargs["timeout"] == 300
        assert "shell" not in kwargs
        assert kwargs["stderr_path"].parent == tmp_path / "logs"


def test_k_one_is_full_byte_identity_and_duplicate_ids_not_deduplicated(full, tmp_path):
    item = request(full, tmp_path, 1)
    result = SupportSelectionWrapper(runner=FakeRunner()).run(item)
    assert gzip.decompress(result.vcf_path.read_bytes()) == gzip.decompress(item.artifact.vcf_path.read_bytes())
    assert result.evidence_table.read_bytes() == item.artifact.evidence_table.read_bytes()
    assert result.selected_records == 3


def test_header_only_selection_is_valid(full, tmp_path):
    artifact, reference = full
    bgzf(artifact.vcf_path, HEADER + record(10, "100"))
    artifact.evidence_table.unlink()
    native = SVHarmonizationRequest("S1", reference, artifact.sources, tmp_path / "work", artifact.vcf_path, artifact.evidence_table, overwrite=True)
    write_evidence_table(native, artifact.vcf_path)
    result = SupportSelectionWrapper(runner=FakeRunner()).run(request((artifact, reference), tmp_path, 3))
    assert result.selected_records == 0
    assert gzip.decompress(result.vcf_path.read_bytes()) == HEADER
    assert len(result.evidence_table.read_bytes().splitlines()) == 1
    assert result.validate()["counts"]["source_records"] == 1


def test_dry_run_no_writes_or_completed_claim(full, tmp_path):
    item = request(full, tmp_path)
    runner = FakeRunner()
    result = SupportSelectionWrapper(runner=runner).run(item, dry_run=True)
    assert result.status == "planned" and result.to_dict()["executed"] is False
    assert not item.output_directory.exists()
    assert all(kwargs.get("dry_run") for _, kwargs in runner.calls)
    with pytest.raises(OutputValidationError, match="Planned"):
        result.validate()


@pytest.mark.parametrize("old,new", [
    (b"SUPP_VEC=110", b"SUPP_VEC=11"),
    (b"SUPP_VEC=110", b"SUPP_VEC=11x"),
    (b"SUPP_VEC=110", b"SUPP_VEC=000"),
    (b"SUPP_VEC=110;", b""),
    (b"SUPP=2", b"SUPP=3"),
    (b"SUPP=2", b"SUPP=2;SUPP=2"),
    (b"SUPP_VEC=110", b"SUPP_VEC=110;SUPP_VEC=110"),
])
def test_malformed_support_refuses_seal(full, tmp_path, old, new):
    artifact, _ = full
    content = gzip.decompress(artifact.vcf_path.read_bytes()).replace(old, new)
    bgzf(artifact.vcf_path, content)
    item = request(full, tmp_path)
    with pytest.raises(OutputValidationError):
        SupportSelectionWrapper(runner=FakeRunner()).run(item)
    assert not item.paths["manifest"].exists()
    assert not item.paths["vcf"].exists()


def test_absent_supp_scalar_is_allowed_when_vector_and_evidence_agree(full, tmp_path):
    artifact, _ = full
    content = gzip.decompress(artifact.vcf_path.read_bytes()).replace(b";SUPP=1", b"").replace(b";SUPP=2", b"").replace(b";SUPP=3", b"")
    bgzf(artifact.vcf_path, content)
    result = SupportSelectionWrapper(runner=FakeRunner()).run(request(full, tmp_path))
    assert result.selected_records == 2


@pytest.mark.parametrize("mode", ["short", "extra", "reordered", "count", "callers"])
def test_evidence_contract_mismatch_refuses_seal(full, tmp_path, mode):
    artifact, _ = full
    lines = artifact.evidence_table.read_bytes().splitlines(keepends=True)
    if mode == "short": lines = lines[:-1]
    elif mode == "extra": lines.append(lines[-1])
    elif mode == "reordered": lines[1], lines[2] = lines[2], lines[1]
    else:
        values = lines[2].rstrip(b"\n").split(b"\t")
        values[7 if mode == "count" else 10] = b"9" if mode == "count" else b"pav"
        lines[2] = b"\t".join(values) + b"\n"
    artifact.evidence_table.write_bytes(b"".join(lines))
    item = request(full, tmp_path)
    with pytest.raises(OutputValidationError):
        SupportSelectionWrapper(runner=FakeRunner()).run(item)
    assert not item.paths["manifest"].exists()


@pytest.mark.parametrize("step", [("bgzip", "-c"), ("tabix", "-p"), ("tabix", "-l")])
def test_native_failure_never_seals_or_publishes_selected_outputs(full, tmp_path, step):
    item = request(full, tmp_path)
    original = item.artifact.vcf_path.read_bytes()
    with pytest.raises(CommandExecutionError):
        SupportSelectionWrapper(runner=FakeRunner(fail=step)).run(item)
    assert item.artifact.vcf_path.read_bytes() == original
    assert item.output_directory.exists()
    assert not item.paths["manifest"].exists()
    assert not item.paths["vcf"].exists()


def test_refuses_existing_bundle_even_when_empty(full, tmp_path):
    item = request(full, tmp_path)
    item.output_directory.mkdir(parents=True)
    with pytest.raises(OutputValidationError, match="fresh bundle"):
        SupportSelectionWrapper(runner=FakeRunner()).run(item)
    with pytest.raises(OutputValidationError, match="fresh bundle"):
        request(full, tmp_path)


@pytest.mark.parametrize("role", ["vcf", "evidence", "manifest", "source"])
def test_validate_detects_tampering(full, tmp_path, role):
    item = request(full, tmp_path)
    result = SupportSelectionWrapper(runner=FakeRunner()).run(item)
    path = item.artifact.vcf_path if role == "source" else item.paths[role]
    path.write_bytes(path.read_bytes() + b"tamper")
    with pytest.raises((InputValidationError, OutputValidationError)):
        result.validate()


def test_vcf44_rejects_old_indexer_before_creating_bundle(full, tmp_path):
    item = request(full, tmp_path)
    with pytest.raises(ToolVersionError):
        SupportSelectionWrapper(runner=FakeRunner(version="1.21")).run(item)
    assert not item.output_directory.exists()


def test_bnd_selected_as_native_record_without_claiming_mate_closure(full, tmp_path):
    artifact, reference = full
    bgzf(artifact.vcf_path, HEADER + record(10, "110", kind="BND", alt="A]chr1:90]"))
    artifact.evidence_table.unlink()
    write_evidence_table(SVHarmonizationRequest("S1", reference, artifact.sources, tmp_path / "work", artifact.vcf_path, artifact.evidence_table, overwrite=True), artifact.vcf_path)
    result = SupportSelectionWrapper(runner=FakeRunner()).run(request(full, tmp_path))
    assert result.validate()["counts"]["selected_bnd_records"] == 1
    assert result.validate()["policy"]["bnd_mate_completeness"] == "NOT_ASSESSED"
