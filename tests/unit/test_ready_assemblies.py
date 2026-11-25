"""Read-only, bounded and explicit external assembly reuse contracts."""
from dataclasses import replace
import hashlib
import os
from pathlib import Path

import pytest

from hifivar.assembly import AssemblyRole
from hifivar.exceptions import InputValidationError
from hifivar.ready_assemblies import (
    _scan_fasta, load_ready_assembly_manifest, protect_ready_assembly_inputs,
    validate_ready_assembly,
)

H1, H2 = AssemblyRole.HAPLOTYPE1, AssemblyRole.HAPLOTYPE2


def fasta(tmp_path, name="input.fa", payload=b">c1\nACGTNacgtryswkmbdhvn\n"):
    path = tmp_path / name
    path.write_bytes(payload)
    return path


def validate(path, role=H1, sample="S1"):
    return validate_ready_assembly(sample, role, path, hashlib.sha256(path.read_bytes()).hexdigest())


def manifest(tmp_path, rows, header="sample_id\trole\tfasta\tsha256"):
    path = tmp_path / "ready.tsv"
    path.write_text(header + "\n" + "\n".join(rows) + "\n", encoding="utf-8")
    return path


def row(path, sample="S1", role="haplotype1"):
    return f"{sample}\t{role}\t{path.name}\t{hashlib.sha256(path.read_bytes()).hexdigest()}"


def test_external_validation_is_read_only_and_does_not_invent_hifiasm(tmp_path):
    path = fasta(tmp_path)
    before = path.read_bytes(), path.stat().st_mtime_ns
    item = validate(path)
    assert (path.read_bytes(), path.stat().st_mtime_ns) == before
    assert set(tmp_path.iterdir()) == {path}
    data = item.to_dict()
    assert data["execution_kind"] == "VALIDATION_ONLY"
    assert data["role_source"] == "USER_DECLARED_NOT_PHASE_TRUTH"
    assert data["hifiasm_version"] is data["source_gfa"] is data["reference_independent"] is None
    assert data["source_version_source"] == "UNKNOWN"
    item.verify()


@pytest.mark.parametrize("payload", [
    b"", b"ACGT\n", b">\nACGT\n", b">a\n", b">a\n>b\nACG\n",
    b">a\nACG\n>b\n", b">a\nACG\n>a duplicate\nACG\n", b">a\nAC GU\n",
    b">a\nACG*\n", b">a\nACGU\n", b">a\nACG\n\n", b">a\nACG\r",
    b">a\nAC\rGT\n", b">a\x00\nACG\n", b">\xff\nACG\n",
    b"\xef\xbb\xbf>a\nACG\n", b">a\nACG\n>  \nACG\n", b"> a\nACG\n",
])
def test_full_fasta_rejects_invalid_records_including_late_records(tmp_path, payload):
    with pytest.raises(InputValidationError):
        validate(fasta(tmp_path, payload=payload))


@pytest.mark.parametrize("chunk", [1, 2, 3, 4, 7, 11, 1024])
@pytest.mark.parametrize("payload", [b">a description\r\nACGTRY\r\n>b\r\nACGT\r\n",
                                      b">a\nACGTRY\n>b\nACGT"])
def test_chunk_boundaries_crlf_and_missing_terminal_newline(tmp_path, chunk, payload):
    path = fasta(tmp_path, payload=payload)
    digest, size, records, bases = _scan_fasta(path, chunk_size=chunk)
    assert (digest, size, records, bases) == (hashlib.sha256(payload).hexdigest(), len(payload), 2, 10)


def test_unwrapped_sequence_crosses_default_chunk_boundary(tmp_path):
    sequence = b"ACGT" * (1024 * 1024)
    path = fasta(tmp_path, payload=b">long\n" + sequence + b"\n")
    assert validate(path).base_count == len(sequence)


@pytest.mark.parametrize("chunk", [0, -1, True, 9 * 1024 * 1024])
def test_unbounded_or_invalid_stream_chunk_is_rejected(tmp_path, chunk):
    with pytest.raises(InputValidationError):
        _scan_fasta(fasta(tmp_path), chunk_size=chunk)


@pytest.mark.parametrize("sample", ["", "../S1", "S1\tX", ".hidden", "S 1", None])
def test_sample_id_is_not_normalized(tmp_path, sample):
    with pytest.raises(InputValidationError):
        validate(fasta(tmp_path), sample=sample)


@pytest.mark.parametrize("role", [AssemblyRole.PRIMARY, "haplotype1", None])
def test_explicit_haplotype_roles_only(tmp_path, role):
    with pytest.raises(InputValidationError):
        validate(fasta(tmp_path), role=role)


@pytest.mark.parametrize("digest", ["", "abc", "x" * 64, None, "0" * 64])
def test_full_expected_sha_is_mandatory(tmp_path, digest):
    with pytest.raises(InputValidationError):
        validate_ready_assembly("S1", H1, fasta(tmp_path), digest)


@pytest.mark.parametrize("name", ["reads.fastq", "assembly.fa.gz", "input.bam"])
def test_no_implicit_format_conversion(tmp_path, name):
    with pytest.raises(InputValidationError):
        validate(fasta(tmp_path, name))


def test_embedded_nul_manifest_path_is_a_clear_input_error(tmp_path):
    path = fasta(tmp_path)
    line = row(path).replace("input.fa", "invalid\x00.fa")
    with pytest.raises(InputValidationError, match="Cannot read"):
        load_ready_assembly_manifest(manifest(tmp_path, [line]), expected_samples=("S1",))


def test_declared_source_is_not_a_probed_version(tmp_path):
    path = fasta(tmp_path)
    item = validate_ready_assembly("S1", H1, path, hashlib.sha256(path.read_bytes()).hexdigest().upper(),
                                   source_tool="external_assembler", source_version="declared-1")
    assert item.to_dict()["source_version_source"] == "USER_DECLARED"
    assert item.to_dict()["hifiasm_version"] is None


def test_changed_fasta_rejected_even_if_same_size(tmp_path):
    path = fasta(tmp_path, payload=b">a\nACGT\n")
    item = validate(path)
    path.write_bytes(b">a\nTGCA\n")
    with pytest.raises(InputValidationError, match="differs"):
        item.verify()


def test_manifest_preserves_declared_order_and_has_deterministic_caller_roles(tmp_path, monkeypatch):
    a, b, c = fasta(tmp_path, "a.fa"), fasta(tmp_path, "b.fa"), fasta(tmp_path, "c.fa")
    path = manifest(tmp_path, [row(a, "S2", "haplotype2"), row(b, "S1"), row(c, "S2")])
    monkeypatch.chdir(tmp_path.parent)
    result = load_ready_assembly_manifest(path, expected_samples=("S1", "S2"))
    assert result.declared_sample_order == ("S2", "S1")
    assert [item.path for item in result.for_sample("S2")] == [c.resolve(), a.resolve()]
    assert [item.role for item in result.assemblies] == [H2, H1, H1]
    assert result.to_dict()["new_hifiasm_execution"] is False
    assert result.to_dict()["execution_kind"] == "VALIDATION_ONLY"
    for item in result.assemblies:
        item.verify()
    with pytest.raises(InputValidationError):
        result.for_sample("absent")


def test_utf8_bom_manifest_supported(tmp_path):
    path = fasta(tmp_path)
    sheet = manifest(tmp_path, [row(path)])
    sheet.write_bytes(b"\xef\xbb\xbf" + sheet.read_bytes())
    assert load_ready_assembly_manifest(sheet, expected_samples=("S1",)).assemblies[0].path == path.resolve()


@pytest.mark.parametrize("kind", ["missing", "extra", "duplicate_role", "duplicate_header", "unknown_column",
                                  "missing_column", "bad_role", "empty_path", "malformed", "empty_sample"])
def test_manifest_negative_contracts(tmp_path, kind):
    a, b = fasta(tmp_path, "a.fa"), fasta(tmp_path, "b.fa")
    rows, header = [row(a)], "sample_id\trole\tfasta\tsha256"
    expected = ("S1", "S2") if kind == "missing" else ("S1",)
    if kind == "extra": rows.append(row(b, "S2"))
    if kind == "duplicate_role": rows.append(row(b))
    if kind == "duplicate_header": header += "\tsha256"
    if kind == "unknown_column": header += "\tsecret"; rows[0] += "\tvalue"
    if kind == "missing_column": header = "sample_id\trole\tfasta"
    if kind == "bad_role": rows[0] = row(a, role="primary")
    if kind == "empty_path": rows[0] = rows[0].replace("\ta.fa\t", "\t\t")
    if kind == "malformed": rows[0] += "\textra"
    if kind == "empty_sample": rows[0] = row(a, sample="")
    with pytest.raises(InputValidationError):
        load_ready_assembly_manifest(manifest(tmp_path, rows, header), expected_samples=expected)


@pytest.mark.parametrize("expected", ["S1", (), ("S1", "S1"), ("S 1",), None])
def test_expected_sample_set_is_strict(tmp_path, expected):
    path = fasta(tmp_path)
    with pytest.raises(InputValidationError):
        load_ready_assembly_manifest(manifest(tmp_path, [row(path)]), expected_samples=expected)


@pytest.mark.parametrize("alias", ["same_path", "hardlink"])
def test_roles_and_samples_must_be_physical_distinct(tmp_path, alias):
    a = fasta(tmp_path, "a.fa")
    b = a
    if alias == "hardlink":
        b = tmp_path / "b.fa"
        os.link(a, b)
    with pytest.raises(InputValidationError, match="alias"):
        load_ready_assembly_manifest(manifest(tmp_path, [row(a), row(b, role="haplotype2")]),
                                     expected_samples=("S1",))


def test_manifest_mutation_rejected_without_changing_original_fasta(tmp_path):
    path = fasta(tmp_path)
    sheet = manifest(tmp_path, [row(path)])
    item = load_ready_assembly_manifest(sheet, expected_samples=("S1",)).assemblies[0]
    sheet.write_bytes(sheet.read_bytes() + b"\n")
    with pytest.raises(InputValidationError, match="manifest changed"):
        item.verify()


@pytest.mark.parametrize("change", [{"role": H2}, {"sample_id": "S2"},
                                    {"source_tool": "invented"}, {"source_version": "invented"}])
def test_manifest_role_and_source_cannot_be_relabelled_after_validation(tmp_path, change):
    path = fasta(tmp_path)
    item = load_ready_assembly_manifest(manifest(tmp_path, [row(path)]), expected_samples=("S1",)).assemblies[0]
    with pytest.raises(InputValidationError, match="manifest row"):
        replace(item, **change).verify()


@pytest.mark.parametrize("alias", ["same_path", "hardlink", "manifest", "workdir"])
def test_caller_output_log_and_workdir_must_not_alias_inputs(tmp_path, alias):
    path = fasta(tmp_path)
    sheet = manifest(tmp_path, [row(path)])
    item = load_ready_assembly_manifest(sheet, expected_samples=("S1",)).assemblies[0]
    target = path
    if alias == "hardlink":
        target = tmp_path / "output.log"; os.link(path, target)
    if alias == "manifest": target = sheet
    with pytest.raises(InputValidationError):
        protect_ready_assembly_inputs((item,), outputs=(target,) if alias != "workdir" else (),
                                      work_directory=tmp_path if alias == "workdir" else tmp_path / "new_work")


def test_forged_counts_cannot_bypass_full_readback(tmp_path):
    item = validate(fasta(tmp_path))
    with pytest.raises(InputValidationError):
        replace(item, base_count=999).verify()


def test_changed_file_during_full_validation_is_rejected(tmp_path, monkeypatch):
    path = fasta(tmp_path)
    original = os.fstat
    calls = 0
    def changing(fd):
        nonlocal calls
        calls += 1
        if calls == 2:
            path.write_bytes(b">changed\nAAAA\n")
        return original(fd)
    monkeypatch.setattr(os, "fstat", changing)
    with pytest.raises(InputValidationError, match="changed"):
        validate(path)
