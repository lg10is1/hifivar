"""Offline review CLI contracts; synthetic files are not biological validation."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest
import yaml

from hifivar.cli import main
from hifivar.config import load_config
from hifivar.exceptions import ConfigurationError, InputValidationError, OutputValidationError
from hifivar.review_cli import import_review_decisions, prepare_review

DEFAULT = Path(__file__).parents[2] / "src/hifivar/resources/configs/default.yaml"
FIELDS = ["review_id", "sample", "variant_id", "status", "reviewer", "reviewed_at", "comment"]


@pytest.fixture
def review_case(tmp_path):
    # Native IGV reference paths require whitespace-free ancestors; output paths
    # with spaces and Unicode remain separate positive coverage below.
    root = tmp_path / "复核_inputs"
    root.mkdir()
    reference = root / "reference.fa"
    reference.write_text(">chr1\n" + "A" * 1000 + "\n>chr2\n" + "C" * 1000 + "\n", encoding="utf-8")
    Path(f"{reference}.fai").write_text("chr1\t1000\t6\t1000\t1001\nchr2\t1000\t1013\t1000\t1001\n", encoding="utf-8")
    bam = root / "S1.bam"
    bam.write_bytes(b"synthetic alignment placeholder")
    Path(f"{bam}.bai").write_bytes(b"synthetic index placeholder")
    vcf = root / "source.vcf"
    vcf.write_text("##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n##contig=<ID=chr2,length=1000>\n"
                   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                   "chr1\t100\tv1\tA\tT\t.\tPASS\t.\tGT\t0/1\n"
                   "chr1\t200\tbnd1\tA\tA[chr2:500[\t.\tPASS\tSVTYPE=BND\tGT\t0/1\n", encoding="utf-8")
    sheet = root / "samples.tsv"
    sheet.write_text("sample_id\tinput\tinput_type\nS1\tS1.bam\tbam\n", encoding="utf-8")
    selection = root / "selection.tsv"
    selection.write_text(
        "review_id\tsample\tvariant_id\tvariant_type\tcontig\tstart\tend\tsource_vcf\tsource_caller\tevidence_class\tmate_contig\tmate_position\n"
        "R1\tS1\tv1\tSNV\tchr1\t100\t100\tsource.vcf\tdeepvariant\tsmall_variant\t\t\n"
        "R2\tS1\tbnd1\tBND\tchr1\t200\t200\tsource.vcf\texplicit\texplicit\tchr2\t500\n",
        encoding="utf-8")
    config_path = root / "analysis.yaml"
    config_path.write_text(yaml.safe_dump({"reference": {"fasta": "reference.fa", "build": "GRCh38"},
                                          "samples": {"sheet": "samples.tsv"},
                                          "review": {"selection_file": "selection.tsv", "flank_bp": 10,
                                                     "igv_executable": "tools/igv.sh"}}), encoding="utf-8")
    config = load_config(default_config=DEFAULT, user_config=config_path)
    return config, root, config_path, tmp_path / "review prepared"


def hashes(root):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in root.iterdir() if p.is_file()}


def decisions(path, rows):
    with path.open("wt", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            value = {key: "" for key in FIELDS}
            value.update({"review_id": "R1", "sample": "S1", "variant_id": "v1", "status": "SUPPORT",
                          "reviewer": "synthetic-reviewer", "reviewed_at": "2026-10-14T10:00:00+08:00",
                          "comment": "Synthetic fixture observation, not an actual review."})
            value.update(row)
            writer.writerow(value)


def test_prepare_writes_pending_bundle_without_external_execution(review_case, monkeypatch):
    config, root, _, output = review_case
    before = hashes(root)
    monkeypatch.setattr("hifivar.command.CommandRunner.run", lambda *a, **k: pytest.fail("no command execution"))
    outcome = prepare_review(config, output)
    assert outcome["status"] == "REVIEW_PREPARED"
    assert outcome["targets"] == 2 and outcome["expected_screenshots"] == 3
    assert not outcome["screenshots_generated"] and outcome["manual_decisions_completed"] == 0
    assert not list((output / "screenshots").iterdir())
    data = json.loads((output / "review_manifest.json").read_text(encoding="utf-8"))
    assert [row["status"] for row in data["results"]] == ["NOT_REVIEWED", "NOT_REVIEWED"]
    receipt = json.loads((output / "review_bundle.json").read_text(encoding="utf-8"))
    assert receipt["raw_input_hash_policy"] == "STAT_ONLY_NOT_FULL_CONTENT_HASH"
    assert receipt["reference_build"] == "GRCh38"
    assert receipt["result"]["igv_command"][0] == str((root / "tools/igv.sh").absolute())
    batch = (output / "review.igv.batch").read_text(encoding="utf-8")
    assert "goto chr1:90-110" in batch and "goto chr2:490-510" in batch
    assert batch.count("snapshot R2.locus") == 2
    assert hashes(root) == before


def test_prepare_dry_run_is_zero_write_even_with_log_file(review_case, capsys):
    _, root, config_path, output = review_case
    log = root / "unowned.log"
    before = hashes(root)
    assert main(["--config", str(config_path), "--log-file", str(log), "--dry-run",
                 "review", "prepare", "--output-dir", str(output)]) == 0
    assert json.loads(capsys.readouterr().out)["status"] == "REVIEW_PLAN"
    assert not output.exists() and not log.exists() and hashes(root) == before


def test_existing_output_refused_and_preserved(review_case):
    config, _, _, output = review_case
    prepare_review(config, output)
    before = hashes(output)
    with pytest.raises(OutputValidationError, match="already exists"):
        prepare_review(config, output)
    assert hashes(output) == before


@pytest.mark.parametrize("status", ["SUPPORT", "NOT_SUPPORT", "UNCERTAIN"])
def test_import_named_decisions_preserves_missing_and_raw_records(review_case, tmp_path, status):
    config, root, _, output = review_case
    prepare_review(config, output)
    before = hashes(output); inputs_before = hashes(root)
    table = tmp_path / "reviewer.tsv"
    decisions(table, [{"status": status}])
    dest = tmp_path / "adjudicated"
    result = import_review_decisions(output, table, dest)
    assert result["reviewed_targets"] == 1 and result["missing_decision_ids"] == ["R2"]
    data = json.loads((dest / "review_manifest.adjudicated.json").read_text(encoding="utf-8"))
    assert [r["status"] for r in data["results"]] == [status, "NOT_REVIEWED"]
    for row in data["results"]:
        assert row["status_is_truth"] is False and row["reviewer_identity_verified"] is False
    assert data["results"][0]["reviewer"] == "synthetic-reviewer"
    assert data["results"][1]["reviewer"] == ""
    original = json.loads((output / "review_manifest.json").read_text(encoding="utf-8"))
    for old, new in zip(original["results"], data["results"], strict=True):
        for key in ("sample", "variant_id", "locus", "source_vcf", "source_bam", "screenshots"):
            assert old[key] == new[key]
    assert hashes(output) == before and hashes(root) == inputs_before


@pytest.mark.parametrize("change", [
    {"review_id": "unknown"}, {"sample": "S2"}, {"variant_id": "wrong"},
    {"status": "PASS"}, {"status": "REJECT"}, {"reviewer": ""}, {"comment": ""},
    {"reviewed_at": ""}, {"reviewed_at": "2026-10-14"},
    {"reviewed_at": "2026-10-14T10:00:00"}, {"reviewed_at": "not-a-date"},
    {"comment": "two\nlines"}, {"status": "NOT_REVIEWED"},
])
def test_invalid_decisions_fail_before_output_creation(review_case, tmp_path, change):
    config, _, _, output = review_case
    prepare_review(config, output)
    table = tmp_path / "invalid.tsv"
    decisions(table, [change])
    dest = tmp_path / "invalid-output"
    with pytest.raises(InputValidationError):
        import_review_decisions(output, table, dest)
    assert not dest.exists()


def test_duplicate_decision_fails(review_case, tmp_path):
    config, _, _, output = review_case
    prepare_review(config, output)
    table = tmp_path / "duplicate.tsv"
    decisions(table, [{}, {}])
    with pytest.raises(InputValidationError, match="duplicate"):
        import_review_decisions(output, table, tmp_path / "out")


@pytest.mark.parametrize("name", ["review_manifest.json", "review.igv.batch", "effective_config.yaml"])
def test_changed_sealed_bundle_refused(review_case, tmp_path, name):
    config, _, _, output = review_case
    prepare_review(config, output)
    with (output / name).open("at", encoding="utf-8") as stream:
        stream.write("changed")
    table = tmp_path / "reviewer.tsv"
    decisions(table, [{}])
    with pytest.raises(InputValidationError, match="changed"):
        import_review_decisions(output, table, tmp_path / "imported")


def test_template_is_editable_and_can_be_partially_imported(review_case, tmp_path):
    config, _, _, output = review_case
    prepare_review(config, output)
    template = output / "manual_decisions.tsv"
    decisions(template, [{}])
    assert import_review_decisions(output, template, tmp_path / "from-template")["reviewed_targets"] == 1


def test_import_dry_run_and_existing_destination_preservation(review_case, tmp_path):
    config, _, _, output = review_case
    prepare_review(config, output)
    table = tmp_path / "decisions.tsv"
    decisions(table, [{}])
    dest = tmp_path / "imported"
    assert import_review_decisions(output, table, dest, dry_run=True)["status"] == "DECISIONS_PLAN"
    assert not dest.exists()
    import_review_decisions(output, table, dest)
    before = hashes(dest)
    with pytest.raises(OutputValidationError, match="already exists"):
        import_review_decisions(output, table, dest)
    assert hashes(dest) == before


@pytest.mark.parametrize("kind", ["duplicate-header", "extra-field", "short-row"])
def test_malformed_selection_refused(review_case, kind):
    config, root, _, output = review_case
    table = root / "selection.tsv"
    lines = table.read_text(encoding="utf-8").splitlines()
    if kind == "duplicate-header": lines[0] += "\treview_id"
    if kind == "extra-field": lines[1] += "\textra"
    if kind == "short-row": lines[1] = "R1\tS1"
    table.write_text("\n".join(lines) + "\n", encoding="utf-8")
    with pytest.raises(InputValidationError): prepare_review(config, output)
    assert not output.exists()


def test_selection_source_mate_identity_and_missing_index(review_case):
    config, root, _, output = review_case
    table = root / "selection.tsv"
    original = table.read_text(encoding="utf-8")
    table.write_text(original.replace("chr2\t500", "chr2\t1001"), encoding="utf-8")
    with pytest.raises(InputValidationError, match="mate differs from the exact selected ALT"):
        prepare_review(config, output)
    table.write_text(original, encoding="utf-8")
    (root / "S1.bam.bai").unlink()
    with pytest.raises(InputValidationError): prepare_review(config, output)
    assert not output.exists()


def test_empty_selection_and_empty_decisions_are_not_review_completion(review_case, tmp_path):
    config, root, _, output = review_case
    table = root / "selection.tsv"
    table.write_text(table.read_text(encoding="utf-8").splitlines()[0] + "\n", encoding="utf-8")
    assert prepare_review(config, output)["targets"] == 0
    decision_table = tmp_path / "empty.tsv"
    decisions(decision_table, [])
    result = import_review_decisions(output, decision_table, tmp_path / "empty-import")
    assert result["reviewed_targets"] == 0 and result["total_targets"] == 0


def test_write_failure_has_no_completion_receipt_and_keeps_sources(review_case, monkeypatch):
    config, root, _, output = review_case
    before = hashes(root)
    def fail(*args, **kwargs): raise OutputValidationError("synthetic write failure")
    monkeypatch.setattr("hifivar.review_cli.write_yaml_atomic", fail)
    with pytest.raises(OutputValidationError, match="write failure"):
        prepare_review(config, output)
    assert output.exists() and not (output / "review_bundle.json").exists()
    assert hashes(root) == before


def test_cli_prepare_and_import(review_case, tmp_path, capsys):
    _, _, config_path, output = review_case
    assert main(["--config", str(config_path), "review", "prepare", "--output-dir", str(output)]) == 0
    assert json.loads(capsys.readouterr().out)["status"] == "REVIEW_PREPARED"
    table = tmp_path / "manual.tsv"
    decisions(table, [{}])
    assert main(["review", "import-decisions", "--bundle-dir", str(output), "--decisions", str(table),
                 "--output-dir", str(tmp_path / "manual-import")]) == 0
    assert json.loads(capsys.readouterr().out)["reviewed_targets"] == 1


def test_missing_config_or_bundle_is_concise_cli_error(tmp_path, capsys):
    assert main(["review", "prepare", "--output-dir", str(tmp_path / "out")]) == 1
    assert "selection_file" in capsys.readouterr().err
    assert main(["review", "import-decisions", "--bundle-dir", str(tmp_path / "missing"),
                 "--decisions", str(tmp_path / "missing.tsv"), "--output-dir", str(tmp_path / "out")]) == 1
    assert "ERROR:" in capsys.readouterr().err


@pytest.mark.parametrize("change", [{"schema_version": True}, {"schema_version": 2},
                                   {"sealed_artifacts": {"../outside": "bad"}}])
def test_invalid_receipt_schema_or_path_list_refused(review_case, tmp_path, change):
    config, _, _, output = review_case
    prepare_review(config, output)
    receipt_path = output / "review_bundle.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8")); receipt.update(change)
    receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
    table = tmp_path / "manual.tsv"; decisions(table, [{}])
    with pytest.raises(InputValidationError):
        import_review_decisions(output, table, tmp_path / "invalid-receipt")


def test_import_cannot_write_inside_original_bundle(review_case, tmp_path):
    config, _, _, output = review_case
    prepare_review(config, output)
    table = tmp_path / "manual.tsv"; decisions(table, [{}])
    with pytest.raises(OutputValidationError, match="outside"):
        import_review_decisions(output, table, output / "nested-import")
    assert not (output / "nested-import").exists()


def test_forged_truth_semantics_rejected_even_with_updated_local_checksum(review_case, tmp_path):
    config, _, _, output = review_case
    prepare_review(config, output)
    manifest_path = output / "review_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["results"][0]["status_is_truth"] = True
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    receipt_path = output / "review_bundle.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    receipt["sealed_artifacts"][manifest_path.name] = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
    table = tmp_path / "manual.tsv"; decisions(table, [{}])
    with pytest.raises(InputValidationError, match="non-truth"):
        import_review_decisions(output, table, tmp_path / "no-truth")
