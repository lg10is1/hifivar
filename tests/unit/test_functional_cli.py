"""Privacy, local allele validation and reachable CLI tests for opt-in AlphaGenome."""
from pathlib import Path
import json
from types import SimpleNamespace
import pytest
import yaml
from hifivar.cli import build_parser, main
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.functional_cli import run_alphagenome_plan, _validate_ref_alleles

def fixture_plan(tmp_path, *, consent=False):
    ref = tmp_path / "ref.fa"
    ref.write_bytes(b">chr21\n" + b"A" * 10000 + b"C" + b"A" * 9999 + b"\n")
    Path(str(ref) + ".fai").write_text("chr21\t20000\t7\t20000\t20001\n", encoding="ascii")
    annotation = tmp_path / "annotation.tsv"
    annotation.write_text("variant\tannotation\nv1\tselected\n", encoding="utf-8")
    selection = tmp_path / "selected.tsv"
    selection.write_text(
        "sample\tsource_variant_id\tvariant_category\tcontig\tposition\treference_bases\talternate_bases\tsource_annotation\tselection_reason\n"
        "HG002\tv1\tsmall\tchr21\t10001\tC\tT\tannotation.tsv\texplicit_research_selection\n",
        encoding="utf-8")
    data = {"schema_version": 1, "selection_file": "selected.tsv",
            "reference": {"fasta": "ref.fa", "build": "GRCh38"},
            "alphagenome": {"model_version": "ALL_FOLDS", "requested_modalities": ["RNA_SEQ"],
                            "ontology_terms": ["UBERON:0000955"], "sequence_length": 16384},
            "data_sharing_consent": consent}
    plan = tmp_path / "plan.yaml"
    plan.write_text(yaml.safe_dump(data), encoding="utf-8")
    return plan, data, annotation

class FakeBackend:
    def __init__(self):
        self.calls = 0
    def score(self, request, *, allow_network=False):
        assert allow_network is True
        self.calls += 1
        assert request.selected_variants[0].source_variant_id == "v1"
        return SimpleNamespace(to_dict=lambda: {
            "score_rows": [{"source_variant_id": "v1", "modality": "RNA_SEQ",
                            "ontology_curie": "UBERON:0000955", "raw_score": 0.1}],
            "service_test": "SYNTHETIC_TRANSPORT_NOT_LIVE_API"})

def test_dry_run_no_backend_no_output(tmp_path):
    plan, _, _ = fixture_plan(tmp_path)
    backend = FakeBackend()
    out = tmp_path / "planned"
    result = run_alphagenome_plan(plan, out, dry_run=True, backend=backend)
    assert result["status"] == "PLANNED"
    assert result["network_execution"] == "NOT_RUN_DRY_RUN"
    assert result["reference"]["local_selected_ref_check"][0]["status"] == "MATCH"
    assert result["reference"]["whole_fasta_sha256"] == "NOT_COMPUTED"
    assert backend.calls == 0 and not out.exists()

@pytest.mark.parametrize("consent,allow", [(False, False), (True, False), (False, True)])
def test_two_explicit_network_gates(tmp_path, consent, allow):
    plan, _, _ = fixture_plan(tmp_path, consent=consent)
    backend = FakeBackend()
    out = tmp_path / "out"
    with pytest.raises(InputValidationError, match="allow-network"):
        run_alphagenome_plan(plan, out, allow_network=allow, backend=backend)
    assert backend.calls == 0 and not out.exists()

def test_ref_mismatch_blocks_network(tmp_path):
    plan, _, _ = fixture_plan(tmp_path, consent=True)
    selected = tmp_path / "selected.tsv"
    selected.write_text(selected.read_text().replace("\tC\tT\t", "\tA\tT\t"))
    backend = FakeBackend()
    with pytest.raises(InputValidationError, match="REF does not match"):
        run_alphagenome_plan(plan, tmp_path / "out", allow_network=True, backend=backend)
    assert backend.calls == 0

def test_literal_key_not_allowed_in_plan(tmp_path):
    plan, data, _ = fixture_plan(tmp_path)
    data["alphagenome"]["api_key"] = "synthetic-not-a-real-credential"
    plan.write_text(yaml.safe_dump(data))
    with pytest.raises(InputValidationError, match="documented fields"):
        run_alphagenome_plan(plan, tmp_path / "out", dry_run=True)

def test_success_keeps_input_bytes_and_rich_scores(tmp_path):
    plan, _, annotation = fixture_plan(tmp_path, consent=True)
    old = {p: p.read_bytes() for p in tmp_path.iterdir() if p.is_file()}
    backend = FakeBackend()
    out = tmp_path / "out"
    result = run_alphagenome_plan(plan, out, allow_network=True, backend=backend)
    saved = json.loads((out / "alphagenome.result.json").read_text())
    assert result == saved and saved["status"] == "COMPLETE"
    assert saved["result"]["score_rows"][0]["ontology_curie"] == "UBERON:0000955"
    assert saved["scientific_policy"]["functional_impact_is_call_confidence"] is False
    assert all(p.read_bytes() == value for p, value in old.items())
    with pytest.raises(OutputValidationError, match="must be new"):
        run_alphagenome_plan(plan, out, allow_network=True, backend=backend)
    assert backend.calls == 1

def test_provider_exception_never_serializes_text(tmp_path):
    plan, _, _ = fixture_plan(tmp_path, consent=True)
    sentinel = "synthetic credential sentinel"
    class Broken:
        def score(self, request, **kwargs):
            raise RuntimeError(sentinel)
    out = tmp_path / "failed"
    with pytest.raises(InputValidationError) as error:
        run_alphagenome_plan(plan, out, allow_network=True, backend=Broken())
    saved = (out / "alphagenome.result.json").read_text()
    assert sentinel not in saved and sentinel not in str(error.value)
    assert json.loads(saved)["status"] == "FAILED"
    assert error.value.__suppress_context__ is True

def test_input_mutation_invalidates_result(tmp_path):
    plan, _, annotation = fixture_plan(tmp_path, consent=True)
    class Mutator(FakeBackend):
        def score(self, request, **kwargs):
            annotation.write_text("changed")
            return super().score(request, **kwargs)
    out = tmp_path / "failed"
    with pytest.raises(InputValidationError, match="sanitized failure"):
        run_alphagenome_plan(plan, out, allow_network=True, backend=Mutator())
    result = json.loads((out / "alphagenome.result.json").read_text())
    assert result["status"] == "FAILED"
    assert result["source_immutability"] == "FAILED"
    assert "result" not in result

@pytest.mark.parametrize("newline", [b"\n", b"\r\n"])
def test_indexed_ref_multiline_and_crlf(tmp_path, newline):
    header = b">chr21" + newline
    ref = tmp_path / "ref.fa"
    ref.write_bytes(header + b"ACGT" + newline + b"ACGT" + newline)
    fai = tmp_path / "ref.fa.fai"
    fai.write_text(f"chr21\t8\t{len(header)}\t4\t{4+len(newline)}\n")
    v = SimpleNamespace(contig="chr21", position=3, reference_bases="GTAC",
                        sample_id="HG002", source_variant_id="v1")
    assert _validate_ref_alleles(ref, fai, (v,))[0]["status"] == "MATCH"

def test_new_cli_parser_and_dispatch(tmp_path, monkeypatch, capsys):
    from hifivar import functional_cli
    seen = {}
    def fake(plan, out, **kwargs):
        seen.update(kwargs)
        return {"status": "PLANNED", "selected_variant_count": 1, "network_execution": "NOT_RUN_DRY_RUN"}
    monkeypatch.setattr(functional_cli, "run_alphagenome_plan", fake)
    assert main(["--dry-run", "functional", "alphagenome", "--plan", str(tmp_path/"p.yaml"),
                 "--output-dir", str(tmp_path/"out")]) == 0
    assert seen == {"dry_run": True, "allow_network": False}
    assert "PLANNED" in capsys.readouterr().out
    args = build_parser().parse_args(["benchmark", "sv-depth", "--manifest", "m.yaml",
        "--output-dir", "out", "--unit-id", "depth5_r1", "--caller", "sniffles2", "--resume"])
    assert args.unit_ids == ["depth5_r1"] and args.callers == ["sniffles2"] and args.resume


def test_model_context_end_blocks_network(tmp_path):
    plan, _, _ = fixture_plan(tmp_path, consent=True)
    selected = tmp_path / "selected.tsv"
    selected.write_text(selected.read_text().replace("10001\tC\tT", "19000\tA\tT"))
    backend = FakeBackend()
    out = tmp_path / "out"
    with pytest.raises(InputValidationError, match="full model context"):
        run_alphagenome_plan(plan, out, allow_network=True, backend=backend)
    assert backend.calls == 0 and not out.exists()

@pytest.mark.parametrize("nested", [False, True])
def test_duplicate_authorization_and_settings_rejected(tmp_path, nested):
    plan, _, _ = fixture_plan(tmp_path, consent=False)
    text = plan.read_text()
    if nested:
        text = text.replace("model_version: ALL_FOLDS", "model_version: ALL_FOLDS\n  model_version: FOLD_0")
    else:
        text += "\ndata_sharing_consent: true\n"
    plan.write_text(text)
    backend = FakeBackend()
    with pytest.raises(InputValidationError, match="duplicate"):
        run_alphagenome_plan(plan, tmp_path/"out", allow_network=True, backend=backend)
    assert backend.calls == 0


def test_plan_changed_after_parse_is_not_bound_to_new_bytes(tmp_path, monkeypatch):
    from hifivar import functional_cli as mod
    plan, _, _ = fixture_plan(tmp_path, consent=True)
    original_path = mod._path
    changed = False
    def mutate(value, base, name):
        nonlocal changed
        if not changed:
            changed = True
            plan.write_text(plan.read_text().replace("ALL_FOLDS", "FOLD_0"))
        return original_path(value, base, name)
    monkeypatch.setattr(mod, "_path", mutate)
    backend = FakeBackend()
    with pytest.raises(InputValidationError, match="changed"):
        run_alphagenome_plan(plan, tmp_path/"out", allow_network=True, backend=backend)
    assert backend.calls == 0 and not (tmp_path/"out").exists()

def test_selection_changed_after_parse_blocks_dispatch(tmp_path, monkeypatch):
    from hifivar import functional_cli as mod
    plan, _, _ = fixture_plan(tmp_path, consent=True)
    original = mod.read_functional_selection
    def mutate(path):
        values = original(path)
        path.write_text(path.read_text().replace("explicit_research_selection", "changed_reason"))
        return values
    monkeypatch.setattr(mod, "read_functional_selection", mutate)
    backend = FakeBackend()
    with pytest.raises(InputValidationError, match="changed"):
        run_alphagenome_plan(plan, tmp_path/"out", allow_network=True, backend=backend)
    assert backend.calls == 0 and not (tmp_path/"out").exists()



def _write_ambiguous_selection(path, case):
    """Mutate only tiny synthetic TSV; the sentinel is never a real credential."""
    import csv
    import io
    fields = ["sample", "source_variant_id", "variant_category", "contig", "position",
              "reference_bases", "alternate_bases", "source_annotation", "selection_reason"]
    sentinel = "SYNTHETIC_ROW_VALUE_DO_NOT_ECHO_08ad3"
    values = ["HG002", "v1", "small", "chr21", "10001", "C", "T",
              "annotation.tsv", "explicit_research_selection"]
    if case == "duplicate_header":
        fields += ["alternate_bases"];values += ["G"]
    elif case == "extra_field":
        values += [sentinel]
    elif case == "short_row":
        values = values[:7]
    elif case == "bad_position":
        values[4] = sentinel
    elif case == "bad_category":
        values[2] = sentinel
    elif case == "missing_annotation":
        values[7] = sentinel
    elif case == "empty_annotation":
        values[7] = ""
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(fields)
    if case == "unclosed_quote":
        stream.write("\t".join(values[:-1]) + '\t"' + sentinel)
    else:
        writer.writerow(values)
    path.write_text(stream.getvalue(), encoding="utf-8", newline="")
    return sentinel


@pytest.mark.parametrize("case", [
    "duplicate_header", "extra_field", "short_row", "unclosed_quote",
    "bad_position", "bad_category", "missing_annotation", "empty_annotation",
])
def test_selection_tsv_ambiguity_and_bad_fields_are_controlled_without_row_echo(tmp_path, case):
    import traceback
    from hifivar.functional import read_functional_selection
    fixture_plan(tmp_path)
    selected = tmp_path / "selected.tsv"
    sentinel = _write_ambiguous_selection(selected, case)
    with pytest.raises(InputValidationError) as caught:
        read_functional_selection(selected)
    assert sentinel not in str(caught.value)
    assert sentinel not in "".join(traceback.format_exception(caught.value))
    if caught.value.__context__ is not None:
        assert caught.value.__suppress_context__ is True


@pytest.mark.parametrize("case", [
    "duplicate_header", "extra_field", "short_row", "unclosed_quote",
    "bad_position", "bad_category", "missing_annotation", "empty_annotation",
])
def test_selection_tsv_rejected_before_backend_credentials_sdk_or_output(tmp_path, monkeypatch, capsys, case):
    from hifivar import alphagenome_backend as adapter
    plan, _, _ = fixture_plan(tmp_path, consent=True)
    sentinel = _write_ambiguous_selection(tmp_path / "selected.tsv", case)
    calls = []
    def forbidden(*args, **kwargs):
        calls.append("FORBIDDEN_BACKEND_OR_SDK")
        raise AssertionError("Malformed selection must stop before backend/credential/SDK access")
    # The local adapter supplies the settings schema; optional SDK/auth factories
    # and backend construction must remain completely unused.
    monkeypatch.setattr(adapter.AlphaGenomeBackend, "__init__", forbidden)
    monkeypatch.setattr(adapter.AlphaGenomeBackend, "score", forbidden)
    monkeypatch.setattr(adapter, "_load_sdk", forbidden)
    monkeypatch.setattr(adapter, "_create_client", forbidden)
    output = tmp_path / "must_not_exist"
    code = main(["functional", "alphagenome", "--plan", str(plan),
                 "--output-dir", str(output), "--allow-network"])
    captured = capsys.readouterr()
    assert code != 0 and calls == [] and not output.exists()
    assert sentinel not in captured.out + captured.err


@pytest.mark.parametrize("newline", ["\n", "\r\n"])
def test_selection_tsv_keeps_bom_reordered_columns_and_quoted_fields(tmp_path, newline):
    import csv
    import io
    from hifivar.functional import read_functional_selection
    plan, _, _ = fixture_plan(tmp_path)
    fields = ["selection_reason", "source_annotation", "alternate_bases", "reference_bases",
              "position", "contig", "variant_category", "source_variant_id", "sample"]
    reason = 'SYNTHETIC explicit\tselection with "quoted text"'
    values = [reason, "annotation.tsv", "T", "C", "10001", "chr21", "small", "v1", "HG002"]
    text = io.StringIO(newline="")
    writer = csv.writer(text, delimiter="\t", lineterminator=newline)
    writer.writerow(fields);writer.writerow(values)
    selected = tmp_path / "selected.tsv"
    selected.write_bytes(b"\xef\xbb\xbf" + text.getvalue().encode("utf-8"))
    parsed = read_functional_selection(selected)
    assert len(parsed) == 1 and parsed[0].selection_reason == reason
    assert parsed[0].reference_bases == "C" and parsed[0].alternate_bases == "T"
    backend = FakeBackend();output = tmp_path / "no_output"
    report = run_alphagenome_plan(plan, output, dry_run=True, backend=backend)
    assert report["status"] == "PLANNED" and backend.calls == 0 and not output.exists()



def test_selection_annotation_expanduser_failure_is_controlled_and_private(tmp_path, monkeypatch):
    import traceback
    from hifivar.functional import read_functional_selection
    plan, _, _ = fixture_plan(tmp_path, consent=True)
    sentinel = "SYNTHETIC_HOME_LOOKUP_VALUE_DO_NOT_ECHO_42bb"
    original = Path.expanduser
    def unavailable_home(path):
        if path.name == "annotation.tsv":
            raise RuntimeError(sentinel)
        return original(path)
    monkeypatch.setattr(Path, "expanduser", unavailable_home)
    with pytest.raises(InputValidationError) as caught:
        read_functional_selection(tmp_path / "selected.tsv")
    assert sentinel not in "".join(traceback.format_exception(caught.value))
    assert caught.value.__suppress_context__ is True
    backend = FakeBackend();output = tmp_path / "no_home_output"
    with pytest.raises(InputValidationError):
        run_alphagenome_plan(plan, output, allow_network=True, backend=backend)
    assert backend.calls == 0 and not output.exists()
