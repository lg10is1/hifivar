"""Source identity is provenance, not normalization, filtering or confidence."""
from dataclasses import replace
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
import struct
import zlib

import pytest

from hifivar.assembly_sv import SVEvidenceSource
from hifivar.command import CommandResult
from hifivar.config import load_config
from hifivar.exceptions import ConfigurationError, InputValidationError, OutputValidationError
from hifivar.harmonization import EvidenceRunStatus, SVEvidenceSourceArtifact, SVHarmonizationRequest, write_evidence_table
from hifivar.jasmine import JasmineWrapper
from hifivar.phase9 import Phase9Settings, phase9_identity_paths
from hifivar.reference import ReferenceGenome
from hifivar.sv_input_identity import POLICY, POLICY_SHA256, SourceIdentityPolicy, SourceIdentityWrapper, verify_source_view

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "src/hifivar/resources/configs/default.yaml"


def bgzf(path, data):
    obj = zlib.compressobj(level=6, wbits=-15)
    packed = obj.compress(data) + obj.flush()
    header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", 18 + len(packed) + 8 - 1)
    path.write_bytes(header + packed + struct.pack("<II", zlib.crc32(data), len(data)))


def header(sample="S1", ending="\n"):
    return ending.join([
        "##fileformat=VCFv4.3", "##contig=<ID=chr1,length=1000>",
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
        '##INFO=<ID=IDLIST,Number=.,Type=String,Description="Members">',
        '##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Sources">',
        '##INFO=<ID=SUPP,Number=1,Type=Integer,Description="Count">',
        '##INFO=<ID=MATEID,Number=.,Type=String,Description="Mate">',
        '##INFO=<ID=PARID,Number=1,Type=String,Description="Parent">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}", "",
    ])


def record(identifier, position=10, *, extra="", gt="0|1", ending="\n"):
    return f"chr1\t{position}\t{identifier}\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END={position+60}{extra}\tGT\t{gt}{ending}"


def policy():
    return SourceIdentityPolicy(POLICY["policy_id"], 1, POLICY_SHA256, POLICY["mode"])


@pytest.fixture
def inputs(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + "A" * 1000 + "\n")
    Path(str(fasta) + ".fai").write_text("chr1\t1000\t6\t1000\t1001\n")
    ref = ReferenceGenome.from_fasta(fasta, build="GRCh38")
    def source(caller="pav", records=None, sample="S1", ending="\n"):
        path = tmp_path / f"S1.{caller}.sv.vcf.gz"
        bgzf(path, (header(sample, ending) + (records if records is not None else record("DUPID") + record("DUPID", 20, gt="1/1"))).encode())
        bgzf(Path(str(path) + ".tbi"), b"TBI\x01")
        return SVEvidenceSourceArtifact("S1", SVEvidenceSource.ASSEMBLY if caller in {"pav", "svim_asm"} else SVEvidenceSource.READ,
            caller, path, Path(str(path) + ".tbi"), EvidenceRunStatus.COMPLETED,
            ("haplotype1", "haplotype2") if caller in {"pav", "svim_asm"} else ())
    return ref, source


class NativeFixtureRunner:
    """Models command contracts; does not claim actual native indexing."""
    def __init__(self, tmp_path, *, merged_transform=None):
        self.tools = {}
        for name in ("jasmine", "bgzip", "tabix"):
            p = tmp_path / name
            p.write_bytes(b"#!/bin/sh\n# executable fixture, never run as a shell\n")
            self.tools[name] = p
        self.calls = []
        self.merged_transform = merged_transform

    def require_executable(self, executable, **kwargs):
        return self.tools[Path(executable).name]

    def run(self, args, **kwargs):
        args = tuple(map(str, args))
        self.calls.append((args, kwargs))
        name, output = Path(args[0]).name, ""
        if kwargs.get("dry_run"):
            return CommandResult(args, 0, "", "", 0, None, False)
        if "--version" in args:
            output = "Jasmine 1.1.5" if name == "jasmine" else f"{name} 1.24"
        elif name == "bgzip":
            bgzf(Path(kwargs["stdout_path"]), Path(args[-1]).read_bytes())
        elif name == "tabix" and args[1:3] == ("-p", "vcf"):
            bgzf(Path(args[-1] + ".tbi"), b"TBI\x01")
        elif name == "tabix" and args[1] == "-l":
            output = "chr1\n"
        elif name == "tabix":
            with gzip.open(args[1], "rt") as reader:
                output = "".join(r for r in reader if not r.startswith("#"))
        elif name == "jasmine":
            listing = Path(next(a.split("=", 1)[1] for a in args if a.startswith("file_list=")))
            paths = [Path(x) for x in listing.read_text().splitlines()]
            identifiers = []
            for p in paths:
                identifiers.extend(r.split("\t")[2] for r in p.read_text().splitlines() if not r.startswith("#"))
            merged = record("J1", extra=";IDLIST=" + ",".join(identifiers) + ";SUPP_VEC=" + "1" * len(paths) + f";SUPP={len(paths)}")
            if self.merged_transform:
                merged = self.merged_transform(merged, identifiers)
            Path(next(a.split("=", 1)[1] for a in args if a.startswith("out_file="))).write_text(header() + merged)
        return CommandResult(args, 0, output, "", .01, None, True)


def prepare(tmp_path, inputs, sources=None):
    ref, create = inputs
    sources = sources or (create(),)
    runner = NativeFixtureRunner(tmp_path)
    result = SourceIdentityWrapper(runner=runner).run(sources, sample_id="S1", reference=ref,
        output_directory=tmp_path / "identity", policy=policy())
    return result, sources, runner


@pytest.mark.parametrize("ending", ["\n", "\r\n"])
def test_duplicate_different_gt_preserves_all_non_id_bytes_and_eol(tmp_path, inputs, ending):
    ref, create = inputs
    source = create(records=record("DUPID", ending=ending) + record("DUPID", 20, gt="1/1", ending=ending), ending=ending)
    native = source.vcf_path.read_bytes()
    result, _, runner = prepare(tmp_path, inputs, (source,))
    assert source.vcf_path.read_bytes() == native
    rows = [json.loads(x) for x in result.ledgers[0].read_text().splitlines()]
    assert [r["original_record_id"] for r in rows] == ["DUPID", "DUPID"]
    assert len({r["analysis_record_id"] for r in rows}) == 2
    assert all(r["analysis_record_id"].startswith("HFVID1_pav_") for r in rows)
    assert verify_source_view(source.vcf_path, result.views[0], result.ledgers[0])["records"] == 2
    assert result.validate()["source_summaries"] == [{"records": 2, "changed_ID_records": 2}]
    assert all(kw.get("inherit_env") is False and "timeout" in kw for _, kw in runner.calls)


def test_unique_ids_unchanged_and_missing_dot_gets_stable_key(tmp_path, inputs):
    _, create = inputs
    source = create(records=record("unique") + record(".", 20))
    result, _, _ = prepare(tmp_path, inputs, (source,))
    rows = [json.loads(x) for x in result.ledgers[0].read_text().splitlines()]
    assert rows[0]["analysis_record_id"] == "unique" and not rows[0]["ID_changed"]
    assert rows[1]["ID_changed"] and rows[1]["original_record_id"] == "."
    assert rows[1]["source_record_key"].endswith(":record:2")


def test_cross_caller_collisions_are_aliased_without_merging_records(tmp_path, inputs):
    _, create = inputs
    a, b = create("sawfish", record("same")), create("pav", record("same", 11))
    result, _, _ = prepare(tmp_path, inputs, (a, b))
    rows = [json.loads(p.read_text()) for p in result.ledgers]
    assert {r["original_record_id"] for r in rows} == {"same"}
    assert len({r["analysis_record_id"] for r in rows}) == 2
    assert result.source_order == ("sawfish", "pav")


def test_identical_content_distinct_callers_keep_distinct_record_keys(tmp_path, inputs):
    _, create = inputs
    a, b = create("sawfish", record("same")), create("pav", record("same"))
    assert a.vcf_path.read_bytes() == b.vcf_path.read_bytes()
    result, _, _ = prepare(tmp_path, inputs, (a,b))
    rows = [json.loads(p.read_text()) for p in result.ledgers]
    assert len({x["source_record_key"] for x in rows}) == 2
    assert all(x["source_record_key"].startswith("caller:" + x["caller"] + ":sha256:") for x in rows)


@pytest.mark.parametrize("link", ["MATEID", "PARID"])
def test_ambiguous_mate_parent_stops_and_preserves_failed_receipt(tmp_path, inputs, link):
    _, create = inputs
    source = create(records=record("same") + record("same", 20) + record("third", 30, extra=f";{link}=same"))
    with pytest.raises(InputValidationError, match="Ambiguous"):
        prepare(tmp_path, inputs, (source,))
    assert not (tmp_path / "identity/manifest.json").exists()
    assert json.loads((tmp_path / "identity/FAILED.json").read_text())["success_manifest_published"] is False


def test_unique_mate_ids_not_rewritten(tmp_path, inputs):
    _, create = inputs
    source = create(records=record("a", extra=";MATEID=b") + record("b", 20, extra=";MATEID=a"))
    result, _, _ = prepare(tmp_path, inputs, (source,))
    assert result.validate()["source_summaries"][0]["changed_ID_records"] == 0


@pytest.mark.parametrize("identifier", ["HFVID1_reserved", "a,b", "a;b", "bad id", ""])
def test_illegal_or_reserved_ids_are_not_silently_repaired(tmp_path, inputs, identifier):
    _, create = inputs
    with pytest.raises(InputValidationError):
        prepare(tmp_path, inputs, (create(records=record(identifier)),))


def test_wrong_sample_is_not_relabeled(tmp_path, inputs):
    _, create = inputs
    with pytest.raises((InputValidationError, OutputValidationError)):
        prepare(tmp_path, inputs, (create(sample="OTHER"),))


def test_same_physical_source_cannot_masquerade_as_second_caller(tmp_path, inputs):
    _, create = inputs
    a = create("sawfish", record("a"))
    with pytest.raises(InputValidationError, match="same physical"):
        prepare(tmp_path, inputs, (a, replace(a, caller="pav")))


def test_existing_output_untouched_and_dry_run_has_no_writes(tmp_path, inputs):
    ref, create = inputs
    source = create()
    runner = NativeFixtureRunner(tmp_path)
    wrapper = SourceIdentityWrapper(runner=runner)
    directory = tmp_path / "identity"
    planned = wrapper.run((source,), sample_id="S1", reference=ref, output_directory=directory, policy=policy(), dry_run=True)
    assert planned["executed"] is False and not directory.exists() and runner.calls == []
    directory.mkdir()
    (directory / "old").write_bytes(b"PRESERVE")
    with pytest.raises(OutputValidationError, match="fresh absolute"):
        wrapper.run((source,), sample_id="S1", reference=ref, output_directory=directory, policy=policy())
    assert (directory / "old").read_bytes() == b"PRESERVE" and runner.calls == []


@pytest.mark.parametrize("role", ["view", "index", "ledger", "database", "canonical", "reference"])
def test_tampering_cannot_pass_identity_validation(tmp_path, inputs, role):
    result, sources, _ = prepare(tmp_path, inputs)
    path = {"view": result.views[0], "index": Path(str(result.views[0]) + ".tbi"),
            "ledger": result.ledgers[0], "database": result.database_path,
            "canonical": sources[0].vcf_path, "reference": inputs[0].fasta}[role]
    with path.open("ab") as writer:
        writer.write(b"TAMPER")
    with pytest.raises(OutputValidationError):
        result.validate()


def test_round_trip_rejects_actual_non_id_edit(tmp_path, inputs):
    result, sources, _ = prepare(tmp_path, inputs)
    view = result.views[0]
    changed = gzip.decompress(view.read_bytes()).replace(b"0|1", b"1|1")
    bgzf(view, changed)
    with pytest.raises(OutputValidationError):
        verify_source_view(sources[0].vcf_path, view, result.ledgers[0])


@pytest.mark.parametrize("bad", [None, True, "yes", {"enabled": 1}, {"enabled": False, "unknown": 1},
    {"enabled": False, "policy_id": "x"}, {"enabled": True}, {**policy().to_dict(), "policy_sha256": "wrong"},
    {**policy().to_dict(), "policy_version": True}, {**policy().to_dict(), "timeout_seconds": 0}])
def test_strict_policy_contract(bad):
    with pytest.raises(ConfigurationError):
        SourceIdentityPolicy.from_config(bad)


def test_configuration_default_off_and_enabled_requires_parent(tmp_path):
    import yaml
    assert Phase9Settings.from_config(load_config(DEFAULT)).input_identity_policy is None
    user = tmp_path / "config.yaml"
    user.write_text(yaml.safe_dump({"sv": {"harmonization": {"enabled": False, "input_identity": policy().to_dict()}}}))
    with pytest.raises(ConfigurationError, match="requires"):
        load_config(DEFAULT, user_config=user)
    user.write_text(yaml.safe_dump({"sv": {"harmonization": {"enabled": True, "input_identity": policy().to_dict()}}}))
    assert Phase9Settings.from_config(load_config(DEFAULT, user_config=user)).input_identity_policy == policy()


def test_jasmine_integration_retains_canonical_sources_and_original_ids(tmp_path, inputs):
    ref, create = inputs
    sources = (create("sawfish", record("same")), create("pav", record("same", 11) + record("same", 20)))
    original = {str(s.vcf_path): s.vcf_path.read_bytes() for s in sources}
    request = SVHarmonizationRequest("S1", ref, sources, tmp_path / "work", tmp_path / "out/S1.harmonized.sv.vcf.gz",
                                    tmp_path / "out/S1.sv.evidence.tsv", input_identity_policy=policy())
    runner = NativeFixtureRunner(tmp_path)
    result = JasmineWrapper(runner=runner).run(request)
    assert result.source_identity is not None and result.artifact.sources == sources
    assert result.artifact.source_identity_manifest == result.source_identity.manifest_path
    with request.evidence_table.open() as reader:
        row = next(csv.DictReader(reader, delimiter="\t"))
    assert row["source_record_ids"] == "same,same,same"
    assert len(set(json.loads(row["source_record_keys"]))) == 3
    assert row["source_record_callers"] == "sawfish,pav,pav"
    assert row["evidence_class"] == "READ_AND_ASSEMBLY" and row["support_count"] == "2"
    assert row["source_files"] == ",".join(str(s.vcf_path) for s in sources)
    assert row["representation_files"] == ",".join(map(str, result.source_identity.views))
    assert all(s.vcf_path.read_bytes() == original[str(s.vcf_path)] for s in sources)
    assert result.to_dict()["source_identity"]["canonical_sources_unchanged"] is True


@pytest.mark.parametrize("corruption", ["omit", "duplicate", "unknown", "wrong_vector"])
def test_evidence_never_seals_missing_duplicated_unknown_or_misattributed_members(tmp_path, inputs, corruption):
    ref, create = inputs
    sources = (create("sawfish", record("a")), create("pav", record("b", 20)))
    def transform(merged, ids):
        if corruption == "omit": return merged.replace("," + ids[-1], "")
        if corruption == "duplicate": return merged.replace(ids[0], ids[0] + "," + ids[0])
        if corruption == "unknown": return merged.replace(ids[0], "unknown")
        return merged.replace("SUPP_VEC=11", "SUPP_VEC=10")
    runner = NativeFixtureRunner(tmp_path, merged_transform=transform)
    request = SVHarmonizationRequest("S1", ref, sources, tmp_path / "work", tmp_path / "out/S1.harmonized.sv.vcf.gz",
                                    tmp_path / "out/S1.sv.evidence.tsv", input_identity_policy=policy())
    with pytest.raises(OutputValidationError):
        JasmineWrapper(runner=runner).run(request)
    assert not request.evidence_table.exists()


def test_identity_paths_are_deterministic_and_reject_escape(tmp_path):
    result = phase9_identity_paths(tmp_path, "S1", ("sawfish", "pav"))
    assert result["pav_ledger"] == tmp_path / "S1/source-identity-v1/01.pav.analysis.records.jsonl"
    with pytest.raises(InputValidationError): phase9_identity_paths(tmp_path, "../escape", ("pav",))
    with pytest.raises(InputValidationError): phase9_identity_paths(tmp_path, "S1", ("pav", "pav"))


def test_deterministic_fresh_runs_preserve_order_and_exact_ledger(tmp_path, inputs):
    ref, create = inputs
    sources = (create(),)
    runner = NativeFixtureRunner(tmp_path)
    wrapper = SourceIdentityWrapper(runner=runner)
    a = wrapper.run(sources, sample_id="S1", reference=ref, output_directory=tmp_path / "a", policy=policy())
    b = wrapper.run(sources, sample_id="S1", reference=ref, output_directory=tmp_path / "b", policy=policy())
    assert a.ledgers[0].read_bytes() == b.ledgers[0].read_bytes()
    assert gzip.decompress(a.views[0].read_bytes()) == gzip.decompress(b.views[0].read_bytes())


@pytest.mark.parametrize("position", [-1, 1002])
def test_out_of_reference_coordinates_are_rejected_not_clamped(tmp_path, inputs, position):
    _, create = inputs
    source = create(records=record("v1", position))
    with pytest.raises(InputValidationError, match="bounds"):
        prepare(tmp_path, inputs, (source,))


@pytest.mark.parametrize("position", [0, 1001])
def test_virtual_telomere_coordinates_are_retained(tmp_path, inputs, position):
    _, create = inputs
    source = create(records=record("v1", position))
    result, _, _ = prepare(tmp_path, inputs, (source,))
    assert verify_source_view(source.vcf_path, result.views[0], result.ledgers[0])["records"] == 1


def test_header_reference_length_mismatch_is_not_renamed(tmp_path, inputs):
    _, create = inputs
    source = create()
    bgzf(source.vcf_path, gzip.decompress(source.vcf_path.read_bytes()).replace(b"length=1000", b"length=999"))
    with pytest.raises(InputValidationError, match="name/length"):
        prepare(tmp_path, inputs, (source,))


@pytest.mark.parametrize("fault", ["missing_index", "bad_index", "old_version", "nonexecuted", "command_failure"])
def test_finalization_faults_never_publish_success_manifest(tmp_path, inputs, fault):
    from hifivar.exceptions import CommandExecutionError
    ref, create = inputs
    source = create()
    class FaultRunner(NativeFixtureRunner):
        def run(self, args, **kwargs):
            if fault == "command_failure" and Path(args[0]).name == "bgzip" and "-c" in args:
                raise CommandExecutionError("synthetic bgzip failed")
            result = super().run(args, **kwargs)
            if fault == "old_version" and "--version" in args:
                return replace(result, stdout="tabix 1.21")
            if fault == "nonexecuted":
                return replace(result, executed=False)
            if Path(args[0]).name == "tabix" and list(args[1:3]) == ["-p", "vcf"]:
                index = Path(str(args[-1]) + ".tbi")
                if fault == "missing_index": index.unlink()
                if fault == "bad_index": index.write_bytes(b"invalid index")
            return result
    with pytest.raises((OutputValidationError, CommandExecutionError)):
        SourceIdentityWrapper(runner=FaultRunner(tmp_path)).run((source,), sample_id="S1", reference=ref,
            output_directory=tmp_path / "identity", policy=policy())
    assert not (tmp_path / "identity/manifest.json").exists()
    assert (tmp_path / "identity/FAILED.json").is_file()


def test_source_change_between_inventory_and_materialization_is_rejected(tmp_path, inputs, monkeypatch):
    import hifivar.sv_input_identity as module
    ref, create = inputs
    source = create()
    original = module._derive
    def altered(db, path, *args, **kwargs):
        bgzf(path, gzip.decompress(path.read_bytes()).replace(b"0|1", b"1|1"))
        return original(db, path, *args, **kwargs)
    monkeypatch.setattr(module, "_derive", altered)
    with pytest.raises(OutputValidationError, match="changed between"):
        prepare(tmp_path, inputs, (source,))
    assert not (tmp_path / "identity/manifest.json").exists()


def test_no_calls_keeps_zero_records_distinct_from_failed_source(tmp_path, inputs):
    _, create = inputs
    empty = replace(create(records=""), status=EvidenceRunStatus.NO_CALLS)
    result, _, _ = prepare(tmp_path, inputs, (empty,))
    assert result.validate()["source_summaries"] == [{"records": 0, "changed_ID_records": 0}]
    with pytest.raises(InputValidationError):
        SourceIdentityWrapper().run((replace(empty, status=EvidenceRunStatus.FAILED),),
            sample_id="S1", reference=inputs[0], output_directory=tmp_path / "failed", policy=policy(), dry_run=True)


def test_phase9_explicit_policy_reaches_wrapper_and_preexisting_bundle_blocks_all_tools(tmp_path, inputs):
    from hifivar.context import AnalysisContext
    from hifivar.phase9 import run_phase9
    from hifivar.sample import InputDataset, Sample
    from hifivar.sample_sheet import SampleRecord
    from unittest.mock import Mock
    ref, create = inputs
    source = create()
    bam = tmp_path / "input.bam"
    bam.write_bytes(b"SYNTHETIC")
    context = AnalysisContext(ref, (SampleRecord(Sample("S1", InputDataset.from_files([bam]))),), {})
    cfg = load_config(DEFAULT).to_dict()
    cfg["sv"]["harmonization"].update(enabled=True, input_identity=policy().to_dict())
    wrapped = Mock(wraps=JasmineWrapper(runner=NativeFixtureRunner(tmp_path)))
    compare = Mock()
    compare.run.return_value.to_dict.return_value = {"status": "planned"}
    report = run_phase9(context, {"S1": (source,)}, output_directory=tmp_path / "out", work_directory=tmp_path / "work",
                        config=cfg, jasmine_wrapper=wrapped, truvari_wrapper=compare, dry_run=True)
    assert wrapped.run.call_args.args[0].input_identity_policy == policy()
    assert report.to_dict()["settings"]["input_identity"]["policy_sha256"] == POLICY_SHA256
    assert not (tmp_path / "out").exists()
    protected = tmp_path / "out/S1/source-identity-v1/protected"
    protected.parent.mkdir(parents=True)
    protected.write_bytes(b"KEEP")
    wrapped.reset_mock()
    with pytest.raises(OutputValidationError, match="fresh per-sample"):
        run_phase9(context, {"S1": (source,)}, output_directory=tmp_path / "out", work_directory=tmp_path / "work",
                   config=cfg, jasmine_wrapper=wrapped, truvari_wrapper=compare)
    wrapped.run.assert_not_called()
    assert protected.read_bytes() == b"KEEP"
