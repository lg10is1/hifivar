"""Generated intent is not a called file; receipts bind explicit downstream sources.

These unit fixtures use tiny structural BGZF/TBI bytes and real isolated Python
children, not real variant callers or an independently readable tabix database.
They do not establish annotation, genotyping or calling accuracy.
"""

from __future__ import annotations

import base64
import copy
import csv
import hashlib
import io
import json
import os
import struct
import sys
import zlib
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

import hifivar.attempt_receipts as receipts
import hifivar.launcher as launcher
import hifivar.package_resources as resources
from hifivar.command import CommandRunner
from hifivar.downstream_handoff import (
    materialize_downstream, prepare_downstream_intent, producer_config,
)
from hifivar.exceptions import InputValidationError, OutputValidationError, ReferenceError
from hifivar.reference import ReferenceGenome


DEFAULT = Path(__file__).parents[2] / "src/hifivar/resources/configs/default.yaml"


def _bgzf(data: bytes) -> bytes:
    compressor = zlib.compressobj(wbits=-15)
    compressed = compressor.compress(data) + compressor.flush()
    size = 18 + len(compressed) + 8
    header = (b"\x1f\x8b\x08\x04" + b"\0" * 4 + b"\0\xff\x06\0BC\x02\0"
              + struct.pack("<H", size - 1))
    return header + compressed + struct.pack("<II", zlib.crc32(data), len(data))


def _vcf(sample: str, *, category: str = "sv", records: bool = True,
         build: str = "GRCh38", ordinary: bool = False) -> bytes:
    declarations = ('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    if category == "sv":
        declarations += '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
        record = "chr1\t100\tv1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL\tGT\t0/1\n"
    elif category == "tr":
        declarations += ('##INFO=<ID=TRID,Number=1,Type=String,Description="Locus">\n'
                         '##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motifs">\n'
                         '##INFO=<ID=STRUC,Number=1,Type=String,Description="Structure">\n')
        record = "chr1\t100\tv1\tA\t<TR>\t.\tPASS\tTRID=L1;MOTIFS=A;STRUC=(A)n\tGT\t0/1\n"
    else:
        declarations += ('##ALT=<ID=NON_REF,Description="Any non-reference allele">\n'
                         '##INFO=<ID=END,Number=1,Type=Integer,Description="Block end">\n'
                         '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">\n')
        record = ("chr1\t100\tv1\tA\tC\t.\tPASS\t.\tGT:PL\t0/1:90,0,90\n" if ordinary
                  else "chr1\t100\t.\tA\t<NON_REF>\t.\t.\tEND=150\tGT:PL\t0/0:0,30,300\n")
    text = (f"##fileformat=VCFv4.3\n##reference={build}\n##contig=<ID=chr1,length=1000>\n"
            + declarations + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
            + (record if records else ""))
    return _bgzf(text.encode("utf-8"))


def _row(sample: str, *, track: str = "sv", producer: str = "sv.sawfish", role: str = "vcf") -> dict:
    return {"sample_id": sample, "track": track, "origin": "GENERATED",
            "producer": producer, "output_role": role, "reference_build": "GRCh38"}


@pytest.fixture
def handoff_fixture(tmp_path: Path, monkeypatch):
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\n" + "A" * 1000 + "\n", encoding="utf-8")
    Path(f"{fasta}.fai").write_text("chr1\t1000\t6\t1000\t1001\n", encoding="utf-8")
    reference = ReferenceGenome.from_fasta(fasta, build="GRCh38")
    context = SimpleNamespace(reference=reference, sample_ids=("S2", "S1"), n_samples=2)
    config = yaml.safe_load(DEFAULT.read_text(encoding="utf-8"))
    config["sv"]["enabled"] = True
    config["small"]["enabled"] = True
    config["cohort"].update(enabled=True, cohort_id="C1")
    config["cohort"]["sv"]["enabled"] = True
    config["cohort"]["small_variants"].update(threads=1, memory_gb=1)
    config["runtime"]["threads"] = 4
    config["workflow"]["memory_mb"] = 64000
    root = tmp_path / "run"
    runtime = {"payload_sha256": "b" * 64, "hifivar_version": "synthetic-unit-runtime"}
    monkeypatch.setattr(launcher, "runtime_identity", lambda _: runtime)
    monkeypatch.setattr(resources, "installed_workflow_root", lambda: root)
    proxy = SimpleNamespace(environ={}, fdopen=os.fdopen, fsync=os.fsync, link=os.link)
    monkeypatch.setattr(receipts, "os", proxy)

    def freeze(document=None, *, dry_run=False):
        root.mkdir()
        source = tmp_path / "source_plan.yaml"
        document = document or {"schema_version": 1, "cohort": [_row("S1"), _row("S2")]}
        source.write_text(yaml.safe_dump(document), encoding="utf-8")
        config["workflow"]["downstream_source_plan"] = str(source)
        controls, inputs, intent = prepare_downstream_intent(config, context, tmp_path, root)
        for name, text in controls.items():
            (root / name).write_text(text, encoding="utf-8")
        (root / "effective.yaml").write_text(yaml.safe_dump(config), encoding="utf-8")
        sheet = "sample_id\tinput\nS2\t/external/S2.bam\nS1\t/external/S1.bam\n"
        (root / "samples.tsv").write_text(sheet, encoding="utf-8")
        manifest = {
            "identity": {"schema": 2, "run_dir": str(root), "runtime": runtime,
                         "sample_sheet": sheet, "config": config},
            "control_sha256": {name: receipts.artifact_identity(root / name)["sha256"]
                               for name in ("effective.yaml", "samples.tsv", "downstream_plan.yaml")},
        }
        (root / "launch_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        attempt = root / "attempts" / "attempt_a"
        attempt.mkdir(parents=True)
        (attempt / "status.json").write_text(
            json.dumps({"status": "PRODUCER_DAG_RUNNING", "dry_run": dry_run}), encoding="utf-8",
        )
        context_path = receipts.create_attempt_context(root, attempt, dry_run=dry_run)
        proxy.environ = {receipts.ATTEMPT_CONTEXT_ENV: str(context_path)}
        return SimpleNamespace(root=root, attempt=attempt, config=config, context=context,
                               source=source, inputs=inputs, intent=intent, proxy=proxy)

    freeze.config = config
    freeze.context = context
    freeze.root = root
    return freeze


def _seal(frozen, sample: str, *, producer="sv.sawfish", empty=False,
          wrong_header_sample=None, header_build="GRCh38", receipt_build="GRCh38",
          ordinary_gvcf=False, version="2.2.1", catalog_id=None, output_path=None):
    track, tool = {
        "small.deepvariant": ("small", "deepvariant"),
        "tr.trgt": ("tr", "trgt"),
    }.get(producer, (producer, producer.split(".")[1]))
    output = output_path or (frozen.root / "results" / sample / (producer.replace(".", "_") + ".vcf.gz"))
    outputs = {"vcf": output, "index": Path(f"{output}.tbi")}
    if track == "small":
        outputs.update(gvcf=output.with_name("deepvariant.g.vcf.gz"),
                       gvcf_index=output.with_name("deepvariant.g.vcf.gz.tbi"))
    if track == "tr":
        outputs.update(spanning_bam=output.with_name("spanning.bam"),
                       spanning_bam_index=output.with_name("spanning.bam.bai"))
    payloads = {}
    for role, path in outputs.items():
        if role.endswith("index") or role == "index":
            payload = _bgzf(b"TBI\x01synthetic-structural-index")
        elif role == "spanning_bam":
            payload = b"synthetic-unused-spanning-bam-fixture"
        else:
            category = track if track in {"small", "tr"} else "sv"
            payload = _vcf(wrong_header_sample or sample, category=category,
                           records=not empty, build=header_build,
                           ordinary=ordinary_gvcf if role == "gvcf" else track == "small")
        payloads[str(path)] = base64.b64encode(payload).decode("ascii")
    with receipts.attempt_step(track, (sample,), outputs=outputs) as receipt:
        CommandRunner().run([
            sys.executable, "-I", "-B", "-c",
            "import json,base64,sys;from pathlib import Path;"
            "d=json.loads(sys.argv[1]);"
            "[(Path(p).parent.mkdir(parents=True,exist_ok=True),"
            "Path(p).write_bytes(base64.b64decode(b))) for p,b in d.items()]",
            json.dumps(payloads),
        ], inherit_env=False, env={"SystemRoot": r"C:\Windows"} if os.name == "nt" else {}, timeout=20)
        receipt.complete(
            outputs=outputs, tool_versions={tool: version},
            validation={"validated": True, "scope": "SYNTHETIC_UNIT_FIXTURE_NOT_SCIENTIFIC_VALIDATION"},
            metadata={"reference_build": receipt_build, "catalog_id": catalog_id,
                      "tool_version_sources": {tool: "UNIT_COMMAND_RECORDED"}},
        )
    return outputs


def _table(path: Path) -> list[dict]:
    return list(csv.DictReader(io.StringIO(path.read_text(encoding="utf-8")), delimiter="\t"))


def test_generated_handoff_requires_selected_attempt_receipts_and_preserves_sources(handoff_fixture):
    frozen = handoff_fixture()
    outputs = [_seal(frozen, sample) for sample in ("S1", "S2")]
    watched = [frozen.source, frozen.root / "downstream_plan.yaml", frozen.root / "effective.yaml",
               frozen.context.reference.fasta, frozen.context.reference.fai]
    watched += [path for result in outputs for path in result.values()]
    before = {str(path): receipts.artifact_identity(path) for path in watched}
    path, record = materialize_downstream(frozen.config, frozen.context, frozen.root, frozen.attempt)
    rows = _table(frozen.attempt / "downstream/cohort.inputs.tsv")
    assert [row["sample"] for row in rows] == ["S2", "S1"]
    assert all(row["state"] == "CALLED" and row["source_tool"] == "sawfish" for row in rows)
    assert all(row["source_version"] == "2.2.1" for row in rows)
    assert len(record["source_receipts"]) == 2
    assert all(link["producer_attempt_id"] == frozen.attempt.name for link in record["source_receipts"])
    assert all(link["version_source"] == "UNIT_COMMAND_RECORDED" for link in record["source_receipts"])
    assert all(row["version_source"] == "SELECTED_ATTEMPT_PRODUCER_RECEIPT"
               for row in record["source_validation"]["cohort"])
    assert record["accuracy"] == "NOT_ASSESSED" and not record["implicit_caller_choice"]
    assert not record["ordinary_vcf_substituted_for_gvcf"] and not record["cross_sample_sv_clustering"]
    effective = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert effective["cohort"]["enabled"] and not effective["small"]["enabled"]
    assert effective["workflow"]["downstream_source_plan"] is None
    assert before == {str(path): receipts.artifact_identity(path) for path in watched}


def test_absent_receipt_does_not_turn_existing_vcf_into_called_source(handoff_fixture):
    frozen = handoff_fixture()
    output = frozen.root / "results" / "S1.sv.vcf.gz"
    output.parent.mkdir()
    output.write_bytes(_vcf("S1"))
    Path(f"{output}.tbi").write_bytes(_bgzf(b"TBI\x01synthetic-index"))
    with pytest.raises(InputValidationError, match="execution receipt"):
        materialize_downstream(frozen.config, frozen.context, frozen.root, frozen.attempt)
    assert not (frozen.attempt / "downstream").exists()


@pytest.mark.parametrize("error", ["receipt_build", "header_build", "header_sample"])
def test_receipt_and_vcf_reference_and_sample_must_all_match(handoff_fixture, error):
    frozen = handoff_fixture()
    updates = {"receipt_build": "GRCh37"} if error == "receipt_build" else (
        {"header_build": "GRCh37"} if error == "header_build" else {"wrong_header_sample": "Other"})
    _seal(frozen, "S1", **updates)
    _seal(frozen, "S2")
    with pytest.raises((InputValidationError, OutputValidationError, ReferenceError)):
        materialize_downstream(frozen.config, frozen.context, frozen.root, frozen.attempt)


def test_missing_sample_is_not_silently_dropped_from_generated_cohort(handoff_fixture):
    with pytest.raises(InputValidationError, match="every declared sample"):
        handoff_fixture({"schema_version": 1, "cohort": [_row("S1")]})
    assert not (handoff_fixture.root / "downstream_plan.yaml").exists()


def test_missing_explicit_caller_selector_is_not_auto_chosen(handoff_fixture):
    rows = [_row("S1"), _row("S2")]
    del rows[0]["producer"]
    with pytest.raises(InputValidationError, match="producer"):
        handoff_fixture({"schema_version": 1, "cohort": rows})


def test_requested_caller_cannot_be_replaced_by_other_caller_receipt(handoff_fixture):
    frozen = handoff_fixture()
    for sample in ("S1", "S2"):
        _seal(frozen, sample, producer="sv.sniffles2")
    with pytest.raises(InputValidationError, match="sv.sawfish"):
        materialize_downstream(frozen.config, frozen.context, frozen.root, frozen.attempt)


def test_valid_header_only_source_is_no_calls_not_missing_or_homref(handoff_fixture):
    frozen = handoff_fixture()
    _seal(frozen, "S1", empty=True)
    _seal(frozen, "S2")
    _, record = materialize_downstream(frozen.config, frozen.context, frozen.root, frozen.attempt)
    rows = _table(frozen.attempt / "downstream/cohort.inputs.tsv")
    states = {row["sample"]: row["state"] for row in rows}
    assert states == {"S2": "CALLED", "S1": "NO_CALLS"}
    assert record["accuracy"] == "NOT_ASSESSED"


@pytest.mark.parametrize("ordinary", [False, True])
def test_cohort_small_requires_real_reference_block_not_vcf_label(handoff_fixture, ordinary):
    handoff_fixture.config["cohort"]["sv"]["enabled"] = False
    handoff_fixture.config["cohort"]["small_variants"]["enabled"] = True
    rows = [_row(sample, track="small_variants", producer="small.deepvariant", role="gvcf")
            for sample in ("S1", "S2")]
    frozen = handoff_fixture({"schema_version": 1, "cohort": rows})
    for sample in ("S1", "S2"):
        _seal(frozen, sample, producer="small.deepvariant", ordinary_gvcf=ordinary, version="1.10.0")
    if ordinary:
        with pytest.raises(InputValidationError, match="GVCF_STRUCTURE_UNVERIFIED"):
            materialize_downstream(frozen.config, frozen.context, frozen.root, frozen.attempt)
    else:
        _, record = materialize_downstream(frozen.config, frozen.context, frozen.root, frozen.attempt)
        assert all(link["output_role"] == "gvcf" for link in record["source_receipts"])
        assert all(row["source_path"].endswith(".g.vcf.gz") for row in _table(frozen.attempt / "downstream/cohort.inputs.tsv"))


def test_ordinary_vcf_role_rejected_in_cohort_small_plan(handoff_fixture):
    handoff_fixture.config["cohort"]["sv"]["enabled"] = False
    handoff_fixture.config["cohort"]["small_variants"]["enabled"] = True
    rows = [_row(sample, track="small_variants", producer="small.deepvariant", role="vcf")
            for sample in ("S1", "S2")]
    with pytest.raises(InputValidationError, match="gVCF"):
        handoff_fixture({"schema_version": 1, "cohort": rows})


def test_dry_run_plan_cannot_materialize_generated_sources(handoff_fixture):
    frozen = handoff_fixture(dry_run=True)
    with pytest.raises(InputValidationError, match="execution receipt"):
        materialize_downstream(frozen.config, frozen.context, frozen.root, frozen.attempt)
    assert not (frozen.attempt / "downstream").exists()


def test_materialization_never_overwrites_an_existing_destination(handoff_fixture):
    frozen = handoff_fixture()
    for sample in ("S1", "S2"):
        _seal(frozen, sample)
    materialize_downstream(frozen.config, frozen.context, frozen.root, frozen.attempt)
    target = frozen.attempt / "downstream/materialization.json"
    before = target.read_bytes()
    with pytest.raises(OutputValidationError, match="already exists"):
        materialize_downstream(frozen.config, frozen.context, frozen.root, frozen.attempt)
    assert target.read_bytes() == before


def test_producer_configuration_only_disables_consumers_without_mutating_original(handoff_fixture):
    original = copy.deepcopy(handoff_fixture.config)
    result = producer_config(original)
    assert original == handoff_fixture.config
    assert not result["annotation"]["enabled"] and not result["cohort"]["enabled"]
    assert result["small"] == original["small"] and result["sv"] == original["sv"]
    assert result["workflow"]["downstream_source_plan"] is None


def test_annotation_preserves_explicit_category_and_receipt_version_lineage(handoff_fixture, tmp_path):
    config = handoff_fixture.config
    config["cohort"]["enabled"] = False
    config["cohort"]["sv"]["enabled"] = False
    regions = tmp_path / "regions.bed"
    regions.write_text("chr1\t90\t170\tGENE1\n", encoding="utf-8")
    config["annotation"].update(enabled=True, overlap_enabled=True,
                                 gene_bed=str(regions), gene_version="unit-database-1")
    row = _row("S1")
    row["variant_ids"] = ["v1"]
    frozen = handoff_fixture({"schema_version": 1, "annotation": [row]})
    output = _seal(frozen, "S1", version="2.2.1-unit")
    _, record = materialize_downstream(config, frozen.context, frozen.root, frozen.attempt)
    rows = _table(frozen.attempt / "downstream/annotation.inputs.tsv")
    assert len(rows) == 1
    assert rows[0]["sample"] == "S1" and rows[0]["variant_category"] == "sv"
    assert rows[0]["source_variant_ids"] == "v1" and rows[0]["source_vcf"] == str(output["vcf"])
    assert record["source_receipts"][0]["source_version"] == "2.2.1-unit"
    assert record["source_validation"]["annotation"][0]["source_origin"] == "GENERATED"


@pytest.mark.parametrize("mismatch", [False, True])
def test_generated_tr_requires_same_catalog_identity_in_plan_and_receipt(handoff_fixture, tmp_path, mismatch):
    config = handoff_fixture.config
    config["cohort"]["sv"]["enabled"] = False
    config["cohort"]["tr"]["enabled"] = True
    catalog = tmp_path / "catalog.bed"
    catalog.write_text("chr1\t90\t170\tID=L1;MOTIFS=A;STRUC=(A)n\n", encoding="utf-8")
    config["tr"].update(enabled=True, catalog=str(catalog), catalog_reference_build="GRCh38")
    identity = "sha256:" + hashlib.sha256(catalog.read_bytes()).hexdigest()
    rows = [_row(sample, track="tr", producer="tr.trgt") for sample in ("S1", "S2")]
    frozen = handoff_fixture({"schema_version": 1, "cohort": rows})
    for sample in ("S1", "S2"):
        _seal(frozen, sample, producer="tr.trgt", catalog_id="sha256:" + "f" * 64 if mismatch else identity)
    if mismatch:
        with pytest.raises(InputValidationError, match="catalog checksum"):
            materialize_downstream(config, frozen.context, frozen.root, frozen.attempt)
        assert not (frozen.attempt / "downstream").exists()
    else:
        _, record = materialize_downstream(config, frozen.context, frozen.root, frozen.attempt)
        assert all(row["catalog_id"] == identity for row in _table(frozen.attempt / "downstream/cohort.inputs.tsv"))
        assert len(record["source_receipts"]) == 2
