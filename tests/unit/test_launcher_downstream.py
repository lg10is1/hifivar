"""Ready external annotation/cohort routing; no calling accuracy is simulated."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import struct
import zlib
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from hifivar import launcher
from hifivar.config import HiFiVarConfig, load_config
from hifivar.exceptions import ConfigurationError, InputValidationError, OutputValidationError, ReferenceError, WorkflowError


DEFAULT = Path(__file__).parents[2] / "src/hifivar/resources/configs/default.yaml"
COHORT_FIELDS = ("sample", "track", "state", "source_path", "index_path", "source_tool", "source_version", "reference_build", "catalog_id")


def bgzf(data: bytes) -> bytes:
    """Tiny structural BGZF fixture, not an independently readable tabix database."""
    compressor = zlib.compressobj(wbits=-15)
    payload = compressor.compress(data) + compressor.flush()
    size = 18 + len(payload) + 8
    header = b"\x1f\x8b\x08\x04" + b"\0" * 4 + b"\0\xff\x06\0BC\x02\0" + struct.pack("<H", size - 1)
    return header + payload + struct.pack("<II", zlib.crc32(data), len(data))


def vcf(path, sample="S1", category="sv", *, contig="chr1", length=1000, build="GRCh38", compressed=True):
    markers = {
        "small": ('##ALT=<ID=NON_REF,Description="Any non-reference allele">\n'
                  '##INFO=<ID=END,Number=1,Type=Integer,Description="Reference-block end">\n'
                  '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Genotype likelihoods">\n'),
        "sv": '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n',
        "tr": '##INFO=<ID=TRID,Number=1,Type=String,Description="Locus">\n##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motif">\n##INFO=<ID=STRUC,Number=1,Type=String,Description="Structure">\n',
    }
    record = (f"{contig}\t100\t.\tA\t<NON_REF>\t.\t.\tEND=150\tGT:PL\t0/0:0,30,300\n"
              if category == "small" else f"{contig}\t100\tv1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=150\tGT\t0/1\n")
    text = (f"##fileformat=VCFv4.3\n##reference={build}\n##contig=<ID={contig},length={length}>\n" + markers[category]
            + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
            + record)
    path.write_bytes(bgzf(text.encode()) if compressed else text.encode())
    if compressed:
        Path(f"{path}.tbi").write_bytes(bgzf(b"TBI\x01" + b"test-structural-index"))
    return path


def write_cohort(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COHORT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def row(sample, track="sv", state="CALLED", **updates):
    result = dict(zip(COHORT_FIELDS, (sample, track, state, f"{sample}.{track}.vcf.gz", f"{sample}.{track}.vcf.gz.tbi", "source-caller", "1.0", "GRCh38", "CAT1" if track == "tr" else "")))
    result.update(updates)
    return result


@pytest.fixture
def ready(tmp_path, monkeypatch):
    (tmp_path / "ref.fa").write_text(">chr1\n" + "A" * 1000 + "\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t1000\t6\t1000\t1001\n", encoding="utf-8")
    for sample in ("S2", "S1"):
        (tmp_path / f"{sample}.bam").write_bytes(b"metadata-context-only-not-a-real-bam")
        for category in ("small", "sv", "tr"):
            vcf(tmp_path / f"{sample}.{category}.vcf.gz", sample, category)
    (tmp_path / "samples.tsv").write_text("sample_id\tinput\nS2\tS2.bam\nS1\tS1.bam\n", encoding="utf-8")
    (tmp_path / "regions.bed").write_text("chr1\t90\t170\tGENE1\n", encoding="utf-8")
    (tmp_path / "annotation.tsv").write_text("sample\tvariant_category\tsource_vcf\tsource_tool\tsource_variant_ids\nS1\tsv\tS1.sv.vcf.gz\tcaller\tv1\n", encoding="utf-8")
    write_cohort(tmp_path / "cohort.tsv", [row("S1"), row("S2")])
    user = tmp_path / "user.yaml"
    user.write_text(yaml.safe_dump({"reference": {"fasta": "ref.fa", "build": "GRCh38"}, "samples": {"sheet": "samples.tsv"},
        "runtime": {"threads": 4}, "workflow": {"memory_mb": 4096},
        "annotation": {"enabled": True, "input_manifest": "annotation.tsv", "overlap_enabled": True,
                       "gene_bed": "regions.bed", "gene_version": "1.0", "threads": 2, "memory_mb": 1000}}), encoding="utf-8")
    workflow = tmp_path / "workflow"
    workflow.mkdir()
    (workflow / "Snakefile").write_text("rule all:\n    input: []\n", encoding="utf-8")
    monkeypatch.setattr(launcher.platform, "system", lambda: "Linux")
    monkeypatch.setattr(launcher, "installed_workflow_root", lambda: workflow)
    return load_config(DEFAULT, user_config=user)


def change(config, **updates):
    values = config.to_dict()
    for section, changes in updates.items():
        values[section].update(changes)
    return HiFiVarConfig(values, config.sources)


def cohort_config(config, track="sv", **updates):
    values = config.to_dict()
    values["annotation"]["enabled"] = False
    values["cohort"].update(enabled=True, cohort_id="C1", input_manifest="cohort.tsv")
    values["cohort"][track].update(enabled=True, threads=1)
    if track != "small_variants":
        values["cohort"][track]["memory_mb"] = 1000
    if track == "small_variants":
        values["cohort"][track]["memory_gb"] = 1
    values["cohort"][track].update(updates)
    return HiFiVarConfig(values, config.sources)


class Scheduler:
    """Only simulate scheduler boundary, never scientific tool outputs."""
    def __init__(self, fail=False):
        self.calls, self.fail = [], fail

    def run(self, command, **kwargs):
        self.calls.append(list(command))
        assert command[0] != "samtools", "Downstream-only run must not probe unused primary BAMs"
        out = "9.26.1\n" if "--version" in command else ""
        code = 1 if self.fail and "--version" not in command else 0
        for key in ("stdout_path", "stderr_path"):
            if kwargs.get(key):
                Path(kwargs[key]).write_text(out if key == "stdout_path" else "", encoding="utf-8")
        return SimpleNamespace(stdout=out, returncode=code, duration_seconds=0.01)


def test_annotation_only_owned_manifest_hashes_and_no_unused_bam_probe(ready, tmp_path):
    original = ready.to_dict()
    source_bytes = (tmp_path / "annotation.tsv").read_bytes()
    effective, context, _, inputs = launcher.prepare_launch(ready, tmp_path / "run")
    assert context.sample_ids == ("S2", "S1")
    assert effective["annotation"]["input_manifest"] == str(tmp_path / "run/annotation.inputs.tsv")
    assert effective["annotation"]["gene_bed"] == str(tmp_path / "regions.bed")
    assert not any(path.suffix == ".bai" for path in inputs)
    assert ready.to_dict() == original
    runner, root = Scheduler(), tmp_path / "run"
    assert launcher.WorkflowLauncher(runner).run(ready, run_dir=root, dry_run=True) == "DRY_RUN_PASS"
    assert (tmp_path / "annotation.tsv").read_bytes() == source_bytes
    copied = list(csv.DictReader(io.StringIO((root / "annotation.inputs.tsv").read_text()), delimiter="\t"))
    assert copied[0]["source_vcf"] == str(tmp_path / "S1.sv.vcf.gz")
    manifest = json.loads((root / "launch_manifest.json").read_text())
    assert "annotation.inputs.tsv" in manifest["control_sha256"]
    actual = {item["path"]: item["sha256"] for item in manifest["identity"]["inputs"]}
    for source in ("annotation.tsv", "regions.bed", "S1.sv.vcf.gz", "S1.sv.vcf.gz.tbi"):
        assert actual[str(tmp_path / source)] == hashlib.sha256((tmp_path / source).read_bytes()).hexdigest()
    audit = json.loads(next((root / "attempts").glob("*/downstream_input_audit.json")).read_text())
    assert audit["accuracy"] == "NOT_ASSESSED"
    assert audit["backend_failure_isolation"] == "NOT_IMPLEMENTED_BY_NATIVE_DAG"
    assert audit["functional_backend"] == "SEPARATE_EXPLICIT_API_CLI_NOT_THIS_DAG"
    assert audit["annotation"][0]["tool_annotation_scope"] == "WHOLE_SOURCE_VCF_NOT_FILTERED_BY_ID_LIST"
    assert not (root / "results/annotation").exists()


def test_cohort_reorders_by_exact_names_and_never_modifies_original(ready, tmp_path):
    source = (tmp_path / "cohort.tsv").read_bytes()
    config, root = cohort_config(ready), tmp_path / "run"
    assert launcher.WorkflowLauncher(Scheduler()).run(config, run_dir=root, dry_run=True) == "DRY_RUN_PASS"
    copied = list(csv.DictReader(io.StringIO((root / "cohort.inputs.tsv").read_text()), delimiter="\t"))
    assert [item["sample"] for item in copied] == ["S2", "S1"]
    assert all(item["source_path"] == str(tmp_path / f"{item['sample']}.sv.vcf.gz") for item in copied)
    assert (tmp_path / "cohort.tsv").read_bytes() == source
    audit = json.loads(next((root / "attempts").glob("*/downstream_input_audit.json")).read_text())
    assert audit["cohort_sv_cross_sample_clustering"] is False
    assert all(item["version_source"] == "USER_DECLARED_NOT_REPROBED" for item in audit["cohort"])


@pytest.mark.parametrize("track,category", [("small_variants", "small"), ("tr", "tr")])
def test_ready_gvcf_and_tr_tracks_are_routed(ready, tmp_path, track, category):
    write_cohort(tmp_path / "cohort.tsv", [row(sample, track, source_path=f"{sample}.{category}.vcf.gz", index_path=f"{sample}.{category}.vcf.gz.tbi") for sample in ("S1", "S2")])
    config = cohort_config(ready, track)
    effective, _, _, inputs = launcher.prepare_launch(config, tmp_path / "run")
    assert effective["cohort"][track]["enabled"]
    assert len([path for path in inputs if str(path).endswith(f".{category}.vcf.gz")]) == 2


@pytest.mark.parametrize("target", ["annotation.tsv", "S1.sv.vcf.gz", "S1.sv.vcf.gz.tbi", "regions.bed", "annotation.inputs.tsv"])
def test_resume_rejects_changed_sources_or_owned_manifest(ready, tmp_path, target):
    root = tmp_path / "run"
    obj = launcher.WorkflowLauncher(Scheduler())
    obj.run(ready, run_dir=root, dry_run=True)
    path = root / target if target == "annotation.inputs.tsv" else tmp_path / target
    with path.open("ab") as handle:
        handle.write(b"\nchanged")
    with pytest.raises((InputValidationError, OutputValidationError), match="RESUME|VCF|index|BED|row|sample|header"):
        obj.run(ready, run_dir=root, resume=True)


@pytest.mark.parametrize("updates,match", [
    ({"sample": "S3"}, "unknown sample"),
    ({"sample": "UnknownSampleName"}, "unknown sample"),
    ({"source_vcf": "not_yet_generated.vcf.gz"}, "exist|missing|file|read"),
    ({"source_variant_ids": ""}, "explicit source_variant_ids"),
    ({"variant_category": "small"}, "no enabled applicable"),
])
def test_annotation_invalid_ready_selection_fails_before_run(ready, tmp_path, updates, match):
    item = {"sample": "S1", "variant_category": "sv", "source_vcf": "S1.sv.vcf.gz", "source_tool": "caller", "source_variant_ids": "v1"}
    item.update(updates)
    fields = tuple(item)
    (tmp_path / "annotation.tsv").write_text(launcher._manifest_text(fields, [item]), encoding="utf-8")
    with pytest.raises(InputValidationError, match=match):
        launcher.prepare_launch(ready, tmp_path / "run")
    assert not (tmp_path / "run").exists()


def test_annotation_empty_manifest_fails(ready, tmp_path):
    (tmp_path / "annotation.tsv").write_text("sample\tvariant_category\tsource_vcf\tsource_tool\tsource_variant_ids\n", encoding="utf-8")
    with pytest.raises((InputValidationError, ReferenceError)):
        launcher.prepare_launch(ready, tmp_path / "run")


@pytest.mark.parametrize("kwargs", [{"sample": "S2"}, {"sample": "UnknownSampleName"}, {"contig": "1"}, {"length": 999}, {"build": "GRCh37"}])
def test_header_sample_contig_length_and_build_fail_closed(ready, tmp_path, kwargs):
    vcf(tmp_path / "S1.sv.vcf.gz", **kwargs)
    with pytest.raises((InputValidationError, ReferenceError)):
        launcher.prepare_launch(ready, tmp_path / "run")


def test_conventional_index_is_required_not_only_csi(ready, tmp_path):
    (tmp_path / "S1.sv.vcf.gz.tbi").rename(tmp_path / "S1.sv.vcf.gz.csi")
    with pytest.raises((InputValidationError, OutputValidationError)):
        launcher.prepare_launch(ready, tmp_path / "run")


def test_plain_annotation_vcf_is_supported_without_fake_index(ready, tmp_path):
    vcf(tmp_path / "plain.vcf", compressed=False)
    (tmp_path / "annotation.tsv").write_text("sample\tvariant_category\tsource_vcf\tsource_tool\tsource_variant_ids\nS1\tsv\tplain.vcf\tcaller\tv1\n", encoding="utf-8")
    _, _, _, inputs = launcher.prepare_launch(ready, tmp_path / "run")
    assert tmp_path / "plain.vcf" in inputs
    assert tmp_path / "plain.vcf.tbi" not in inputs


@pytest.mark.parametrize("content", ["1\t0\t10\n", "chr1\t0\t1001\n", "chr1\t10\t10\n", "chr1\t-1\t10\n", "# empty\n"])
def test_overlap_bed_reference_bounds_and_nonempty_are_checked(ready, tmp_path, content):
    (tmp_path / "regions.bed").write_text(content, encoding="utf-8")
    with pytest.raises(InputValidationError):
        launcher.prepare_launch(ready, tmp_path / "run")


def test_overlap_bed_need_not_be_sorted(ready, tmp_path):
    (tmp_path / "regions.bed").write_text("chr1\t100\t200\nchr1\t0\t50\n", encoding="utf-8")
    launcher.prepare_launch(ready, tmp_path / "run")


@pytest.mark.parametrize("rows,match", [
    ([row("S1")], "every declared sample"),
    ([row("S1"), row("S2"), row("S3")], "unknown sample"),
    ([row("S1"), row("S1"), row("S2")], "duplicate"),
    ([row("S1", "unexpected"), row("S2")], "unknown sample or track"),
    ([row("S1", state="PENDING"), row("S2")], "Unknown cohort sample state"),
    ([row("S1", source_path="future.vcf.gz", index_path="future.vcf.gz.tbi"), row("S2")], "exist|missing|file|read"),
    ([row("S1", reference_build="GRCh37"), row("S2")], "build|reference"),
    ([row("S1", source_version=""), row("S2")], "source version"),
    ([row("S1", state="FAILED"), row("S2")], "Noncallable"),
    ([row(s, state="NOT_RUN", source_path="", index_path="") for s in ("S1", "S2")], "no ready callable"),
])
def test_cohort_invalid_manifest_never_creates_run(ready, tmp_path, rows, match):
    write_cohort(tmp_path / "cohort.tsv", rows)
    with pytest.raises(InputValidationError, match=match):
        launcher.prepare_launch(cohort_config(ready), tmp_path / "run")
    assert not (tmp_path / "run").exists()


def test_noncallable_sv_sample_state_preserved_not_turned_homref(ready, tmp_path):
    write_cohort(tmp_path / "cohort.tsv", [row("S1"), row("S2", state="FAILED", source_path="", index_path="")])
    config = cohort_config(ready)
    root = tmp_path / "run"
    launcher.WorkflowLauncher(Scheduler()).run(config, run_dir=root, dry_run=True)
    rows = list(csv.DictReader(io.StringIO((root / "cohort.inputs.tsv").read_text()), delimiter="\t"))
    assert rows[0]["sample"] == "S2" and rows[0]["state"] == "FAILED" and rows[0]["source_path"] == ""


def test_glnexus_never_drops_noncallable_sample(ready, tmp_path):
    write_cohort(tmp_path / "cohort.tsv", [row("S1", "small_variants", source_path="S1.small.vcf.gz", index_path="S1.small.vcf.gz.tbi"), row("S2", "small_variants", state="FAILED", source_path="", index_path="")])
    with pytest.raises(InputValidationError, match="every cohort sample"):
        launcher.prepare_launch(cohort_config(ready, "small_variants"), tmp_path / "run")


@pytest.mark.parametrize("track,category", [("small_variants", "sv"), ("tr", "sv"), ("sv", "small")])
def test_track_specific_header_contract(ready, tmp_path, track, category):
    write_cohort(tmp_path / "cohort.tsv", [row(s, track, source_path=f"{s}.{category}.vcf.gz", index_path=f"{s}.{category}.vcf.gz.tbi") for s in ("S1", "S2")])
    with pytest.raises(InputValidationError, match="required declarations"):
        launcher.prepare_launch(cohort_config(ready, track), tmp_path / "run")


def test_tr_must_use_identical_explicit_catalog(ready, tmp_path):
    write_cohort(tmp_path / "cohort.tsv", [row("S1", "tr"), row("S2", "tr", catalog_id="OTHER")])
    with pytest.raises(InputValidationError, match="identical explicit catalog_id"):
        launcher.prepare_launch(cohort_config(ready, "tr"), tmp_path / "run")


@pytest.mark.parametrize("track,updates,match", [("sv", {"threads": 2}, "threads must be 1"), ("sv", {"memory_mb": 5000}, "memory_mb exceeds"), ("small_variants", {"memory_gb": 5}, "memory_mb exceeds"), ("small_variants", {"threads": 5}, "threads exceeds")])
def test_downstream_resource_budgets(ready, tmp_path, track, updates, match):
    with pytest.raises(ConfigurationError, match=match):
        launcher.prepare_launch(cohort_config(ready, track, **updates), tmp_path / "run")


def test_no_empty_cohort_track_can_report_success(ready, tmp_path):
    config = cohort_config(ready)
    value = config.to_dict()
    value["cohort"]["sv"]["enabled"] = False
    with pytest.raises(ConfigurationError, match="at least one enabled cohort track"):
        launcher.prepare_launch(HiFiVarConfig(value, config.sources), tmp_path / "run")


@pytest.mark.parametrize("section", ["annotation", "cohort"])
def test_downstream_overwrite_not_allowed(ready, tmp_path, section):
    config = ready if section == "annotation" else cohort_config(ready)
    with pytest.raises(ConfigurationError, match="overwrite=false"):
        launcher.prepare_launch(change(config, **{section: {"overwrite": True}}), tmp_path / "run")


def test_alphagenome_requires_separate_explicit_cli_not_main_dag(ready):
    value = ready.to_dict()
    value["annotation"]["functional_enabled"] = True
    with pytest.raises(ConfigurationError, match="hifivar functional alphagenome.*explicit plan and user credential.*not an automatic main-run route"):
        launcher._enabled_tracks(value)


def test_failed_scheduler_cannot_become_success_or_independent_backend_isolation(ready, tmp_path):
    root = tmp_path / "run"
    with pytest.raises(WorkflowError, match="No automatic retry"):
        launcher.WorkflowLauncher(Scheduler(fail=True)).run(ready, run_dir=root)
    status = json.loads(next((root / "attempts").glob("*/status.json")).read_text())
    assert status["status"] == "FAILED" and status["returncode"] == 1
    assert not (root / "results/annotation").exists()


@pytest.mark.parametrize('invalid', [None, 'missing_index', 'wrong_sample', 'wrong_contig_length'])
def test_ready_harmonization_uses_complete_analysis_context(ready, invalid):
    cfg = ready
    root = Path(cfg.sources['user']).parent
    data = cfg.to_dict()
    data['annotation']['enabled'] = False
    data['sv']['enabled'] = False
    external = {}
    for sample in ('S2', 'S1'):
        path = root / f'{sample}.pav.assembly.sv.vcf.gz'
        vcf(path, sample='DIFFERENT' if invalid == 'wrong_sample' and sample == 'S1' else sample,
            length=999 if invalid == 'wrong_contig_length' and sample == 'S1' else 1000)
        external[sample] = {'pav': str(path)}
    if invalid == 'missing_index':
        Path(str(root / 'S1.pav.assembly.sv.vcf.gz') + '.tbi').unlink()
    data['sv']['harmonization'].update(enabled=True, threads=2, memory_mb=1000, input_vcfs=external)
    config = HiFiVarConfig(data, cfg.sources)
    destination = root / 'harmonization_run'
    if invalid is not None:
        with pytest.raises((InputValidationError, ReferenceError)):
            launcher.prepare_launch(config, destination)
    else:
        effective, context, _, inputs = launcher.prepare_launch(config, destination)
        assert context.sample_ids == ('S2', 'S1')
        for sample in context.sample_ids:
            selected = Path(effective['sv']['harmonization']['input_vcfs'][sample]['pav'])
            assert selected in inputs and Path(str(selected) + '.tbi') in inputs
    assert not destination.exists()
