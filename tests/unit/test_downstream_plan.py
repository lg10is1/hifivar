"""Pure source-intent parsing; no caller execution or VCF-content certification."""

from __future__ import annotations

import hashlib
from copy import deepcopy
from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest
import yaml

from hifivar.command import CommandRunner
from hifivar.config import load_config, validate_config
from hifivar.downstream_plan import read_downstream_plan
from hifivar.exceptions import InputValidationError


DEFAULT = Path(__file__).parents[2] / "src/hifivar/resources/configs/default.yaml"


def configuration():
    config = load_config(DEFAULT).to_dict()
    config["annotation"].update(enabled=True, input_manifest="provided-by-launcher.tsv", vep_enabled=True,
                                vep_cache_directory="cache", vep_cache_version="116")
    config["small"]["enabled"] = True
    config["sv"]["enabled"] = True
    config["sv"]["harmonization"]["enabled"] = True
    config["tr"]["enabled"] = True
    config["phasing"]["enabled"] = True
    config["assembly"]["enabled"] = True
    config["assembly_sv"]["enabled"] = True
    return config


def generated(sample="S1", track="sv", producer="sv.sawfish", role="vcf", **updates):
    row = dict(sample_id=sample, track=track, origin="GENERATED", producer=producer, output_role=role)
    row.update(updates)
    return row


def external(tmp_path, sample="S1", track="sv", **updates):
    source = tmp_path / f"{sample}.{track}.vcf.gz"
    index = tmp_path / f"{sample}.{track}.vcf.gz.tbi"
    source.write_bytes(b"Existing path fixture, deliberately not a validated VCF")
    index.write_bytes(b"Existing index fixture, deliberately not a validated index")
    row = dict(sample_id=sample, track=track, origin="EXTERNAL", path=source.name,
               index_path=index.name, source_tool="sawfish", source_version="2.2.1", reference_build="GRCh38")
    row.update(updates)
    return row


def plan(tmp_path, *, annotation=None, cohort=None, **extra):
    document = {"schema_version": 1, "annotation": annotation or [], "cohort": cohort or []}
    document.update(extra)
    path = tmp_path / "source_plan.yaml"
    path.write_text(yaml.safe_dump(document, allow_unicode=True, sort_keys=False), encoding="utf-8")
    return path


def read(path, config=None, sample_ids=("S1", "S2"), reference_build="GRCh38"):
    return read_downstream_plan(path, sample_ids=sample_ids, reference_build=reference_build,
                                config=config if config is not None else configuration())


def cohort_enabled(config, *tracks):
    config["cohort"].update(enabled=True, cohort_id="C1", input_manifest="materialized-later.tsv")
    for track in tracks:
        config["cohort"][track]["enabled"] = True
    config["cohort"]["small_variants"]["memory_gb"] = 1
    return config


def test_generated_plan_is_frozen_pure_and_has_exact_input_hash(tmp_path, monkeypatch):
    config = configuration()
    original = deepcopy(config)
    source = plan(tmp_path, annotation=[generated(variant_ids=["v1", "v2"])])
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    monkeypatch.setattr(CommandRunner, "run", lambda *a, **k: pytest.fail("Plan parsing must not run commands"))
    result = read(source, config)
    assert result.schema_version == 1 and result.path == source.resolve()
    assert result.sha256 == hashlib.sha256(before[source.name]).hexdigest()
    assert result.annotation[0].producer == "sv.sawfish"
    assert result.annotation[0].path is None and result.annotation[0].source_version is None
    assert result.annotation[0].variant_ids == ("v1", "v2")
    assert isinstance(result.annotation, tuple) and result.cohort == ()
    with pytest.raises(FrozenInstanceError):
        result.schema_version = 2
    with pytest.raises(FrozenInstanceError):
        result.annotation[0].state = "FAILED"
    assert config == original
    assert {p.name: p.read_bytes() for p in tmp_path.iterdir()} == before


def test_compatible_real_config_and_empty_optional_groups(tmp_path):
    config = configuration()
    config["tr"].update(catalog="catalog.bed", catalog_reference_build="GRCh38")
    validate_config(config)
    path = tmp_path / "empty.yaml"
    path.write_text("schema_version: 1\n", encoding="utf-8")
    result = read(path, config)
    assert result.annotation == result.cohort == ()


@pytest.mark.parametrize("track,producer,role", [
    ("small", "small.deepvariant", "vcf"),
    ("small", "phasing.hiphase", "vcf"),
    ("sv", "sv.sawfish", "vcf"), ("sv", "sv.sniffles2", "vcf"),
    ("sv", "sv.pbsv", "vcf"), ("sv", "sv.cutesv", "vcf"),
    ("sv", "sv.harmonization", "vcf"),
    ("sv", "assembly_sv.pav", "vcf"), ("sv", "assembly_sv.svim_asm", "vcf"),
])
def test_all_supported_non_tr_generated_annotation_selectors(tmp_path, track, producer, role):
    result = read(plan(tmp_path, annotation=[generated(track=track, producer=producer, role=role)]))
    assert result.annotation[0].producer == producer


def test_generated_gvcf_covers_all_cohort_samples_without_future_files(tmp_path):
    config = cohort_enabled(configuration(), "small_variants")
    result = read(plan(tmp_path, cohort=[generated(s, "small_variants", "small.deepvariant", "gvcf")
                                       for s in ("S2", "S1")]), config)
    assert [row.sample_id for row in result.cohort] == ["S2", "S1"]
    assert all(row.path is None and row.output_role == "gvcf" for row in result.cohort)


@pytest.mark.parametrize("track,producer,role", [
    ("small", "sv.sawfish", "vcf"), ("sv", "small.deepvariant", "vcf"),
    ("small", "small.deepvariant", "gvcf"), ("sv", "sv.sawfish", "gvcf"),
    ("sv", "unknown.caller", "vcf"), ("sv", "sv.sawfish", "bam"),
])
def test_invalid_generated_category_and_role_fail(tmp_path, track, producer, role):
    with pytest.raises(InputValidationError):
        read(plan(tmp_path, annotation=[generated(track=track, producer=producer, role=role)]))


@pytest.mark.parametrize("producer,disabled", [
    ("small.deepvariant", ("small",)), ("phasing.hiphase", ("phasing",)),
    ("phasing.hiphase", ("small",)), ("sv.sawfish", ("sv",)),
    ("sv.sawfish", ("sv", "sawfish")), ("sv.harmonization", ("sv", "harmonization")),
    ("assembly_sv.pav", ("assembly_sv", "pav")), ("assembly_sv.svim_asm", ("assembly",)),
])
def test_generated_producers_require_enabled_config(tmp_path, producer, disabled):
    config = configuration()
    section = config[disabled[0]]
    if len(disabled) == 2:
        section = section[disabled[1]]
    section["enabled"] = False
    track = "small" if producer.startswith(("small.", "phasing.")) else "sv"
    with pytest.raises(InputValidationError, match="disabled|enabled"):
        read(plan(tmp_path, annotation=[generated(track=track, producer=producer)]), config)


@pytest.mark.parametrize("producer", ["assembly_sv.pav", "assembly_sv.svim_asm"])
@pytest.mark.parametrize("destination", ["annotation", "cohort"])
def test_generated_ready_assembly_caller_output_is_valid_source_intent(tmp_path, monkeypatch, producer, destination):
    config = configuration()
    config["assembly"]["enabled"] = False
    # Parsing intent is pure. The launcher must validate this real manifest
    # and each FASTA later; the plan does not certify future caller outputs.
    config["assembly_sv"]["input_manifest"] = "ready.tsv"
    if destination == "cohort":
        cohort_enabled(config, "sv")
    before = deepcopy(config)
    monkeypatch.setattr(CommandRunner, "run", lambda *a, **k: pytest.fail("Intent parsing cannot execute a caller"))
    path = plan(tmp_path, **{destination: [generated(producer=producer)]})
    result = read(path, config, sample_ids=("S1",))
    row = getattr(result, destination)[0]
    assert row.producer == producer and row.origin == "GENERATED"
    assert row.path is None and row.source_version is None
    assert config == before and not (tmp_path / "ready.tsv").exists()


@pytest.mark.parametrize("producer", ["assembly_sv.pav", "assembly_sv.svim_asm"])
@pytest.mark.parametrize("destination", ["annotation", "cohort"])
@pytest.mark.parametrize("new_assembly,ready_manifest", [(False, None), (True, "ready.tsv")])
def test_generated_assembly_source_requires_exactly_one_input_mode(tmp_path, producer, destination, new_assembly, ready_manifest):
    config = configuration()
    config["assembly"]["enabled"] = new_assembly
    config["assembly_sv"]["input_manifest"] = ready_manifest
    if destination == "cohort":
        cohort_enabled(config, "sv")
    with pytest.raises(InputValidationError, match="assembly|input_manifest"):
        read(plan(tmp_path, **{destination: [generated(producer=producer)]}), config, sample_ids=("S1",))


@pytest.mark.parametrize("producer", ["assembly_sv.pav", "assembly_sv.svim_asm"])
@pytest.mark.parametrize("value", ["", " ", "bad\tpath", False, [], Path("ready.tsv")])
def test_generated_ready_source_rejects_malformed_manifest_intent(tmp_path, producer, value):
    config = configuration()
    config["assembly"]["enabled"] = False
    config["assembly_sv"]["input_manifest"] = value
    with pytest.raises(InputValidationError, match="input_manifest"):
        read(plan(tmp_path, annotation=[generated(producer=producer)]), config)


@pytest.mark.parametrize("producer", ["assembly_sv.pav", "assembly_sv.svim_asm"])
@pytest.mark.parametrize("disable_section", [True, False])
def test_ready_source_intent_does_not_enable_a_disabled_caller(tmp_path, producer, disable_section):
    config = configuration()
    config["assembly"]["enabled"] = False
    config["assembly_sv"]["input_manifest"] = "ready.tsv"
    if disable_section:
        config["assembly_sv"]["enabled"] = False
    else:
        config["assembly_sv"][producer.split(".")[1]]["enabled"] = False
    with pytest.raises(InputValidationError, match="disabled"):
        read(plan(tmp_path, annotation=[generated(producer=producer)]), config)


def test_external_paths_are_absolute_without_claiming_vcf_validation(tmp_path):
    row = external(tmp_path)
    source = plan(tmp_path, annotation=[row])
    result = read(source)
    assert result.annotation[0].path == (tmp_path / row["path"]).resolve()
    assert result.annotation[0].index_path == (tmp_path / row["index_path"]).resolve()
    assert result.annotation[0].producer is None
    assert result.annotation[0].source_tool == "sawfish"
    # These fixtures are not real VCF/TBI payloads: readiness is a later gate.
    assert result.annotation[0].path.read_bytes().startswith(b"Existing path fixture")


@pytest.mark.parametrize("updates", [
    {"path": "future.vcf.gz"}, {"index_path": "missing.tbi"},
    {"path": None}, {"index_path": None}, {"source_tool": ""},
    {"source_version": None}, {"reference_build": "GRCh37"},
    {"producer": "sv.sawfish"}, {"output_role": "vcf"},
])
def test_invalid_external_metadata_never_becomes_ready(tmp_path, updates):
    row = external(tmp_path, **updates)
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    source = plan(tmp_path, annotation=[row])
    with pytest.raises(InputValidationError):
        read(source)
    assert all((tmp_path / key).read_bytes() == value for key, value in before.items())


def test_external_build_is_required_not_inferred(tmp_path):
    row = external(tmp_path)
    row.pop("reference_build")
    with pytest.raises(InputValidationError, match="explicit reference_build"):
        read(plan(tmp_path, annotation=[row]))


def test_external_source_index_alias_is_rejected(tmp_path):
    row = external(tmp_path)
    row["index_path"] = row["path"]
    with pytest.raises(InputValidationError, match="distinct"):
        read(plan(tmp_path, annotation=[row]))


def test_external_hardlinked_index_is_rejected_without_mutation(tmp_path):
    row = external(tmp_path)
    alias = tmp_path / "hardlinked.tbi"
    alias.hardlink_to(tmp_path / row["path"])
    original = alias.read_bytes()
    row["index_path"] = alias.name
    with pytest.raises(InputValidationError, match="distinct"):
        read(plan(tmp_path, annotation=[row]))
    assert alias.read_bytes() == original


@pytest.mark.parametrize("field,value", [
    ("path", "future.vcf.gz"), ("index_path", "future.tbi"),
    ("source_tool", "guessed-tool"), ("source_version", "1.0"),
    ("state", "FAILED"),
])
def test_generated_intent_cannot_fabricate_receipts(tmp_path, field, value):
    with pytest.raises(InputValidationError):
        read(plan(tmp_path, annotation=[generated(**{field: value})]))


@pytest.mark.parametrize("state", ["NOT_RUN", "FAILED", "DISABLED", "MISSING_INPUT", "NOT_OBSERVED"])
def test_noncallable_cohort_states_preserved_without_fake_files(tmp_path, state):
    config = cohort_enabled(configuration(), "sv")
    row = dict(sample_id="S2", track="sv", origin="EXTERNAL", state=state)
    result = read(plan(tmp_path, cohort=[generated(), row]), config)
    assert result.cohort[1].state == state
    assert result.cohort[1].path is None and result.cohort[1].source_version is None


def test_no_calls_is_callable_and_requires_actual_external_artifacts(tmp_path):
    config = cohort_enabled(configuration(), "sv")
    row = external(tmp_path, "S2", state="NO_CALLS")
    result = read(plan(tmp_path, cohort=[generated(), row]), config)
    assert result.cohort[1].state == "NO_CALLS" and result.cohort[1].path is not None
    row.pop("path")
    with pytest.raises(InputValidationError):
        read(plan(tmp_path, cohort=[generated(), row]), config)


@pytest.mark.parametrize("state", ["PASS", "PENDING", "NO_VARIANTS", "failed", True, 1, None])
def test_unknown_or_wrong_typed_states_fail(tmp_path, state):
    config = cohort_enabled(configuration(), "sv")
    with pytest.raises(InputValidationError):
        read(plan(tmp_path, cohort=[generated(), dict(sample_id="S2", track="sv", origin="EXTERNAL", state=state)]), config)


@pytest.mark.parametrize("field,value", [("path", "future.vcf.gz"), ("producer", "sv.sawfish"),
                                         ("source_tool", "caller"), ("source_version", "1.0"),
                                         ("variant_ids", ["v1"])])
def test_noncallable_cohort_cannot_carry_future_paths_or_guessed_metadata(tmp_path, field, value):
    config = cohort_enabled(configuration(), "sv")
    row = dict(sample_id="S2", track="sv", origin="EXTERNAL", state="FAILED", **{field: value})
    with pytest.raises(InputValidationError, match="Noncallable"):
        read(plan(tmp_path, cohort=[generated(), row]), config)


@pytest.mark.parametrize("rows", [[], [generated()], [generated("S1"), generated("S1")],
                                  [generated("S1"), generated("S3")]])
def test_cohort_never_silently_omits_or_duplicates_samples(tmp_path, rows):
    config = cohort_enabled(configuration(), "sv")
    with pytest.raises(InputValidationError, match="cover|duplicate|unknown sample"):
        read(plan(tmp_path, cohort=rows), config)


def test_every_enabled_cohort_track_needs_all_samples(tmp_path):
    config = cohort_enabled(configuration(), "sv", "small_variants")
    rows = [generated(s) for s in ("S1", "S2")]
    with pytest.raises(InputValidationError, match="Every enabled"):
        read(plan(tmp_path, cohort=rows), config)


def test_duplicate_annotation_pair_rejected_even_for_distinct_producers(tmp_path):
    with pytest.raises(InputValidationError, match="duplicate"):
        read(plan(tmp_path, annotation=[generated(), generated(producer="sv.pbsv")]))


def test_generated_tr_catalog_is_content_identity_and_immutable(tmp_path):
    catalog = tmp_path / "catalog.bed"
    catalog.write_bytes(b"chr1\t1\t20\tID=L1;MOTIFS=AC;STRUC=(AC)n\n")
    config = cohort_enabled(configuration(), "tr")
    config["tr"].update(catalog=str(catalog), catalog_reference_build="GRCh38")
    before = catalog.read_bytes()
    digest = hashlib.sha256(before).hexdigest()
    rows = [generated(s, "tr", "tr.trgt") for s in ("S2", "S1")]
    rows[0]["catalog_id"] = digest.upper()
    result = read(plan(tmp_path, annotation=[generated(track="tr", producer="tr.trgt")], cohort=rows), config)
    assert all(row.catalog_id == "sha256:" + digest for row in (*result.annotation, *result.cohort))
    assert catalog.read_bytes() == before


def test_generated_tr_declared_hash_must_match_configured_catalog(tmp_path):
    catalog = tmp_path / "catalog.bed"
    catalog.write_bytes(b"reference catalog\n")
    config = configuration()
    config["tr"]["catalog"] = str(catalog)
    with pytest.raises(InputValidationError, match="differs"):
        read(plan(tmp_path, annotation=[generated(track="tr", producer="tr.trgt", catalog_id="a" * 64)]), config)


def test_catalog_relative_to_retained_user_config_not_different_plan_directory(tmp_path):
    resources = tmp_path / "resources"
    resources.mkdir()
    catalog = resources / "catalog.bed"
    catalog.write_bytes(b"reference catalog\n")
    user = tmp_path / "user.yaml"
    user.write_text(yaml.safe_dump({
        "annotation": {"enabled": True, "input_manifest": "ready.tsv", "vep_enabled": True,
                       "vep_cache_directory": "cache", "vep_cache_version": "116"},
        "tr": {"enabled": True, "catalog": "resources/catalog.bed"},
    }), encoding="utf-8")
    config = load_config(DEFAULT, user_config=user)
    plans = tmp_path / "plans"
    plans.mkdir()
    result = read(plan(plans, annotation=[generated(track="tr", producer="tr.trgt")]), config)
    assert result.annotation[0].catalog_id == "sha256:" + hashlib.sha256(catalog.read_bytes()).hexdigest()


def test_external_tr_can_use_declared_hash_without_new_tr_calling(tmp_path):
    config = cohort_enabled(configuration(), "tr")
    config["tr"]["enabled"] = False
    digest = "a" * 64
    rows = [external(tmp_path, s, "tr", catalog_id=("sha256:" if s == "S1" else "") + digest)
            for s in ("S1", "S2")]
    result = read(plan(tmp_path, cohort=rows), config)
    assert {row.catalog_id for row in result.cohort} == {"sha256:" + digest}


@pytest.mark.parametrize("catalog_id", [None, "CAT1", "catalog.bed", "sha256:1234", 123, "b" * 64])
def test_external_tr_catalog_identity_strict(tmp_path, catalog_id):
    config = cohort_enabled(configuration(), "tr")
    catalog = tmp_path / "catalog.bed"
    catalog.write_bytes(b"reference catalog\n")
    config["tr"]["catalog"] = str(catalog)
    actual = hashlib.sha256(catalog.read_bytes()).hexdigest()
    rows = [external(tmp_path, "S1", "tr", catalog_id=actual),
            external(tmp_path, "S2", "tr", catalog_id=catalog_id)]
    with pytest.raises(InputValidationError, match="catalog"):
        read(plan(tmp_path, cohort=rows), config)


def test_inconsistent_external_tr_catalogs_fail_without_config_catalog(tmp_path):
    config = cohort_enabled(configuration(), "tr")
    rows = [external(tmp_path, "S1", "tr", catalog_id="a" * 64),
            external(tmp_path, "S2", "tr", catalog_id="b" * 64)]
    with pytest.raises(InputValidationError, match="same explicit catalog"):
        read(plan(tmp_path, cohort=rows), config)


@pytest.mark.parametrize("updates", [{"catalog": None}, {"catalog": "missing.bed"},
                                    {"catalog_reference_build": "GRCh37"}, {"enabled": False}])
def test_generated_tr_requires_configured_catalog_and_build(tmp_path, updates):
    catalog = tmp_path / "catalog.bed"
    catalog.write_bytes(b"reference catalog\n")
    config = configuration()
    config["tr"].update(catalog=str(catalog), catalog_reference_build="GRCh38")
    config["tr"].update(updates)
    with pytest.raises(InputValidationError):
        read(plan(tmp_path, annotation=[generated(track="tr", producer="tr.trgt")]), config)


@pytest.mark.parametrize("document", [
    {"schema_version": True}, {"schema_version": "1"}, {"schema_version": 2},
    {"schema_version": 1, "annotation": None}, {"schema_version": 1, "annotation": {}},
    {"schema_version": 1, "cohort": "no"}, {"schema_version": 1, "unexpected": []},
    {"annotation": []}, [], "text",
])
def test_strict_top_level_schema(tmp_path, document):
    path = tmp_path / "invalid.yaml"
    path.write_text(yaml.safe_dump(document), encoding="utf-8")
    with pytest.raises(InputValidationError):
        read(path)


@pytest.mark.parametrize("updates", [
    {"unexpected": "value"}, {"sample_id": "S3"}, {"sample_id": "bad\tname"},
    {"sample_id": 1}, {"track": "small_variants"}, {"origin": "AUTO"},
    {"producer": "sv.sawfish\nnext"}, {"producer": " sv.sawfish"},
    {"variant_ids": "v1"}, {"variant_ids": ["v1", "v1"]},
    {"variant_ids": ["v1;v2"]}, {"variant_ids": ["v1\rnext"]},
    {"variant_ids": ["v1\u200b"]}, {"variant_ids": ["v1\u2028next"]},
    {"variant_ids": ["v1\u2029next"]}, {"catalog_id": "a" * 64},
])
def test_strict_row_fields_and_control_injection(tmp_path, updates):
    with pytest.raises(InputValidationError):
        read(plan(tmp_path, annotation=[generated(**updates)]))


@pytest.mark.parametrize("field", ["path", "index_path", "source_tool", "source_version", "reference_build"])
def test_external_control_injection_is_rejected(tmp_path, field):
    row = external(tmp_path)
    row[field] = "bad\nvalue"
    with pytest.raises(InputValidationError, match="control"):
        read(plan(tmp_path, annotation=[row]))


@pytest.mark.parametrize("content", [
    "schema_version: 1\nschema_version: 1\n",
    "schema_version: 1\nannotation: &rows []\ncohort: *rows\n",
    "schema_version: 1\n<<: {annotation: []}\n",
    "schema_version: 1\nannotation:\n- sample_id: S1\n  sample_id: S2\n",
    "schema_version: 1\nannotation: [\n",
    "schema_version: 1\n1: ignored\n",
])
def test_duplicate_yaml_keys_aliases_merge_and_malformed_yaml_fail(tmp_path, content):
    path = tmp_path / "invalid.yaml"
    path.write_text(content, encoding="utf-8")
    with pytest.raises(InputValidationError):
        read(path)


def test_bounded_plan_input_bytes_lines_and_rows(tmp_path):
    path = tmp_path / "large.yaml"
    path.write_bytes(b"#" + b"x" * (1024 * 1024))
    with pytest.raises(InputValidationError, match="1 MiB"):
        read(path)
    path.write_text("schema_version: 1\n" + "# comment\n" * 4096, encoding="utf-8")
    with pytest.raises(InputValidationError, match="4096-line"):
        read(path)
    path.write_text("schema_version: 1\nannotation: [" + ",".join(["{}"] * 4097) + "]\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="4096-row"):
        read(path)


@pytest.mark.parametrize("samples", [(), "S1", ("S1", "S1"), ("bad\nname",)])
def test_invalid_declared_sample_sets_fail(tmp_path, samples):
    with pytest.raises(InputValidationError):
        read(plan(tmp_path), sample_ids=samples)


@pytest.mark.parametrize("group", ["annotation", "cohort"])
def test_disabled_destinations_cannot_silently_discard_plan_rows(tmp_path, group):
    config = configuration()
    config[group]["enabled"] = False
    with pytest.raises(InputValidationError, match="disabled"):
        read(plan(tmp_path, **{group: [generated()]}), config)


def test_wrong_typed_config_enabled_flag_not_coerced(tmp_path):
    config = configuration()
    config["sv"]["sawfish"]["enabled"] = "yes"
    with pytest.raises(InputValidationError, match="booleans"):
        read(plan(tmp_path, annotation=[generated()]), config)


def test_unused_configured_catalog_does_not_probe_unrelated_file(tmp_path):
    config = configuration()
    config["tr"]["catalog"] = "not-needed-for-this-sv-plan.bed"
    assert read(plan(tmp_path, annotation=[generated()]), config).annotation[0].track == "sv"


def test_utf8_bom_unicode_paths_and_build_alias_are_explicit(tmp_path):
    source = tmp_path / "变异.vcf.gz"
    index = tmp_path / "变异.vcf.gz.tbi"
    source.write_bytes(b"file")
    index.write_bytes(b"index")
    row = dict(sample_id="S1", track="sv", origin="EXTERNAL", path=source.name,
               index_path=str(index.resolve()), source_tool="caller", source_version="1.0", reference_build="hg38")
    path = plan(tmp_path, annotation=[row])
    path.write_bytes(b"\xef\xbb\xbf" + path.read_bytes())
    result = read(path)
    assert result.annotation[0].path == source.resolve()
    assert result.annotation[0].reference_build == "GRCh38"
    assert result.sha256 == hashlib.sha256(path.read_bytes()).hexdigest()


def test_generated_harmonization_native_and_derived_are_explicit_even_without_read_calling(tmp_path):
    from hifivar.breakends import BndPolicy, POLICY, POLICY_SHA256
    config = configuration()
    config['sv']['enabled'] = False
    config['sv']['harmonization']['bnd_representation'] = BndPolicy(POLICY['policy_id'],1,POLICY_SHA256,POLICY['mode']).to_dict()
    for role in ('vcf','bnd_vcf'):
        result = read(plan(tmp_path, annotation=[generated(producer='sv.harmonization',role=role)]), config)
        assert result.annotation[0].output_role == role
        assert result.annotation[0].bnd_bundle_manifest is None


def test_generated_derived_role_requires_enabled_exact_policy(tmp_path):
    config = configuration()
    with pytest.raises(InputValidationError, match='fixed BND'):
        read(plan(tmp_path, annotation=[generated(producer='sv.harmonization',role='bnd_vcf')]),config)


@pytest.mark.parametrize('producer',['sv.sawfish','assembly_sv.pav','small.deepvariant'])
def test_derived_role_not_accepted_from_unrelated_producer(tmp_path,producer):
    with pytest.raises(InputValidationError):
        read(plan(tmp_path, annotation=[generated(producer=producer,role='bnd_vcf')]))


def test_derived_role_not_substituted_into_cohort(tmp_path):
    from hifivar.breakends import BndPolicy, POLICY, POLICY_SHA256
    config=cohort_enabled(configuration(),'sv')
    config['sv']['harmonization']['bnd_representation']=BndPolicy(POLICY['policy_id'],1,POLICY_SHA256,POLICY['mode']).to_dict()
    rows=[generated(s,producer='sv.harmonization',role='bnd_vcf') for s in ('S1','S2')]
    with pytest.raises(InputValidationError,match='not cohort'):
        read(plan(tmp_path,cohort=rows),config)


@pytest.mark.parametrize('group,track,origin',[('annotation','small','EXTERNAL'),('cohort','sv','EXTERNAL'),('annotation','sv','GENERATED')])
def test_bundle_path_is_not_fabricated_for_generated_or_other_categories(tmp_path,group,track,origin):
    manifest=tmp_path/'manifest.json';manifest.write_text('{}',encoding='utf-8')
    config=cohort_enabled(configuration(),'sv') if group=='cohort' else configuration()
    row=generated(track=track,bnd_bundle_manifest=str(manifest)) if origin=='GENERATED' else external(tmp_path,track=track,bnd_bundle_manifest=str(manifest))
    with pytest.raises(InputValidationError):
        read(plan(tmp_path,**{group:[row]}),config)


def test_external_derived_plan_retains_explicit_bundle_but_does_not_certify_contents(tmp_path):
    manifest=tmp_path/'manifest.json';manifest.write_text('{}',encoding='utf-8')
    row=external(tmp_path,bnd_bundle_manifest=manifest.name)
    result=read(plan(tmp_path,annotation=[row]))
    assert result.annotation[0].bnd_bundle_manifest==manifest.resolve()
    assert result.annotation[0].producer is None
