"""No-network AlphaGenome adapter regressions."""
from __future__ import annotations

from dataclasses import replace
import json
from types import SimpleNamespace

import pytest

from hifivar.alphagenome_backend import (
    AlphaGenomeBackend, AlphaGenomeRequest, AlphaGenomeSettings,
    DEFAULT_CREDENTIAL_ENV, OFFICIAL_ENDPOINT, SDK_VERSION,
    _create_client, _deadline_interceptor,
)
from hifivar.annotation import VariantCategory
from hifivar.exceptions import CommandExecutionError, ConfigurationError, InputValidationError
from hifivar.functional import PrioritizedVariant


@pytest.fixture
def selected(tmp_path):
    annotation = tmp_path / "selected.annotation.tsv"
    annotation.write_text("source\n", encoding="utf-8")
    return PrioritizedVariant("S1", "v1", VariantCategory.SMALL, "chr1", 1000000,
                              "A", "C", annotation, "explicit_research_selection")


class Frame:
    def __init__(self, rows):
        self.rows = rows

    def __len__(self):
        return len(self.rows)

    def to_dict(self, *, orient):
        assert orient == "records"
        return self.rows


class Interval:
    def __init__(self, start, width=1, chromosome="chr1", strand="."):
        self.start, self.width = start, width
        self.end, self.chromosome, self.strand = start + width, chromosome, strand

    def resize(self, width):
        return Interval(self.start - (width - self.width) // 2, width, self.chromosome, self.strand)

    def __str__(self):
        return f"{self.chromosome}:{self.start}-{self.end}"


class Variant:
    def __init__(self, *, chromosome, position, reference_bases, alternate_bases):
        self.chromosome, self.position = chromosome, position
        self.reference_bases, self.alternate_bases = reference_bases, alternate_bases
        self.reference_interval = Interval(position - 1, len(reference_bases), chromosome)

    def __str__(self):
        return f"{self.chromosome}:{self.position}:{self.reference_bases}>{self.alternate_bases}"


def apparatus(rows=None, *, error=None, version=SDK_VERSION):
    state = {"calls": [], "closed": 0, "credential_names": []}
    if rows is None:
        rows = [
            {"output_type": "RNA_SEQ", "raw_score": 0.2, "quantile_score": 0.7,
             "track_name": "brain_plus", "track_strand": "+", "ontology_curie": "UBERON:0000955",
             "gene_id": "ENSG001", "variant_scorer": "GeneMaskLFCScorer"},
            {"output_type": "RNA_SEQ", "raw_score": -0.3, "track_name": "brain_minus",
             "track_strand": "-", "ontology_curie": "UBERON:0000955",
             "gene_id": "ENSG002", "variant_scorer": "GeneMaskLFCScorer"},
        ]

    def score(**kwargs):
        state["calls"].append(kwargs)
        if error is not None:
            raise error
        return [object() for _ in kwargs["variant_scorers"]]

    def tidy(scores, **kwargs):
        assert kwargs == {"match_gene_strand": False, "include_extended_metadata": True}
        return Frame(rows)

    sdk = SimpleNamespace(
        version=version, genome=SimpleNamespace(Variant=Variant, Interval=Interval),
        dna_client=SimpleNamespace(Organism=SimpleNamespace(HOMO_SAPIENS="human")),
        variant_scorers=SimpleNamespace(RECOMMENDED_VARIANT_SCORERS={"RNA_SEQ": "rna", "DNASE": "dnase"}, tidy_scores=tidy),
    )

    def factory(binding, settings, credential):
        state["passed_secret"] = credential
        return SimpleNamespace(score_variant=score), lambda: state.update(closed=state["closed"] + 1)

    def credential(name):
        state["credential_names"].append(name)
        return "synthetic-test-key"

    backend = AlphaGenomeBackend(sdk_loader=lambda: sdk, credential_reader=credential, client_factory=factory)
    return backend, state


@pytest.mark.parametrize("kwargs", [
    {"requested_modalities": None}, {"ontology_terms": None}, {"requested_modalities": ({},)}, {"ontology_terms": ({},)},
    {"reference_build": "GRCh37"}, {"model_version": "latest"}, {"sequence_length": 2048},
    {"max_variants": 21}, {"max_variants": True}, {"max_rows_per_variant": 100001},
    {"max_allele_length": 10001}, {"requested_modalities": ()},
    {"requested_modalities": ("RNA_SEQ", "RNA_SEQ")}, {"requested_modalities": ("made_up",)},
    {"ontology_terms": ("brain",)}, {"credential_env": "literal-key-value"},
    {"rpc_timeout_seconds": 0}, {"rpc_timeout_seconds": float("inf")},
    {"connect_timeout_seconds": True},
])
def test_settings_reject_invalid_contracts(kwargs):
    with pytest.raises(InputValidationError):
        AlphaGenomeSettings(**kwargs)


@pytest.mark.parametrize("changes", [
    {"alternate_bases": "<DEL>"}, {"alternate_bases": "A,C"},
    {"alternate_bases": "A]chr2:300]"}, {"alternate_bases": "*"},
    {"alternate_bases": "."}, {"alternate_bases": "N"}, {"alternate_bases": "A"},
    {"contig": "1"}, {"contig": "chrM"}, {"position": 1},
    {"variant_category": VariantCategory.TR},
    {"variant_category": VariantCategory.SV, "reference_bases": "AC", "alternate_bases": "GT"},
])
def test_request_rejects_unsupported_before_network(selected, changes):
    with pytest.raises(InputValidationError):
        AlphaGenomeRequest((replace(selected, **changes),), AlphaGenomeSettings())


def test_request_duplicate_and_limit(selected):
    with pytest.raises(InputValidationError):
        AlphaGenomeRequest((selected, selected), AlphaGenomeSettings())
    with pytest.raises(InputValidationError):
        AlphaGenomeRequest((selected, replace(selected, source_variant_id="v2")), AlphaGenomeSettings(max_variants=1))


def test_sequence_resolved_sv_insertion_allowed(selected):
    request = AlphaGenomeRequest((replace(selected, variant_category=VariantCategory.SV, alternate_bases="A" + "C" * 100),), AlphaGenomeSettings())
    assert request.selected_variants[0].variant_category is VariantCategory.SV


def test_network_flag_prevents_sdk_and_credential_reads(selected):
    def forbidden(*args):
        raise AssertionError("must not be called")
    backend = AlphaGenomeBackend(sdk_loader=forbidden, credential_reader=forbidden)
    with pytest.raises(ConfigurationError, match="explicit"):
        backend.score(AlphaGenomeRequest((selected,), AlphaGenomeSettings()))


def test_legacy_variable_not_discovered(selected):
    # Do not mutate/read the host environment, even to back up a prior value.
    synthetic_environment = {"ALPHAGENOME_API_KEY": "synthetic-legacy-key"}
    looked_up = []
    def credential_reader(name):
        looked_up.append(name)
        return synthetic_environment.get(name)
    backend = AlphaGenomeBackend(sdk_loader=lambda: pytest.fail("no SDK load"),
                                 credential_reader=credential_reader)
    with pytest.raises(ConfigurationError, match="missing"):
        backend.score(AlphaGenomeRequest((selected,), AlphaGenomeSettings()), allow_network=True)
    assert looked_up == [DEFAULT_CREDENTIAL_ENV]


def test_preserves_tracks_tissues_scores_and_provenance(selected):
    backend, state = apparatus()
    request = AlphaGenomeRequest((selected,), AlphaGenomeSettings(model_version="FOLD_1"))
    value = backend.score(request, allow_network=True).to_dict()
    rows = value["variant_scores"][0]["score_rows"]
    assert len(rows) == 2
    assert [r["gene_id"] for r in rows] == ["ENSG001", "ENSG002"]
    assert [r["track_strand"] for r in rows] == ["+", "-"]
    assert [r["raw_score"] for r in rows] == [0.2, -0.3]
    assert rows[0]["quantile_score"] == 0.7
    assert value["scope"]["genomic_variant_data_transmitted_to_provider"] is True
    assert value["request"]["settings"]["model_version"] == "FOLD_1"
    assert value["scientific_policy"]["functional_impact_is_call_confidence"] is False
    assert value["runtime_policy"]["whole_batch_hard_deadline"] is False
    assert state["calls"][0]["merge_stranded_gene_tracks"] is False
    assert state["credential_names"] == [DEFAULT_CREDENTIAL_ENV]
    assert state["closed"] == 1
    assert "synthetic-test-key" not in json.dumps(value)
    rows[0]["raw_score"] = 999
    assert backend.score(request, allow_network=True).to_dict()["variant_scores"][0]["score_rows"][0]["raw_score"] == 0.2


def test_ontology_filter_is_explicit_after_scoring(selected):
    backend, state = apparatus()
    request = AlphaGenomeRequest((selected,), AlphaGenomeSettings(ontology_terms=("UBERON:99999",)))
    value = backend.score(request, allow_network=True).to_dict()
    assert value["variant_scores"][0]["status"] == "NO_MATCHING_SCORE_ROWS"
    assert value["variant_scores"][0]["returned_rows_before_ontology_filter"] == 2
    assert value["variant_scores"][0]["score_rows"] == []
    assert "ontology_terms" not in state["calls"][0]


@pytest.mark.parametrize("rows", [
    [{"output_type": "RNA_SEQ", "raw_score": float("nan")}],
    [{"output_type": "RNA_SEQ", "raw_score": True}],
    [{"output_type": "RNA_SEQ", "raw_score": "0.2"}],
    [{"output_type": "DNASE", "raw_score": 1.0}],
])
def test_bad_provider_scores_fail_closed(selected, rows):
    backend, state = apparatus(rows)
    with pytest.raises(CommandExecutionError):
        backend.score(AlphaGenomeRequest((selected,), AlphaGenomeSettings()), allow_network=True)
    assert state["closed"] == 1


def test_row_bound_before_materialization(selected):
    backend, state = apparatus()
    with pytest.raises(CommandExecutionError):
        backend.score(AlphaGenomeRequest((selected,), AlphaGenomeSettings(max_rows_per_variant=1)), allow_network=True)
    assert state["closed"] == 1


def test_provider_error_and_chain_do_not_expose_key(selected):
    backend, state = apparatus(error=RuntimeError("synthetic-test-key provider details"))
    with pytest.raises(CommandExecutionError) as caught:
        backend.score(AlphaGenomeRequest((selected,), AlphaGenomeSettings()), allow_network=True)
    assert "synthetic-test-key" not in str(caught.value)
    assert caught.value.__cause__ is None
    assert caught.value.__context__ is None
    assert state["closed"] == 1


def test_initialization_failure_sanitized(selected):
    backend, _ = apparatus()
    def bad_factory(*args):
        raise RuntimeError("synthetic-test-key connection metadata")
    backend._client_factory = bad_factory
    with pytest.raises(CommandExecutionError) as caught:
        backend.score(AlphaGenomeRequest((selected,), AlphaGenomeSettings()), allow_network=True)
    assert "synthetic-test-key" not in str(caught.value)
    assert caught.value.__context__ is None


def test_provider_metadata_redacted(selected):
    backend, _ = apparatus([{"output_type": "RNA_SEQ", "raw_score": 1.0,
                            "track_name": "synthetic-test-key"}])
    out = backend.score(AlphaGenomeRequest((selected,), AlphaGenomeSettings()), allow_network=True).to_dict()
    assert out["variant_scores"][0]["score_rows"][0]["track_name"] == "[REDACTED]"
    assert "synthetic-test-key" not in json.dumps(out)


def test_sdk_version_guard(selected):
    backend, state = apparatus(version="0.9.0")
    with pytest.raises(CommandExecutionError):
        backend.score(AlphaGenomeRequest((selected,), AlphaGenomeSettings()), allow_network=True)
    assert not state["calls"]


class DetailsBase:
    pass


class InterceptorBase:
    pass


def fake_grpc():
    state = {}
    channel = SimpleNamespace(close=lambda: state.update(closed=True))
    def secure_channel(address, credentials, *, options):
        state.update(address=address, credentials=credentials, options=options)
        return channel
    def ready(ch):
        assert ch is channel
        return SimpleNamespace(result=lambda *, timeout: state.update(connect_timeout=timeout))
    def intercept(ch, interceptor):
        state["interceptor"] = interceptor
        return "intercepted-channel"
    return SimpleNamespace(ClientCallDetails=DetailsBase,
                           StreamStreamClientInterceptor=InterceptorBase,
                           secure_channel=secure_channel, ssl_channel_credentials=lambda: "TLS",
                           channel_ready_future=ready, intercept_channel=intercept), state


@pytest.mark.parametrize("existing,expected", [(None, 7), (2, 2), (20, 7)])
def test_deadline_interceptor_public_api_shape(existing, expected):
    grpc, _ = fake_grpc()
    interceptor = _deadline_interceptor(grpc, 7)
    detail = SimpleNamespace(method="/ScoreVariant", timeout=existing,
                             metadata=(("field", "value"),), credentials=None,
                             wait_for_ready=False, compression=None)
    seen = {}
    result = interceptor.intercept_stream_stream(
        lambda d, it: seen.update(details=d, iterator=it) or "result", detail, iter([1]))
    assert result == "result"
    assert seen["details"].timeout == expected
    assert seen["details"].metadata == detail.metadata
    assert list(seen["iterator"]) == [1]


def test_secure_client_public_constructor_contract():
    grpc, state = fake_grpc()
    def client(**kwargs):
        state["client_args"] = kwargs
        return "client"
    sdk = SimpleNamespace(grpc=grpc, dna_client=SimpleNamespace(
        DnaClient=client, ModelVersion={"ALL_FOLDS": "selected-enum"}))
    actual, close = _create_client(sdk, AlphaGenomeSettings(), "synthetic-test-key")
    assert actual == "client"
    assert state["address"] == OFFICIAL_ENDPOINT
    assert state["credentials"] == "TLS"
    assert state["connect_timeout"] == 15
    assert state["client_args"]["channel"] == "intercepted-channel"
    assert state["client_args"]["model_version"] == "selected-enum"
    assert state["client_args"]["metadata"] == (("x-goog-api-key", "synthetic-test-key"),)
    close()
    assert state["closed"] is True


def test_deadline_interceptor_real_grpc_type_when_available():
    grpc = pytest.importorskip("grpc")
    assert isinstance(_deadline_interceptor(grpc, 7), grpc.StreamStreamClientInterceptor)


@pytest.mark.parametrize("ref,alt", [("A", "CTTT"), ("AC", "G"), ("AC", "GT" + "A" * 60), ("A", "AC")])
def test_sv_requires_pure_anchored_50bp_indel(selected, ref, alt):
    with pytest.raises(InputValidationError):
        AlphaGenomeRequest((replace(selected, variant_category=VariantCategory.SV,
                                    reference_bases=ref, alternate_bases=alt),),
                           AlphaGenomeSettings())


def test_credential_reader_failure_suppressed(selected):
    def failure(name):
        raise RuntimeError("synthetic-test-key")
    backend = AlphaGenomeBackend(credential_reader=failure)
    with pytest.raises(ConfigurationError) as caught:
        backend.score(AlphaGenomeRequest((selected,), AlphaGenomeSettings()), allow_network=True)
    assert "synthetic-test-key" not in str(caught.value)
    assert caught.value.__context__ is None


def test_settings_serialization_contains_only_credential_reference():
    payload = AlphaGenomeSettings().to_dict()
    encoded = json.dumps(payload, allow_nan=False)
    assert payload["credential_env"] == DEFAULT_CREDENTIAL_ENV
    assert "synthetic-test-key" not in encoded
    assert not {"api_key", "token", "authorization", "credential_value"}.intersection(payload)
    assert payload["model_version"] == "ALL_FOLDS"
    assert payload["requested_modalities"] == ["RNA_SEQ"]

def test_realistic_sdk_genome_objects_preserved_as_structured_json(selected):
    variant = Variant(chromosome="chr1", position=1000000, reference_bases="A", alternate_bases="C")
    interval = variant.reference_interval.resize(1048576)
    backend, _ = apparatus([{"output_type": "RNA_SEQ", "raw_score": 0.2,
                            "variant_id": variant, "scored_interval": interval,
                            "track_name": "brain"}])
    out = backend.score(AlphaGenomeRequest((selected,), AlphaGenomeSettings()), allow_network=True).to_dict()
    row = out["variant_scores"][0]["score_rows"][0]
    assert row["variant_id"]["position"] == 1000000
    assert row["variant_id"]["reference_bases"] == "A"
    assert row["variant_id"]["coordinate_system"] == "1_based_vcf"
    assert row["variant_id"]["identifier"] == str(variant)
    assert row["scored_interval"]["start"] == interval.start
    assert row["scored_interval"]["end"] == interval.end
    assert row["scored_interval"]["coordinate_system"] == "0_based_half_open"
    json.dumps(out, allow_nan=False)


@pytest.mark.parametrize("column,value", [
    ("variant_id", Variant(chromosome="chr1", position=1000001, reference_bases="A", alternate_bases="C")),
    ("scored_interval", Interval(1, 1048576, "chr2")),
    ("unknown", object()),
])
def test_mismatched_or_unknown_metadata_objects_fail_closed(selected, column, value):
    backend, _ = apparatus([{"output_type": "RNA_SEQ", "raw_score": 0.2, column: value}])
    with pytest.raises(CommandExecutionError):
        backend.score(AlphaGenomeRequest((selected,), AlphaGenomeSettings()), allow_network=True)
