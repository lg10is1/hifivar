"""Opt-in AlphaGenome v1 scoring via a user's own API key.

No credentials are discovered automatically, written to disk, or forwarded to a
custom endpoint. This backend returns track-level scores, not clinical labels.
The CLI caller must verify REF alleles against its declared local GRCh38 FASTA.
"""

from __future__ import annotations

from collections import namedtuple
import copy
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib import import_module, metadata
import math
import os
import re
from typing import Any, Callable, Mapping

from hifivar.annotation import VariantCategory
from hifivar.exceptions import (
    CommandExecutionError, ConfigurationError, InputValidationError,
    OutputValidationError, ToolNotFoundError, ToolVersionError,
)
from hifivar.functional import PrioritizedVariant

SDK_VERSION = "0.10.0"
OFFICIAL_ENDPOINT = "dns:///gdmscience.googleapis.com:443"
DEFAULT_CREDENTIAL_ENV = "HIFIVAR_ALPHAGENOME_API_KEY"
SUPPORTED_SEQUENCE_LENGTHS = (16384, 131072, 524288, 1048576)
MODEL_VERSIONS = ("ALL_FOLDS", "FOLD_0", "FOLD_1", "FOLD_2", "FOLD_3")
SUPPORTED_MODALITIES = ("ATAC", "CONTACT_MAPS", "DNASE", "CHIP_TF", "CHIP_HISTONE", "CAGE", "PROCAP", "RNA_SEQ", "SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS")
_ALLOWED_CONTIGS = frozenset([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])
_TOKEN = re.compile(r"^[A-Z][A-Z0-9_]*$")
_ONTOLOGY = re.compile(r"^(?:UBERON|CL|CLO|EFO):[0-9]+$")
_DNA = re.compile(r"^[ACGT]+$")


def _bound_integer(value: object, label: str, lower: int, upper: int) -> int:
    if type(value) is not int or not lower <= value <= upper:
        raise InputValidationError(f"{label} must be an integer from {lower} to {upper}.")
    return value


@dataclass(frozen=True, slots=True)
class AlphaGenomeSettings:
    reference_build: str = "GRCh38"
    model_version: str = "ALL_FOLDS"
    requested_modalities: tuple[str, ...] = ("RNA_SEQ",)
    ontology_terms: tuple[str, ...] = ()
    sequence_length: int = 1048576
    max_variants: int = 20
    max_rows_per_variant: int = 100000
    max_allele_length: int = 10000
    credential_env: str = DEFAULT_CREDENTIAL_ENV
    connect_timeout_seconds: float = 15.0
    rpc_timeout_seconds: float = 60.0

    def __post_init__(self) -> None:
        if self.reference_build != "GRCh38":
            raise InputValidationError("AlphaGenome adapter accepts declared human GRCh38 only.")
        if self.model_version not in MODEL_VERSIONS:
            raise InputValidationError("AlphaGenome model_version must name an explicit supported checkpoint.")
        for name, bound in (("max_variants", 20), ("max_rows_per_variant", 100000), ("max_allele_length", 10000)):
            _bound_integer(getattr(self, name), name, 1, bound)
        if type(self.sequence_length) is not int or self.sequence_length not in SUPPORTED_SEQUENCE_LENGTHS:
            raise InputValidationError("AlphaGenome sequence_length must be one of the pinned SDK lengths.")
        if not isinstance(self.requested_modalities, (tuple, list)):
            raise InputValidationError("requested_modalities must be an explicit sequence.")
        modalities = tuple(self.requested_modalities)
        if any(not isinstance(v, str) for v in modalities):
            raise InputValidationError("AlphaGenome modality names must be strings.")
        if not modalities or len(modalities) > 20 or len(set(modalities)) != len(modalities):
            raise InputValidationError("Choose 1 to 20 distinct AlphaGenome scorer names.")
        if any(not isinstance(v, str) or v not in SUPPORTED_MODALITIES for v in modalities):
            raise InputValidationError("AlphaGenome modality/scorer names must be explicit uppercase identifiers.")
        if not isinstance(self.ontology_terms, (tuple, list)):
            raise InputValidationError("ontology_terms must be an explicit sequence.")
        terms = tuple(self.ontology_terms)
        if any(not isinstance(v, str) for v in terms):
            raise InputValidationError("Ontology terms must be strings.")
        if len(terms) > 20 or len(set(terms)) != len(terms) or any(not isinstance(v, str) or not _ONTOLOGY.fullmatch(v) for v in terms):
            raise InputValidationError("Ontology filters must be up to 20 distinct UBERON/CL/CLO/EFO CURIEs.")
        if not isinstance(self.credential_env, str) or not _TOKEN.fullmatch(self.credential_env):
            raise InputValidationError("credential_env must be an environment-variable name, never the API key.")
        for name, upper in (("connect_timeout_seconds", 60), ("rpc_timeout_seconds", 300)):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 1 <= value <= upper:
                raise InputValidationError(f"{name} must be finite and between 1 and {upper}.")
        object.__setattr__(self, "requested_modalities", modalities)
        object.__setattr__(self, "ontology_terms", terms)

    def to_dict(self) -> dict[str, Any]:
        return {name: list(value) if isinstance(value, tuple) else value
                for name in self.__dataclass_fields__
                for value in (getattr(self, name),)}


@dataclass(frozen=True, slots=True)
class AlphaGenomeRequest:
    selected_variants: tuple[PrioritizedVariant, ...]
    settings: AlphaGenomeSettings

    def __post_init__(self) -> None:
        if not isinstance(self.settings, AlphaGenomeSettings):
            raise InputValidationError("AlphaGenome request requires validated settings.")
        if not isinstance(self.selected_variants, (tuple, list)):
            raise InputValidationError("AlphaGenome selection must be an explicit sequence.")
        variants = tuple(self.selected_variants)
        if not variants or len(variants) > self.settings.max_variants:
            raise InputValidationError("AlphaGenome selection must be non-empty and within max_variants.")
        keys: set[tuple[str, str]] = set()
        for variant in variants:
            if not isinstance(variant, PrioritizedVariant):
                raise InputValidationError("AlphaGenome requires explicit PrioritizedVariant records.")
            key = variant.sample_id, variant.source_variant_id
            if key in keys:
                raise InputValidationError("AlphaGenome explicit selection contains duplicate source identities.")
            keys.add(key)
            if variant.contig not in _ALLOWED_CONTIGS:
                raise InputValidationError("AlphaGenome adapter supports chr1-chr22/chrX/chrY only; contigs are not renamed.")
            ref, alt = variant.reference_bases, variant.alternate_bases
            if not _DNA.fullmatch(ref) or not _DNA.fullmatch(alt) or ref == alt:
                raise InputValidationError("AlphaGenome requires explicit unequal A/C/G/T alleles; symbolic, BND, missing and multi-ALT records are unsupported.")
            if max(len(ref), len(alt)) > self.settings.max_allele_length:
                raise InputValidationError("Selected allele exceeds the configured engineering allele-length bound.")
            if variant.variant_category == VariantCategory.TR:
                raise InputValidationError("TR genotypes require explicit allele reconstruction and are not supported by this adapter.")
            if variant.variant_category == VariantCategory.SV:
                anchored_indel = ((len(ref) == 1 and alt.startswith(ref)) or
                                  (len(alt) == 1 and ref.startswith(alt)))
                if not anchored_indel or abs(len(ref) - len(alt)) < 50:
                    raise InputValidationError("SV selection requires an anchored sequence-resolved INS/DEL of at least 50bp; INV/BND/CPX are unsupported.")
            if max(len(ref), len(alt)) >= self.settings.sequence_length:
                raise InputValidationError("The selected allele must fit inside the model context.")
            if variant.position - 1 - (self.settings.sequence_length - len(ref)) // 2 < 0:
                raise InputValidationError("The model interval would start before chromosome zero; padded inference is not enabled.")
        object.__setattr__(self, "selected_variants", variants)

    def to_dict(self) -> dict[str, Any]:
        return {"settings": self.settings.to_dict(),
                "selected_variants": [v.to_dict() for v in self.selected_variants],
                "selection_is_explicit": True}


@dataclass(frozen=True, slots=True)
class AlphaGenomeResult:
    """JSON-safe response with no retained authentication value."""

    payload: Mapping[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return copy.deepcopy(dict(self.payload))


def _load_sdk() -> Any:
    failed = False
    sdk = None
    try:
        sdk = type("AlphaGenomeSDK", (), {
            "version": metadata.version("alphagenome"),
            "genome": import_module("alphagenome.data.genome"),
            "dna_client": import_module("alphagenome.models.v1.dna_client"),
            "variant_scorers": import_module("alphagenome.models.v1.variant_scorers"),
            "grpc": import_module("grpc"),
        })
    except Exception:
        failed = True
    if failed:
        raise ToolNotFoundError("Install the optional AlphaGenome backend dependency in this environment.") from None
    if sdk.version != SDK_VERSION:
        raise ToolVersionError(f"This adapter requires alphagenome=={SDK_VERSION}.") from None
    return sdk


def _deadline_interceptor(grpc: Any, seconds: float) -> Any:
    class CallDetails(namedtuple("CallDetails", "method timeout metadata credentials wait_for_ready compression"), grpc.ClientCallDetails):
        pass

    class Deadline(grpc.StreamStreamClientInterceptor):
        def intercept_stream_stream(self, continuation: Any, details: Any, requests: Any) -> Any:
            timeout = seconds if details.timeout is None else min(details.timeout, seconds)
            bounded = CallDetails(details.method, timeout, details.metadata,
                                  details.credentials, details.wait_for_ready, details.compression)
            return continuation(bounded, requests)

    return Deadline()


def _create_client(sdk: Any, settings: AlphaGenomeSettings, secret: str) -> tuple[Any, Callable[[], None]]:
    grpc = sdk.grpc
    channel = grpc.secure_channel(
        OFFICIAL_ENDPOINT, grpc.ssl_channel_credentials(),
        options=(("grpc.max_send_message_length", -1), ("grpc.max_receive_message_length", -1)))
    try:
        grpc.channel_ready_future(channel).result(timeout=settings.connect_timeout_seconds)
        bounded = grpc.intercept_channel(channel, _deadline_interceptor(grpc, settings.rpc_timeout_seconds))
        client = sdk.dna_client.DnaClient(
            channel=bounded, model_version=sdk.dna_client.ModelVersion[settings.model_version],
            metadata=(("x-goog-api-key", secret),))
        return client, channel.close
    except Exception:
        channel.close()
        raise


def _json_value(value: Any, secret: str, *, genome_module: Any = None) -> Any:
    # Only trusted, exact SDK data types receive a string representation.
    # Arbitrary provider objects never reach str/repr as a serialization fallback.
    interval_type = getattr(genome_module, "Interval", None)
    variant_type = getattr(genome_module, "Variant", None)
    if isinstance(interval_type, type) and type(value) is interval_type:
        return _json_value({
            "type": "genomic_interval", "chromosome": value.chromosome,
            "start": value.start, "end": value.end, "strand": value.strand,
            "coordinate_system": "0_based_half_open", "identifier": str(value),
        }, secret)
    if isinstance(variant_type, type) and type(value) is variant_type:
        return _json_value({
            "type": "genomic_variant", "chromosome": value.chromosome,
            "position": value.position, "reference_bases": value.reference_bases,
            "alternate_bases": value.alternate_bases,
            "coordinate_system": "1_based_vcf", "identifier": str(value),
        }, secret)
    if value is None or type(value) in (bool, int):
        return value
    if isinstance(value, str):
        return value.replace(secret, "[REDACTED]")
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if hasattr(value, "item"):
        return _json_value(value.item(), secret, genome_module=genome_module)
    if type(value) in (list, tuple):
        return [_json_value(x, secret, genome_module=genome_module) for x in value]
    if isinstance(value, Mapping):
        return {str(k).replace(secret, "[REDACTED]"): _json_value(v, secret, genome_module=genome_module) for k, v in value.items()}
    raise OutputValidationError("AlphaGenome score metadata contains an unsupported non-JSON value.")


class AlphaGenomeBackend:
    """Sequential, explicit-selection scorer; it neither writes outputs nor ranks calls."""

    def __init__(self, *, sdk_loader: Callable[[], Any] | None = None,
                 credential_reader: Callable[[str], str | None] | None = None,
                 client_factory: Callable[[Any, AlphaGenomeSettings, str], tuple[Any, Callable[[], None]]] | None = None) -> None:
        self._sdk_loader = sdk_loader or _load_sdk
        self._credential_reader = credential_reader or os.environ.get
        self._client_factory = client_factory or _create_client

    def score(self, request: AlphaGenomeRequest, *, allow_network: bool = False) -> AlphaGenomeResult:
        if not isinstance(request, AlphaGenomeRequest):
            raise InputValidationError("AlphaGenome score requires a validated request.")
        if allow_network is not True:
            raise ConfigurationError("AlphaGenome network execution requires explicit allow_network=True.")
        secret = None
        credential_failure = False
        try:
            secret = self._credential_reader(request.settings.credential_env)
        except Exception:
            credential_failure = True
        if credential_failure:
            raise ConfigurationError("Unable to read the explicitly configured credential variable.") from None
        if not isinstance(secret, str) or not 8 <= len(secret) <= 4096 or secret != secret.strip() or any(c in secret for c in "\r\n"):
            raise ConfigurationError("The explicitly configured AlphaGenome credential variable is missing or invalid.")
        failed = False
        output: list[Mapping[str, Any]] = []
        closer: Callable[[], None] | None = None
        started = datetime.now(timezone.utc).isoformat()
        try:
            sdk = self._sdk_loader()
            if sdk.version != SDK_VERSION:
                raise ToolVersionError(f"This adapter requires alphagenome=={SDK_VERSION}.")
            configs = sdk.variant_scorers.RECOMMENDED_VARIANT_SCORERS
            if any(name not in configs for name in request.settings.requested_modalities):
                raise InputValidationError("A requested scorer is unsupported by the pinned AlphaGenome SDK.")
            scorers = [configs[name] for name in request.settings.requested_modalities]
            client, closer = self._client_factory(sdk, request.settings, secret)
            for selected in request.selected_variants:
                variant = sdk.genome.Variant(
                    chromosome=selected.contig, position=selected.position,
                    reference_bases=selected.reference_bases, alternate_bases=selected.alternate_bases)
                interval = variant.reference_interval.resize(request.settings.sequence_length)
                if interval.start < 0 or interval.width != request.settings.sequence_length:
                    raise InputValidationError("AlphaGenome generated context does not meet the declared interval contract.")
                scores = client.score_variant(
                    interval=interval, variant=variant, variant_scorers=scorers,
                    organism=sdk.dna_client.Organism.HOMO_SAPIENS,
                    merge_stranded_gene_tracks=False)
                if not isinstance(scores, (list, tuple)) or len(scores) != len(scorers):
                    raise OutputValidationError("AlphaGenome scorer result cardinality is invalid.")
                frame = sdk.variant_scorers.tidy_scores(scores, match_gene_strand=False, include_extended_metadata=True)
                if frame is None or len(frame) > request.settings.max_rows_per_variant:
                    raise OutputValidationError("AlphaGenome returned no score table or exceeded the output-row bound.")
                rows = []
                for raw in frame.to_dict(orient="records"):
                    variant_type = getattr(sdk.genome, "Variant", None)
                    interval_type = getattr(sdk.genome, "Interval", None)
                    row_variant = raw.get("variant_id")
                    if isinstance(variant_type, type) and type(row_variant) is variant_type:
                        if any(getattr(row_variant, field) != getattr(variant, field)
                               for field in ("chromosome", "position", "reference_bases", "alternate_bases")):
                            raise OutputValidationError("AlphaGenome returned metadata for a different variant.")
                    row_interval = raw.get("scored_interval")
                    if isinstance(interval_type, type) and type(row_interval) is interval_type:
                        if any(getattr(row_interval, field) != getattr(interval, field)
                               for field in ("chromosome", "start", "end")):
                            raise OutputValidationError("AlphaGenome scored-interval metadata does not match the request.")
                    row = _json_value(raw, secret, genome_module=sdk.genome)
                    if not isinstance(row.get("output_type"), str) or row["output_type"] not in request.settings.requested_modalities:
                        raise OutputValidationError("AlphaGenome score table returned an unrequested or missing output type.")
                    score = row.get("raw_score")
                    if isinstance(score, bool) or not isinstance(score, (int, float)) or not math.isfinite(score):
                        raise OutputValidationError("AlphaGenome raw scores must be finite numeric values.")
                    if request.settings.ontology_terms and row.get("ontology_curie") not in request.settings.ontology_terms:
                        continue
                    rows.append(row)
                output.append({"variant": _json_value(selected.to_dict(), secret),
                               "status": "SCORED" if rows else "NO_MATCHING_SCORE_ROWS",
                               "returned_rows_before_ontology_filter": len(frame),
                               "score_rows": rows})
        except Exception:
            failed = True
        finally:
            if closer is not None:
                try:
                    closer()
                except Exception:
                    failed = True
        if failed:
            raise CommandExecutionError(
                "AlphaGenome request failed. Check account authorization, quota, network and pinned SDK compatibility; provider exception text was suppressed.") from None
        payload = {
            "schema_version": "1.0", "backend": "alphagenome_v1",
            "sdk_version": SDK_VERSION, "request": request.to_dict(),
            "variant_scores": output, "started_utc": started,
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "endpoint": OFFICIAL_ENDPOINT,
            "authentication": {"source": "EXPLICIT_ENVIRONMENT_VARIABLE",
                               "variable_name": request.settings.credential_env,
                               "credential_value_stored": False},
            "scientific_policy": {"functional_impact_is_call_confidence": False,
                                  "functional_impact_is_truth": False,
                                  "clinical_pathogenicity_classification": False,
                                  "large_sv_effect_accuracy_established": False},
            "runtime_policy": {"connect_timeout_seconds": request.settings.connect_timeout_seconds,
                               "per_rpc_timeout_seconds": request.settings.rpc_timeout_seconds,
                               "sdk_max_transient_attempts": 5, "additional_hifivar_retries": 0,
                               "whole_batch_hard_deadline": False},
            "scope": {"organism": "HOMO_SAPIENS", "model_reference": "GRCh38.p13.genome.fa",
                      "local_reference_match": "CALLER_MUST_VERIFY_REF_BEFORE_NETWORK",
                      "ontology_selection": "POST_SCORING_RESULT_FILTER_NOT_REMOTE_REQUEST_FILTER",
                      "genomic_variant_data_transmitted_to_provider": True}}
        return AlphaGenomeResult(_json_value(payload, secret))



__all__ = ["AlphaGenomeSettings", "AlphaGenomeRequest", "AlphaGenomeResult",
           "AlphaGenomeBackend", "SDK_VERSION", "OFFICIAL_ENDPOINT",
           "DEFAULT_CREDENTIAL_ENV", "SUPPORTED_SEQUENCE_LENGTHS", "MODEL_VERSIONS", "SUPPORTED_MODALITIES"]
