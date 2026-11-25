"""Sequence-resolved PCR/Sanger candidates, separate from experimental validation.

Coordinates entering this module are genomic 1-based closed. Primer3 offsets
are 0-based; its RIGHT coordinate is the rightmost (5-prime) template base.
Only an explicitly selected biallelic SNV or anchored INS/DEL is supported.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import math
import re
from typing import Mapping

from hifivar.exceptions import InputValidationError, OutputValidationError

_DNA = re.compile(r"[ACGT]+\Z")
_TEMPLATE_DNA = re.compile(r"[ACGTN]+\Z")


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise InputValidationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _number(value: object, name: str, minimum: float, maximum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or not minimum <= value <= maximum:
        raise InputValidationError(f"{name} must be a finite number in [{minimum}, {maximum}].")
    return float(value)


def reverse_complement(sequence: str) -> str:
    """Reverse-complement a validated DNA string (including N in templates)."""
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


@dataclass(frozen=True, slots=True)
class VariantTemplate:
    sample: str
    contig: str
    position: int
    ref: str
    alt: str
    sequence_start: int
    reference_sequence: str
    source_variant_id: str
    source_record_sha256: str | None = None

    def __post_init__(self) -> None:
        for name in ("sample", "contig", "source_variant_id"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip() or any(c in value for c in "\t\r\n\x00"):
                raise InputValidationError(f"{name} must be nonempty and contain no control separators.")
        _integer(self.position, "position", 1, 2**63 - 1)
        _integer(self.sequence_start, "sequence_start", 1, self.position)
        for name in ("ref", "alt"):
            value = getattr(self, name)
            if not isinstance(value, str) or len(value) > 100_000 or not _DNA.fullmatch(value):
                raise InputValidationError(f"{name} must be explicit uppercase A/C/G/T sequence; symbolic/complex SVs are unsupported.")
        if not isinstance(self.reference_sequence, str) or not 2 <= len(self.reference_sequence) <= 100_000 or not _TEMPLATE_DNA.fullmatch(self.reference_sequence):
            raise InputValidationError("reference_sequence must contain 2..100000 uppercase A/C/G/T/N bases.")
        if self.ref == self.alt:
            raise InputValidationError("REF and ALT must differ.")
        if not (len(self.ref) == len(self.alt) == 1 or
                (len(self.ref) == 1 < len(self.alt) and self.alt.startswith(self.ref)) or
                (len(self.alt) == 1 < len(self.ref) and self.ref.startswith(self.alt))):
            raise InputValidationError("Only biallelic SNV and simple anchored INS/DEL are supported; resolve/normalize complex alleles explicitly.")
        if self.offset + len(self.ref) >= len(self.reference_sequence) or self.offset == 0:
            raise InputValidationError("Local template must contain flanking sequence on both sides of REF.")
        if self.reference_sequence[self.offset:self.offset + len(self.ref)] != self.ref:
            raise InputValidationError("REF allele does not match the local reference template at its genomic coordinate.")
        if self.source_record_sha256 is not None and not re.fullmatch(r"[0-9a-f]{64}", self.source_record_sha256):
            raise InputValidationError("source_record_sha256 must be a lowercase SHA256 hex string.")

    @property
    def offset(self) -> int:
        return self.position - self.sequence_start

    @property
    def variant_type(self) -> str:
        if len(self.ref) == len(self.alt):
            return "SNV"
        return "INS" if len(self.alt) > len(self.ref) else "DEL"

    @property
    def alternate_sequence(self) -> str:
        return self.reference_sequence[:self.offset] + self.alt + self.reference_sequence[self.offset + len(self.ref):]

    def identity(self) -> dict[str, object]:
        return {"sample": self.sample, "contig": self.contig, "position": self.position,
                "ref": self.ref, "alt": self.alt, "source_variant_id": self.source_variant_id,
                "source_record_sha256": self.source_record_sha256, "sequence_start": self.sequence_start,
                "reference_template_sha256": hashlib.sha256(self.reference_sequence.encode("ascii")).hexdigest(),
                "reference_template_length": len(self.reference_sequence), "variant_type": self.variant_type}


@dataclass(frozen=True, slots=True)
class PrimerDesignSettings:
    requested_pairs: int = 3
    candidate_multiplier: int = 5
    product_size_ranges: tuple[tuple[int, int], ...] = ((200, 800),)
    primer_min_size: int = 18
    primer_opt_size: int = 20
    primer_max_size: int = 25
    primer_min_tm: float = 57.0
    primer_opt_tm: float = 60.0
    primer_max_tm: float = 63.0
    primer_min_gc: float = 20.0
    primer_max_gc: float = 80.0
    max_pair_tm_difference: float = 3.0
    min_pair_shift: int = 5
    timeout_seconds: float = 120.0
    thermodynamic_parameters_path: str | None = None

    def __post_init__(self) -> None:
        _integer(self.requested_pairs, "requested_pairs", 1, 20)
        _integer(self.candidate_multiplier, "candidate_multiplier", 1, 10)
        for name in ("primer_min_size", "primer_opt_size", "primer_max_size"):
            _integer(getattr(self, name), name, 15, 36)
        if not self.primer_min_size <= self.primer_opt_size <= self.primer_max_size:
            raise InputValidationError("Primer lengths must satisfy min <= opt <= max.")
        for name in ("primer_min_tm", "primer_opt_tm", "primer_max_tm"):
            _number(getattr(self, name), name, 0, 100)
        if not self.primer_min_tm <= self.primer_opt_tm <= self.primer_max_tm:
            raise InputValidationError("Primer Tm must satisfy min <= opt <= max.")
        _number(self.primer_min_gc, "primer_min_gc", 0, 100)
        _number(self.primer_max_gc, "primer_max_gc", self.primer_min_gc, 100)
        _number(self.max_pair_tm_difference, "max_pair_tm_difference", 0, 100)
        _integer(self.min_pair_shift, "min_pair_shift", 1, 1000)
        _number(self.timeout_seconds, "timeout_seconds", 0.01, 3600)
        if not isinstance(self.product_size_ranges, (tuple, list)) or not 1 <= len(self.product_size_ranges) <= 10:
            raise InputValidationError("product_size_ranges must contain 1..10 [min,max] pairs.")
        for pair in self.product_size_ranges:
            if not isinstance(pair, (tuple, list)) or len(pair) != 2:
                raise InputValidationError("Every product range must have two integers.")
            _integer(pair[0], "product minimum", 30, 100_000)
            _integer(pair[1], "product maximum", pair[0], 100_000)
        if self.thermodynamic_parameters_path is not None and (not isinstance(self.thermodynamic_parameters_path, str) or not self.thermodynamic_parameters_path.strip() or any(c in self.thermodynamic_parameters_path for c in "\r\n\x00")):
            raise InputValidationError("thermodynamic_parameters_path must be a nonempty, single-line path.")

    def to_dict(self) -> dict[str, object]:
        result = asdict(self)
        result['product_size_ranges'] = [list(pair) for pair in self.product_size_ranges]
        return result


def build_boulder_input(template: VariantTemplate, settings: PrimerDesignSettings) -> str:
    """Build one record, constraining both oligos to unchanged flanks."""
    if not isinstance(template, VariantTemplate) or not isinstance(settings, PrimerDesignSettings):
        raise InputValidationError("Design requires VariantTemplate and PrimerDesignSettings.")
    right_start = template.offset + len(template.ref)
    tags = {
        "SEQUENCE_ID": hashlib.sha256(repr(template.identity()).encode("utf-8")).hexdigest()[:24],
        "SEQUENCE_TEMPLATE": template.reference_sequence,
        "SEQUENCE_TARGET": f"{template.offset},{len(template.ref)}",
        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": f"0,{template.offset},{right_start},{len(template.reference_sequence) - right_start}",
        "PRIMER_TASK": "generic", "PRIMER_FIRST_BASE_INDEX": 0,
        "PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_PICK_INTERNAL_OLIGO": 0, "PRIMER_EXPLAIN_FLAG": 1,
        "PRIMER_NUM_RETURN": settings.requested_pairs * settings.candidate_multiplier,
        "PRIMER_PRODUCT_SIZE_RANGE": " ".join(f"{lo}-{hi}" for lo, hi in settings.product_size_ranges),
        "PRIMER_MIN_SIZE": settings.primer_min_size, "PRIMER_OPT_SIZE": settings.primer_opt_size,
        "PRIMER_MAX_SIZE": settings.primer_max_size,
        "PRIMER_MIN_TM": settings.primer_min_tm, "PRIMER_OPT_TM": settings.primer_opt_tm,
        "PRIMER_MAX_TM": settings.primer_max_tm,
        "PRIMER_MIN_GC": settings.primer_min_gc, "PRIMER_MAX_GC": settings.primer_max_gc,
        "PRIMER_PAIR_MAX_DIFF_TM": settings.max_pair_tm_difference,
        "PRIMER_MAX_NS_ACCEPTED": 0,
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1,
        "P3_FILE_FLAG": 0,
    }
    if settings.thermodynamic_parameters_path is not None:
        tags["PRIMER_THERMODYNAMIC_PARAMETERS_PATH"] = settings.thermodynamic_parameters_path
    return "\n".join(f"{key}={value}" for key, value in tags.items()) + "\n=\n"


def parse_boulder_record(text: str) -> dict[str, str]:
    """Require one complete response and reject duplicate tags or extra records."""
    if not isinstance(text, str) or len(text) > 10_000_000:
        raise OutputValidationError("Primer3 response is missing or exceeds its bounded size.")
    result: dict[str, str] = {}
    ended = False
    for line in text.splitlines():
        if not line:
            continue
        if line == "=":
            if ended:
                raise OutputValidationError("Expected exactly one Primer3 output record.")
            ended = True
            continue
        if ended or "=" not in line:
            raise OutputValidationError("Malformed or multiple Primer3 output records.")
        key, value = line.split("=", 1)
        if not re.fullmatch(r"[A-Z0-9_]+", key) or key in result:
            raise OutputValidationError("Malformed or duplicate Primer3 output tag.")
        result[key] = value
    if not ended:
        raise OutputValidationError("Primer3 output record has no terminating '='.")
    if result.get("PRIMER_ERROR"):
        raise OutputValidationError(f"Primer3 rejected the design: {result['PRIMER_ERROR']}")
    return result


def _output_int(value: str | None, name: str, low: int, high: int) -> int:
    try:
        number = int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError) as error:
        raise OutputValidationError(f"Primer3 {name} is not an integer.") from error
    if not low <= number <= high:
        raise OutputValidationError(f"Primer3 {name} is outside allowed bounds.")
    return number


def _output_float(value: str | None, name: str) -> float:
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError) as error:
        raise OutputValidationError(f"Primer3 {name} is not numeric.") from error
    if not math.isfinite(number):
        raise OutputValidationError(f"Primer3 {name} must be finite.")
    return number


def _position(tags: Mapping[str, str], key: str, length: int) -> tuple[int, int]:
    try:
        raw = tags[key].split(",")
    except KeyError as error:
        raise OutputValidationError(f"Primer3 output is missing {key}.") from error
    if len(raw) != 2:
        raise OutputValidationError(f"Primer3 {key} must be a start,length pair.")
    return _output_int(raw[0], key, 0, length - 1), _output_int(raw[1], key, 1, length)


def evaluate_candidates(template: VariantTemplate, settings: PrimerDesignSettings,
                        tags: Mapping[str, str]) -> dict[str, object]:
    """Validate Primer3 pairs and reconstruct products; never assert specificity."""
    count = _output_int(tags.get("PRIMER_PAIR_NUM_RETURNED"), "PRIMER_PAIR_NUM_RETURNED", 0,
                        settings.requested_pairs * settings.candidate_multiplier)
    accepted: list[dict[str, object]] = []
    rejected: list[dict[str, object]] = []
    seen: list[tuple[int, int]] = []
    for index in range(count):
        prefix = f"PRIMER_PAIR_{index}"
        left_start, left_length = _position(tags, f"PRIMER_LEFT_{index}", len(template.reference_sequence))
        right_end, right_length = _position(tags, f"PRIMER_RIGHT_{index}", len(template.reference_sequence))
        right_start = right_end - right_length + 1
        if (left_start + left_length > template.offset or right_start < template.offset + len(template.ref)
                or right_end >= len(template.reference_sequence) or right_start < 0):
            raise OutputValidationError("Primer3 returned primers outside unchanged flanks or template bounds.")
        left = template.reference_sequence[left_start:left_start + left_length]
        right = reverse_complement(template.reference_sequence[right_start:right_end + 1])
        if not _DNA.fullmatch(left) or not _DNA.fullmatch(right):
            raise OutputValidationError("Primer3 returned ambiguous primer binding sequence.")
        if tags.get(f"PRIMER_LEFT_{index}_SEQUENCE") != left or tags.get(f"PRIMER_RIGHT_{index}_SEQUENCE") != right:
            raise OutputValidationError("Primer3 primer sequence/orientation disagrees with template coordinates.")
        ref_product = template.reference_sequence[left_start:right_end + 1]
        variant_offset = template.offset - left_start
        alt_product = ref_product[:variant_offset] + template.alt + ref_product[variant_offset + len(template.ref):]
        if _output_int(tags.get(prefix + "_PRODUCT_SIZE"), prefix + "_PRODUCT_SIZE", 1, 100_000) != len(ref_product):
            raise OutputValidationError("Primer3 product length disagrees with reconstructed product.")
        reason = None
        if not any(lo <= len(ref_product) <= hi for lo, hi in settings.product_size_ranges):
            reason = "REFERENCE_PRODUCT_OUTSIDE_REQUESTED_RANGE"
        elif not any(lo <= len(alt_product) <= hi for lo, hi in settings.product_size_ranges):
            reason = "ALTERNATE_PRODUCT_OUTSIDE_REQUESTED_RANGE"
        elif any(abs(left_start - old_left) < settings.min_pair_shift and abs(right_end - old_right) < settings.min_pair_shift for old_left, old_right in seen):
            reason = "INSUFFICIENT_POSITIONAL_DIVERSITY"
        left_tm = _output_float(tags.get(f"PRIMER_LEFT_{index}_TM"), "left Tm")
        right_tm = _output_float(tags.get(f"PRIMER_RIGHT_{index}_TM"), "right Tm")
        penalty = _output_float(tags.get(prefix + "_PENALTY"), "pair penalty")
        if not (settings.primer_min_size <= left_length <= settings.primer_max_size and settings.primer_min_size <= right_length <= settings.primer_max_size):
            reason = "PRIMER_LENGTH_OUTSIDE_REQUESTED_RANGE"
        if not (settings.primer_min_tm <= left_tm <= settings.primer_max_tm and settings.primer_min_tm <= right_tm <= settings.primer_max_tm and abs(left_tm - right_tm) <= settings.max_pair_tm_difference):
            reason = "TM_OUTSIDE_REQUESTED_RANGE"
        gc_left = 100 * (left.count("G") + left.count("C")) / len(left)
        gc_right = 100 * (right.count("G") + right.count("C")) / len(right)
        if not all(settings.primer_min_gc <= gc <= settings.primer_max_gc for gc in (gc_left, gc_right)):
            reason = "GC_OUTSIDE_REQUESTED_RANGE"
        if penalty < 0:
            raise OutputValidationError("Primer3 pair penalty must be nonnegative.")
        if reason:
            rejected.append({"primer3_index": index, "reason": reason})
            continue
        seen.append((left_start, right_end))
        if len(accepted) >= settings.requested_pairs:
            continue
        structure = {}
        for tag in (f"PRIMER_LEFT_{index}_SELF_ANY_TH", f"PRIMER_LEFT_{index}_SELF_END_TH",
                    f"PRIMER_LEFT_{index}_HAIRPIN_TH", f"PRIMER_RIGHT_{index}_SELF_ANY_TH",
                    f"PRIMER_RIGHT_{index}_SELF_END_TH", f"PRIMER_RIGHT_{index}_HAIRPIN_TH",
                    prefix + "_COMPL_ANY_TH", prefix + "_COMPL_END_TH"):
            if tag in tags:
                structure[tag] = _output_float(tags[tag], tag)
        accepted.append({
            "rank": len(accepted) + 1, "primer3_index": index, "primer3_pair_penalty": penalty,
            "forward_sequence_5to3": left, "reverse_sequence_5to3": right,
            "forward_start_1based": template.sequence_start + left_start,
            "forward_end_1based": template.sequence_start + left_start + left_length - 1,
            "reverse_start_1based": template.sequence_start + right_start,
            "reverse_end_1based": template.sequence_start + right_end,
            "reverse_binding_strand": "-", "forward_binding_strand": "+",
            "forward_tm_celsius": left_tm, "reverse_tm_celsius": right_tm,
            "forward_gc_percent": gc_left, "reverse_gc_percent": gc_right,
            "reference_product_bp": len(ref_product), "alternate_product_bp": len(alt_product),
            "reference_product_sequence": ref_product, "alternate_product_sequence": alt_product,
            "product_length_delta_bp": len(alt_product) - len(ref_product),
            "thermodynamic_metrics_reported_by_primer3": structure,
            "specificity_status": "NOT_ASSESSED", "sample_nearby_variants_status": "NOT_ASSESSED",
            "experimental_status": "NOT_RUN",
        })
    return {
        "status": "CANDIDATES_GENERATED" if accepted else "NO_CANDIDATES",
        "variant": template.identity(), "settings": settings.to_dict(),
        "assay": "PCR_SANGER" if template.variant_type == "SNV" else "FLANKING_PCR",
        "snv_gel_discriminates_alleles": False if template.variant_type == "SNV" else None,
        "gel_resolution_status": "NOT_APPLICABLE_EQUAL_LENGTH_SNV" if template.variant_type == "SNV" else "NOT_ASSESSED",
        "requested_pair_count": settings.requested_pairs, "primer3_returned_pair_count": count,
        "candidate_pair_count": len(accepted), "requested_count_satisfied": len(accepted) >= settings.requested_pairs,
        "candidates": accepted, "rejected_candidates": rejected,
        "primer3_warning": tags.get("PRIMER_WARNING"),
        "primer3_explanations": {key: value for key, value in tags.items() if key.endswith("_EXPLAIN")},
        "specificity_status": "NOT_ASSESSED", "experimental_status": "NOT_RUN",
        "limitations": ["Candidates require genomic specificity and nearby sample-variant assessment.",
                        "Predicted product lengths do not establish gel resolution, amplification or genotype.",
                        "No allele dropout, primer dimer or sample haplotype absence is experimentally excluded."],
    }


__all__ = ["VariantTemplate", "PrimerDesignSettings", "build_boulder_input", "parse_boulder_record",
           "evaluate_candidates", "reverse_complement"]
