"""Reference-genome metadata shared by future HiFiVar analysis modules."""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass, replace
from pathlib import Path
from types import MappingProxyType
from typing import Final, Mapping

from hifivar import validation
from hifivar.exceptions import InputValidationError, ReferenceError
from hifivar.logging_utils import get_logger


PathInput = str | Path

REFERENCE_BUILD_ALIASES: Final[Mapping[str, str]] = MappingProxyType(
    {
        "grch37": "GRCh37",
        "hg19": "GRCh37",
        "grch38": "GRCh38",
        "hg38": "GRCh38",
        "t2t-chm13": "T2T-CHM13",
        "custom": "custom",
        "unknown": "unknown",
    }
)

_COMPRESSED_FASTA_SUFFIXES = (".fa.gz", ".fasta.gz", ".fna.gz")
_SHA256_PATTERN = re.compile(r"[0-9a-fA-F]{64}")
_LOGGER = get_logger(__name__)


@dataclass(frozen=True, slots=True)
class Contig:
    """One named reference sequence and its positive length in bases."""

    name: str
    length: int

    def __post_init__(self) -> None:
        """Reject ambiguous or unusable contig metadata."""
        if not isinstance(self.name, str) or not self.name:
            raise ReferenceError("Contig name must be a non-empty string.")
        if (
            not isinstance(self.length, int)
            or isinstance(self.length, bool)
            or self.length <= 0
        ):
            raise ReferenceError(
                f"Contig '{self.name}' length must be a positive integer."
            )

    def to_dict(self) -> dict[str, object]:
        """Return JSON/YAML-friendly contig metadata."""
        return {"name": self.name, "length": self.length}


@dataclass(frozen=True, slots=True)
class ReferenceGenome:
    """Immutable metadata for one indexed, uncompressed reference FASTA.

    Use :meth:`from_fasta` when loading a reference from disk. Direct
    construction is intended for already validated metadata, such as a future
    provenance manifest, and still enforces internal metadata invariants.
    """

    fasta: Path
    fai: Path
    build: str | None
    contigs: tuple[Contig, ...]
    sha256: str | None = None

    def __post_init__(self) -> None:
        """Normalize immutable values and enforce metadata invariants."""
        object.__setattr__(self, "fasta", _coerce_model_path(self.fasta, "FASTA"))
        object.__setattr__(self, "fai", _coerce_model_path(self.fai, "FASTA index"))
        object.__setattr__(self, "build", normalize_reference_build(self.build))
        object.__setattr__(self, "contigs", tuple(self.contigs))
        object.__setattr__(self, "sha256", _normalize_sha256(self.sha256))

        if not self.contigs:
            raise ReferenceError("Reference contig collection is empty.")

        names: set[str] = set()
        for contig in self.contigs:
            if not isinstance(contig, Contig):
                raise ReferenceError(
                    "Reference contigs must contain only Contig objects."
                )
            if contig.name in names:
                raise ReferenceError(
                    f"Reference contig collection contains duplicate contig "
                    f"'{contig.name}'."
                )
            names.add(contig.name)

    @classmethod
    def from_fasta(
        cls,
        fasta: PathInput,
        *,
        build: str | None = None,
        compute_checksum: bool = False,
        sha256: str | None = None,
    ) -> ReferenceGenome:
        """Load metadata from a FASTA and its required conventional FAI.

        The FASTA receives the lightweight Phase 0 validation only. Contig
        names and lengths are streamed from ``<fasta>.fai`` in index order.
        The complete FASTA is read only when ``compute_checksum`` is true.

        Raises:
            ReferenceError: If the FASTA, FAI, or reference metadata is invalid.
        """
        fasta_path = _coerce_model_path(fasta, "FASTA")
        if _has_compressed_fasta_suffix(fasta_path):
            raise ReferenceError(
                f"Compressed FASTA is currently unsupported as a primary "
                f"workflow reference: '{fasta_path}'. Use an uncompressed "
                ".fa, .fasta, or .fna reference."
            )

        _LOGGER.debug("Loading reference FASTA: %s", fasta_path)
        try:
            validation.validate_fasta(fasta_path)
        except InputValidationError as error:
            raise ReferenceError(
                f"Invalid reference FASTA '{fasta_path}': {error}"
            ) from error

        fai_path = Path(f"{fasta_path}.fai")
        _LOGGER.debug("Loading reference FASTA index: %s", fai_path)
        try:
            indexed_contigs = validation.read_fai_contigs(fai_path)
        except ReferenceError as error:
            if not fai_path.exists():
                raise ReferenceError(
                    f"FASTA index missing for '{fasta_path}': expected "
                    f"'{fai_path}'. Create FASTA index before using this "
                    "reference."
                ) from error
            raise

        normalized_provided_checksum = _normalize_sha256(sha256)
        checksum = normalized_provided_checksum
        if compute_checksum:
            try:
                checksum = validation.compute_sha256(fasta_path)
            except InputValidationError as error:
                raise ReferenceError(
                    f"Unable to checksum reference FASTA '{fasta_path}': {error}"
                ) from error
            if (
                normalized_provided_checksum is not None
                and checksum != normalized_provided_checksum
            ):
                raise ReferenceError(
                    f"Provided SHA256 does not match reference FASTA "
                    f"'{fasta_path}'."
                )

        reference = cls(
            fasta=fasta_path,
            fai=fai_path,
            build=build,
            contigs=tuple(
                Contig(name=name, length=length)
                for name, length in indexed_contigs.items()
            ),
            sha256=checksum,
        )
        _LOGGER.debug(
            "Reference loaded: fasta=%s, fai=%s, build=%s, contig_count=%d",
            reference.fasta,
            reference.fai,
            reference.build,
            len(reference.contigs),
        )
        return reference

    @property
    def contig_names(self) -> tuple[str, ...]:
        """Return contig names in deterministic FAI order."""
        return tuple(contig.name for contig in self.contigs)

    def get_contig(self, name: str) -> Contig:
        """Return named contig metadata without normalizing the name."""
        for contig in self.contigs:
            if contig.name == name:
                return contig
        raise ReferenceError(f"Reference does not contain contig '{name}'.")

    def validate_contigs(self, query_contigs: Iterable[str]) -> None:
        """Require query contigs to be an exact-name subset of this reference."""
        validation.validate_contig_compatibility(
            self.contig_names,
            query_contigs,
        )

    def with_checksum(self) -> ReferenceGenome:
        """Return a new instance with a freshly computed FASTA SHA256."""
        try:
            checksum = validation.compute_sha256(self.fasta)
        except InputValidationError as error:
            raise ReferenceError(
                f"Unable to checksum reference FASTA '{self.fasta}': {error}"
            ) from error
        return replace(self, sha256=checksum)

    def to_dict(self, *, include_contigs: bool = False) -> dict[str, object]:
        """Return JSON/YAML-friendly reference metadata.

        The default summary omits the potentially large contig list while
        retaining its count. Callers may request full contig metadata.
        """
        metadata: dict[str, object] = {
            "build": self.build,
            "fasta": str(self.fasta),
            "fai": str(self.fai),
            "sha256": self.sha256,
            "contig_count": len(self.contigs),
        }
        if include_contigs:
            metadata["contigs"] = [contig.to_dict() for contig in self.contigs]
        return metadata


def normalize_reference_build(build: str | None) -> str | None:
    """Canonicalize only explicit, unambiguous reference-build aliases."""
    if build is None:
        return None
    if not isinstance(build, str) or not build.strip():
        raise ReferenceError(
            "Reference build must be a non-empty string or None."
        )
    stripped = build.strip()
    return REFERENCE_BUILD_ALIASES.get(stripped.casefold(), stripped)


def _coerce_model_path(path: PathInput, label: str) -> Path:
    """Normalize a model path without resolving or translating platforms."""
    if not isinstance(path, (str, Path)):
        raise ReferenceError(f"{label} path must be a string or Path.")
    if isinstance(path, str) and not path.strip():
        raise ReferenceError(f"{label} path must not be empty.")
    return Path(path).expanduser()


def _normalize_sha256(checksum: str | None) -> str | None:
    """Validate an optional hexadecimal SHA256 value."""
    if checksum is None:
        return None
    if not isinstance(checksum, str) or _SHA256_PATTERN.fullmatch(checksum) is None:
        raise ReferenceError("Reference SHA256 must contain 64 hexadecimal characters.")
    return checksum.lower()


def _has_compressed_fasta_suffix(path: Path) -> bool:
    """Return whether a path uses a recognized gzip FASTA suffix."""
    lowered_name = path.name.lower()
    return any(lowered_name.endswith(suffix) for suffix in _COMPRESSED_FASTA_SUFFIXES)


__all__ = [
    "Contig",
    "REFERENCE_BUILD_ALIASES",
    "ReferenceGenome",
    "normalize_reference_build",
]
