"""Streaming SV evidence model and conservative normalization boundary."""

from __future__ import annotations

import gzip
from contextlib import ExitStack
import json
import sqlite3
import uuid
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from collections.abc import Iterator

from hifivar.assembly_sv import AssemblySvArtifact, SVEvidenceSource
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.reference import ReferenceGenome
from hifivar.sv import StructuralVariantArtifact, validate_sv_vcf
from hifivar.validation import validate_output_file
from hifivar.sv_input_identity import SourceIdentityPolicy, SourceIdentityResult


class EvidenceRunStatus(str, Enum):
    COMPLETED = "completed"
    NO_CALLS = "no_calls"
    NOT_RUN = "not_run"
    FAILED = "failed"
    DISABLED = "disabled"


class EvidenceClass(str, Enum):
    READ_ONLY = "READ_ONLY"
    ASSEMBLY_ONLY = "ASSEMBLY_ONLY"
    READ_AND_ASSEMBLY = "READ_AND_ASSEMBLY"
    UNRESOLVED = "UNRESOLVED"


@dataclass(frozen=True, slots=True)
class SVEvidenceSourceArtifact:
    sample_id: str
    source: SVEvidenceSource
    caller: str
    vcf_path: Path | None
    index_path: Path | None
    status: EvidenceRunStatus
    haplotypes: tuple[str, ...] = ()
    error: str | None = None

    @classmethod
    def from_read(cls, artifact: StructuralVariantArtifact) -> "SVEvidenceSourceArtifact":
        return cls(artifact.sample_id, SVEvidenceSource.READ, artifact.caller.value,
                   artifact.vcf_path, artifact.index_path, EvidenceRunStatus.COMPLETED)

    @classmethod
    def from_assembly(cls, artifact: AssemblySvArtifact) -> "SVEvidenceSourceArtifact":
        return cls(artifact.sample_id, SVEvidenceSource.ASSEMBLY, artifact.caller.value,
                   artifact.vcf_path, artifact.index_path, EvidenceRunStatus.COMPLETED,
                   tuple(item.role.value for item in artifact.assemblies))

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id, "source": self.source.value, "caller": self.caller,
            "vcf_path": str(self.vcf_path) if self.vcf_path else None,
            "index_path": str(self.index_path) if self.index_path else None,
            "status": self.status.value, "haplotypes": list(self.haplotypes), "error": self.error,
            "haplotype_role_status": (
                "NOT_RUN" if self.status not in {EvidenceRunStatus.COMPLETED, EvidenceRunStatus.NO_CALLS}
                else "DECLARED_SOURCE_INPUT_ROLES_NOT_VARIANT_PHASE_ACCURACY" if self.haplotypes
                else "UNVERIFIED_ASSEMBLY_INPUT_ROLES" if self.source is SVEvidenceSource.ASSEMBLY
                else "NOT_APPLICABLE_READ_SOURCE"
            ),
        }


@dataclass(frozen=True, slots=True)
class SVEvidenceCollection:
    sample_id: str
    sources: tuple[SVEvidenceSourceArtifact, ...]

    def __post_init__(self) -> None:
        if not self.sample_id or not self.sources:
            raise InputValidationError("SV evidence collection requires a sample and sources.")
        if any(item.sample_id != self.sample_id for item in self.sources):
            raise InputValidationError("All SV evidence sources must belong to one sample.")
        identities = [(item.source, item.caller) for item in self.sources]
        if len(identities) != len(set(identities)):
            raise InputValidationError("SV evidence collection contains a duplicate caller source.")

    def to_dict(self) -> dict[str, object]:
        return {"sample_id": self.sample_id, "sources": [item.to_dict() for item in self.sources]}
@dataclass(frozen=True, slots=True)
class SVEvidenceRecord:
    sample_id: str
    evidence_source: SVEvidenceSource
    caller: str
    source_vcf: Path
    original_record_id: str
    deterministic_record_id: str
    contig: str
    start: int
    end: int | None
    svtype: str
    native_svtype: str | None
    svlen: int | None
    insertion_sequence_length: int | None
    breakpoint_uncertainty: tuple[str | None, str | None]
    haplotypes: tuple[str, ...]
    native_info: tuple[tuple[str, str | None], ...]
    unresolved: bool

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id, "evidence_source": self.evidence_source.value,
            "caller": self.caller, "source_vcf": str(self.source_vcf),
            "original_record_id": self.original_record_id,
            "deterministic_record_id": self.deterministic_record_id,
            "contig": self.contig, "start": self.start, "end": self.end,
            "svtype": self.svtype, "native_svtype": self.native_svtype,
            "svlen": self.svlen, "insertion_sequence_length": self.insertion_sequence_length,
            "breakpoint_uncertainty": list(self.breakpoint_uncertainty),
            "haplotypes": list(self.haplotypes), "native_info": dict(self.native_info),
            "unresolved": self.unresolved,
        }

@dataclass(frozen=True, slots=True)
class SVHarmonizationRequest:
    sample_id: str
    reference: ReferenceGenome
    sources: tuple[SVEvidenceSourceArtifact, ...]
    work_directory: Path
    output_vcf: Path
    evidence_table: Path
    max_dist: int = 1000
    distance_type: str = "linear"
    overwrite: bool = False
    input_identity_policy: SourceIdentityPolicy | None = None

    def __post_init__(self) -> None:
        if self.input_identity_policy is not None and not isinstance(self.input_identity_policy, SourceIdentityPolicy):
            raise InputValidationError("Harmonization input identity requires an explicit supported policy.")
        if not self.sample_id or not isinstance(self.reference, ReferenceGenome):
            raise InputValidationError("Harmonization requires sample_id and ReferenceGenome.")
        if not self.sources or any(item.sample_id != self.sample_id for item in self.sources):
            raise InputValidationError("Harmonization sources must belong to one sample.")
        completed = [item for item in self.sources if item.status in {EvidenceRunStatus.COMPLETED, EvidenceRunStatus.NO_CALLS}]
        if not completed:
            raise InputValidationError("Harmonization requires at least one completed caller source.")
        if self.distance_type not in {"linear", "nonlinear"}:
            raise InputValidationError("Harmonization distance_type must be linear or nonlinear.")
        if isinstance(self.max_dist, bool) or not isinstance(self.max_dist, int) or self.max_dist <= 0:
            raise InputValidationError("Harmonization max_dist must be positive.")
        output = Path(self.output_vcf)
        if output.name != f"{self.sample_id}.harmonized.sv.vcf.gz":
            raise InputValidationError("Harmonized output must follow {sample}.harmonized.sv.vcf.gz.")
        table = Path(self.evidence_table)
        if table.name != f"{self.sample_id}.sv.evidence.tsv":
            raise InputValidationError("Evidence table must follow {sample}.sv.evidence.tsv.")
        object.__setattr__(self, "work_directory", Path(self.work_directory))
        object.__setattr__(self, "output_vcf", output)
        object.__setattr__(self, "evidence_table", table)
        for path in (output, Path(f"{output}.tbi"), table):
            if path.exists() and not self.overwrite:
                raise OutputValidationError(f"Harmonization output already exists: '{path}'.")

    @property
    def output_index(self) -> Path:
        return Path(f"{self.output_vcf}.tbi")

    @property
    def runnable_sources(self) -> tuple[SVEvidenceSourceArtifact, ...]:
        return tuple(item for item in self.sources if item.status in {EvidenceRunStatus.COMPLETED, EvidenceRunStatus.NO_CALLS})

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id, "reference": self.reference.to_dict(),
            "sources": [item.to_dict() for item in self.sources],
            "work_directory": str(self.work_directory), "output_vcf": str(self.output_vcf),
            "output_index": str(self.output_index), "evidence_table": str(self.evidence_table),
            "max_dist": self.max_dist, "distance_type": self.distance_type,
            "overwrite": self.overwrite,
            **({"input_identity": self.input_identity_policy.to_dict()} if self.input_identity_policy else {}),
        }

@dataclass(frozen=True, slots=True)
class HarmonizedSvArtifact:
    sample_id: str
    reference_fasta: Path
    reference_build: str | None
    vcf_path: Path
    index_path: Path
    evidence_table: Path
    sources: tuple[SVEvidenceSourceArtifact, ...]
    jasmine_version: str
    commands: tuple[tuple[str, ...], ...]
    harmonized: bool = True
    intermediate_files: tuple[Path, ...] = ()
    truth_label: None = None
    source_identity_manifest: Path | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id, "reference_fasta": str(self.reference_fasta),
            "reference_build": self.reference_build, "vcf_path": str(self.vcf_path),
            "index_path": str(self.index_path), "evidence_table": str(self.evidence_table),
            "sources": [item.to_dict() for item in self.sources],
            "jasmine_version": self.jasmine_version,
            "commands": [list(item) for item in self.commands],
            "harmonized": self.harmonized, "truth_label": self.truth_label,
            "intermediate_files": [str(path) for path in self.intermediate_files],
            **({"source_identity_manifest": str(self.source_identity_manifest)} if self.source_identity_manifest else {}),
        }


def iter_sv_evidence(source: SVEvidenceSourceArtifact) -> Iterator[SVEvidenceRecord]:
    """Stream canonical fields while retaining native INFO and representation."""
    if source.status not in {EvidenceRunStatus.COMPLETED, EvidenceRunStatus.NO_CALLS}:
        return
    if source.vcf_path is None:
        raise InputValidationError("Completed evidence source lacks a VCF path.")
    with _open_vcf(source.vcf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 8:
                raise OutputValidationError(f"Malformed SV record in '{source.vcf_path}'.")
            contig, position, native_id, ref, alt, _, _, info_text = fields[:8]
            info = _parse_info(info_text)
            native_type = info.get("SVTYPE")
            known = {"DEL", "INS", "DUP", "INV", "BND", "CNV"}
            canonical = native_type.upper() if native_type and native_type.upper() in known else "UNRESOLVED"
            start = int(position)
            end = _integer(info.get("END"))
            svlen = _integer((info.get("SVLEN") or "").split(",")[0])
            insertion_length = len(alt) - len(ref) if canonical == "INS" and not alt.startswith("<") else None
            record_id = native_id if native_id not in {"", "."} else f"{contig}:{start}:{ref}:{alt}"
            yield SVEvidenceRecord(
                source.sample_id, source.source, source.caller, source.vcf_path,
                record_id, f"{source.caller}:{record_id}", contig, start, end,
                canonical, native_type, svlen, insertion_length,
                (info.get("CIPOS"), info.get("CIEND")), source.haplotypes,
                tuple(info.items()), canonical == "UNRESOLVED",
            )

def write_evidence_table(request: SVHarmonizationRequest, merged_vcf: Path,
                         *, source_identity: SourceIdentityResult | None = None) -> Path:
    """Stream Jasmine records into an evidence table; never assign confidence."""
    request.evidence_table.parent.mkdir(parents=True, exist_ok=True)
    temporary = request.evidence_table.with_name(f".{request.evidence_table.name}.hifivar.tmp")
    if temporary.exists():
        raise OutputValidationError(f"Evidence-table temporary exists: '{temporary}'.")
    source_by_caller = {item.caller: item for item in request.sources}
    header = (
        "harmonized_variant_id\tsample\tsvtype\tcontig\tstart\tend\t"
        "evidence_class\tsupport_count\tread_support_count\tassembly_support_count\t"
        "supporting_callers\tsource_record_ids\tsource_files\thaplotypes\n"
    )
    identity_payload = None
    if source_identity is not None:
        identity_payload = source_identity.validate()
        if source_identity.source_order != tuple(s.caller for s in request.runnable_sources):
            raise OutputValidationError("Evidence identity source order differs from the actual Jasmine list.")
        for source, recorded in zip(request.runnable_sources, identity_payload["canonical_sources"]):
            if (str(source.vcf_path.resolve()) != recorded["vcf"]["path"]
                    or str(source.index_path.resolve()) != recorded["index"]["path"]):
                raise OutputValidationError("Evidence identity refers to different canonical inputs.")
        header = header.rstrip("\n") + "\tanalysis_record_ids\tsource_record_keys\tsource_record_callers\trepresentation_files\tsource_identity_manifest\n"
    with ExitStack() as stack:
        identity_db = emissions = None
        if source_identity is not None:
            identity_db = sqlite3.connect(source_identity.database_path.as_uri() + "?mode=ro", uri=True)
            stack.callback(identity_db.close)
            request.work_directory.mkdir(parents=True, exist_ok=True)
            emission_path = request.work_directory / f"evidence_membership_{uuid.uuid4().hex}.sqlite3"
            emissions = sqlite3.connect(emission_path)
            stack.callback(emissions.close)
            emissions.execute("CREATE TABLE used (analysis_id TEXT PRIMARY KEY)")
        writer = stack.enter_context(temporary.open("x", encoding="utf-8", newline="\n"))
        writer.write(header)
        with _open_vcf(merged_vcf) as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 8:
                    raise OutputValidationError("Malformed harmonized SV record.")
                info = _parse_info(fields[7])
                identifiers = tuple(item for item in (info.get("IDLIST") or fields[2]).split(",") if item)
                support_vector = info.get("SUPP_VEC")
                runnable = request.runnable_sources
                if source_identity is not None and (support_vector is None or not info.get("IDLIST")):
                    raise OutputValidationError("Explicit source identity requires SUPP_VEC and IDLIST, without fallback inference.")
                if support_vector is not None:
                    if len(support_vector) != len(runnable) or any(
                        bit not in {"0", "1"} for bit in support_vector
                    ):
                        raise OutputValidationError("Jasmine SUPP_VEC does not match the input source list.")
                    callers = tuple(
                        source.caller for source, bit in zip(runnable, support_vector) if bit == "1"
                    )
                    if not callers:
                        raise OutputValidationError("Jasmine SUPP_VEC does not identify a supporting source.")
                else:
                    # Legacy Jasmine output may lack SUPP_VEC.  In that case only,
                    # retain the conservative caller-prefixed ID fallback.  Native
                    # caller IDs are not otherwise a reliable membership contract.
                    mapped_callers = tuple(
                        item.split(":", 1)[0]
                        for item in identifiers
                        if ":" in item and item.split(":", 1)[0] in source_by_caller
                    )
                    callers = (
                        tuple(dict.fromkeys(mapped_callers))
                        if len(mapped_callers) == len(identifiers)
                        else ()
                    )
                read = tuple(name for name in callers if source_by_caller[name].source is SVEvidenceSource.READ)
                assembly = tuple(name for name in callers if source_by_caller[name].source is SVEvidenceSource.ASSEMBLY)
                if not callers:
                    evidence_class = EvidenceClass.UNRESOLVED
                elif read and assembly:
                    evidence_class = EvidenceClass.READ_AND_ASSEMBLY
                elif assembly:
                    evidence_class = EvidenceClass.ASSEMBLY_ONLY
                else:
                    evidence_class = EvidenceClass.READ_ONLY
                files = tuple(str(source_by_caller[name].vcf_path) for name in callers)
                haplotypes = tuple(hap for name in callers for hap in source_by_caller[name].haplotypes)
                row = (
                    fields[2], request.sample_id, info.get("SVTYPE") or "UNRESOLVED",
                    fields[0], fields[1], info.get("END") or "", evidence_class.value,
                    str(len(callers)), str(len(read)), str(len(assembly)), ",".join(callers),
                    ",".join(identifiers), ",".join(files), ",".join(haplotypes),
                )
                if source_identity is not None:
                    original_ids, keys, record_callers = [], [], []
                    for identifier in identifiers:
                        found = identity_db.execute("SELECT source,ordinal,original_id FROM records WHERE analysis_id=?", (identifier,)).fetchone()
                        if found is None:
                            raise OutputValidationError("IDLIST contains an unknown analysis record ID.")
                        index, ordinal, original_id = found
                        source = runnable[index]
                        if source.caller not in callers:
                            raise OutputValidationError("IDLIST source membership contradicts SUPP_VEC.")
                        try:
                            emissions.execute("INSERT INTO used VALUES (?)", (identifier,))
                        except sqlite3.IntegrityError as error:
                            raise OutputValidationError("A source record appears more than once in the harmonized output.") from error
                        original_ids.append(original_id)
                        # The key comes from the bound canonical file, never from parsing arbitrary caller IDs.
                        source_payload = identity_payload["canonical_sources"][index]
                        keys.append(f"caller:{source.caller}:sha256:{source_payload['vcf']['sha256']}:record:{ordinal}")
                        record_callers.append(source.caller)
                    if set(record_callers) != set(callers):
                        raise OutputValidationError("IDLIST omits a SUPP_VEC supporting source.")
                    row = (*row[:11], ",".join(original_ids), *row[12:], ",".join(identifiers),
                           json.dumps(keys, separators=(",", ":")), ",".join(record_callers),
                           ",".join(str(source_identity.views[i]) for i, s in enumerate(runnable) if s.caller in callers),
                           str(source_identity.manifest_path))
                writer.write("\t".join(row) + "\n")
        if source_identity is not None:
            expected = identity_db.execute("SELECT count(*) FROM records").fetchone()[0]
            if emissions.execute("SELECT count(*) FROM used").fetchone()[0] != expected:
                raise OutputValidationError("Harmonization omitted a canonical source record; identity evidence cannot be sealed.")
            source_identity.validate()
    temporary.replace(request.evidence_table)
    validate_output_file(request.evidence_table)
    return request.evidence_table


def validate_harmonized_artifact(artifact: HarmonizedSvArtifact, reference: ReferenceGenome) -> HarmonizedSvArtifact:
    validate_sv_vcf(artifact.vcf_path, artifact.index_path, sample_id=artifact.sample_id, reference=reference)
    validate_output_file(artifact.evidence_table)
    if not artifact.sources:
        raise OutputValidationError("Harmonized artifact lacks source provenance.")
    return artifact


def _open_vcf(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if str(path).lower().endswith(".gz") else path.open("r", encoding="utf-8", newline="")


def _parse_info(text: str) -> dict[str, str | None]:
    result: dict[str, str | None] = {}
    if text in {"", "."}:
        return result
    for item in text.split(";"):
        key, separator, value = item.partition("=")
        result[key] = value if separator else None
    return result


def _integer(value: str | None) -> int | None:
    if value in {None, "", "."}:
        return None
    try:
        return int(value)
    except ValueError:
        return None


__all__ = [
    "EvidenceClass", "EvidenceRunStatus", "HarmonizedSvArtifact",
    "SVEvidenceCollection", "SVEvidenceRecord", "SVEvidenceSourceArtifact", "SVHarmonizationRequest",
    "iter_sv_evidence", "validate_harmonized_artifact", "write_evidence_table",
]
