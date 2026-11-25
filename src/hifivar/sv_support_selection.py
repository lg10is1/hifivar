"""Optional source-support views, retaining native harmonization as authority.

This module does not change Jasmine clustering or infer variant confidence.
Only a fresh, separate bundle is published; its manifest seals success last.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import re

from hifivar.command import CommandRunner
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar.harmonization import EvidenceRunStatus, HarmonizedSvArtifact, validate_harmonized_artifact
from hifivar.reference import ReferenceGenome
from hifivar.sample import validate_sample_id
from hifivar.sv import validate_sv_vcf
from hifivar.sv_input_identity import file_identity

_VERSION = re.compile(r"(\d+(?:\.\d+)+)")
_REQUIRED_COLUMNS = ("harmonized_variant_id", "sample", "contig", "start", "support_count", "supporting_callers")
_POLICY = {
    "operator": ">=", "count_from": "SUPP_VEC_POPCOUNT", "source_order": "ACTUAL_JASMINE_INPUT_ORDER",
    "record_representation": "VERBATIM_RETAINED_LINES", "canonical_output_role_unchanged": True,
    "caller_count_is_confidence": False, "harmonization_is_truth": False,
    "bnd_mate_completeness": "NOT_ASSESSED", "downstream_consumer_opt_in_required": True,
}


def _positive_integer(value: object) -> int:
    if type(value) is not int or value < 1:
        raise InputValidationError("minimum_caller_support must be a positive integer, not a boolean.")
    return value


def support_selection_paths(outputroot: str | Path, sample: str, k: int) -> dict[str, Path]:
    validate_sample_id(sample)
    _positive_integer(k)
    directory = Path(outputroot).absolute() / sample / f"support-ge-{k}"
    stem = f"{sample}.support-ge-{k}"
    vcf = directory / f"{stem}.sv.vcf.gz"
    return {"directory": directory, "vcf": vcf, "index": Path(f"{vcf}.tbi"),
            "evidence": directory / f"{stem}.sv.evidence.tsv", "manifest": directory / "manifest.json"}


@dataclass(frozen=True, slots=True)
class SupportSelectionRequest:
    artifact: HarmonizedSvArtifact
    reference: ReferenceGenome
    minimum_caller_support: int
    output_directory: Path

    def __post_init__(self) -> None:
        if not isinstance(self.artifact, HarmonizedSvArtifact) or not isinstance(self.reference, ReferenceGenome):
            raise InputValidationError("Support selection requires a full HarmonizedSvArtifact and ReferenceGenome.")
        validate_sample_id(self.artifact.sample_id)
        k = _positive_integer(self.minimum_caller_support)
        if (Path(self.artifact.reference_fasta).absolute() != self.reference.fasta.absolute()
                or self.artifact.reference_build != self.reference.build):
            raise InputValidationError("Support selection reference differs from native harmonization.")
        sources = self.runnable_sources
        if not sources or len({source.caller for source in sources}) != len(sources):
            raise InputValidationError("Support selection requires distinct actual runnable caller sources.")
        if any(source.sample_id != self.artifact.sample_id for source in sources):
            raise InputValidationError("Support selection sources belong to another sample.")
        if k > len(sources):
            raise InputValidationError("minimum_caller_support exceeds the actual Jasmine input source count.")
        directory = Path(self.output_directory).absolute()
        if directory.name != f"support-ge-{k}" or directory.parent.name != self.artifact.sample_id:
            raise InputValidationError("Support selection bundle must follow {sample}/support-ge-{k}.")
        if directory.exists() or directory.is_symlink():
            raise OutputValidationError("Support selection requires a fresh bundle; overwrite is unsupported.")
        object.__setattr__(self, "output_directory", directory)

    @property
    def runnable_sources(self):
        return tuple(source for source in self.artifact.sources
                     if source.status in {EvidenceRunStatus.COMPLETED, EvidenceRunStatus.NO_CALLS})

    @property
    def paths(self) -> dict[str, Path]:
        return support_selection_paths(self.output_directory.parent.parent, self.artifact.sample_id,
                                       self.minimum_caller_support)

    def to_dict(self) -> dict[str, object]:
        return {"sample_id": self.artifact.sample_id, "minimum_caller_support": self.minimum_caller_support,
                "source_order": [source.caller for source in self.runnable_sources],
                "source_vcf": str(self.artifact.vcf_path), "source_index": str(self.artifact.index_path),
                "source_evidence": str(self.artifact.evidence_table),
                "paths": {key: str(value) for key, value in self.paths.items()}, "policy": dict(_POLICY)}


def _identity(path: Path) -> dict[str, object]:
    value = file_identity(path)
    return {key: value[key] for key in ("path", "sha256", "bytes")}


def _source_identities(request: SupportSelectionRequest) -> dict[str, dict[str, object]]:
    artifact = request.artifact
    result = {"vcf": _identity(artifact.vcf_path), "index": _identity(artifact.index_path),
              "evidence": _identity(artifact.evidence_table)}
    if artifact.source_identity_manifest is not None:
        result["source_identity_manifest"] = _identity(artifact.source_identity_manifest)
    return result


def _parse_support(parts: list[bytes], n: int) -> tuple[int, bytes]:
    found: dict[bytes, bytes] = {}
    for item in parts[7].split(b";"):
        key, separator, value = item.partition(b"=")
        if key in {b"SUPP_VEC", b"SUPP"}:
            if not separator or key in found:
                raise OutputValidationError("Missing value or duplicate SUPP_VEC/SUPP INFO field.")
            found[key] = value
    vector = found.get(b"SUPP_VEC")
    if vector is None or len(vector) != n or any(bit not in b"01" for bit in vector):
        raise OutputValidationError("SUPP_VEC must match the actual Jasmine source order and length.")
    count = vector.count(b"1")
    if count == 0:
        raise OutputValidationError("SUPP_VEC cannot be all zero.")
    if b"SUPP" in found and (not found[b"SUPP"].isdigit() or int(found[b"SUPP"]) != count):
        raise OutputValidationError("SUPP differs from authoritative SUPP_VEC popcount.")
    return count, vector


def _iter_joined(request: SupportSelectionRequest):
    """Yield native VCF lines paired to original evidence rows in record order."""
    with gzip.open(request.artifact.vcf_path, "rb") as native, request.artifact.evidence_table.open("rb") as table:
        header = table.readline()
        columns = header.rstrip(b"\r\n").decode("utf-8").split("\t")
        if len(columns) != len(set(columns)) or any(key not in columns for key in _REQUIRED_COLUMNS):
            raise OutputValidationError("Full evidence table lacks unique required support-selection columns.")
        indices = {name: columns.index(name) for name in _REQUIRED_COLUMNS}
        yield None, header, None, None
        seen_header = False
        seen_record = False
        for line in native:
            if line.startswith(b"#"):
                if seen_record:
                    raise OutputValidationError("VCF header appears after data.")
                if line.startswith(b"#CHROM\t"):
                    if seen_header or line.rstrip(b"\r\n").split(b"\t")[9:] != [request.artifact.sample_id.encode("utf-8")]:
                        raise OutputValidationError("Native VCF sample header differs from support-selection sample.")
                    seen_header = True
                yield line, None, None, None
                continue
            if not seen_header:
                raise OutputValidationError("Native VCF is missing its sample header.")
            seen_record = True
            parts = line.rstrip(b"\r\n").split(b"\t")
            if len(parts) != 10:
                raise OutputValidationError("Support selection requires one-sample native VCF records.")
            count, vector = _parse_support(parts, len(request.runnable_sources))
            row = table.readline()
            if not row:
                raise OutputValidationError("Evidence table is shorter than the full native VCF.")
            values = row.rstrip(b"\r\n").split(b"\t")
            if len(values) != len(columns):
                raise OutputValidationError("Evidence table row has an inconsistent column count.")
            callers = b",".join(source.caller.encode("utf-8") for source, bit in zip(request.runnable_sources, vector) if bit == ord("1"))
            expected = {"harmonized_variant_id": parts[2], "sample": request.artifact.sample_id.encode("utf-8"),
                        "contig": parts[0], "start": parts[1], "support_count": str(count).encode("ascii"),
                        "supporting_callers": callers}
            if any(values[indices[name]] != value for name, value in expected.items()):
                raise OutputValidationError("Full evidence row differs from its native VCF record/source membership.")
            is_bnd = b"[" in parts[4] or b"]" in parts[4] or b"SVTYPE=BND" in parts[7].split(b";")
            yield line, row, count, is_bnd
        if not seen_header or table.readline():
            raise OutputValidationError("Missing VCF header or extra evidence rows.")


def _derive(request: SupportSelectionRequest, plain: Path, evidence: Path) -> dict[str, object]:
    histogram: Counter[int] = Counter()
    selected = total = bnd = 0
    with plain.open("xb") as writer, evidence.open("xb") as table:
        for line, row, count, is_bnd in _iter_joined(request):
            if count is None:
                if line is not None:
                    writer.write(line)
                if row is not None:
                    table.write(row)
                continue
            total += 1
            histogram[count] += 1
            if count >= request.minimum_caller_support:
                selected += 1
                bnd += bool(is_bnd)
                writer.write(line)
                table.write(row)
    return {"source_records": total, "selected_records": selected, "excluded_records": total - selected,
            "source_support_histogram": {str(k): histogram[k] for k in sorted(histogram)},
            "selected_bnd_records": bnd, "bnd_mate_completeness": "NOT_ASSESSED"}


def _readback(request: SupportSelectionRequest, vcf: Path, evidence: Path) -> int:
    count = 0
    with gzip.open(vcf, "rb") as result, evidence.open("rb") as table:
        for native, row, support, _ in _iter_joined(request):
            if support is not None and support < request.minimum_caller_support:
                continue
            if native is not None and next(result, None) != native:
                raise OutputValidationError("Selected VCF header/record bytes or order changed.")
            if row is not None and next(table, None) != row:
                raise OutputValidationError("Selected evidence row bytes or order changed.")
            count += support is not None
        if next(result, None) is not None or next(table, None) is not None:
            raise OutputValidationError("Selected bundle contains extra VCF or evidence records.")
    return count


@dataclass(frozen=True, slots=True)
class SupportSelectionResult:
    request: SupportSelectionRequest
    status: str
    commands: tuple[tuple[str, ...], ...]
    tool_versions: dict[str, str] | None = None
    manifest_path: Path | None = None
    manifest_sha256: str | None = None
    selected_records: int | None = None
    source_records: int | None = None

    @property
    def vcf_path(self) -> Path:
        return self.request.paths["vcf"]

    @property
    def index_path(self) -> Path:
        return self.request.paths["index"]

    @property
    def evidence_table(self) -> Path:
        return self.request.paths["evidence"]

    def to_dict(self) -> dict[str, object]:
        return {"status": self.status, "executed": self.status == "completed", "request": self.request.to_dict(),
                "commands": [list(command) for command in self.commands], "tool_versions": self.tool_versions,
                "manifest_path": str(self.manifest_path) if self.manifest_path else None,
                "manifest_sha256": self.manifest_sha256, "selected_records": self.selected_records,
                "source_records": self.source_records, "canonical_output_role_unchanged": True,
                "caller_count_is_confidence": False}

    def validate(self) -> dict[str, object]:
        if self.status != "completed" or self.manifest_path is None or self.manifest_sha256 is None:
            raise OutputValidationError("Planned support selection cannot be validated as completed.")
        if _identity(self.manifest_path)["sha256"] != self.manifest_sha256:
            raise OutputValidationError("Support-selection manifest changed after publication.")
        payload = json.loads(self.manifest_path.read_text(encoding="utf-8"))
        if (payload.get("schema_version") != "1.0" or payload.get("status") != "completed"
                or payload.get("request") != self.request.to_dict() or payload.get("policy") != _POLICY):
            raise OutputValidationError("Support-selection manifest does not match its bound request/policy.")
        if _source_identities(self.request) != payload.get("source_files"):
            raise OutputValidationError("Full canonical source bytes changed after support selection.")
        for role in ("vcf", "index", "evidence"):
            if _identity(self.request.paths[role]) != payload["output_files"].get(role):
                raise OutputValidationError("Selected output bytes changed after support selection.")
        validate_sv_vcf(self.vcf_path, self.index_path, sample_id=self.request.artifact.sample_id, reference=self.request.reference)
        if _readback(self.request, self.vcf_path, self.evidence_table) != payload["counts"]["selected_records"]:
            raise OutputValidationError("Support-selection record count differs from manifest.")
        return payload


class SupportSelectionWrapper:
    def __init__(self, *, bgzip_executable="bgzip", tabix_executable="tabix", runner: CommandRunner | None = None):
        self.bgzip = bgzip_executable
        self.tabix = tabix_executable
        self.runner = runner or CommandRunner()

    def _commands(self, plain: Path, vcf: Path) -> tuple[tuple[str, ...], ...]:
        return ((self.bgzip, "-c", str(plain)), (self.bgzip, "-t", str(vcf)),
                (self.tabix, "-p", "vcf", str(vcf)), (self.tabix, "-l", str(vcf)))

    def run(self, request: SupportSelectionRequest, *, dry_run=False, stderr_path=None) -> SupportSelectionResult:
        if not isinstance(request, SupportSelectionRequest):
            raise InputValidationError("Support selection requires SupportSelectionRequest.")
        paths = request.paths
        stage = paths["directory"] / ".working"
        plain, vcf = stage / "selected.vcf", stage / paths["vcf"].name
        table = stage / paths["evidence"].name
        commands = self._commands(plain, vcf)
        if dry_run:
            for index, command in enumerate(commands):
                self.runner.run(command, dry_run=True, stdout_path=vcf if index == 0 else None)
            return SupportSelectionResult(request, "planned", commands)
        validate_harmonized_artifact(request.artifact, request.reference)
        before = _source_identities(request)
        versions = {}
        for name, executable in (("bgzip", self.bgzip), ("tabix", self.tabix)):
            result = self.runner.run([executable, "--version"])
            match = _VERSION.search("\n".join(item for item in (result.stdout, result.stderr) if item))
            if match is None:
                raise ToolVersionError(f"Unable to parse {name} version.")
            versions[name] = match[1]
        with gzip.open(request.artifact.vcf_path, "rb") as handle:
            fileformat = handle.readline().rstrip(b"\r\n")
        if fileformat == b"##fileformat=VCFv4.4" and tuple(int(x) for x in versions["tabix"].split(".")) < (1, 23, 1):
            raise ToolVersionError("VCFv4.4 support selection requires HTSlib/tabix >=1.23.1.")
        paths["directory"].parent.mkdir(parents=True, exist_ok=True)
        # The exclusive directory creation is the no-overwrite reservation.
        # A failed bundle has no success manifest and is retained for diagnosis.
        try:
            paths["directory"].mkdir()
        except FileExistsError as error:
            raise OutputValidationError("Support selection requires a fresh bundle; overwrite is unsupported.") from error
        stage.mkdir()
        counts = _derive(request, plain, table)
        records = []
        base_log = Path(stderr_path) if stderr_path is not None else stage / "support_selection.log"
        for index, command in enumerate(commands):
            log = base_log.with_name(f"{base_log.stem}.{index + 1}{base_log.suffix}")
            result = self.runner.run(command, stdout_path=vcf if index == 0 else None,
                                     stderr_path=log, timeout=300)
            records.append({"command": list(command), "returncode": result.returncode,
                            "duration_seconds": result.duration_seconds, "stderr_path": str(log)})
        index_path = Path(f"{vcf}.tbi")
        validate_sv_vcf(vcf, index_path, sample_id=request.artifact.sample_id, reference=request.reference)
        if _readback(request, vcf, table) != counts["selected_records"]:
            raise OutputValidationError("Selected count differs from retained-record readback.")
        if _source_identities(request) != before:
            raise OutputValidationError("Full canonical source bytes changed during support selection.")
        for origin, role in ((vcf, "vcf"), (index_path, "index"), (table, "evidence")):
            origin.rename(paths[role])
        outputs = {role: _identity(paths[role]) for role in ("vcf", "index", "evidence")}
        payload = {"schema_version": "1.0", "status": "completed", "request": request.to_dict(),
                   "policy": dict(_POLICY), "reference": request.reference.to_dict(),
                   "source_files": before, "output_files": outputs, "counts": counts,
                   "tool_versions": versions, "commands": records,
                   "validation": {"retained_vcf_lines_unchanged": True, "retained_evidence_lines_unchanged": True,
                                  "source_files_unchanged": True, "bgzf_integrity": True, "tabix_readable": True,
                                  "empty_selection_is_valid": True}}
        with paths["manifest"].open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(payload, handle, ensure_ascii=True, sort_keys=True, indent=2)
            handle.write("\n")
        result = SupportSelectionResult(request, "completed", commands, versions, paths["manifest"],
                                        str(_identity(paths["manifest"])["sha256"]), counts["selected_records"], counts["source_records"])
        result.validate()
        return result


__all__ = ["SupportSelectionRequest", "SupportSelectionResult", "SupportSelectionWrapper", "support_selection_paths"]
