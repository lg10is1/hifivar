"""Independent Ensembl VEP offline/cache wrapper."""

from __future__ import annotations

import os
import re
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path

from hifivar.annotation import (
    AnnotationArtifact,
    AnnotationDatabase,
    AnnotationInput,
    AnnotationResult,
    AnnotationRunStatus,
    AnnotationSource,
)
from hifivar.command import CommandRunner, _default_stdout_log_path
from hifivar.annotation_paths import protect_annotation_outputs
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar.validation import validate_output_file


_VERSION = re.compile(
    r"(?:ensembl-vep|\bvep\b|version|release)\s*(?:version\s*)?[:= ]\s*v?(\d+(?:\.\d+)*)",
    re.IGNORECASE,
)


# Ensembl VEP's default tabular fields; flag/plugin columns may be appended.
# The wrapper requests --tab without narrowing --fields. Do not interpret a
# nonempty tool log, partial table or HTML statistics file as annotation output.
_TSV_REQUIRED_COLUMNS = frozenset({
    "Uploaded_variation", "Location", "Allele", "Gene", "Feature",
    "Feature_type", "Consequence", "cDNA_position", "CDS_position",
    "Protein_position", "Amino_acids", "Codons", "Existing_variation",
})


def _validate_vep_tsv(path: Path) -> None:
    """Stream the native --tab structure without rewriting or judging effects.

    Header-only output is legal: VEP need not emit one row per input variant.
    Metadata precedes a single column header; missing native values use '-'.
    This is a format gate, not evidence of biological accuracy or completeness.
    """
    validate_output_file(path)
    columns: list[str] | None = None
    required_indexes: tuple[int, ...] = ()
    line_number = 0

    def invalid(reason: str) -> OutputValidationError:
        # Do not echo annotation payloads or tool logs into error messages.
        return OutputValidationError(
            f"VEP_TSV_INVALID: '{path}', line {line_number}: {reason}."
        )

    try:
        with path.open(encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, start=1):
                line = raw.rstrip("\r\n")
                if not line:
                    continue
                if any(ord(char) < 32 and char != "\t" for char in line):
                    raise invalid("unexpected control character")
                if columns is None:
                    if line.startswith("##"):
                        continue
                    if not line.startswith("#"):
                        raise invalid("data appears before the tabular column header")
                    columns = line[1:].split("\t")
                    if any(not field.strip() for field in columns):
                        raise invalid("empty column name")
                    if len(set(columns)) != len(columns):
                        raise invalid("duplicate column name")
                    if not _TSV_REQUIRED_COLUMNS.issubset(columns):
                        raise invalid("required default --tab columns are missing")
                    required_indexes = tuple(
                        i for i, field in enumerate(columns)
                        if field in _TSV_REQUIRED_COLUMNS
                    )
                    continue
                if line.startswith("#"):
                    raise invalid("repeated header or metadata after the column header")
                fields = line.split("\t")
                if len(fields) != len(columns):
                    raise invalid("row column count differs from the header")
                if any(not fields[i].strip() for i in required_indexes):
                    raise invalid("empty native field; missing values must use '-'")
    except (OSError, UnicodeError) as error:
        raise invalid("output is not readable UTF-8 tabular text") from error
    if columns is None:
        raise invalid("tabular column header is missing")


@dataclass(frozen=True, slots=True)
class VepRequest:
    input: AnnotationInput
    cache_directory: Path
    cache_version: str
    species: str
    assembly: str
    output_tsv: Path
    threads: int = 1
    overwrite: bool = False

    def __post_init__(self) -> None:
        cache = Path(self.cache_directory).expanduser()
        if not cache.is_dir():
            raise InputValidationError(f"VEP cache directory is missing: '{cache}'.")
        for name in ("cache_version", "species", "assembly"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise InputValidationError(f"VEP {name} must be non-empty.")
        if self.input.reference.build != self.assembly:
            raise InputValidationError(
                f"VEP assembly '{self.assembly}' does not match reference build "
                f"'{self.input.reference.build}'."
            )
        if not isinstance(self.threads, int) or isinstance(self.threads, bool) or self.threads <= 0:
            raise InputValidationError("VEP threads must be a positive integer.")
        if not isinstance(self.overwrite, bool):
            raise InputValidationError("VEP overwrite must be boolean.")
        object.__setattr__(self, "cache_directory", cache)
        object.__setattr__(self, "output_tsv", Path(self.output_tsv).expanduser())

    def to_dict(self) -> dict[str, object]:
        return {
            "input": self.input.to_dict(),
            "cache_directory": str(self.cache_directory),
            "cache_version": self.cache_version,
            "species": self.species,
            "assembly": self.assembly,
            "reference_fasta": str(self.input.reference.fasta),
            "output_tsv": str(self.output_tsv),
            "threads": self.threads,
            "overwrite": self.overwrite,
            "offline": True,
            "cache_download_performed": False,
        }


class VepWrapper:
    """Run VEP in explicit offline cache mode through CommandRunner."""

    def __init__(self, executable: str = "vep", *, runner: CommandRunner | None = None) -> None:
        if not isinstance(executable, str) or not executable.strip():
            raise InputValidationError("VEP executable must be non-empty.")
        self.executable = executable
        self.runner = runner or CommandRunner()

    def build_command(self, request: VepRequest) -> tuple[str, ...]:
        return (
            self.executable,
            "--input_file", str(request.input.source_vcf.absolute()),
            "--output_file", str(request.output_tsv.absolute()),
            "--format", "vcf",
            "--tab",
            "--offline",
            "--cache",
            "--dir_cache", str(request.cache_directory.absolute()),
            "--cache_version", request.cache_version,
            "--species", request.species,
            "--assembly", request.assembly,
            "--fasta", str(request.input.reference.fasta.absolute()),
            "--fork", str(request.threads),
        )

    def detect_version(self) -> str:
        self.runner.require_executable(self.executable)
        result = self.runner.run((self.executable, "--help"))
        output = "\n".join(value for value in (result.stdout, result.stderr) if value)
        match = _VERSION.search(output)
        if match is None:
            raise ToolVersionError("Unable to parse Ensembl VEP version from --help output.")
        return match.group(1)

    def _protect_paths(
        self, request: VepRequest, stdout_path: Path | None,
        stderr_path: Path | None, *, staged_output: Path | None = None,
    ) -> None:
        """Guard the actual logs, including CommandRunner's implicit stdout."""
        effective_stdout = stdout_path
        if effective_stdout is None and stderr_path is not None:
            if not isinstance(stderr_path, (str, Path)):
                raise InputValidationError("VEP stderr path must be a string or Path.")
            effective_stdout = _default_stdout_log_path(Path(stderr_path).expanduser())
        protect_annotation_outputs(
            outputs={"annotation_tsv": request.output_tsv, "stdout": effective_stdout,
                     "stderr": stderr_path, "staged_tsv": staged_output},
            immutable_files=(request.input.source_vcf, Path(f"{request.input.source_vcf}.tbi"),
                             Path(f"{request.input.source_vcf}.csi"),
                             request.input.reference.fasta, request.input.reference.fai),
            resource_directories=(request.cache_directory,),
            executable=self.executable,
        )

    @staticmethod
    def _check_destination(request: VepRequest) -> None:
        if request.output_tsv.exists() and not request.overwrite:
            raise OutputValidationError(f"VEP output already exists: '{request.output_tsv}'.")
        if request.output_tsv.exists() and request.output_tsv.is_dir():
            raise OutputValidationError(f"VEP output is a directory: '{request.output_tsv}'.")

    def run(
        self,
        request: VepRequest,
        *,
        dry_run: bool = False,
        stdout_path: Path | None = None,
        stderr_path: Path | None = None,
    ) -> AnnotationResult:
        self._protect_paths(request, stdout_path, stderr_path)
        command = self.build_command(request)
        if dry_run:
            self.runner.run(command, dry_run=True)
            return AnnotationResult(
                request.input, AnnotationSource.VEP,
                AnnotationRunStatus.PLANNED, command,
            )
        self._check_destination(request)
        version = self.detect_version()
        self._protect_paths(request, stdout_path, stderr_path)
        self._check_destination(request)
        request.output_tsv.parent.mkdir(parents=True, exist_ok=True)
        # VEP can create statistics sidecars beside its native output. A fresh
        # owned directory prevents them following aliases at the final prefix.
        # Preserve staging/sidecars on failure for diagnosis; never remove an
        # existing user result before a replacement has passed validation.
        stage = Path(tempfile.mkdtemp(prefix=".vep-", dir=request.output_tsv.parent))
        staged_output = stage / request.output_tsv.name
        command = self.build_command(replace(request, output_tsv=staged_output))
        self._protect_paths(request, stdout_path, stderr_path, staged_output=staged_output)
        result = self.runner.run(
            command,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
        )
        self._protect_paths(request, stdout_path, stderr_path, staged_output=staged_output)
        _validate_vep_tsv(staged_output)
        self._check_destination(request)
        try:
            if request.overwrite:
                os.replace(staged_output, request.output_tsv)
            else:
                # Same-filesystem exclusive publication: a result created while
                # VEP was running must not be replaced by a check-then-rename.
                # Keep the staged hardlink as native execution provenance.
                os.link(staged_output, request.output_tsv)
        except OSError as error:
            raise OutputValidationError(
                f"VEP validated TSV could not be atomically published to '{request.output_tsv}': {error}"
            ) from error
        validate_output_file(request.output_tsv)
        database = AnnotationDatabase(
            name=f"ensembl-vep-cache:{request.species}",
            version=request.cache_version,
            path=request.cache_directory,
            reference_build=request.assembly,
        )
        artifact = AnnotationArtifact(
            request.input,
            AnnotationSource.VEP,
            request.output_tsv,
            "tsv",
            version,
            (database,),
            command,
        )
        return AnnotationResult(
            request.input,
            AnnotationSource.VEP,
            AnnotationRunStatus.COMPLETED,
            command,
            version,
            result.duration_seconds,
            artifact,
        )


__all__ = ["VepRequest", "VepWrapper"]
