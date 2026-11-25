"""Lossless Phase 12 SV and TR cohort table builders."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import sqlite3
import stat
import uuid
from pathlib import Path
from typing import Iterator, TextIO

from hifivar.cohort import (
    CohortDefinition,
    CohortSampleInput,
    CohortTrack,
    CohortTrackResult,
    SampleCallState,
    validate_track_inputs,
)
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.logging_utils import get_logger


_LOGGER = get_logger(__name__)


def _file_snapshot(path: Path, identity: tuple[int, int]) -> tuple | None:
    """Read a regular owned file conservatively; never follow a replacement link."""
    try:
        before = path.lstat()
        if not stat.S_ISREG(before.st_mode) or (before.st_dev, before.st_ino) != identity:
            return None
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            opened = os.fstat(handle.fileno())
            if (opened.st_dev, opened.st_ino) != identity:
                return None
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        after = path.lstat()
        signature = lambda value: (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_mode)
        if signature(before) != signature(after):
            return None
        return (*signature(after), digest.hexdigest())
    except OSError:
        return None


class _CohortOutputPair:
    """Validate in owned staging, then exclusively publish each final file.

    Captured publication errors conservatively roll back unchanged owned files.
    Two separate paths are not physically atomic as a pair: process death, power
    loss or arbitrary concurrent replacements remain outside this guarantee.
    Failed staging is retained and a retry receives new, independent paths.
    """

    def __init__(self, targets: tuple[Path, Path], sources: tuple[Path, ...]) -> None:
        self.targets = targets
        self.sources = sources
        self.paths: list[Path] = []
        self.directories: dict[Path, tuple[int, int]] = {}
        self.owned: dict[Path, tuple[int, int]] = {}
        self.sealed: dict[Path, tuple] = {}

    def _preflight(self) -> None:
        _reserve_outputs(*self.targets)
        resolved = [path.resolve() for path in self.targets]
        if resolved[0].is_relative_to(resolved[1]) or resolved[1].is_relative_to(resolved[0]):
            raise InputValidationError("Cohort output paths must be distinct and not ancestors.")
        for target, real_target in zip(self.targets, resolved):
            if target.is_symlink():
                raise OutputValidationError(f"Refusing a linked cohort output: '{target}'.")
            for parent in target.parents:
                if parent.exists() and not parent.is_dir():
                    raise OutputValidationError(f"Cohort output parent is not a directory: '{parent}'.")
            for source in self.sources:
                real_source = source.resolve()
                if real_target.is_relative_to(real_source) or real_source.is_relative_to(real_target):
                    raise InputValidationError(f"Cohort output conflicts with protected input: '{target}'.")

    def _stage(self, parent: Path) -> Path:
        real_parent = parent.resolve()
        for directory in self.directories:
            if directory.parent.resolve() == real_parent:
                return directory
        parent.mkdir(parents=True, exist_ok=True)
        directory = parent / f".hifivar-cohort-{uuid.uuid4().hex}"
        directory.mkdir(mode=0o700)
        metadata = directory.lstat()
        self.directories[directory] = (metadata.st_dev, metadata.st_ino)
        return directory

    def __enter__(self) -> "_CohortOutputPair":
        self._preflight()
        try:
            self.paths = [self._stage(path.parent) / f"output-{index}.tsv"
                          for index, path in enumerate(self.targets)]
        except Exception as error:
            self._retain(error)
            raise
        return self

    def open_text(self, index: int) -> TextIO:
        path = self.paths[index]
        handle = _new_text(path)
        try:
            metadata = os.fstat(handle.fileno())
        except Exception:
            handle.close()
            raise
        self.owned[path] = (metadata.st_dev, metadata.st_ino)
        return handle

    def database(self, requested: Path) -> Path:
        """Treat the requested scratch path as a parent hint, never overwrite it."""
        parent = requested.parent.resolve()
        if any(parent.is_relative_to(path.resolve()) for path in (*self.targets, *self.sources)):
            raise InputValidationError("Cohort scratch parent conflicts with a final output or protected input.")
        path = self._stage(requested.parent) / "scratch.sqlite"
        with path.open("xb") as handle:
            metadata = os.fstat(handle.fileno())
            self.owned[path] = (metadata.st_dev, metadata.st_ino)
        _LOGGER.info("Cohort scratch database: %s", path)
        return path

    def seal(self, path: Path) -> tuple:
        snapshot = _file_snapshot(path, self.owned[path])
        if snapshot is None:
            raise OutputValidationError(f"Cohort staging identity changed: '{path}'.")
        self.sealed[path] = snapshot
        return snapshot

    @staticmethod
    def _remove_unchanged(path: Path, snapshot: tuple) -> bool:
        if _file_snapshot(path, snapshot[:2]) != snapshot:
            return False
        path.unlink()
        return True

    def publish(self) -> None:
        snapshots = [self.seal(path) for path in self.paths]
        attempted: list[tuple[Path, tuple]] = []
        try:
            for staged, target, snapshot in zip(self.paths, self.targets, snapshots):
                # Staging is a sibling, so each link stays on its destination's
                # filesystem. link() is exclusive; never replace a late file.
                # Register before the call: it can create the link and then be
                # interrupted before returning. Rollback still requires the
                # original sealed inode, metadata and full content snapshot.
                attempted.append((target, snapshot))
                os.link(staged, target)
            for target, snapshot in attempted:
                if _file_snapshot(target, snapshot[:2]) != snapshot:
                    raise OutputValidationError(
                        f"Published cohort output differs from validated staging: '{target}'."
                    )
        except BaseException as error:
            preserved = []
            for target, snapshot in reversed(attempted):
                if not os.path.lexists(target):
                    continue
                try:
                    removed = self._remove_unchanged(target, snapshot)
                except OSError:
                    removed = False
                if not removed:
                    preserved.append(str(target))
            detail = f" Changed/unremovable outputs preserved for inspection: {preserved}." if preserved else ""
            if isinstance(error, (OSError, OutputValidationError)):
                raise OutputValidationError(f"Cohort output publication failed: {error}.{detail}") from error
            if preserved:
                _LOGGER.error("COHORT_ROLLBACK_PRESERVED_CHANGED_OUTPUTS: %s", preserved)
                add_note = getattr(error, "add_note", None)
                if add_note is not None:
                    add_note(detail.strip())
            # A captured interrupt receives best-effort rollback, but remains
            # the original interrupt. SIGKILL/power loss cannot enter here.
            raise

    def _retain(self, error: BaseException) -> None:
        if self.directories:
            message = "COHORT_STAGING_RETAINED: directories=" + repr([str(path) for path in self.directories])
            message += "; files=" + repr([str(path) for path in self.owned])
            _LOGGER.error(message)
            add_note = getattr(error, "add_note", None)
            if add_note is not None:
                add_note(message)

    def __exit__(self, error_type, error, traceback) -> bool:
        if error is not None:
            self._retain(error)
            return False
        for path, snapshot in self.sealed.items():
            try:
                removed = self._remove_unchanged(path, snapshot)
            except OSError:
                removed = False
            if not removed:
                _LOGGER.warning("Cohort staging file changed or could not be cleaned; preserved: %s", path)
        for directory, identity in reversed(tuple(self.directories.items())):
            try:
                metadata = directory.lstat()
                if stat.S_ISDIR(metadata.st_mode) and (metadata.st_dev, metadata.st_ino) == identity:
                    directory.rmdir()  # Empty owned directories only; never recurse.
            except OSError:
                _LOGGER.warning("Cohort staging directory retained: %s", directory)
        return False


def _protected_inputs(cohort: CohortDefinition, inputs: tuple[CohortSampleInput, ...]) -> tuple[Path, ...]:
    return (cohort.reference.fasta, cohort.reference.fai, *(
        path for item in inputs for path in (item.source_path, item.index_path) if path is not None
    ))


def build_sv_cohort_tables(
    cohort: CohortDefinition,
    inputs: tuple[CohortSampleInput, ...],
    *,
    site_table: Path,
    sample_matrix: Path,
) -> CohortTrackResult:
    """Write native-representation SV rows without inventing cross-sample clustering."""
    validate_track_inputs(cohort, inputs)
    site_table, sample_matrix = Path(site_table), Path(sample_matrix)
    callable_count = sum(item.callable for item in inputs)
    disambiguated_count = 0
    with _CohortOutputPair((site_table, sample_matrix), _protected_inputs(cohort, inputs)) as pending:
        database = pending.database(site_table.parent / "sv_source_keys.sqlite")
        connection: sqlite3.Connection | None = None
        try:
            connection = sqlite3.connect(database)
            variant_count = _stage_sv_source_records(connection, cohort, inputs)
            with pending.open_text(0) as sites, pending.open_text(1) as matrix:
                sites.write("cohort_site_id\tsource_sample\tsource_variant_id\tcontig\tstart\tend\tsvtype\tref\talt\tsample_support_count\tcallable_sample_count\tsample_support_fraction\tsource_vcf\n")
                matrix.write("cohort_site_id\tsample\tstate\tgenotype\tsource_variant_id\tsource_vcf\n")
                records = connection.execute(
                    "SELECT sample, record_ordinal, record_id, legacy_key, fields, source "
                    "FROM sv_records ORDER BY row_number"
                )
                for sample_id, ordinal, record_id, legacy_key, encoded, source_path in records:
                    fields = json.loads(encoded)
                    count = connection.execute(
                        "SELECT record_count FROM sv_key_counts WHERE legacy_key = ?", (legacy_key,)
                    ).fetchone()[0]
                    if count == 1:
                        cohort_site_id = legacy_key
                    else:
                        # Sample IDs cannot start with '@' or contain ':'. This
                        # reserved namespace is disjoint from EVERY legacy key,
                        # even when a caller ID contains this exact key text.
                        cohort_site_id = f"@hifivar-sv-record-v2:{sample_id}:{ordinal}"
                        disambiguated_count += 1
                    connection.execute("INSERT INTO sv_assigned_keys VALUES (?)", (cohort_site_id,))
                    info = _parse_info(fields[7])
                    svtype = (info.get("SVTYPE") or "UNRESOLVED").upper()
                    end = info.get("END") or ""
                    fraction = 1 / callable_count if callable_count else 0.0
                    sites.write("\t".join((cohort_site_id, sample_id, record_id, fields[0], fields[1], end, svtype, fields[3], fields[4], "1", str(callable_count), f"{fraction:.12g}", source_path)) + "\n")
                    genotype = _genotype(fields)
                    for sample in inputs:
                        if sample.sample_id == sample_id:
                            state = SampleCallState.CALLED
                            cell_gt = genotype
                        elif sample.callable:
                            state = SampleCallState.NOT_OBSERVED
                            cell_gt = ""
                        else:
                            state = sample.state
                            cell_gt = ""
                        matrix.write("\t".join((cohort_site_id, sample.sample_id, state.value, cell_gt, record_id, source_path)) + "\n")
            connection.commit()
        except sqlite3.Error as error:
            raise OutputValidationError(f"Unable to index unique SV cohort source keys: {error}") from error
        finally:
            if connection is not None:
                connection.close()
        pending.seal(database)
        pending.publish()
    return CohortTrackResult(
        CohortTrack.SV, True, SampleCallState.CALLED, "native_harmonized_sv", None,
        (site_table, sample_matrix), (), inputs,
        {"sample_count": len(inputs), "callable_sample_count": callable_count,
         "native_site_count": variant_count, "cross_sample_clustering": False,
         "site_key_policy": "legacy_unique_else_reserved_record_ordinal_v2",
         "disambiguated_record_count": disambiguated_count,
         "site_key_uniqueness_verified": True},
        "Lossless source-native representation; no new cross-sample clustering was applied.",
    )


def _stage_sv_source_records(
    connection: sqlite3.Connection,
    cohort: CohortDefinition,
    inputs: tuple[CohortSampleInput, ...],
) -> int:
    """Read each VCF once; keep records and complete key counts on disk.

    Counts are complete before any key is allocated, so all members of a
    collision group are disambiguated. No in-memory whole-VCF index or second
    source-file read is needed. An ordinal identifies a record in this exact
    source order, not an allele or a shared cross-sample locus.
    """
    connection.execute(
        "CREATE TABLE sv_records (row_number INTEGER PRIMARY KEY, sample TEXT, "
        "record_ordinal INTEGER, record_id TEXT, legacy_key TEXT, fields TEXT, source TEXT)"
    )
    connection.execute("CREATE TABLE sv_key_counts (legacy_key TEXT PRIMARY KEY, record_count INTEGER)")
    connection.execute("CREATE TABLE sv_assigned_keys (cohort_site_id TEXT PRIMARY KEY)")
    count = 0
    for source in inputs:
        if not source.callable:
            continue
        assert source.source_path is not None
        header_sample, contigs = _vcf_header(source.source_path)
        if header_sample != source.sample_id:
            raise OutputValidationError(f"SV VCF sample mismatch for '{source.sample_id}': observed '{header_sample}'.")
        cohort.reference.validate_contigs(contigs)
        for ordinal, fields in enumerate(_vcf_records(source.source_path), start=1):
            cohort.reference.validate_contigs((fields[0],))
            record_id = fields[2] if fields[2] not in {"", "."} else f"{fields[0]}:{fields[1]}:{fields[3]}:{fields[4]}"
            legacy_key = f"{source.sample_id}:{record_id}"
            count += 1
            connection.execute("INSERT INTO sv_records VALUES (?, ?, ?, ?, ?, ?, ?)",
                               (count, source.sample_id, ordinal, record_id, legacy_key,
                                json.dumps(fields, ensure_ascii=False), str(source.source_path)))
            connection.execute(
                "INSERT INTO sv_key_counts VALUES (?, 1) ON CONFLICT(legacy_key) "
                "DO UPDATE SET record_count = record_count + 1", (legacy_key,)
            )
    connection.commit()
    return count


def build_tr_cohort_tables(
    cohort: CohortDefinition,
    inputs: tuple[CohortSampleInput, ...],
    *,
    locus_table: Path,
    sample_matrix: Path,
    scratch_database: Path,
) -> CohortTrackResult:
    """Stream TRGT records through an on-disk index and emit a locus matrix."""
    validate_track_inputs(cohort, inputs)
    catalog_ids = {item.catalog_id for item in inputs if item.callable}
    if None in catalog_ids or len(catalog_ids) != 1:
        raise InputValidationError("Callable TR cohort inputs require one identical explicit catalog_id.")
    locus_table, sample_matrix, scratch_database = Path(locus_table), Path(sample_matrix), Path(scratch_database)
    with _CohortOutputPair((locus_table, sample_matrix), _protected_inputs(cohort, inputs)) as pending:
        database = pending.database(scratch_database)
        connection: sqlite3.Connection | None = None
        primary_error: BaseException | None = None
        try:
            connection = sqlite3.connect(database)
            connection.execute("CREATE TABLE calls (trid TEXT, sample TEXT, contig TEXT, pos INTEGER, end INTEGER, motifs TEXT, genotype TEXT, source TEXT, PRIMARY KEY (trid, sample))")
            connection.execute("CREATE INDEX calls_order ON calls(contig, pos, trid)")
            for source in inputs:
                if not source.callable:
                    continue
                assert source.source_path is not None
                header_sample, contigs = _vcf_header(source.source_path)
                if header_sample != source.sample_id:
                    raise OutputValidationError(f"TR VCF sample mismatch for '{source.sample_id}': observed '{header_sample}'.")
                cohort.reference.validate_contigs(contigs)
                for fields in _vcf_records(source.source_path):
                    cohort.reference.validate_contigs((fields[0],))
                    info = _parse_info(fields[7])
                    trid, position, end = _tr_locus_identity(fields, source.source_path)
                    try:
                        connection.execute(
                            "INSERT INTO calls VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                            (trid, source.sample_id, fields[0], position, end, info.get("MOTIFS") or "", _genotype(fields), str(source.source_path)),
                        )
                    except sqlite3.IntegrityError as error:
                        raise OutputValidationError(f"Duplicate TR locus '{trid}' for sample '{source.sample_id}'.") from error
            connection.commit()
            callable_count = sum(item.callable for item in inputs)
            loci = connection.execute("SELECT trid, MIN(contig), MIN(pos), MIN(end), MIN(motifs), COUNT(DISTINCT sample), COUNT(DISTINCT contig || ':' || pos || ':' || end || ':' || motifs) FROM calls GROUP BY trid ORDER BY MIN(contig), MIN(pos), trid")
            locus_count = 0
            with pending.open_text(0) as sites, pending.open_text(1) as matrix:
                sites.write("trid\tcontig\tstart\tend\tmotifs\tsample_support_count\tcallable_sample_count\tsample_support_fraction\tcatalog_id\n")
                matrix.write("trid\tsample\tstate\tgenotype\tsource_vcf\n")
                for trid, contig, position, end, motifs, support, representations in loci:
                    locus_count += 1
                    if representations != 1:
                        raise OutputValidationError(f"TR locus '{trid}' has inconsistent coordinates or motifs across samples.")
                    sites.write("\t".join((trid, contig, str(position), str(end), motifs, str(support), str(callable_count), f"{support / callable_count:.12g}" if callable_count else "", next(iter(catalog_ids)))) + "\n")
                    observed = {row[0]: (row[1], row[2]) for row in connection.execute("SELECT sample, genotype, source FROM calls WHERE trid = ?", (trid,))}
                    for sample in inputs:
                        if sample.sample_id in observed:
                            genotype, source_path = observed[sample.sample_id]
                            state = SampleCallState.CALLED
                        elif sample.callable:
                            genotype, source_path, state = "", "", SampleCallState.NOT_OBSERVED
                        else:
                            genotype, source_path, state = "", "", sample.state
                        matrix.write("\t".join((trid, sample.sample_id, state.value, genotype, source_path)) + "\n")
            result = CohortTrackResult(
                CohortTrack.TR, True, SampleCallState.CALLED, "trgt_cohort_matrix", None,
                (locus_table, sample_matrix), (), inputs,
                {"sample_count": len(inputs), "callable_sample_count": callable_count, "locus_count": locus_count, "catalog_id": next(iter(catalog_ids))},
            )
        except sqlite3.Error as error:
            primary_error = OutputValidationError(f"Unable to build TR cohort index: {error}")
            raise primary_error from error
        except BaseException as error:
            primary_error = error
            raise
        finally:
            if connection is not None:
                try:
                    connection.close()
                except sqlite3.Error as error:
                    if primary_error is None:
                        raise OutputValidationError(f"Unable to close TR cohort index: {error}") from error
                    detail = f"TR cohort index also failed to close: {error}"
                    _LOGGER.error(detail)
                    add_note = getattr(primary_error, "add_note", None)
                    if add_note is not None:
                        add_note(detail)
        pending.seal(database)
        pending.publish()
    return result


def write_track_result(path: Path, result: CohortTrackResult) -> None:
    path = Path(path)
    _reserve_outputs(path)
    with _new_text(path) as handle:
        json.dump(result.to_dict(), handle, indent=2, sort_keys=True)
        handle.write("\n")


def _vcf_header(path: Path) -> tuple[str | None, tuple[str, ...]]:
    samples: tuple[str, ...] | None = None
    contigs: list[str] = []
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("##contig=<ID="):
                contigs.append(line.split("ID=", 1)[1].split(",", 1)[0].split(">", 1)[0])
            elif line.startswith("#CHROM\t"):
                samples = tuple(line.rstrip("\r\n").split("\t")[9:])
                break
    if samples is None or len(samples) != 1:
        raise OutputValidationError(f"Expected exactly one sample in VCF '{path}', observed {samples!r}.")
    return samples[0], tuple(contigs)


def _vcf_records(path: Path) -> Iterator[list[str]]:
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 10:
                raise OutputValidationError(f"Malformed single-sample VCF record in '{path}'.")
            yield fields


def _genotype(fields: list[str]) -> str:
    keys = fields[8].split(":")
    values = fields[9].split(":")
    return values[keys.index("GT")] if "GT" in keys and keys.index("GT") < len(values) else ""


def _parse_info(text: str) -> dict[str, str | None]:
    result: dict[str, str | None] = {}
    if text not in {"", "."}:
        for item in text.split(";"):
            key, separator, value = item.partition("=")
            result[key] = value if separator else None
    return result


def _safe_int(value: str | None, default: int) -> int:
    try:
        return int(value) if value not in {None, "", "."} else default
    except ValueError:
        return default


def _tr_locus_identity(fields: list[str], source: Path) -> tuple[str, int, int]:
    """Require explicit TRGT identity; never impute coordinates from POS/REF.

    VCF ID may be '.', while INFO/TRID must be one real catalog identifier.
    Coordinates retain the existing VCF/padding convention and the on-disk
    cohort index's positive signed-64-bit integer boundary.
    """
    values: dict[str, list[str | None]] = {"TRID": [], "END": []}
    for item in fields[7].split(";"):
        key, separator, value = item.partition("=")
        if key in values:
            values[key].append(value if separator else None)
    for key, observed in values.items():
        if len(observed) != 1:
            raise OutputValidationError(f"TR VCF requires exactly one INFO/{key} in '{source}'.")
    trid = values["TRID"][0]
    if trid is None or not trid or trid == "." or "," in trid or any(char.isspace() for char in trid):
        raise OutputValidationError(f"TR VCF record lacks a valid INFO/TRID in '{source}'.")
    position = _tr_coordinate(fields[1], "POS", source)
    end = _tr_coordinate(values["END"][0], "END", source)
    if end < position:
        raise OutputValidationError(f"TR VCF END precedes POS in '{source}'.")
    return trid, position, end


def _tr_coordinate(value: str | None, field: str, source: Path) -> int:
    """Parse an explicit coordinate without accepting floats or missing marks."""
    if value is None or not value.isascii() or not value.isdecimal():
        raise OutputValidationError(f"TR VCF {field} must be an explicit positive integer in '{source}'.")
    significant = value.lstrip("0")
    if not significant or len(significant) > 19:
        raise OutputValidationError(f"TR VCF {field} is outside the cohort integer range in '{source}'.")
    number = int(significant)
    if number > 9223372036854775807:
        raise OutputValidationError(f"TR VCF {field} is outside the cohort integer range in '{source}'.")
    return number


def _open_text(path: Path) -> TextIO:
    try:
        return gzip.open(path, "rt", encoding="utf-8", newline="") if str(path).lower().endswith(".gz") else path.open("r", encoding="utf-8", newline="")
    except OSError as error:
        raise InputValidationError(f"Unable to open cohort VCF '{path}': {error}") from error


def _reserve_outputs(*paths: Path) -> None:
    for path in paths:
        if Path(path).exists():
            raise OutputValidationError(f"Refusing to overwrite cohort artifact: '{path}'.")


def _new_text(path: Path) -> TextIO:
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.open("x", encoding="utf-8", newline="\n")


__all__ = ["build_sv_cohort_tables", "build_tr_cohort_tables", "write_track_result"]
