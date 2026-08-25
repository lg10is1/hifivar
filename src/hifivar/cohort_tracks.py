"""Lossless Phase 12 SV and TR cohort table builders."""

from __future__ import annotations

import gzip
import json
import sqlite3
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
    _reserve_outputs(site_table, sample_matrix)
    callable_count = sum(item.callable for item in inputs)
    variant_count = 0
    with _new_text(site_table) as sites, _new_text(sample_matrix) as matrix:
        sites.write("cohort_site_id\tsource_sample\tsource_variant_id\tcontig\tstart\tend\tsvtype\tref\talt\tsample_support_count\tcallable_sample_count\tsample_support_fraction\tsource_vcf\n")
        matrix.write("cohort_site_id\tsample\tstate\tgenotype\tsource_variant_id\tsource_vcf\n")
        for source in inputs:
            if not source.callable:
                continue
            assert source.source_path is not None
            header_sample, contigs = _vcf_header(source.source_path)
            if header_sample != source.sample_id:
                raise OutputValidationError(f"SV VCF sample mismatch for '{source.sample_id}': observed '{header_sample}'.")
            cohort.reference.validate_contigs(contigs)
            for fields in _vcf_records(source.source_path):
                cohort.reference.validate_contigs((fields[0],))
                variant_count += 1
                info = _parse_info(fields[7])
                record_id = fields[2] if fields[2] not in {"", "."} else f"{fields[0]}:{fields[1]}:{fields[3]}:{fields[4]}"
                cohort_site_id = f"{source.sample_id}:{record_id}"
                svtype = (info.get("SVTYPE") or "UNRESOLVED").upper()
                end = info.get("END") or ""
                fraction = 1 / callable_count if callable_count else 0.0
                sites.write("\t".join((cohort_site_id, source.sample_id, record_id, fields[0], fields[1], end, svtype, fields[3], fields[4], "1", str(callable_count), f"{fraction:.12g}", str(source.source_path))) + "\n")
                genotype = _genotype(fields)
                for sample in inputs:
                    if sample.sample_id == source.sample_id:
                        state = SampleCallState.CALLED
                        cell_gt = genotype
                    elif sample.callable:
                        state = SampleCallState.NOT_OBSERVED
                        cell_gt = ""
                    else:
                        state = sample.state
                        cell_gt = ""
                    matrix.write("\t".join((cohort_site_id, sample.sample_id, state.value, cell_gt, record_id, str(source.source_path))) + "\n")
    return CohortTrackResult(
        CohortTrack.SV, True, SampleCallState.CALLED, "native_harmonized_sv", None,
        (site_table, sample_matrix), (), inputs,
        {"sample_count": len(inputs), "callable_sample_count": callable_count, "native_site_count": variant_count, "cross_sample_clustering": False},
        "Lossless source-native representation; no new cross-sample clustering was applied.",
    )


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
    _reserve_outputs(locus_table, sample_matrix, scratch_database)
    scratch_database.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(scratch_database)
    completed = False
    try:
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
                trid = info.get("TRID")
                if not trid:
                    raise OutputValidationError(f"TR VCF record lacks TRID in '{source.source_path}'.")
                end = _safe_int(info.get("END"), _safe_int(fields[1], 0))
                try:
                    connection.execute(
                        "INSERT INTO calls VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                        (trid, source.sample_id, fields[0], int(fields[1]), end, info.get("MOTIFS") or "", _genotype(fields), str(source.source_path)),
                    )
                except sqlite3.IntegrityError as error:
                    raise OutputValidationError(f"Duplicate TR locus '{trid}' for sample '{source.sample_id}'.") from error
        connection.commit()
        callable_count = sum(item.callable for item in inputs)
        loci = connection.execute("SELECT trid, MIN(contig), MIN(pos), MIN(end), MIN(motifs), COUNT(DISTINCT sample), COUNT(DISTINCT contig || ':' || pos || ':' || end || ':' || motifs) FROM calls GROUP BY trid ORDER BY MIN(contig), MIN(pos), trid")
        locus_count = 0
        with _new_text(locus_table) as sites, _new_text(sample_matrix) as matrix:
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
        completed = True
        return result
    finally:
        connection.close()
        if completed and scratch_database.exists():
            scratch_database.unlink()


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
            if len(fields) < 10:
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
