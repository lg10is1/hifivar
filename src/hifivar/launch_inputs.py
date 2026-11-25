"""Read-only samtools preflight for the indexed-BAM batch launcher."""
from __future__ import annotations

import re
import os
from pathlib import Path

from hifivar.command import CommandRunner
from hifivar.exceptions import InputValidationError, ToolVersionError
from hifivar.reference import ReferenceGenome


def unique_bam_index(bam: Path) -> Path:
    """Require one unambiguous, nonempty conventional BAI/CSI; never create it."""
    candidates = {Path(str(bam) + suffix) for suffix in (".bai", ".csi")}
    candidates.update(bam.with_suffix(suffix) for suffix in (".bai", ".csi"))
    found = sorted(path for path in candidates if path.exists())
    if len(found) != 1 or not found[0].is_file() or found[0].stat().st_size == 0:
        raise InputValidationError(f"BAM requires exactly one nonempty BAI/CSI: {bam}")
    return found[0].resolve()


def parse_header_records(header: str) -> dict[str, list[dict[str, str]]]:
    """Parse identity fields without accepting duplicate keys."""
    records: dict[str, list[dict[str, str]]] = {"@HD": [], "@SQ": [], "@RG": []}
    for line in header.splitlines():
        parts = line.split("\t")
        if parts[0] not in records:
            continue
        fields: dict[str, str] = {}
        for field in parts[1:]:
            if ":" not in field:
                raise InputValidationError("Malformed BAM header field.")
            key, value = field.split(":", 1)
            if key in fields:
                raise InputValidationError(f"Duplicate BAM header field: {key}")
            fields[key] = value
        records[parts[0]].append(fields)
    return records


def validate_read_groups(records: dict, sample_id: str, *, pacbio: bool = False) -> None:
    """Require unique read groups and every SM, including on a raw uBAM."""
    ids: set[str] = set()
    if not records["@RG"]:
        raise InputValidationError(f"BAM has no @RG/SM for declared sample {sample_id}.")
    for rg in records["@RG"]:
        group = rg.get("ID", "")
        if not group.strip() or group in ids or rg.get("SM") != sample_id:
            raise InputValidationError(f"BAM @RG IDs must be unique and every SM must equal {sample_id}.")
        if pacbio and rg.get("PL", "").upper() != "PACBIO":
            raise InputValidationError("uBAM read groups must declare PL:PACBIO.")
        ids.add(group)


def validate_bam_header(header: str, sample_id: str, reference: ReferenceGenome) -> None:
    """Check declared sort, exact reference dictionary, and every read group's SM."""
    records = parse_header_records(header)
    if len(records["@HD"]) != 1 or records["@HD"][0].get("SO") != "coordinate":
        raise InputValidationError("Aligned BAM must declare SO:coordinate.")
    observed: dict[str, int] = {}
    for sq in records["@SQ"]:
        name = sq.get("SN", "")
        try:
            length = int(sq.get("LN", ""))
        except ValueError as error:
            raise InputValidationError("Invalid BAM @SQ length.") from error
        if not name or name in observed or length <= 0:
            raise InputValidationError("Empty/duplicate BAM @SQ or invalid length.")
        observed[name] = length
    if observed != {contig.name: contig.length for contig in reference.contigs}:
        raise InputValidationError("REFERENCE_CONTIG_MISMATCH: BAM names/lengths must exactly match FASTA FAI.")
    validate_read_groups(records, sample_id)


class IndexedBamPreflight:
    """Validate without rewriting inputs; all external execution uses CommandRunner."""

    def __init__(self, runner: CommandRunner, executable: str) -> None:
        self.runner = runner
        self.executable = executable
        # Read only the named non-secret variables; never enumerate/log the
        # parent environment or inherit unrelated credentials into probes.
        allowed = ("PATH", "HOME", "USER", "LOGNAME", "LANG", "LC_ALL", "TZ", "CONDA_PREFIX", "CONDA_EXE", "SystemRoot", "SYSTEMROOT")
        self.environment = {key: os.environ[key] for key in allowed if key in os.environ}

    def run(self, bam: Path, index: Path, sample_id: str,
            reference: ReferenceGenome, log_dir: Path | None = None) -> dict[str, object]:
        """Check header, EOF and indexed access; this is not a full BAM record audit."""
        self.runner.require_executable(self.executable)
        steps = []

        def execute(label: str, args: list[str]) -> str:
            stdout = log_dir / f"{sample_id}.{label}.stdout.log" if log_dir else None
            stderr = log_dir / f"{sample_id}.{label}.stderr.log" if log_dir else None
            result = self.runner.run(
                [self.executable, *args], stdout_path=stdout, stderr_path=stderr,
                timeout=60, env=self.environment, inherit_env=False,
            )
            steps.append({"command": list(result.args), "returncode": result.returncode,
                          "duration_seconds": result.duration_seconds,
                          "stdout_path": str(stdout) if stdout else None,
                          "stderr_path": str(stderr) if stderr else None})
            return stdout.read_text(encoding="utf-8", errors="strict") if stdout else (result.stdout or "")

        version = execute("version", ["--version"])
        match = re.search(r"(?m)^samtools\s+(\S+)", version)
        if not match:
            raise ToolVersionError("Unable to identify samtools version for BAM preflight.")
        execute("quickcheck", ["quickcheck", "-v", str(bam)])
        header = execute("header", ["view", "--no-PG", "-H", str(bam)])
        validate_bam_header(header, sample_id, reference)
        region = f"{reference.contigs[0].name}:1-1"
        count = execute("index_probe", ["view", "-c", "-X", str(bam), str(index), region])
        if not count.strip().isdigit():
            raise InputValidationError("BAM indexed-access probe did not return a count.")
        return {"sample_id": sample_id, "samtools_version": match.group(1),
                "header_identity": "PASS", "quickcheck": "PASS", "indexed_access": "PASS",
                "scope": "header_EOF_and_one_indexed_region_not_full_record_scan", "steps": steps}

    def read_header(self, bam: Path) -> str:
        """Read only header text without creating logs, outputs or indexes."""
        self.runner.require_executable(self.executable)
        result = self.runner.run([self.executable, "view", "--no-PG", "-H", str(bam)], timeout=60,
                                 env=self.environment, inherit_env=False)
        return result.stdout or ""
