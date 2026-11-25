"""Native PacBio uBAM alignment, without an implicit FASTQ conversion."""
from __future__ import annotations

from pathlib import Path

from hifivar import validation
from hifivar.command import CommandRunner
from hifivar.exceptions import InputValidationError, OutputValidationError
from hifivar.launch_inputs import parse_header_records, validate_read_groups
from hifivar.pbmm2 import Pbmm2Options, Pbmm2Wrapper
from hifivar.reference import ReferenceGenome


def validate_ubam_header(header: str, sample_id: str) -> None:
    """Require no reference dictionary and consistent PacBio sample metadata."""
    records = parse_header_records(header)
    if records["@SQ"] or len(records["@HD"]) != 1 or records["@HD"][0].get("SO") == "coordinate":
        raise InputValidationError("Native uBAM requires one @HD, no @SQ and non-coordinate sort order.")
    validate_read_groups(records, sample_id, pacbio=True)


class Pbmm2UbamWrapper:
    """Align an explicitly validated uBAM through pbmm2; no SM/tag rewriting."""

    def __init__(self, runner: CommandRunner, executable: str, samtools_executable: str,
                 options: Pbmm2Options) -> None:
        self.runner = runner
        self.executable = executable
        self.samtools = samtools_executable
        self.options = options

    def inspect(self, source: Path, sample_id: str, *, timeout: float) -> dict:
        """Check header, EOF, unmapped primary/unpaired flags and record count."""
        validation.validate_alignment_file(source)
        self.runner.require_executable(self.samtools)
        self.runner.run([self.samtools, "quickcheck", "-u", "-v", str(source)], timeout=timeout)
        header = self.runner.run([self.samtools, "view", "--no-PG", "-H", str(source)], timeout=timeout)
        validate_ubam_header(header.stdout or "", sample_id)
        counts = {}
        for label, flags in (("total", []), ("mapped", ["-F", "4"]),
                             ("paired_secondary_supplementary", ["--rf", "2305"])):
            result = self.runner.run([self.samtools, "view", "-c", *flags, str(source)], timeout=timeout)
            text = (result.stdout or "").strip()
            if not text.isdigit():
                raise InputValidationError("uBAM record scan did not return a count.")
            counts[label] = int(text)
        if counts["total"] == 0 or counts["mapped"] or counts["paired_secondary_supplementary"]:
            raise InputValidationError("uBAM must contain nonempty unpaired primary unmapped reads only.")
        assigned = 0
        for group in parse_header_records(header.stdout or "")["@RG"]:
            result = self.runner.run([self.samtools, "view", "-c", "-d", f"RG:{group['ID']}", str(source)], timeout=timeout)
            if not (result.stdout or "").strip().isdigit():
                raise InputValidationError("uBAM RG scan did not return a count.")
            assigned += int(result.stdout.strip())
        if assigned != counts["total"]:
            raise InputValidationError("uBAM records contain missing or undeclared RG tags.")
        counts["assigned_to_declared_read_groups"] = assigned
        return {"sample_id": sample_id, "record_counts": counts,
                "scope": "header_EOF_flags_and_count_not_read_quality_certification"}

    def run(self, source: Path, sample_id: str, reference: ReferenceGenome, output: Path,
            *, threads: int, timeout: float, dry_run: bool = False) -> dict:
        """Run once at a new destination; caller separately validates/indexes BAM."""
        validation.validate_fasta(reference.fasta, require_fai=True)
        inspection = self.inspect(source, sample_id, timeout=timeout)
        if output.exists():
            raise OutputValidationError(f"Native uBAM alignment output already exists: {output}")
        if type(threads) is not int or threads <= 0:
            raise InputValidationError("pbmm2 mapping threads must be positive.")
        if self.options.sort_threads is None:
            raise InputValidationError("Native uBAM alignment requires explicit sort_threads.")
        command = [self.executable, "align", str(reference.fasta), str(source), str(output),
                   "--preset", self.options.preset, "--sort", "--bam-index", "NONE",
                   "-j", str(threads), "-J", str(self.options.sort_threads),
                   "--log-level", self.options.log_level]
        if dry_run:
            return {"status": "PLANNED", "command": command, "input": inspection}
        version = Pbmm2Wrapper(runner=self.runner, executable=self.executable).detect_version()
        output.parent.mkdir(parents=True, exist_ok=True)
        result = self.runner.run(command, timeout=timeout, stderr_path=output.parent / "pbmm2.stderr.log")
        if not result.executed or result.returncode != 0:
            raise OutputValidationError("Native uBAM pbmm2 did not execute successfully.")
        validation.validate_output_file(output)
        return {"status": "COMPLETED", "command": command, "tool_version": version,
                "duration_seconds": result.duration_seconds, "returncode": result.returncode,
                "input": inspection, "conversion": "NONE_NATIVE_UBAM"}
