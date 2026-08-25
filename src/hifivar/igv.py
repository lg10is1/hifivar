"""IGV Desktop batch generation and execution through CommandRunner."""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Sequence

from hifivar.command import CommandRunner, format_command
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar.review import ReviewEvidence, ReviewTarget
from hifivar.validation import validate_output_file


_VERSION = re.compile(r"\d+(?:\.\d+)+")


class IgvRunStatus(str, Enum):
    PLANNED = "planned"
    COMPLETED = "completed"


@dataclass(frozen=True, slots=True)
class IgvBatchPlan:
    """Deterministic IGV batch artifact and shell-free launcher command."""

    targets: tuple[ReviewTarget, ...]
    batch_path: Path
    snapshot_directory: Path
    batch_text: str
    args: tuple[str, ...]
    display_command: str

    @property
    def screenshots(self) -> tuple[Path, ...]:
        return tuple(path for target in self.targets for path in target.screenshot_paths)

    def to_dict(self) -> dict[str, object]:
        return {
            "targets": [target.review_id for target in self.targets],
            "batch_path": str(self.batch_path),
            "snapshot_directory": str(self.snapshot_directory),
            "batch_text": self.batch_text,
            "args": list(self.args),
            "display_command": self.display_command,
            "shell": False,
            "screenshots": [str(path) for path in self.screenshots],
        }


@dataclass(frozen=True, slots=True)
class IgvRunResult:
    plan: IgvBatchPlan
    status: IgvRunStatus
    tool_version: str | None
    duration_seconds: float
    evidence: tuple[ReviewEvidence, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "status": self.status.value,
            "tool_version": self.tool_version,
            "duration_seconds": self.duration_seconds,
            "plan": self.plan.to_dict(),
            "evidence": [item.to_dict() for item in self.evidence],
        }


class IgvWrapper:
    """Execute IGV Desktop `--batch`; never simulate GUI interaction in Python."""

    def __init__(self, executable: str = "igv.sh", *, runner: CommandRunner | None = None) -> None:
        if not isinstance(executable, str) or not executable.strip():
            raise InputValidationError("IGV executable must be a non-empty string.")
        self.executable = executable
        self.runner = runner or CommandRunner()

    def plan(
        self,
        targets: Sequence[ReviewTarget],
        *,
        batch_path: str | Path,
        snapshot_directory: str | Path,
    ) -> IgvBatchPlan:
        ordered = tuple(targets)
        if any(not isinstance(target, ReviewTarget) for target in ordered):
            raise InputValidationError("IGV targets must contain only ReviewTarget values.")
        ids = [target.review_id for target in ordered]
        if len(ids) != len(set(ids)):
            raise InputValidationError("IGV targets contain duplicate review IDs.")
        snapshots = Path(snapshot_directory).expanduser()
        for target in ordered:
            if target.output_directory / "screenshots" != snapshots:
                raise InputValidationError(
                    "IGV snapshot_directory must match every target output directory."
                )
        batch = Path(batch_path).expanduser()
        lines: list[str] = []
        for target in ordered:
            lines.extend((
                "new",
                f"genome {_batch_arg(target.reference_fasta)}",
                f"snapshotDirectory {_batch_arg(snapshots)}",
                f"load {_batch_arg(target.alignment_path)}",
                f"load {_batch_arg(target.source_vcf)}",
            ))
            for locus, screenshot in zip(target.loci, target.screenshot_paths, strict=True):
                lines.append(f"goto {locus.igv_locus}")
                lines.append(f"snapshot {screenshot.name}")
        lines.append("exit")
        text = "\n".join(lines) + "\n"
        args = (self.executable, "--batch", str(batch))
        return IgvBatchPlan(ordered, batch, snapshots, text, args, format_command(args))

    def detect_version(self) -> str:
        self.runner.require_executable(self.executable)
        result = self.runner.run([self.executable, "--version"])
        output = "\n".join(value for value in (result.stdout, result.stderr) if value)
        match = _VERSION.search(output)
        if match is None:
            raise ToolVersionError("Unable to parse IGV Desktop version from --version output.")
        return match.group(0)

    def run(
        self,
        plan: IgvBatchPlan,
        *,
        dry_run: bool = False,
        overwrite: bool = False,
        stdout_path: Path | None = None,
        stderr_path: Path | None = None,
    ) -> IgvRunResult:
        if not isinstance(plan, IgvBatchPlan):
            raise InputValidationError("IGV run requires IgvBatchPlan.")
        if not isinstance(overwrite, bool):
            raise InputValidationError("IGV overwrite must be boolean.")
        if dry_run:
            self.runner.run(plan.args, dry_run=True)
            return IgvRunResult(
                plan, IgvRunStatus.PLANNED, None, 0.0,
                _evidence(plan, generated=False),
            )

        _protect_outputs(plan, overwrite=overwrite)
        if not plan.targets:
            _write_batch(plan, overwrite=overwrite)
            plan.snapshot_directory.mkdir(parents=True, exist_ok=True)
            return IgvRunResult(
                plan, IgvRunStatus.COMPLETED, None, 0.0,
                (),
            )
        version = self.detect_version()
        _write_batch(plan, overwrite=overwrite)
        plan.snapshot_directory.mkdir(parents=True, exist_ok=True)
        result = self.runner.run(
            plan.args,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
        )
        for screenshot in plan.screenshots:
            validate_output_file(screenshot)
        return IgvRunResult(
            plan,
            IgvRunStatus.COMPLETED,
            version,
            result.duration_seconds,
            _evidence(plan, generated=True),
        )


def _evidence(plan: IgvBatchPlan, *, generated: bool) -> tuple[ReviewEvidence, ...]:
    return tuple(
        ReviewEvidence(
            target=target,
            batch_path=plan.batch_path,
            batch_command=plan.args,
            screenshots=target.screenshot_paths,
            generated=generated,
        )
        for target in plan.targets
    )


def _batch_arg(path: Path) -> str:
    text = str(path.absolute())
    if any(character in text for character in ('"', "\r", "\n")):
        raise InputValidationError(f"IGV batch path contains unsupported characters: '{path}'.")
    return f'"{text}"'


def _protect_outputs(plan: IgvBatchPlan, *, overwrite: bool) -> None:
    owned = (plan.batch_path, *plan.screenshots)
    for path in owned:
        if path.exists() and path.is_dir():
            raise OutputValidationError(f"IGV owned output is a directory: '{path}'.")
        if path.exists() and not overwrite:
            raise OutputValidationError(f"IGV output already exists: '{path}'.")
    if overwrite:
        for screenshot in plan.screenshots:
            screenshot.unlink(missing_ok=True)


def _write_batch(plan: IgvBatchPlan, *, overwrite: bool) -> None:
    temporary = plan.batch_path.with_name(f".{plan.batch_path.name}.hifivar.tmp")
    if temporary.exists():
        if not overwrite:
            raise OutputValidationError(f"IGV batch temporary file exists: '{temporary}'.")
        temporary.unlink()
    try:
        plan.batch_path.parent.mkdir(parents=True, exist_ok=True)
        temporary.write_text(plan.batch_text, encoding="utf-8", newline="\n")
        os.replace(temporary, plan.batch_path)
    except OSError as error:
        temporary.unlink(missing_ok=True)
        raise OutputValidationError(f"Unable to write IGV batch '{plan.batch_path}': {error}") from error


__all__ = ["IgvBatchPlan", "IgvRunResult", "IgvRunStatus", "IgvWrapper"]
