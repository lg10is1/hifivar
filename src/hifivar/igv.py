"""IGV Desktop batch generation and execution through CommandRunner."""

from __future__ import annotations

import os
import math
import re
import shutil
import stat
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Sequence

from hifivar.command import CommandRunner, _default_stdout_log_path, format_command
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar.review import ReviewEvidence, ReviewTarget
from hifivar.junction_consumers import revalidate_junction_record
from hifivar.snapshot_validation import require_png_decoder, validate_png_snapshots
from hifivar.snapshot_session import validate_igv_session


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
    max_panel_height: int | None = None
    sleep_interval_ms: int | None = None
    verify_session: bool = False

    @property
    def screenshots(self) -> tuple[Path, ...]:
        return tuple(path for target in self.targets for path in target.screenshot_paths)

    @property
    def sessions(self) -> tuple[Path, ...]:
        return tuple(path.with_suffix('.session.xml') for path in self.screenshots) if self.verify_session else ()

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
            "verify_session": self.verify_session,
            "sessions": [str(path) for path in self.sessions],
            "execution_cwd": str(self.snapshot_directory.absolute()) if self.verify_session else None,
            "display_settings": {"max_panel_height": self.max_panel_height,
                                 "sleep_interval_ms": self.sleep_interval_ms},
        }


@dataclass(frozen=True, slots=True)
class IgvRunResult:
    plan: IgvBatchPlan
    status: IgvRunStatus
    tool_version: str | None
    duration_seconds: float
    evidence: tuple[ReviewEvidence, ...]
    snapshot_validation: tuple[dict[str, object], ...] = ()
    session_validation: tuple[dict[str, object], ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "status": self.status.value,
            "tool_version": self.tool_version,
            "duration_seconds": self.duration_seconds,
            "plan": self.plan.to_dict(),
            "evidence": [item.to_dict() for item in self.evidence],
            "snapshot_validation": {
                "expected_count": len(self.plan.screenshots),
                "validated_count": len(self.snapshot_validation),
                "images": list(self.snapshot_validation),
                "human_review": "NOT_REVIEWED",
                "track_loading": "NOT_ASSESSED",
            },
            "session_validation": {"enabled": self.plan.verify_session,
                                   "expected_count": len(self.plan.sessions),
                                   "validated_count": len(self.session_validation),
                                   "sessions": list(self.session_validation)},
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
        max_panel_height: int | None = None,
        sleep_interval_ms: int | None = None,
        verify_session: bool = False,
    ) -> IgvBatchPlan:
        ordered = tuple(targets)
        if not isinstance(verify_session, bool):
            raise InputValidationError('IGV verify_session must be boolean.')
        if any(not isinstance(target, ReviewTarget) for target in ordered):
            raise InputValidationError("IGV targets must contain only ReviewTarget values.")
        for target in ordered:
            if target._junction_record:
                revalidate_junction_record(target._junction_record)
        ids = [target.review_id for target in ordered]
        if len(ids) != len(set(ids)):
            raise InputValidationError("IGV targets contain duplicate review IDs.")
        for name, value, lower, upper in (("max_panel_height", max_panel_height, 64, 50000),
                                          ("sleep_interval_ms", sleep_interval_ms, 0, 60000)):
            if value is not None and (isinstance(value, bool) or not isinstance(value, int) or not lower <= value <= upper):
                raise InputValidationError(f"IGV {name} must be an integer in [{lower}, {upper}].")
        snapshots = Path(snapshot_directory).expanduser()
        for target in ordered:
            if target.output_directory / "screenshots" != snapshots:
                raise InputValidationError(
                    "IGV snapshot_directory must match every target output directory."
                )
        batch = Path(batch_path).expanduser()
        lines: list[str] = []
        references = {str(target.reference_fasta.absolute()) for target in ordered}
        for reference in references:
            if any(c.isspace() or ord(c) < 32 or c in {'"', "'"} for c in reference):
                raise InputValidationError('IGV_REFERENCE_PATH_UNSUPPORTED: native IGV 2.19.8 reference paths must not contain whitespace, controls or quotes; no input was moved or substituted.')
        # BatchRunner inspects the first command to decide whether to load its
        # default genome. Select the explicit local reference before any UI
        # setting, avoiding an unrelated default-genome network initialization.
        if ordered:
            lines.append(f"genome {_genome_batch_arg(ordered[0].reference_fasta)}")
        if max_panel_height is not None:
            lines.append(f"maxPanelHeight {max_panel_height}")
        if sleep_interval_ms is not None:
            lines.append(f"setSleepInterval {sleep_interval_ms}")
        for target in ordered:
            lines.append('new')
            lines.append(f"genome {_genome_batch_arg(target.reference_fasta)}")
            lines.extend((
                f"snapshotDirectory {_batch_arg(snapshots)}",
                f"load {_batch_arg(target.alignment_path)}",
                f"load {_batch_arg(target.source_vcf)}",
            ))
            for locus, screenshot in zip(target.loci, target.screenshot_paths, strict=True):
                lines.append(f"goto {locus.igv_locus}")
                if verify_session:
                    # saveSession tests the .xml suffix before stripping quotes.
                    # A safe relative filename and owned cwd also support spaces.
                    lines.append(f"saveSession ./{screenshot.with_suffix('.session.xml').name}")
                lines.append(f"snapshot {screenshot.name}")
        lines.append("exit")
        text = "\n".join(lines) + "\n"
        program = self.executable
        if verify_session and ('/' in program or '\\' in program):
            program = str(Path(program).expanduser().absolute())
        args = (program, "--batch", str(batch.absolute() if verify_session else batch))
        return IgvBatchPlan(ordered, batch, snapshots, text, args, format_command(args),
                            max_panel_height, sleep_interval_ms, verify_session)

    def detect_version(self) -> str:
        self.runner.require_executable(self.executable)
        result = self.runner.run([self.executable, "--version"], timeout=30)
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
        timeout_seconds: float = 14400,
    ) -> IgvRunResult:
        if not isinstance(plan, IgvBatchPlan):
            raise InputValidationError("IGV run requires IgvBatchPlan.")
        if not isinstance(overwrite, bool):
            raise InputValidationError("IGV overwrite must be boolean.")
        if (isinstance(timeout_seconds, bool) or not isinstance(timeout_seconds, (int, float))
                or not math.isfinite(timeout_seconds) or timeout_seconds <= 0):
            raise InputValidationError("IGV timeout_seconds must be finite and positive.")
        _protect_write_paths(
            plan, executable=self.executable, stdout_path=stdout_path, stderr_path=stderr_path,
        )
        if dry_run:
            self.runner.run(plan.args, dry_run=True)
            return IgvRunResult(
                plan, IgvRunStatus.PLANNED, None, 0.0,
                _evidence(plan, generated=False),
            )

        for target in plan.targets:
            if target._junction_record:
                revalidate_junction_record(target._junction_record)
        _protect_outputs(plan, overwrite=overwrite)
        if not plan.targets:
            _write_batch(plan, overwrite=overwrite)
            plan.snapshot_directory.mkdir(parents=True, exist_ok=True)
            return IgvRunResult(
                plan, IgvRunStatus.COMPLETED, None, 0.0,
                (),
            )
        require_png_decoder()
        version = self.detect_version()
        # A version adapter must not retarget a log or snapshot into an input.
        _protect_write_paths(
            plan, executable=self.executable, stdout_path=stdout_path, stderr_path=stderr_path,
        )
        _protect_outputs(plan, overwrite=overwrite)
        if overwrite:
            for screenshot in (*plan.screenshots, *plan.sessions):
                screenshot.unlink(missing_ok=True)
        _write_batch(plan, overwrite=overwrite)
        plan.snapshot_directory.mkdir(parents=True, exist_ok=True)
        result = self.runner.run(
            plan.args,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            timeout=timeout_seconds,
            cwd=plan.snapshot_directory.absolute() if plan.verify_session else None,
        )
        snapshot_validation = validate_png_snapshots(plan.screenshots)
        session_validation = []
        if plan.verify_session:
            for target in plan.targets:
                for locus, screenshot in zip(target.loci, target.screenshot_paths, strict=True):
                    session_validation.append(validate_igv_session(screenshot.with_suffix('.session.xml'),
                        reference=target.reference_fasta, alignment=target.alignment_path,
                        source_vcf=target.source_vcf, contig=locus.contig,
                        start=locus.window_start, end=locus.window_end))
        for target in plan.targets:
            if target._junction_record:
                revalidate_junction_record(target._junction_record)
        return IgvRunResult(
            plan,
            IgvRunStatus.COMPLETED,
            version,
            result.duration_seconds,
            _evidence(plan, generated=True),
            snapshot_validation,
            tuple(session_validation),
        )


def _evidence(plan: IgvBatchPlan, *, generated: bool) -> tuple[ReviewEvidence, ...]:
    return tuple(
        ReviewEvidence(
            target=target,
            batch_path=plan.batch_path,
            batch_command=plan.args,
            screenshots=target.screenshot_paths,
            generated=generated and bool(target.screenshot_paths),
        )
        for target in plan.targets
    )


def _batch_arg(path: Path) -> str:
    text = str(path.absolute())
    if any(character in text for character in ('"', "\r", "\n")):
        raise InputValidationError(f"IGV batch path contains unsupported characters: '{path}'.")
    return f'"{text}"'


def _genome_batch_arg(path: Path) -> str:
    """IGV 2.19.x genome tokens retain quotes instead of stripping them.

    Use an unquoted local path. Whitespace paths fail before any execution;
    no implicit copy, symlink or genome substitution is made.
    """
    text = str(path.absolute())
    if any(c.isspace() or c in {'"', "'"} or ord(c) < 32 for c in text):
        raise InputValidationError('IGV genome batch requires a reference path without whitespace or quotes; no reference was moved or substituted.')
    return text


def _protect_outputs(plan: IgvBatchPlan, *, overwrite: bool) -> None:
    owned = (plan.batch_path, *plan.screenshots, *plan.sessions)
    for path in owned:
        if path.exists() and path.is_dir():
            raise OutputValidationError(f"IGV owned output is a directory: '{path}'.")
        if path.exists() and not overwrite:
            raise OutputValidationError(f"IGV output already exists: '{path}'.")


def _guard_path(path: str | Path) -> tuple[Path, tuple[int, int] | None, bool]:
    """Resolve one finite write/resource path without creating it or reading data."""
    if not isinstance(path, (str, Path)) or not str(path).strip():
        raise OutputValidationError("IGV_OUTPUT_PROTECTION: paths must be nonempty.")
    try:
        value = Path(path).expanduser().absolute()
        resolved = value.resolve()
        try:
            info = value.stat()
        except FileNotFoundError:
            if value.is_symlink():
                raise OutputValidationError(f"IGV_OUTPUT_PROTECTION: dangling symlink '{value}'.")
            return resolved, None, False
        return resolved, (info.st_dev, info.st_ino), stat.S_ISDIR(info.st_mode)
    except (OSError, ValueError, RuntimeError) as error:
        raise OutputValidationError(f"IGV_OUTPUT_PROTECTION: cannot inspect '{path}': {error}") from error


def _protect_write_paths(
    plan: IgvBatchPlan, *, executable: str,
    stdout_path: str | Path | None, stderr_path: str | Path | None,
    additional_outputs: Sequence[tuple[str, str | Path]] = (),
) -> None:
    """Reject direct, resolved and hardlink aliases before any write or command.

    Check the actual implicit stdout destination used by CommandRunner, batch
    staging path, every screenshot and the finite selected BND bundle members.
    This is a pre-execution guard, not a recursive resource scan or OS sandbox.
    """
    outputs: list[tuple[str, str | Path]] = [
        ("batch", plan.batch_path),
        ("batch_temporary", plan.batch_path.with_name(f".{plan.batch_path.name}.hifivar.tmp")),
        *((f"screenshot_{i}", path) for i, path in enumerate(plan.screenshots)),
        *((f"session_{i}", path) for i, path in enumerate(plan.sessions)),
        *additional_outputs,
    ]
    if stdout_path is not None:
        outputs.append(("stdout", stdout_path))
    elif stderr_path is not None:
        _guard_path(stderr_path)
        outputs.append(("implicit_stdout", _default_stdout_log_path(Path(stderr_path).expanduser())))
    if stderr_path is not None:
        outputs.append(("stderr", stderr_path))
    checked = [(label, path, *_guard_path(path)) for label, path in outputs]
    for i, (label, path, resolved, identity, directory) in enumerate(checked):
        if directory:
            raise OutputValidationError(f"IGV_OUTPUT_PROTECTION: output '{label}' is a directory: '{path}'.")
        for other_label, other_path, other_resolved, other_identity, _ in checked[:i]:
            if (resolved == other_resolved or resolved in other_resolved.parents
                    or other_resolved in resolved.parents
                    or (identity is not None and identity == other_identity)):
                raise OutputValidationError(
                    f"IGV_OUTPUT_PROTECTION: output aliases collide: '{label}' '{path}' and '{other_label}' '{other_path}'."
                )
    resources: set[Path] = set()
    for target in plan.targets:
        resources.update((target.source_vcf, target.alignment_path, target.reference_fasta,
                          Path(str(target.reference_fasta) + ".fai")))
        for path, suffixes in ((target.source_vcf, (".tbi", ".csi")),
                               (target.alignment_path, (".bai", ".csi", ".crai"))):
            for suffix in suffixes:
                for index in (Path(str(path) + suffix), path.with_suffix(suffix)):
                    if index.exists() or index.is_symlink():
                        resources.add(index)
        if target.trgt_visualization_path is not None:
            resources.add(target.trgt_visualization_path)
        if target.bnd_bundle_manifest is not None:
            resources.add(Path(target.bnd_bundle_manifest))
        if target._junction_record is not None and target._junction_record.bundle_binding is not None:
            binding = target._junction_record.bundle_binding
            try:
                resources.add(Path(binding["manifest"]["path"]))
                resources.update(Path(item["path"]) for item in binding["artifacts"].values())
            except (KeyError, TypeError, AttributeError) as error:
                raise InputValidationError("IGV output guard requires a complete selected BND bundle binding.") from error
    # Discovery only, including a symlink's actual target. Never execute a tool
    # to decide whether output paths are safe.
    try:
        program = shutil.which(executable)
        supplied = Path(executable).expanduser()
        if program is not None:
            resources.add(Path(program))
        elif supplied.exists() or supplied.is_symlink() or supplied.is_absolute():
            resources.add(supplied)
    except (OSError, ValueError, RuntimeError) as error:
        raise OutputValidationError(f"IGV_OUTPUT_PROTECTION: cannot inspect executable: {error}") from error
    for resource in sorted(resources, key=str):
        resolved, identity, _ = _guard_path(resource)
        for label, path, target, target_identity, _ in checked:
            if (target == resolved or target in resolved.parents or resolved in target.parents
                    or (identity is not None and identity == target_identity)):
                raise OutputValidationError(
                    f"IGV_OUTPUT_PROTECTION: output '{label}' '{path}' aliases or overlaps immutable input '{resource}'."
                )


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
