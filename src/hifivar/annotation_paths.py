"""Read-only, fail-closed output guards for annotation adapters.

No directory or file is created by this module. Resource identity is checked
using resolved paths and filesystem identities, not expensive database hashes.
This is pre-execution collision protection, not an operating-system sandbox.
"""

from __future__ import annotations

import os
from pathlib import Path
import shutil
from typing import Mapping, Sequence

from hifivar.exceptions import OutputValidationError


def _absolute(path: str | Path) -> Path:
    if not isinstance(path, (str, Path)) or not str(path).strip():
        raise OutputValidationError("ANNOTATION_OUTPUT_PROTECTION: nonempty paths are required.")
    try:
        return Path(path).expanduser().absolute()
    except (OSError, ValueError, RuntimeError) as error:
        raise OutputValidationError(f"ANNOTATION_OUTPUT_PROTECTION: invalid path: {error}") from error


def _resolved(path: Path) -> Path:
    try:
        return path.resolve()
    except (OSError, ValueError, RuntimeError) as error:
        raise OutputValidationError(f"ANNOTATION_OUTPUT_PROTECTION: cannot resolve '{path}': {error}") from error


def _identity(path: Path) -> tuple[int, int] | None:
    try:
        value = path.stat()
        return value.st_dev, value.st_ino
    except FileNotFoundError:
        if path.is_symlink():
            raise OutputValidationError(f"ANNOTATION_OUTPUT_PROTECTION: dangling symlink '{path}'.")
        return None
    except OSError as error:
        raise OutputValidationError(f"ANNOTATION_OUTPUT_PROTECTION: cannot inspect '{path}': {error}") from error


def protect_annotation_outputs(
    *, outputs: Mapping[str, str | Path | None], immutable_files: Sequence[str | Path],
    resource_directories: Sequence[str | Path], output_prefix: str | Path | None = None,
    executable: str | Path | None = None, protect_perl_helpers: bool = False,
) -> None:
    """Reject output/log aliases before any execution, mkdir, or unlink.

    An output may be a new sibling of input files in an owned run directory.
    It may not equal/contain an input, lie within/contain a resource directory,
    or hardlink/symlink to a protected resource. Nested resource-directory links
    are followed read-only with an inode visited set, including links outside
    the apparent cache/database path. Unreadable resources fail closed.

    Existing ANNOVAR prefix-derived helper files are also checked: native
    adapters can write them even though only two final products are collected.
    They are not treated as immutable merely because they live in the run.
    """
    targets: list[tuple[str, Path, Path, tuple[int, int] | None]] = []
    for label, supplied in outputs.items():
        if supplied is None:
            continue
        path = _absolute(supplied)
        if path.is_dir():
            raise OutputValidationError(f"ANNOTATION_OUTPUT_PROTECTION: output/log '{label}' is a directory: '{path}'.")
        targets.append((label, path, _resolved(path), _identity(path)))
    if output_prefix is not None:
        prefix = _absolute(output_prefix)
        if prefix.parent.exists():
            try:
                helpers = [path for path in prefix.parent.iterdir() if path.name.startswith(prefix.name + ".")]
            except OSError as error:
                raise OutputValidationError(f"ANNOTATION_OUTPUT_PROTECTION: cannot inspect output prefix: {error}") from error
            declared = {item[1] for item in targets}
            for path in helpers:
                if path not in declared:
                    targets.append(("annovar_prefix_helper", path, _resolved(path), _identity(path)))
    for i, (label, path, resolved, identity) in enumerate(targets):
        for other_label, other_path, other_resolved, other_identity in targets[:i]:
            if (resolved == other_resolved or resolved in other_resolved.parents
                    or other_resolved in resolved.parents
                    or (identity is not None and identity == other_identity)):
                raise OutputValidationError(
                    f"ANNOTATION_OUTPUT_PROTECTION: output/log aliases collide: '{label}' '{path}' and '{other_label}' '{other_path}'."
                )

    def check_resource(path: Path, *, directory: bool) -> None:
        resolved = _resolved(path)
        identity = _identity(path)
        for label, target, resolved_target, target_identity in targets:
            if (resolved_target == resolved or resolved_target in resolved.parents
                    or (directory and resolved in resolved_target.parents)
                    or (identity is not None and identity == target_identity)):
                raise OutputValidationError(
                    f"ANNOTATION_OUTPUT_PROTECTION: output/log '{label}' '{target}' aliases or overlaps immutable resource '{path}'."
                )

    program_files: list[Path] = []
    if executable is not None:
        try:
            text = str(executable)
            # Discovery only: never require, execute, or install a tool here.
            found = shutil.which(text)
            supplied_program = _absolute(executable)
            if found is not None:
                program_files.append(_resolved(Path(found)))
            elif supplied_program.exists() or supplied_program.is_symlink() or "/" in text or "\\" in text:
                program_files.append(_resolved(supplied_program))
            if protect_perl_helpers:
                for program in tuple(program_files):
                    if program.parent.is_dir():
                        program_files.extend(path for path in program.parent.iterdir() if path.name.endswith(".pl"))
        except (OSError, ValueError, RuntimeError) as error:
            raise OutputValidationError(f"ANNOTATION_OUTPUT_PROTECTION: cannot inspect program/helper identities: {error}") from error
    for supplied in (*immutable_files, *program_files):
        check_resource(_absolute(supplied), directory=False)
    roots = tuple(_absolute(path) for path in resource_directories)
    # Fast structural checks run before traversing potentially large caches.
    for root in roots:
        if not root.is_dir():
            raise OutputValidationError(f"ANNOTATION_OUTPUT_PROTECTION: resource directory is unavailable: '{root}'.")
        check_resource(root, directory=True)
    pending = list(roots)
    visited: set[tuple[int, int]] = set()
    while pending:
        directory = pending.pop()
        check_resource(directory, directory=True)
        identity = _identity(directory)
        if identity in visited:
            continue
        if identity is None:
            raise OutputValidationError(f"ANNOTATION_OUTPUT_PROTECTION: resource disappeared: '{directory}'.")
        visited.add(identity)
        try:
            with os.scandir(directory) as entries:
                for entry in entries:
                    child = Path(entry.path)
                    is_directory = entry.is_dir(follow_symlinks=True)
                    check_resource(child, directory=is_directory)
                    if is_directory:
                        pending.append(child)
        except OSError as error:
            raise OutputValidationError(f"ANNOTATION_OUTPUT_PROTECTION: cannot inspect resource tree '{directory}': {error}") from error


__all__ = ["protect_annotation_outputs"]
