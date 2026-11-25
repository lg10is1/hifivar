"""Bounded native Primer3 adapter; no implicit sequence upload or shell."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re

from hifivar.command import CommandRunner
from hifivar.exceptions import InputValidationError, OutputValidationError, ToolVersionError
from hifivar.primer_design import (PrimerDesignSettings, VariantTemplate, build_boulder_input,
                                  evaluate_candidates, parse_boulder_record)


class Primer3Wrapper:
    """Generate computational candidates in a fresh, caller-owned unit directory."""

    def __init__(self, executable: str = "primer3_core", *, runner: CommandRunner | None = None) -> None:
        if not isinstance(executable, str) or not executable.strip() or any(c in executable for c in "\r\n\x00"):
            raise InputValidationError("Primer3 executable must be a nonempty path/name.")
        self.executable = executable
        self.runner = runner or CommandRunner()

    def detect_version(self, *, timeout_seconds: float = 30.0) -> str:
        self.runner.require_executable(self.executable)
        result = self.runner.run([self.executable, "--about"], timeout=timeout_seconds)
        text = "\n".join(part for part in (result.stdout, result.stderr) if part)
        match = re.search(r"(?:primer3(?:_core)?\s*(?:release\s*)?)?(\d+\.\d+(?:\.\d+)?)", text, re.IGNORECASE)
        if match is None:
            raise ToolVersionError("Unable to determine native Primer3 version from --about.")
        return match.group(1)

    def design(self, template: VariantTemplate, output_dir: str | Path,
               settings: PrimerDesignSettings | None = None, *, dry_run: bool = False) -> dict[str, object]:
        settings = settings or PrimerDesignSettings()
        boulder = build_boulder_input(template, settings)
        if not isinstance(output_dir, (str, Path)) or not str(output_dir).strip():
            raise InputValidationError("Primer3 output_dir must be a nonempty path.")
        if not isinstance(dry_run, bool):
            raise InputValidationError("dry_run must be boolean.")
        output = Path(output_dir).expanduser().absolute()
        if output.exists() or output.is_symlink():
            raise OutputValidationError(f"Primer3 output directory already exists: '{output}'.")
        input_path = output / "primer3.input.txt"
        stdout_path = output / "primer3.output.txt"
        stderr_path = output / "primer3.stderr.log"
        command = [self.executable, "--strict_tags", "--io_version=4", str(input_path)]
        if dry_run:
            return {"status": "PLANNED", "variant": template.identity(), "settings": settings.to_dict(),
                    "command": command, "shell": False, "output_directory": str(output),
                    "specificity_status": "NOT_ASSESSED", "experimental_status": "NOT_RUN"}
        if settings.thermodynamic_parameters_path is not None and not Path(settings.thermodynamic_parameters_path).expanduser().is_dir():
            raise InputValidationError("Primer3 thermodynamic parameter directory does not exist.")
        version = self.detect_version(timeout_seconds=min(settings.timeout_seconds, 30))
        # Recheck after an adapter/tool probe: do not follow an output planted by it.
        if output.exists() or output.is_symlink():
            raise OutputValidationError("Primer3 output appeared during preflight; refusing overwrite.")
        output.mkdir(parents=True, exist_ok=False)
        input_path.write_text(boulder, encoding="utf-8", newline="\n")
        result = self.runner.run(command, cwd=output, timeout=settings.timeout_seconds,
                                 stdout_path=stdout_path, stderr_path=stderr_path)
        if not result.executed or result.returncode != 0:
            raise OutputValidationError("Primer3 command did not complete with executed=True and exit code zero.")
        if not stdout_path.is_file() or stdout_path.stat().st_size == 0 or stdout_path.stat().st_size > 10_000_000:
            raise OutputValidationError("Primer3 output is empty/missing or larger than 10 MB.")
        try:
            response = parse_boulder_record(stdout_path.read_text(encoding="utf-8"))
        except UnicodeDecodeError as error:
            raise OutputValidationError("Primer3 output is not valid UTF-8 Boulder text.") from error
        expected_id = next(line.split("=", 1)[1] for line in boulder.splitlines() if line.startswith("SEQUENCE_ID="))
        if response.get("SEQUENCE_ID") != expected_id:
            raise OutputValidationError("Primer3 response sequence identity differs from the requested template.")
        if response.get("SEQUENCE_TEMPLATE") != template.reference_sequence:
            raise OutputValidationError("Primer3 response does not echo the requested reference template.")
        outcome = evaluate_candidates(template, settings, response)
        outcome["provenance"] = {"tool": "primer3_core", "tool_version": version,
                                 "command": command, "returncode": result.returncode,
                                 "duration_seconds": result.duration_seconds,
                                 "stdout_path": str(stdout_path), "stderr_path": str(stderr_path),
                                 "boulder_input_sha256": hashlib.sha256(boulder.encode("utf-8")).hexdigest(),
                                 "output_sha256": hashlib.sha256(stdout_path.read_bytes()).hexdigest(),
                                 "shell": False}
        (output / "primer_design.json").write_text(json.dumps(outcome, ensure_ascii=False, indent=2), encoding="utf-8")
        return outcome


__all__ = ["Primer3Wrapper"]
