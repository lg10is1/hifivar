"""BAM/CRAM header parsing via samtools for alignment input classification.

Uses samtools view -H through CommandRunner to extract header metadata
without requiring pysam or other Python HTSlib bindings.
"""

from __future__ import annotations

import re
from pathlib import Path

from hifivar.command import CommandRunner
from hifivar.exceptions import (
    InputValidationError,
    ToolVersionError,
)
from hifivar.logging_utils import get_logger


_LOGGER = get_logger(__name__)

_HD_SORT_PATTERN = re.compile(r"SO:(\S+)")
_SQ_PATTERN = re.compile(r"SN:([^\t]+)")
_RG_ID_PATTERN = re.compile(r"ID:([^\t]+)")
_RG_SM_PATTERN = re.compile(r"SM:([^\t]+)")
_VERSION_PATTERN = re.compile(
    r"samtools\s+([0-9]+(?:\.[0-9]+){1,2}(?:[-+._A-Za-z0-9]*)?)",
    re.IGNORECASE,
)


class SamtoolsHeaderParser:
    """Parse BAM/CRAM headers using samtools view -H via CommandRunner."""

    def __init__(
        self,
        *,
        runner: CommandRunner | None = None,
        executable: str | Path = "samtools",
    ) -> None:
        self.runner = runner or CommandRunner()
        self.executable = str(executable)

    def detect_version(self) -> str:
        """Require samtools and return its parsed version string."""
        self.runner.require_executable(self.executable)
        result = self.runner.run([self.executable, "--version"])
        output = "\n".join(
            text for text in (result.stdout, result.stderr) if text is not None
        )
        match = _VERSION_PATTERN.search(output)
        if match is None:
            raise ToolVersionError(
                "Unable to parse samtools version from `samtools --version` output."
            )
        version = match.group(1)
        _LOGGER.debug("Detected samtools version: %s", version)
        return version

    def parse_header(self, path: Path) -> AlignmentHeaderMetadata:
        """Parse BAM/CRAM header and return structured metadata.

        Args:
            path: Path to BAM or CRAM file.

        Returns:
            Parsed header metadata.

        Raises:
            InputValidationError: If header parsing fails.
        """
        from hifivar.alignment_input import (
            AlignmentHeaderMetadata,
            SortOrder,
        )

        input_path = Path(path)
        if not input_path.is_file():
            raise InputValidationError(
                f"BAM/CRAM header input is not a file: '{input_path}'."
            )

        result = self.runner.run(
            [self.executable, "view", "-H", str(input_path)],
            check=False,
        )

        if result.returncode != 0:
            stderr_text = result.stderr or ""
            raise InputValidationError(
                f"samtools view -H failed for '{input_path}' "
                f"(return code {result.returncode}): {stderr_text.strip()}"
            )

        header_text = result.stdout or ""
        return _parse_header_text(header_text)


def _parse_header_text(header_text: str) -> AlignmentHeaderMetadata:
    """Parse samtools header text into AlignmentHeaderMetadata."""
    from hifivar.alignment_input import (
        AlignmentHeaderMetadata,
        SortOrder,
    )

    has_hd_header = False
    sort_order = SortOrder.UNKNOWN
    has_sq_lines = False
    sample_names: list[str] = []
    read_group_ids: list[str] = []
    reference_contigs: list[str] = []

    for line in header_text.splitlines():
        line = line.strip()
        if not line:
            continue

        if line.startswith("@HD"):
            has_hd_header = True
            sort_match = _HD_SORT_PATTERN.search(line)
            if sort_match:
                sort_value = sort_match.group(1).lower()
                try:
                    sort_order = SortOrder(sort_value)
                except ValueError:
                    sort_order = SortOrder.UNKNOWN

        elif line.startswith("@SQ"):
            has_sq_lines = True
            sq_match = _SQ_PATTERN.search(line)
            if sq_match:
                reference_contigs.append(sq_match.group(1))

        elif line.startswith("@RG"):
            rg_id_match = _RG_ID_PATTERN.search(line)
            if rg_id_match:
                read_group_ids.append(rg_id_match.group(1))
            rg_sm_match = _RG_SM_PATTERN.search(line)
            if rg_sm_match:
                sample_names.append(rg_sm_match.group(1))

    return AlignmentHeaderMetadata(
        sort_order=sort_order,
        has_sq_lines=has_sq_lines,
        sample_names=tuple(sample_names),
        read_group_ids=tuple(read_group_ids),
        has_hd_header=has_hd_header,
        reference_contigs=tuple(reference_contigs),
    )


__all__ = ["SamtoolsHeaderParser"]
