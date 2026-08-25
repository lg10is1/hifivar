"""Write one deterministic Phase 0 marker and its rule log."""

from pathlib import Path
from typing import Any


def main(workflow: Any) -> None:
    """Create parent directories and write the declared marker and log."""
    marker_path = Path(str(workflow.output[0]))
    log_path = Path(str(workflow.log[0]))
    marker_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    marker_path.write_text(
        str(workflow.params.marker_text),
        encoding="utf-8",
    )
    log_path.write_text(
        f"rule={workflow.params.rule_name}\nstatus=success\n",
        encoding="utf-8",
    )


main(snakemake)  # type: ignore[name-defined]
