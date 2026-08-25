"""Aggregate independent Phase 12 track provenance without merging variants."""

import json
import hashlib
import os
from datetime import datetime, timezone
from pathlib import Path

import yaml

from hifivar import __version__
from hifivar.context import AnalysisContext

config = snakemake.config  # type: ignore[name-defined]
context = AnalysisContext.from_config(config)
tracks = [json.loads(Path(str(path)).read_text(encoding="utf-8")) for path in snakemake.input]  # type: ignore[name-defined]
payload = {
    "schema_version": 1,
    "hifivar_version": __version__,
    "git_commit": os.environ.get("HIFIVAR_GIT_COMMIT"),
    "created_at": datetime.now(timezone.utc).isoformat(),
    "effective_config": config,
    "cohort": {
        "cohort_id": str(config["cohort"]["cohort_id"]),
        "sample_ids": list(context.sample_ids),
        "sample_order_sha256": hashlib.sha256(("\n".join(context.sample_ids) + "\n").encode("utf-8")).hexdigest(),
        "reference": context.reference.to_dict(include_contigs=True),
    },
    "tracks": tracks,
}
Path(str(snakemake.output.json)).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")  # type: ignore[name-defined]
Path(str(snakemake.output.yaml)).write_text(yaml.safe_dump(payload, allow_unicode=True, sort_keys=False), encoding="utf-8")  # type: ignore[name-defined]
