"""Fail-closed figure selection for a checksum-pinned historical benchmark."""

import hashlib
import json
from pathlib import Path


def select_cmrg_callers(source: Path, policy: Path) -> tuple[list[str], dict]:
    """Return permitted caller labels without changing source data or metrics.

    A decision applies only to the exact archived TSV. A changed input needs
    renewed review rather than silently inheriting an old identity disposition.
    Inclusion is historical, qualified evidence, not a new identity assertion.
    """
    import csv

    payload = source.read_bytes()
    digest = hashlib.sha256(payload).hexdigest()
    qualification = json.loads(policy.read_text(encoding="utf-8"))
    if digest != qualification["source_sha256"]:
        raise ValueError("CMRG source changed: renew figure identity review")
    rows = list(csv.DictReader(payload.decode("utf-8").splitlines(), delimiter="\t"))
    callers = [row["caller"] for row in rows]
    dispositions = qualification["callers"]
    if len(callers) != len(set(callers)) or set(callers) != set(dispositions):
        raise ValueError("CMRG caller membership does not match reviewed policy")
    allowed = []
    excluded = {}
    for caller in callers:
        decision = dispositions[caller]
        if decision["figure_disposition"] == "INCLUDE_QUALIFIED_HISTORICAL":
            allowed.append(caller)
        elif decision["figure_disposition"] == "EXCLUDE_UNRESOLVED_SAMPLE_IDENTITY":
            excluded[caller] = decision["reason"]
        else:
            raise ValueError(f"Unreviewed CMRG figure disposition: {caller}")
    if not allowed:
        raise ValueError("No qualified CMRG callers remain for plotting")
    return allowed, {
        "source_sha256": digest,
        "policy_sha256": hashlib.sha256(policy.read_bytes()).hexdigest(),
        "included_callers": allowed,
        "excluded_callers": excluded,
        "scope": qualification["scope"],
        "interpretation": qualification["interpretation"],
    }
