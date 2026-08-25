# Phase 14 — Reporting, packaging, release and CI hardening

Status: implementation complete; independent Linux release validation pending.

## Final reporting

`hifivar.report` defines `FinalStatus`, `ReportArtifact`, `ToolRecord`,
`TrackReport`, and `FinalRunReport`. Tracks retain `COMPLETE`, `PARTIAL`,
`FAILED`, `NOT_RUN`, and `DISABLED`; disabled tracks do not become failures and
enabled missing/partial tracks do not become complete. The report records run,
project/schema/Git identity, reference, samples/cohort, enabled phases, tools,
artifacts, QC, benchmark, warnings, effective config and provenance.

JSON/YAML are machine-readable source records. Markdown and self-contained HTML
provide offline human summaries for sample/cohort, small variants, SV, TR,
phasing, assembly, manual review, annotation, cohort and benchmark sections.
The report performs no clinical interpretation and does not alter artifacts.

## Artifact and reproducibility bundle

`hifivar.bundle` creates a directory with `reports/`, `manifests/`, `configs/`,
`provenance/`, and explicitly selected `results/`. Selection is opt-in through
`ReportArtifact.selected_for_bundle` or explicit `BundleItem` objects. Large
BAM/CRAM/SAM/FASTQ/GFA/FASTA inputs become small pointer JSON files by default;
they are copied only with explicit `include_large=True`.

The reproducibility record contains redacted effective config, software
versions, redacted commands, environment summary, Git SHA, reference metadata
and optional sample sheet. Semantic secret keys and caller-declared redaction
values never enter stored config/command artifacts. Checksums are calculated
only for explicitly configured bundle items, not arbitrary WGS resources.

## Packaging and CI

Wheel/sdist now include packaged configuration plus the modular Snakemake
Snakefile, rules, bridges and environment YAML via standard setuptools data
files. `installed_workflow_root()` locates the installed workflow without a
source-tree or private-path assumption. A clean virtual-environment test builds
and installs the wheel outside the repository, verifies imports/resources, and
runs CLI help/version/config validation.

Lightweight GitHub Actions CI covers Python 3.10/3.12 Linux and Python 3.12
Windows, compileall, build, pytest, CLI/package tests and mock/dry-run Snakemake.
Heavy real tools and WGS data remain a separate Linux/HPC release gate.

## Release and security boundaries

- Version remains `0.0.1.dev0`; `0.1.0` is only a future candidate.
- No PyPI/GitHub/container publication, commit, tag or push is automatic.
- Wheel tests reject known local private paths and verify necessary resources.
- Audit evidence, credentials, keys, SIFs, logs, caches and raw genomic data
  are not package inputs.
- `RELEASE_CHECKLIST.md` requires human review and authorization.

## Known limitations

- Phase 13 hap.py/Truvari real benchmark is pending Linux validation.
- VEP remains environment-blocked; ANNOVAR has real evidence.
- Final reports summarize supplied phase records; they do not rediscover every
  artifact from an arbitrary directory or infer scientific success from file existence.
- HTML is deliberately static and contains no interactive/external assets.

Phase 15 is not started.
