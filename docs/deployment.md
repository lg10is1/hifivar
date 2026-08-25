# Validated deployment matrix

This matrix records evidence available in the repository as of Phase 14. A
wrapper/mock PASS is not represented as a real-tool PASS. Site paths, images,
databases and credentials stay outside the Python package.

| Component | Version/backend | Evidence status | Production contract |
|---|---|---|---|
| PAV | 2.4.6, Apptainer | REAL PASS | `apptainer exec --bind <root> <image> snakemake`; host copy of the same image's `/opt/pav` workflow is required and refreshed on image upgrade. |
| IGV | 2.19.6 under Xvfb | REAL PASS | IGV batch mode through `IgvWrapper`; headless display is site-provisioned. IGV may attempt optional default-track network access. |
| GLnexus | 1.4.1 + bcftools 1.21 | REAL PASS | Pinned `workflow/envs/glnexus.yaml`; GLnexus/bcftools wrapper and Snakemake cohort track validated. |
| ANNOVAR | 2020Jun08 + explicit hg38 database | REAL PASS | Native Perl adapter; database root/version provisioned outside HiFiVar; no automatic download. |
| VEP | offline/cache boundary | ENVIRONMENT BLOCKED | Wrapper/mock audited, but no real executable/cache was available. Do not claim real support until revalidated. |
| hap.py | compatibility target 0.3.15 | REAL PASS | Native/managed Linux environment; explicit truth VCF, confident BED, reference and resource versions. Empty `--version` output requires an explicit trustworthy configured version. |
| Truvari comparison | 5.4.0 native environment | VERSION/SMOKE PASS | Phase 9 version detection and comparison boundary audited. |
| Truvari Phase 13 benchmark | 5.4.0 compatibility target | REAL PASS | Explicit truth/confident regions and thresholds; Phase 13 Linux benchmark validation completed. |

## HiFiVar core distributions

| Distribution | Identity | Status |
|---|---|---|
| Python wheel/sdist | `0.1.0rc1` | Windows clean-install PASS; Linux release-candidate wheel PASS before this distribution delta. |
| Conda source environment/recipe | repository `environment.yml` / `conda-recipe` | `LINUX_REAL_VERIFICATION: NOT_RUN` |
| Docker core | repository `Dockerfile`, `python:3.12-slim-bookworm` base | `LINUX_REAL_VERIFICATION: NOT_RUN` |
| Apptainer core | repository `containers/hifivar.def`, `python:3.12-slim-bookworm` base | `LINUX_REAL_VERIFICATION: NOT_RUN` |

The core distribution statuses do not change any caller-level validation. They
must be revalidated independently after the public repository URL and immutable
container digests are assigned.

## Container boundary

PAV is the only production container deployment asserted here, using the
validated Apptainer pattern. This document does not introduce or validate an
alternative PAV Docker/native backend. Container images, launchers, extracted
PAV workflow trees and checksums are deployment artifacts—not Python package
source and not wheel contents.

## Revalidation policy

Any executable, database/cache, environment YAML, image digest, launcher,
reference build, truth set or confident-region change invalidates the matching
real-tool evidence until a tiny Linux run is repeated. Preserve commands,
stdout/stderr, effective config, checksums and output validation evidence.
