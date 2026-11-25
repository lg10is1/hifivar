# Validated deployment matrix

[简体中文](zh_CN/deployment.md)

## Current development context / 当前开发部署范围

Review date: 21 October 2026. See [capabilities and evidence](capabilities.md)
for current source-plan routes. The RC4 table below is a **historical deployment
snapshot**, not proof of absence in later deployments and not a list of bundled
tools. Later identified development cases include native VEP 116.3 with an
explicit offline cache, IGV 2.19.8 with Java 21/Xvfb, DeepVariant 1.10.0 through
Apptainer and Truvari 5.4.0 with a separately pinned compatible HTS stack.

Native annotation and producer-to-benchmark/review cases were independently
accepted within their recorded regions and package identities. This does not
clear a separate ready-input VEP deployment failure, incomplete latest ready-PAV
route, or establish every backend combination. Use the exact interpreter, tool
version, image digest, cache/database identity and runtime receipt for the actual
route. A newer executable elsewhere on PATH is not the configured deployment.

后续VEP真实运行通过不追认旧失败，也不代表已自动部署全部依赖。IGV/Xvfb、容器、
注释数据库和合法授权仍由站点准备；核心Python包不包含这些资源。以下RC4时期表格
保留历史事实，其中“blocked/NOT_RUN”只适用于当时该部署，不作为现在的全局缺失结论。

## RC4 historical deployment snapshot

This matrix records evidence available through the 0.1.0rc4 validation cycle.
A wrapper/mock PASS is not represented as a real-tool PASS. Site paths, images,
databases and credentials stay outside the Python package.

| Component | Version/backend | Evidence status | Production contract |
|---|---|---|---|
| PAV | 2.4.6, Apptainer | REAL CALLER + SV-ONLY ADAPTER PASS | `apptainer exec --bind <root> <image> snakemake`; host copy of the same image's `/opt/pav` workflow is required and refreshed on image upgrade. PAV's mixed root VCF is retained immutably; the Linux-validated SV-only handoff uses host bgzip/tabix 1.23.1 and PAV's version-locked VARTYPE boundary. |
| IGV | 2.19.6 under Xvfb | REAL PASS | IGV batch mode through `IgvWrapper`; headless display is site-provisioned. IGV may attempt optional default-track network access. |
| GLnexus | 1.4.1 + bcftools 1.21 | REAL PASS | Pinned `workflow/envs/glnexus.yaml`; GLnexus/bcftools wrapper and Snakemake cohort track validated. `memory_gb` is an explicit run-specific scheduler/CLI cap; a three-sample WGS run peaked near 153 GB and used 192–200 GB planning, which is evidence for that workload rather than a universal default. |
| Jasmine finalization | Jasmine 1.1.5 + tabix/HTSlib 1.23.1 | VCFv4.4 COMPATIBILITY CONTRACT | HiFiVar preserves the original header, restores only source-defined INFO/FORMAT metadata, and requires a compatible tabix/HTSlib deployment for VCFv4.4 rather than downgrading the header. PAV input must use the separate SV-only artifact, not the mixed root VCF. |
| ANNOVAR | 2020Jun08 + explicit hg38 database | REAL PASS | Native Perl adapter; database root/version provisioned outside HiFiVar; no automatic download. |
| VEP | offline/cache boundary | ENVIRONMENT BLOCKED | Wrapper/mock audited, but no real executable/cache was available. Do not claim real support until revalidated. |
| hap.py | compatibility target 0.3.15 | REAL PASS | Native/managed Linux environment; explicit truth VCF, confident BED, reference and resource versions. Empty `--version` output requires an explicit trustworthy configured version. |
| Truvari comparison | 5.4.0 native environment | VERSION/SMOKE PASS | Phase 9 version detection and comparison boundary audited. |
| Truvari Phase 13 benchmark | 5.4.0 compatibility target | REAL PASS | Explicit truth/confident regions and thresholds; Phase 13 Linux benchmark validation completed. |

## HiFiVar core distributions

| Distribution | Identity | Status |
|---|---|---|
| Python wheel/sdist | `0.1.0rc4` | Windows and independent Linux clean-install, metadata, packaged-resource, content, and security audits PASS. |
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
