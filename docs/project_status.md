# HiFiVar Project Status

[简体中文](zh_CN/project_status.md)

Status date: 2026-10-30

Previous independently validated public baseline: `0.1.0rc4` (pre-release)

Current source/package candidate: `0.1.0rc5`
Release availability: check [GitHub Releases](https://github.com/lg10is1/hifivar/releases)
for the actual tag, assets and release status; source presence does not establish publication.
Production platform: Linux/HPC  
License: Apache-2.0

## Current boundary

See the bilingual [capability and evidence matrix](capabilities.md) for current
routes. Dated test/build results below identify their own payloads, not every
subsequent package sharing the same version.

HiFiVar implements the Phase 0–15 modular contracts for configuration,
provenance, QC, alignment, small variants, read- and assembly-based structural
variants, tandem repeats, phasing, harmonization, manual review, annotation,
cohort analysis, benchmarking, and reporting. External bioinformatics tools
remain independently provisioned and versioned.

The packaged 0.1.0rc4 Snakemake small/SV/TR branches consume an existing
indexed BAM or CRAM. The Phase 2 Python API supports pbmm2 alignment planning
and execution, but there is not yet a unified `hifivar run` command or packaged
FASTQ/uBAM-to-results alignment rule in that release. The post-0.1.0rc4
development line also contains an independently tested minimap2 wrapper,
FASTQ/BAM input classification, alignment validation, and a standalone
Snakemake alignment rule. Those development changes are not part of 0.1.0rc4
and are not in that published release. The current development-source
`hifivar run` connects selected per-sample calling branches to indexed BAMs.
Opt-in sequential FASTQ/native-uBAM alignment preparation then hands validated
BAMs to that calling DAG in the same invocation. This is a bounded two-stage
launcher, not all-branch execution or a unified raw-input Snakemake DAG. It is
included in the [identified 15 October development payload](rc5_candidate.md#validated-development-payload--已验收开发包),
not the initial security-only RC5 build; see [raw-input guide](raw_input_run.md).

## Validation summary

- Windows regression: 917 passed, 2 skipped, 0 failed on 2026-09-06; the two
  skips require Windows symbolic-link privileges.
- Independent Linux RC4 package, clean-install, workflow-resource, security,
  PAV SV-only, and release validation: PASS.
- The identified 15 October RC5 development payload passed Linux regression
  (1274 tests), 36 Snakemake integration tests, source/artifact audits and isolated
  wheel/sdist installation checks. Windows regression passed 1270 tests with
  four platform-specific skips. These dated results do not identify a later
  rebuilt package without its own manifest and checksums; the Linux installation
  reused an existing read-only dependency pool, not a fresh all-tool Conda solve.
- Real Linux executions include pbmm2, DeepVariant, the four read-SV callers,
  TRGT, HiPhase, hifiasm, PAV, SVIM-asm, Jasmine, GLnexus, IGV, ANNOVAR,
  hap.py, and Truvari within their documented scopes.
- Later identified candidate cases include scoped native ANNOVAR and VEP
  execution with explicit source receipts. This does not clear a separate
  ready-input VEP deployment failure or validate every annotation route.
  The later AlphaGenome v1 entry implements explicit selected-candidate scoring
  with the user's own environment credential and remote-data consent. SDK/offline
  contracts have identified acceptance; authenticated API execution remains
  NOT_RUN and is not implied by that software acceptance.
- Phase 15 low-depth implementation has passed local and Linux engineering
  validation, including nested raw-FASTQ sampling and a tiny end-to-end run.
  The small-variant four-fraction, three-execution baseline matrix and paired difficult-region
  evaluations have completed independent audits (12 baseline runs, 24 regional
  evaluations). These are one-HG002 study results, not universal minimum-depth
  or sequencing-cost claims. Historical executed-source manifests are distinct
  from the newly prepared candidate's identity.
- RC5 adds hap.py environment isolation and success/failure runinfo permission
  protection. This is not an operating-system sandbox or a substitute for
  deployment-specific credential management.

## Dedicated SV depth and functional entries

The newer RC5 source provides [SV depth studies](sv_depth_study.md) and
[AlphaGenome selected-candidate scoring](alphagenome.md) as separate commands.
The SV study reuses declared BAMs with measured coverage, runs explicitly selected
read-SV callers and compares fixed-protocol INS/DEL query views. Bounded regional
native cases have been accepted under their own payload identities; the planned
12-unit by four-caller whole-genome depth matrix has **not run**. Completed
small-variant depth results cannot substitute for that SV experiment.

AlphaGenome requires the user's own service authorization and explicit data-sharing
consent. Offline validation is distinct from a successful authenticated service
request. The candidate-selection TSV rejects duplicate column names and rows
with missing or surplus fields before scoring; legitimate input records retain
their declared alleles. This input-contract change needs its own identified
candidate and does not alter old release assets or establish functional accuracy.

新增SV深度入口与AlphaGenome入口由 RC5 提供。小变异深度研究、SV区域工程
验证和完整SV全基因组矩阵是不同证据；后者尚未执行。AlphaGenome软件已实现、离线
接口已验收，真实API调用仍未完成。依赖工具和数据需用户配置，不能把本地源码描述
理解为历史公开包已经具有这些功能。

## Frozen public benchmark scopes

The HG002 small-variant benchmark uses PacBio HiFi CCS 11 kb data, GRCh38,
DeepVariant 1.10.0, hap.py 0.3.15, and chr1–chr22 confident regions. The read-SV
benchmark uses GIAB CMRG GRCh38 v1.00 and covers challenging medically relevant
genes only. CMRG values are not whole-genome SV performance. See
`benchmarking.md` and `known_limitations.md` before reusing any metric.

## Known product gaps

1. The development launcher routes explicitly configured assembly/assembly-SV,
   annotation/cohort, benchmark and review stages, plus per-sample harmonization
   from generated or declared ready sources. These are staged DAGs, not one
   all-module WGS acceptance. Current ready-PAV native route closure and other
   backend combinations still require their own successful exact-payload evidence.
2. Native uBAM routing requires pbmm2 and exact sample identity. CRAM and
   implicit uBAM-to-FASTQ conversion are intentionally not exposed.
3. Alignment preparation is sequential and separate from the calling DAG;
   distributed Slurm execution and a combined raw-input DAG remain pending.
4. Cross-sample SV clustering is intentionally unsupported.
5. Native VEP success is scoped to identified executed cases; a later failed
   route remains failed. AlphaGenome's opt-in software route is implemented;
   authenticated service execution and functional prediction accuracy remain
   unverified. The main calling DAG does not send candidates to the API.
6. Phase 15 study labels 5x/10x/15x/20x denote intended depths, not measured
   coverage. Analysis must use observed depth and identify nested technical
   replicates. Full-input repetitions are not independent biological samples.
   No general minimum depth or optimized caller profile is established, and
   Phase 15 is not part of published RC4. PCR validation remains future work.
7. Core Docker/Apptainer and Conda distribution claims remain separate from
   tool-specific real validation and must not be overstated.

## Documentation policy

Phase documents retain point-in-time design history. This page and the current
README describe the latest supported boundary. Historical audit failures remain
preserved in the private evidence repository and may be superseded only by an
explicit later validation record.

## Recommended next milestones

1. Approve the exact source/version, documentation and candidate-installation
   identity before publishing; keep
   unreviewed manuscripts, human data and sensitive provenance private.
2. Independently review and package the development launcher; distinguish tiny
   engineering acceptance from whole-genome or accuracy validation.
3. Integrate the completed Phase 15 study with transparent actual-depth,
   regional and biological-sample limitations; complete reviewer sign-off
   separately from screenshot generation. Additional biological datasets are
   needed for broader generalization, not to relabel existing technical runs.
4. Add an optional IGV-to-PCR evidence branch with at least three in-silico
   primer candidates, explicit specificity checks, predicted allele products,
   and independent wet-lab result records.
5. Extend validation to independent GIAB samples before making generalized
   performance claims.

## Current closeout disposition

The later installed development candidate has a separate HG002 regional case:
fresh DeepVariant, Sniffles2, TRGT and HiPhase producers followed by same-attempt
small/SV benchmark and TR technical review; reporting was invoked separately.
Native receipts, selected query/index bytes and output metrics were independently
read back. This closes that bounded GENERATED route only. It does not establish
TR/phase accuracy, a second biological sample, combined all-module WGS execution
or public release approval. Earlier execution records below remain separate.

The development source has archived real Linux acceptance for a two-sample,
seven-tool regional calling run, plus separate tiny raw-input routing tests.
These have different datasets and scopes and must not be combined into a claim
of raw-input, multi-sample whole-genome validation. Candidate documents now
describe the launcher separately from the earlier RC5 security-only build;
every candidate still requires an exact source/artifact manifest.

The first regional test completed a 17-job small/read-SV/TRGT/HiPhase DAG. A
separate, later two-sample test completed a fresh 13-job read-SV calling and
within-sample harmonization DAG, including eight caller jobs and two
harmonization jobs. Neither is an all-module WGS test or a new accuracy
benchmark. They must remain separate execution records.

The 24-locus human IGV adjudication is deferred. Review manifests and screenshots
remain available; no human decisions or manually validated locus results are
inferred. Manuscript preparation can proceed using the engineering and matched
truth-based results. Release approval still requires maintainer review of the
final source, documentation and package identity, followed by explicit release
approval and public-download validation.

The historical Sawfish benchmark has an unresolved output-sample identity
qualification; see [benchmarking](benchmarking.md). Its retained numeric row
must not be used as an unqualified accuracy comparison.

The development source additionally exposes an offline [review CLI](review_cli.md)
for IGV batch/pending-manifest preparation and explicit decision import to a
separate directory. It does not execute IGV, authenticate a reviewer, complete
the deferred 24-locus adjudication or rewrite caller outputs. This is a new
candidate payload included in the identified 15 October build. Its
dedicated checks remain separate from historical RC1–RC4 and earlier RC5
package validation; human adjudication is still not complete.

An additional opt-in [calling preflight](preflight.md) reports fresh-launcher
contracts, configured dependencies and explicitly requested native version/BAM
metadata probes. It creates no run directory or analysis artifacts, does not
start containers, and leaves unmeasured items NOT_CHECKED. The original doctor
is still foundation-only. Its offline-review/preflight payload has separate
Linux regression, package and bounded real-input acceptance; this does not
certify all workflows or complete human review. Later changes still require
their own source/artifact identity and validation.

The latest [within-sample read-SV integration](sv_harmonization_run.md) routes
generated caller outputs to separate per-sample Jasmine/evidence/source-comparison
artifacts. It never pools samples, adds a confidence rule or changes the frozen
benchmark data. Historical RC1–RC4 and earlier RC5 assets remain unchanged.
