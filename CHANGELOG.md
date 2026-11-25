# Changelog

[简体中文](CHANGELOG.zh-CN.md)

## RC5 distribution / RC5 发行

The fixed RC5 release ties the Git tag, wheel, sdist, checksums and bilingual user documentation to one source. Runtime coverage is described per route in the deployment matrix.

RC5 使用同一 tag 记录源码、wheel、源码包、校验值与中英文说明。各分析路线的验证范围见部署文档。


- Development: route explicit per-sample ready SV sources through the launcher; bind actual VCF/index inputs to frozen external identities or current-attempt producer receipts. Include JSON policy bytes in runtime identity. This is workflow/provenance hardening, not new caller or accuracy evidence.

All notable changes to HiFiVar will be documented in this file.

## Development increment — 2026-11-03 (unreleased)

- Add explicit local `review snapshots` execution for existing external source
  plans, native BAM preflight, scientific close/overview/breakpoint windows,
  a local image index and compatible independent decision import. All IGV
  execution paths require complete PNG structure/CRC/decode checks with the
  optional review dependency. Human decisions remain pending after capture.
- Correct IGV 2.19.8's quote-retaining `genome` batch token and require per-image
  saved-session reference/locus/resource checks on the direct snapshot entry.
  Valid PNG files alone cannot certify correct capture. Whitespace reference
  paths remain explicitly unsupported in the tested 2.19.8 deployment and fail
  before execution. Whitespace output paths use safe session filenames and the
  recorded owned working directory without substituting a genome.
- Add `primers design --plan` for exact biallelic SNV PCR/Sanger and
  sequence-resolved simple INS/DEL flanking candidates. Record source allele,
  predicted REF/ALT products and genotype-dependent length classes; fewer than
  three candidates remain explicit. Genome-wide specificity, nearby sample
  variation and wet-lab validation are not assessed by this first increment.
- Existing callers, benchmark results and historical release assets are unchanged.

## Development increment — 2026-10-30 (unreleased)

- Add an explicit four-read-caller SV depth/support post-processing entry. Validate
  complete parent-unit receipts, merge once per unit, derive support ≥1–4 from the
  same union and independently benchmark each query against the fixed original
  truth protocol. Preserve caller/union files, observed-depth labels, both TP
  sides and failed/empty-refinement boundaries. No caller algorithm change,
  new native depth result or universal support-threshold recommendation.

- Reject ambiguous AlphaGenome candidate TSVs: duplicate column names and missing
  or surplus data fields fail before scoring. Preserve valid field ordering and
  quoting; report invalid rows without echoing their contents.
- Clarify current SV-depth and AlphaGenome documentation: software implementation,
  bounded regional acceptance, full WGS depth studies and authenticated API
  execution remain separate evidence. Earlier candidates and release assets are unchanged.

## Development increment — 2026-10-29 (unreleased)

- Handle Truvari 5.4.0 perfect-match zero-refinement exits only with verified command, parameter, partition, summary and log evidence. Record NOT_NEEDED with native bench metrics; do not synthesize refined files or relabel historical failed attempts. Preserve strict rejection of partial or unproven output bundles.

- Add a dedicated manifest-driven read-SV depth study entry that reuses existing depth/replicate BAMs, preserves per-caller failures and compares INS/DEL query views using a fixed original-truth protocol. No SV depth result or cost claim is implied by implementation.
- Add opt-in official AlphaGenome 0.10.0 scoring with user-owned environment credentials, local REF/full-context validation, explicit remote-data consent and structured modality/tissue/gene score output. No legacy credential or automatic whole-VCF invocation.
- Preserve complete receipts/failed attempts, reject ambiguous YAML authorization and protect all historical release/manuscript artifacts. These changes do not replace RC1–RC4 assets.


## [0.1.0rc5] - source/package candidate

This candidate includes the post-RC4 Phase 15 experiment framework, standalone
alignment-input development, bounded calling/review entries and subprocess
privacy hardening listed below. It is not a new calling algorithm, a universal
minimum-depth claim or a sequencing-cost result. RC1-RC4 artifacts remain unchanged.

### TR visualization consumer

- Development: accept only the exact TRGT version reports `5.1.0` and the
  observed deployment build `5.1.0-ec66463`; reject other suffixes or conflicting
  reports. Preserve raw version logs and record both normalized `tool_version`
  and complete `tool_version_reported`. This does not retroactively accept an
  earlier failed native attempt; new package/native acceptance remains separate.
- Development: add the explicit `review tr-plot --plan` consumer for one existing
  TRGT 5.1.0 locus. Strict YAML binds seven input files and the parsed plan;
  only fresh-directory allele/motifs SVG is supported. Record separate version
  and native logs, reject active/external SVG structures, verify receipt/output
  identities and retain failures without overwrite, cleanup or automatic retry.
  This is not TR genotyping, parental-phase validation, manual review or new TR
  accuracy evidence, and is not retroactively included in historical releases.

### Annotation output safety

- Development: validate VEP's default tabular header and stream row structure
  before publishing staged output. Reject nonempty error text, malformed tables
  and unreadable UTF-8 without replacing previous TSVs. Preserve legal
  header-only results, extra columns and native bytes. This is format validation,
  not new annotation execution or annotation-accuracy evidence.

- Development: preserve existing VEP TSVs until a newly executed, nonempty
  staged result is ready. Publish with atomic replacement only when overwrite
  was explicitly requested; otherwise use exclusive same-filesystem publication
  and reject late-created outputs. Keep native sidecars/failed staging inside
  the owned output-parent directory. Protect CommandRunner's actual implicit
  stdout companion and recheck final/staged/resource aliases after version and
  native execution. This does not change VEP annotations or repair an unrelated
  Perl/PATH deployment failure. Exact-package native annotation acceptance
  remains separate; see [VEP output safety](docs/vep_output_safety.md).

### Output and catalog validation

- Development: require explicit, unambiguous TRGT locus IDs and POS/END before
  cohort indexing; reject missing, malformed or reversed endpoints without
  coordinate imputation. Preserve native missing VCF IDs, encoded TRID text and
  genotype values. Type SQLite lifecycle failures while retaining the primary
  exception or interrupt and existing successful tracks. This is cohort input
  validation, not new repeat calling, clinical interpretation or accuracy evidence.

- Development: disambiguate repeated SV cohort source keys with reserved,
  deterministic source-record ordinals after disk-backed key counting. Preserve
  unique legacy keys, TSV columns, source IDs/fallbacks, GT and missing states;
  keep duplicate records rather than merging them. Record key-policy and
  uniqueness metrics. Reject extra columns in one-sample SV/TR VCF records.
  This is not cross-sample clustering, caller modification or new benchmark
  evidence; prior packages and scientific results are unchanged.
- Development: validate SV/TR cohort tables in owned staging before exclusive
  publication, verify published bytes and conservatively roll back unchanged
  owned files on captured errors. Preserve failed diagnostics and use a unique
  logged TR scratch database. Table schemas, native source keys, genotypes and
  missing-state semantics are unchanged; this is not physically atomic
  two-file publication or a fix for repeated native SV identifiers.
- Development: reject TRGT output/log paths that alias protected inputs,
  indexes, references, catalogs or executables before probing or running tools.
  Reject catalog BED ends beyond the declared reference contig length without
  truncating the catalog. Calling parameters and variant interpretation do not
  change.
- Development: preflight the full release-bundle write set, including large-file
  pointer names, before writing reports. Reserve generated metadata namespaces,
  reject nonportable/escaping paths and file aliases, and protect input bytes
  even with explicit overwrite. This does not claim concurrent-writer isolation
  or whole-directory transactional publication. See
  [output path safety](docs/output_path_safety.md).

### Documentation reconciliation

- Add a complete bilingual indexed-BAM batch guide and editable config/sample/
  generated-source examples for small variants, four read-SV callers, TRGT,
  joint phasing and three cohort tracks. Reconcile native joint execution wording
  with the separately identified regional cases, without extending it to WGS,
  independent biological samples, accuracy or historical assets. No runtime
  code/resource/algorithm change; rebuilt artifacts need their own identities.
- Add a bilingual capability/evidence matrix and correct outdated launcher,
  benchmark and VEP absence claims in current documentation. Preserve dated
  payload/release boundaries and separate scoped native execution from accuracy,
  human review and complete all-route acceptance. No caller or matching-policy
  change; documentation rebuilds require newly measured artifact identities.
- Include the Chinese README explicitly in the source distribution so its
  bilingual navigation does not point to a missing root document. Keep private
  handoffs and test work directories excluded; no runtime payload change.

### Launcher and review additions

- Development: add explicit optional joint HiPhase inputs for one pbsv/Sawfish SV
  track and/or TRGT alongside small variants. Preserve separate phased VCFs,
  default small-only commands, actual producer dependencies, per-output indexing
  and receipt/source identity bindings. Stream ordered record keys and unordered
  GT-allele preservation checks in joint mode. Global realignment is configured
  explicitly or selected by the documented small/joint `auto` behavior. This is
  not new calling, phase-accuracy evidence or a retrospective release change;
  exact-package native joint acceptance remains a separate gate.

- Development: add opt-in benchmark.source_plan to the owned launcher for
  existing small/read-SV/assembly-SV/TR comparison engines. Freeze truth/regions,
  resolve one declared truth sample's native query per enabled category from
  external identities or selected-attempt producers, and seal commands, controls,
  metric artifacts and report bindings. Reject gVCF/derived-BND substitution,
  changed inputs and missing receipts. Reporting stays explicit; no new calling,
  matching policy, universal accuracy, or historical release change.
  Declare native Truvari directories as directory outputs so the scheduler
  does not pre-create them; preserve its strict no-overwrite contract and seal
  each required native file in the execution receipt.

- Development: route explicit review.source_plan through a final review-only
  consumer DAG, after producers and any annotation/cohort consumer. Bind actual
  BAM preflight/alignment handoff and selected-attempt source receipts; seal IGV
  commands, pending manifests, provenance and every expected endpoint PNG. Reports
  distinguish technical execution from NOT_REVIEWED human status. Default-disabled,
  standalone selection mode and historical assets are preserved; report remains
  a separate explicit CLI command. No new calling or accuracy policy.

- Add a bounded, read-only review source library for explicit external records
  or selected-attempt producer receipts. Preserve exact SITE/raw-ID, independent
  TRID and complete BND bundle bindings without running IGV or rewriting calls.
  This library-only checkpoint preceded the later launcher execution bridge
  listed above; see docs/review_sources.md.

- Reject IGV batch, staging, screenshot and stdout/stderr paths that alias an
  input/reference/index, selected BND bundle member, optional TR visualization
  or executable. Check resolved and hardlink identities plus CommandRunner's
  implicit stdout log before any write/version command and recheck after the
  version probe. Preserve old screenshots when version detection fails. Also
  protect Phase 10 JSON/YAML/TSV manifests and its fixed TSV staging path,
  with a recheck after IGV before manifest serialization. This
  is write-path protection, not a sandbox, calling change or human validation.

- Add a library-level, validation-only contract for explicitly declared existing
  haplotype FASTAs. Strict sample/role/SHA256 manifests, full streaming FASTA
  checks and caller output/log alias guards support standalone PAV/SVIM-asm
  requests without inventing a new hifiasm execution. Generated and external
  assemblies cannot be mixed in one request. Opt-in `assembly_sv.input_manifest`
  routes this mode through the launcher/DAG with an attempt-bound validation-only
  input record, new caller receipts and a separate report track; hifiasm stays
  disabled. Preserve the generated FASTQ route, source declarations and native
  caller versions without inventing GFA paths or producer executions. Generated
  consumer stages do not reactivate the ready producer. See `docs/ready_assemblies.md`.

- Allow explicit generated annotation/cohort selectors for a new PAV/SVIM-asm
  output when its inputs are ready FASTAs and hifiasm is disabled. Reject absent
  or conflicting assembly input modes and retain enabled-caller, ready-input
  validation and current-attempt output receipt requirements. This fixes source
  intent routing; it does not add calling, joint SV genotyping or accuracy claims.

- Extend opt-in native preflight to enabled hifiasm, SVIM-asm (including its
  alignment/finalization dependencies), VEP and GLnexus/bcftools through existing
  wrapper version methods. Preserve command-sourced/version-only evidence,
  disabled branches, bounded credential-free probes and per-group failures.
  PAV/ANNOVAR declarations remain required NOT_CHECKED, not native versions.
  No calling, assembly, database download or scientific decision changes.

- Support explicit `TRID:<locus>` overlap selection for the independent TRGT
  category when raw VCF ID is `.`. Preserve raw IDs, record hashes and ordinal
  provenance; annotate only the repeat interval, excluding its VCF left anchor.
  Reject ambiguous/missing TRIDs and changed sources without rewriting calls.
  Ordinary VCF-ID and BND selection remain separate. This is annotation routing,
  not a new repeat caller, accuracy result or clinical interpretation.

- Bind alignment preparation to the exact launcher attempt and freeze the
  actual calling YAML, ordered BAM sample sheet and selected BAM/index hashes.
  Distinguish new alignment from this-attempt validation of external or prior
  BAM reuse. Legacy completion records do not acquire a retrospective producer
  identity. This is execution provenance, not a new calling or accuracy claim.

- Disable independently enabled harmonization and producer-only representation
  policies in the consumer configuration. Preserve original frozen producer
  policies and existing-bundle overwrite guards; do not rerun integration when
  a downstream stage consumes its same-attempt ready output.

- Development: allow explicitly selected same-attempt harmonization `bnd_vcf` annotation roles with their complete native/ledger bundle. Bind every member and copied canonical native output to producer receipts; keep `vcf` separate and reject implicit fallback, missing bundle, disabled policy or cross-attempt substitution. External source plans may explicitly freeze a ready BND bundle. Standalone harmonization does not require read calling enabled. This is provenance/workflow integration, not BND accuracy or a new caller.


- Development: check standalone harmonization dependencies independently of read
  calling; retain required NOT_CHECKED for uncovered native probes. Stop guessing
  assembly haplotype roles from filenames; bind generated roles/FASTA identities
  to the current assembly receipt and retain unknown roles for unbound sources.
  Source-input role metadata is not per-variant phasing accuracy.

- Development: route explicit per-sample ready SV sources through the launcher; bind actual VCF/index inputs to frozen external identities or current-attempt producer receipts. Include JSON policy bytes in runtime identity. This is workflow/provenance hardening, not new caller or accuracy evidence.

- Add a default-disabled, policy-pinned SV source identity view to Phase 9.
  Duplicate/missing IDs receive reversible analysis keys without changing any
  non-ID record field or canonical caller file. Evidence distinguishes original
  IDs, unique source keys, actual analysis IDs and representation files; optional
  workflow outputs and attempt receipts bind the retained identity bundle.
  Ambiguous MATEID/PARID references fail rather than being guessed. This is an
  unpublished provenance increment, not a caller, normalization, confidence
  method or a completed new installed/native acceptance claim.

- Wire the explicit BND representation policy into Phase 9 configuration and
  per-sample workflow outputs, default-disabled. Enabled runs preserve native
  concordance inputs and bind the ten-file reversible bundle to attempt
  receipts. Dry-runs record plans only; pre-existing bundle directories fail
  before analysis even when native overwrite is requested. This unpublished
  source increment does not redirect annotation/benchmark consumers, validate
  BND accuracy, or update historical release assets.

- Add an explicit, separate BND representation API and versioned policy resource
  for Jasmine 1.1.5 VCFv4.4 outputs. It preserves native binary copies, all
  non-BND record bytes and every protected allele/genotype/source field; only
  BND END/SVLEN tokens are externalized into a mandatory reversible record
  ledger. A two-endpoint junction table keeps local/mate coordinates distinct.
  This is not automatically selected by the launcher, annotation or benchmark,
  and is not an accuracy or confidence claim. Native acceptance is a separate
  gate, not established merely by unit tests.

- Register producer rules only when their own module is enabled. Ready sources
  in the workflow output namespace must not reactivate a disabled caller during
  consumer DAG construction. The read-SV wildcard also excludes disabled
  individual callers. External-only harmonization keeps its independent switch.
  This corrects a native generated-source handoff failure; fresh installed
  end-to-end acceptance of the corrected increment remains required.

- Unpublished downstream-source development separates explicit source plans,
  producer receipts and a consumer-only second DAG. Consumer receipts bind the
  actual VCF/index selection, stage, manifests, configuration and materialization
  identity. Missing receipts are not inferred from old outputs. Terminal failed
  attempts can deliver failed-artifact diagnostics without marking those units
  complete. Cross-attempt receipt reuse and new installed Linux/native acceptance
  remain open; this does not alter historical release assets or calling methods.

- Route opt-in generated read-SV calls through the existing per-sample Phase 9
  Jasmine/evidence/Truvari workflow in `hifivar run`. Validate its scheduler
  budgets and dependencies, and restrict each wildcard bridge to its own sample.
  External/assembly inputs and cross-sample SV clustering remain outside this
  launcher. Default-disabled behavior, caller algorithms and historical assets
  are unchanged; support counts remain descriptive evidence, not confidence.
  Opt-in preflight additionally blocks an indexer below the existing shared
  VCFv4.4 compatibility floor before callers execute; it does not upgrade tools.

- Add opt-in `doctor --preflight` for fresh calling-launcher metadata checks
  and explicitly requested bounded native version/BAM probes. FAIL and
  NOT_CHECKED remain distinct; no analysis, input repair, installation,
  container pull/start, automatic resource tuning or accuracy claim occurs.
  The no-option foundation doctor and historical package payloads are unchanged.

- Add standalone offline `hifivar review prepare` and `review import-decisions`
  entries for the existing Phase 10 module. Preparation leaves all targets
  NOT_REVIEWED; user-supplied decisions retain reviewer, timestamp and comments
  in separate, non-truth evidence. No IGV/caller process is started, no raw
  calls are rewritten, and no manual review is claimed independently verified.
  This extension is absent from RC1–RC4 and earlier RC5 candidate packages.

- Add optional, explicit `small.region` for bounded DeepVariant calling through
  the wrapper, Phase 3 API and batch workflow. Exact reference interval checks,
  scope provenance and resume identity are retained; the default command is
  unchanged. This development-only option is not a Phase 15 setting, a new
  caller, an accuracy claim, or a replacement of earlier release assets.

### Documentation and validation scope

- Distinguish the initial security-only RC5 package from the separately identified
  15 October development payload and any later documentation rebuild. Record
  source manifests and artifact SHA256 rather than relying on a shared version
  string. Available release tags/assets are listed on GitHub Releases; candidate
  preparation alone does not establish publication or pass the final release gate.
- Document separate successful regional launcher tests: a 17-job, two-sample
  small/read-SV/TRGT/HiPhase run and a later 13-job, two-sample read-SV calling plus
  within-sample harmonization run. These are engineering acceptance records, not
  a combined experiment, an accuracy comparison or all-module WGS validation.
- Reconcile offline-review and preflight availability with their completed
  bounded Linux/package checks. No human review or unsupported branch is promoted
  to validated status by these documentation updates.

### Added

- Add Phase 15 deterministic low-depth experiment models, explicit
  profile/depth/replicate planning, preregistered per-class accuracy thresholds,
  conservative replicated assessments, nested per-replicate sampling across
  depths, and JSON/YAML result serialization.
- Add a CommandRunner-only samtools aligned-BAM subsampling/index boundary plus
  coordinate-sort and quickcheck validation, and an optional config-driven
  per-replicate Snakemake materialization branch. This introduces no universal
  minimum-depth claim, caller change, or Phase 16 behavior.
- Add a multi-file, streaming raw HiFi FASTQ mode that hashes PacBio movie/ZMW
  molecule identities once per replicate, fans records into exact nested depth
  subsets, preserves distinct replicate seeds, and keeps all source reads
  immutable. The completed single-HG002 baseline and regional studies retain
  their own executed-source, input and observed-depth identities.
- Add an opt-in Phase 15 raw FASTQ-to-benchmark branch that reuses the existing
  pbmm2, samtools, DeepVariant and hap.py wrappers, isolates each replicate's
  temporary storage, measures observed coverage and records benchmark
  observations without making a minimum-depth or cost-reduction claim.

### Fixed

- Recover only used missing Jasmine FILTER declarations verbatim from explicit
  source headers. Unknown or conflicting missing definitions fail closed;
  raw output, body fields and source calls are preserved. This header-only
  correction does not resolve the separately identified upstream BND END
  representation issue, establish warning-free output or validate a full merged
  callset against truth.

- Isolate hap.py version and benchmark subprocess environments using an explicit
  runtime allowlist; restrict newly written runinfo permissions after success
  and failure. This is not an operating-system sandbox; local permissions and
  environment isolation do not invalidate copied credentials.

- Parse `samtools coverage` from either CommandRunner's in-memory stdout or its
  persisted sibling stdout log, so Phase 15 coverage measurement remains
  compatible with tool-specific stderr logging.

### Scope clarification

- The source also contains the standalone minimap2 wrapper, explicit FASTQ/BAM
  input routing, read-group escaping, transactional BAM/index protection and
  structured command-chain provenance. The bounded calling launcher now supports
  opt-in sequential FASTQ/native-uBAM preparation; CRAM, a unified raw-input DAG
  and all-branch routing remain deferred.
- Phase 15 baseline and regional experiments have completed within their
  documented single-HG002 scope. Candidate packaging does not retroactively
  identify the historical executed source as this new wheel.

## [0.1.0rc4] - 2026-09-02

This release candidate contains only the independently Linux-validated
Jasmine/HTSlib compatibility correction and PAV 2.4.6 SV-only assembly-evidence
handoff. It adds no caller and no new confidence or biological-truth rule.

### Fixed

- Restore missing Jasmine INFO and FORMAT declarations only from immutable
  source VCF headers before indexing. This preserves typed BND fields such as
  `CIPOS` for downstream Truvari parsing and rejects identifiers that have no
  authoritative definition instead of inventing metadata.
- Preserve Jasmine's original VCF version declaration. VCFv4.4 output now
  requires the Linux-validated tabix/HTSlib 1.23.1 compatibility floor rather
  than silently rewriting the header for HTSlib 1.21.
- Preserve PAV 2.4.6's native mixed SNV/indel/SV VCF as immutable raw evidence
  while producing the assembly-SV artifact with PAV's own VARTYPE boundary:
  INV records and INS/DEL records with absolute `SVLEN >= 50`. The streaming
  derived output retains PASS and non-PASS records and all native fields before
  explicit bgzip/tabix finalization.

### Known limitations

- The PAV SV-only selection contract is validated only for PAV 2.4.6; other
  PAV releases remain unsupported until independently verified.
- VCFv4.4 finalization requires tabix/HTSlib 1.23.1 or newer under the
  validated deployment contract.

## [0.1.0rc3] - 2026-08-28

This release candidate contains only the independently Linux-validated Phase 9
evidence-source attribution correction and its regression coverage. It adds no
caller, cross-sample clustering, scientific confidence rule, or VCF rewrite.

### Fixed
- Derive per-sample Jasmine evidence membership from the validated `SUPP_VEC`
  input-order contract instead of guessing caller identity from heterogeneous
  native `IDLIST` formats. This corrects caller/source counts and paths in the
  derived SV evidence TSV while leaving caller and harmonized VCFs unchanged.

## [0.1.0rc2] - 2026-08-27

This release candidate contains deployment/runtime hardening and output
validation fixes only. It adds no caller and changes no scientific decision
rule.

### Added

- Adopt Apache License 2.0 for the public source distribution.
- Add a Conda/Mamba source environment and local noarch Conda recipe for the
  HiFiVar core plus Snakemake.
- Add non-root Docker and Apptainer core definitions without bundling external
  callers, references or licensed databases.
- Add task-oriented installation, quick-start, Linux/HPC, container, output and
  troubleshooting documentation plus a minimal existing-BAM example.
- Add contribution, security and citation metadata for public distribution.

### Changed

- Bound the optional workflow dependency to tested Snakemake major versions 8
  and 9.
- Isolate each DeepVariant sample in its own data-disk temporary directory and
  pass the same `TMPDIR` through native, Docker, and Apptainer execution.
- Serialize DeepVariant samples by default through the global
  `deepvariant_slots` resource; sites may raise `small.max_concurrent_samples`
  only after validating storage, file-descriptor, CPU, and memory limits.
- Require an explicit GLnexus `cohort.small_variants.memory_gb` when that track
  is enabled, and forward the configured cap to both Snakemake and GLnexus.
- Validate GLnexus cohort identity by exact sample-set equality while retaining
  both declared and physical VCF sample order in provenance. Sample-level QC is
  keyed by header sample name, and the raw cohort VCF is never reordered.
- Document the `0.1.0rc2` execution boundary: packaged small/SV/TR workflow
  branches consume indexed BAM/CRAM, while unified FASTQ-to-calling execution
  is not yet exposed through the CLI/DAG.

## [0.1.0rc1] - 2026-08-24

### Added

- Add Phase 14 five-state final run reports with redacted JSON/YAML provenance
  and static offline Markdown/HTML summaries across all analysis categories.
- Add explicit release/reproducibility bundle creation, selected-result copying,
  large primary-data pointers, optional checksums, command redaction, and stable
  reports/manifests/configs/provenance directory conventions.
- Package the modular Snakemake workflow in wheel/sdist data files and add an
  installed-resource locator plus out-of-tree clean-install tests.
- Add lightweight Linux/Windows CI, installation/deployment documentation, a
  validated-versus-pending tool matrix, and a human-authorized release checklist.

- Add Phase 13 benchmark contracts with pinned truth/region provenance,
  explicit PASS/PARTIAL/NOT_RUN/UNSUPPORTED states, atomic JSON/YAML/TSV
  manifests, and immutable-query scientific safeguards.
- Add a CommandRunner-only hap.py adapter with deterministic dry-run commands,
  explicit engine/resources/stratification inputs, and column-name-based
  SNP/INDEL precision/recall/F1 parsing.
- Extend the existing Truvari adapter compatibly with explicit benchmark
  thresholds, confident-region input, official summary parsing, and streaming
  SVTYPE/user-defined size-bin summaries that exclude unsafe BND/complex
  length classification.
- Add exact-catalog TR genotype/allele comparison and independent optional
  Phase 13 Snakemake benchmark tracks for small, read-SV, assembly-SV, and TR.

### Fixed

- Adapt the PAV 2.4.6 workflow boundary to run its default DAG, retain the
  upstream `{sample}.vcf.gz`/TBI pair, and atomically bridge those validated
  files to HiFiVar's deterministic Phase 8 output paths.
- Run text-only Jasmine launchers without a shebang through `bash`, stream
  BGZF inputs into ordered uncompressed work files for Jasmine 1.1.5, and
  retain all original inputs.
- Preserve Jasmine raw output, restore missing FORMAT declarations only from
  source VCF headers, and deterministically external-sort records by declared
  contig order and position before BGZF/tabix finalization.
- Detect Truvari 5.4 with its supported `truvari version` subcommand while
  retaining a legacy `--version` fallback.
- Accept DeepVariant 1.10 gVCFs that omit the legacy `NON_REF` ALT declaration
  by validating native DeepVariant version, RefCall, and MIN_DP/MED_DP header
  markers instead.
- Remove the unsupported `--threads` option from `bcftools sort` in the
  Phase 5 TRGT finalization chain while retaining supported threading for
  `bcftools index` and samtools.
- Reject partial TRGT results when the caller reports a zero-exit
  `[ERROR] - Locus processing:` diagnostic on stderr, preserving the raw
  products and command log for investigation.
- Pre-create writable Docker/Apptainer bind sources before building the
  DeepVariant container command.
- Persist stdout to a sibling `*.stdout.log` whenever a tool caller supplies
  only `stderr_path`, preserving DeepVariant stage diagnostics without changing
  binary stdout destinations.
- Reject real DeepVariant launch when the POSIX open-file soft limit is below
  4096, with an actionable 65536 recommendation and shard guidance.
- Quarantine completed DeepVariant VCF/gVCF/index products when validation
  fails so workflow failure cleanup cannot erase diagnostic artifacts.

### Added

- Phase 12 cohort/order/reference/state/provenance contracts that keep small
  variants, SVs, and TRs independent and never equate a missing record with
  homozygous reference.
- A CommandRunner-only GLnexus 1.4.1/bcftools adapter with deterministic gVCF
  order, BCF/CSI and BGZF-VCF/TBI outputs, exact sample validation, streaming
  QC, dry-run, overwrite protection, and failure propagation.
- Lossless Phase 9 source-native SV cohort tables, catalog-consistent TRGT
  locus/sample matrices, independent optional Snakemake tracks, and a pinned
  GLnexus Linux environment plus validation handoff.

- Phase 11 tool-neutral small-variant/SV/TR annotation artifacts with explicit
  source-variant, reference-build, database-version, command, and tool-version
  provenance while retaining raw caller VCFs unchanged.
- Independent CommandRunner-only ANNOVAR and Ensembl VEP adapters with
  deterministic commands, offline/local database boundaries, dry-run, output
  validation, overwrite protection, and failure propagation.
- Explicit SV/TR gene, exon, regulatory, repeat, and segmental-duplication BED
  overlap evidence that preserves original source breakpoints and SV types.
- An explicit-selection-only AlphaGenome-compatible backend protocol that
  retains model/modality/source provenance and states that functional impact is
  neither call confidence nor truth.
- Optional Phase 11 Snakemake annotation rules plus mock end-to-end and dry-run
  coverage; annotation failure cannot modify upstream calling artifacts.
- Phase 10 review contracts for explicitly selected small variants, structural
  variants, tandem repeats, phased variants, and harmonized evidence without
  changing source VCF or alignment artifacts.
- Configurable variant-centred review loci, including anchor-based insertion
  windows and two-locus BND handling, with stable screenshot paths and a
  separate TRGT-visualization metadata boundary.
- CommandRunner-only IGV batch planning/execution, dry-run previews, output and
  overwrite validation, and JSON/YAML/TSV manual-review manifests whose statuses
  never imply truth, pathogenicity, or clinical classification.
- An optional, config-driven Phase 10 Snakemake branch and fake-IGV end-to-end
  coverage; variant calling remains independent of manual review.
- Initial repository skeleton.
- Python package metadata and base exception hierarchy.
- Shared HiFiVar console and UTF-8 file logging infrastructure.
- Layered YAML configuration, validation, presets, and effective-config output.
- Safe external command execution with dry-run, redaction, and file output.
- Argparse CLI with packaged config resources, config actions, and doctor.
- Lightweight, streaming input validation and SHA256 checksum infrastructure.
- Modular Snakemake smoke DAG with effective-config and resource conventions.
- Phase 0 end-to-end integration, packaging, wheel-resource, and smoke checks.
- Phase 1 reference, sample/input, sample-sheet, pedigree, and cohort data models.
- Phase 1 run-level `AnalysisContext` cross-validation and config construction.
- Versioned, redacted, atomic JSON/YAML `RunManifest` provenance snapshots with
  optional input checksums.
- Lightweight input QC models, deterministic run-level aggregation, and atomic
  JSON/YAML QC reports for FASTQ/BAM/CRAM filesystem metadata.
- Tool-neutral, immutable alignment requests, deterministic multi-sample
  planning, and a minimal backend command-construction protocol.
- Mixed FASTQ/BAM/CRAM alignment plans, resource specifications, explicit reuse,
  no-overwrite policy, and dry-run command/result models.
- Dedicated pbmm2 wrapper with versioning, sorted BAM output, multi-FASTQ FOFN,
  read-group metadata, output validation, and CommandRunner-only execution.
- Alignment artifacts, explicit sort-order provenance, safe BAI/CSI/CRAI index
  planning, a minimal samtools index wrapper, and lightweight alignment QC.
- Complete Phase 2 orchestration and atomic JSON/YAML provenance for input QC,
  alignment/reuse, indexing, and alignment QC.
- DeepVariant PACBIO request/result models and native, Docker, and Apptainer
  runtime isolation with CommandRunner-only execution.
- Indexed BAM/CRAM single-sample calling with deterministic small VCF/gVCF/TBI
  outputs, BGZF/header/sample/contig validation, and strict overwrite handling.
- Config-driven modular DeepVariant Snakemake rule and narrow wrapper bridge.
- Complete Phase 3 AnalysisContext/Phase 2 alignment handoff and atomic
  JSON/YAML command, version, runtime, and artifact provenance.
- Independent Sawfish, Sniffles2, pbsv, and cuteSV wrappers with deterministic
  resources, dry-run commands, version detection, and explicit input contracts.
- Separate per-caller structural-variant artifacts with lightweight
  VCF/BGZF/TBI validation, modular Snakemake targets, and Phase 4 provenance.
- Default cuteSV genotype generation and ownership-marked work-directory
  replacement that refuses to remove directories not created by HiFiVar.
- Platform-aware path identity and a shared semantic secret-redaction policy
  for CLI output and run manifests.
- Phase 5 reference-specific `TandemRepeatCatalog`, TRGT request/result models,
  deterministic shell-free genotype commands, version detection, and strict
  indexed-BAM input handling.
- Official TRGT post-processing path using bcftools/samtools to sort and index
  `sample.tr.vcf.gz` plus `sample.tr.spanning.bam`, with lightweight TRID,
  MOTIFS, STRUC, sample, contig, BGZF, and TBI validation.
- Config-driven modular Phase 5 Snakemake rule and atomic JSON/YAML provenance
  recording catalog, karyotype, commands, resources, tool versions, runtimes,
  and final TR artifacts.
- Phase 6 indexed-alignment/small-variant phasing contracts, a
  CommandRunner-only HiPhase wrapper, deterministic phased VCF/TBI validation,
  modular Snakemake integration, and atomic JSON/YAML run provenance.
- Explicit Phase 6 sample/reference compatibility, index, overwrite, dry-run,
  command, tool-version, runtime, PS-format, sample, and contig checks.
- Phase 7 FASTQ-only hifiasm request/result models with ordered multi-file input,
  deterministic primary and haplotype GFA/FASTA outputs, and strict overwrite
  protection.
- Streaming atomic GFA-to-FASTA conversion that preserves raw hifiasm evidence,
  plus modular Phase 7 Snakemake integration and atomic JSON/YAML provenance.
- Phase 8 independent PAV workflow-adapter and SVIM-asm wrapper branches with
  explicit haplotype inputs, deterministic artifacts, validation, dry-run,
  failure propagation, and provenance.
- Phase 8 modular Snakemake integration preserves PAV and SVIM-asm raw evidence
  without harmonization.
- Phase 9 streaming read/assembly SV evidence model, conservative native-field
  boundary, Jasmine primary clustering, and Truvari comparison-only support.
- Phase 9 deterministic harmonized VCF, evidence table, partial-caller statuses,
  scientific safeguards, Snakemake integration, and run provenance.

### Unpublished junction-aware consumer hardening

- Bind explicitly selected BND annotation/review to actual source records, both ALT endpoints and source identities; never treat END as an event span.
- Retain single/virtual endpoint statuses and unresolved mate QC, without guessed links or confidence claims.
- Accept explicit ready BND-derived annotation/review only with complete policy/reference/ledger round-trip validation and attempt-bound bundle inputs.
- Preserve ordinary SV/TR overlap, raw VCF bytes and human NOT_REVIEWED states. No new caller or scientific decision rule.

### Configurable caller-support selection (development source)

- Add `sv.harmonization.minimum_caller_support` (default 1), bounded by actual
  per-sample caller sources. Higher thresholds produce a fresh, separate bundle
  with exact retained VCF/evidence rows, source order and checksums.
- Preserve complete native Jasmine outputs, BND roles and concordance inputs;
  selection does not change clustering or calibrate correctness probabilities.
  Historical release assets and earlier candidates remain unchanged.
