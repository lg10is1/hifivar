# Changelog

All notable changes to HiFiVar will be documented in this file.

## Unreleased

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
- Document the `0.1.0rc1` execution boundary: packaged small/SV/TR workflow
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
