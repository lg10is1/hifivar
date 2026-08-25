# Phase 2: QC and Alignment Foundation

**Phase 2 status: COMPLETE**

## Phase 2.1: QC Framework / Input QC Data Model — COMPLETE

The implementation was re-audited against the complete Phase 2 scope before
later alignment work. Its public API and lightweight/no-external-tool boundary
remain unchanged.

- Frozen `QCMetric`, `QCIssue`, `QCResult`, and `RunQCReport` models
- Explicit PASS, WARN, FAIL, and NOT_CHECKED status vocabulary
- Deterministic FAIL > WARN > PASS > NOT_CHECKED aggregation
- Lightweight FASTQ file-count, size, path, and compression metadata
- Lightweight BAM/CRAM file metadata and conventional index-presence checks
- Missing BAM/CRAM index represented as `ALIGNMENT_INDEX_MISSING` warning
- Missing/unreadable inputs retained as validation exceptions
- Ordered single- and multi-sample `AnalysisContext` QC
- Reference build, contig count, and checksum-availability run metadata
- Atomic UTF-8 JSON/YAML reports with overwrite protection and Unicode paths
- No implicit input/reference checksums, full FASTQ scans, or external tools

## Phase 2.2: Alignment Interface — COMPLETE

- Frozen tool-neutral `AlignmentRequest` model
- Frozen `AlignmentPlan`, `AlignmentResult`, and shell-free
  `AlignmentCommandPlan` models
- Scheduler-neutral `AlignmentResources` for threads, memory, and runtime
- Explicit `AlignmentTool` values for planned pbmm2 and minimap2 backends
- Explicit BAM/CRAM `AlignmentOutputFormat` with suffix validation
- FASTQ mapped to ALIGN plans; existing BAM/CRAM mapped to explicit REUSE plans
- Positive thread-count and input/reference/output collision validation
- Existing output refusal unless `overwrite=True` is explicit
- Deterministic, ordered mixed-input planning from `AnalysisContext`
- Side-effect-free output planning with Unicode path support
- Minimal `AlignmentBackend` protocol returning `list[str]` commands
- Generic backend command preview with `shell=False` reproducibility metadata
- Standard JSON/YAML-friendly request serialization
- Unit and Phase 1-to-2.2 integration coverage

## Phase 2.3: pbmm2 Integration — COMPLETE

- Dedicated `Pbmm2Wrapper`; no pbmm2 command exists in CLI/report modules
- Executable detection and parsed `pbmm2 --version` provenance
- Official FASTQ/FASTQ-FOFN to sorted BAM command shape with HiFi preset
- Stable sample/read-group metadata and explicit thread allocation
- pbmm2 automatic BAM indexing disabled so Phase 2.4 owns indexing explicitly
- All execution delegated to `CommandRunner` with `shell=False`
- Dry-run works without a locally installed pbmm2 and creates no files
- Optional timeout, redaction values, and tool-specific stderr log forwarding
- Input/reference revalidation, exit-code propagation, and expected BAM validation
- Deterministic, atomic ordered FOFN creation for multiple FASTQ inputs
- Default no-overwrite behavior, including output races after planning;
  explicit `overwrite=True` replaces only the exact requested output
- Alignment YAML settings for tool, format, resources, overwrite, preset, and logs
- Mock/fake-runner unit tests and effective-config dry-run integration coverage

## Phase 2.4: Alignment Post-processing/QC — COMPLETE

- Frozen `AlignmentArtifact` provenance for generated and existing alignments
- Explicit GENERATED/EXISTING source and COORDINATE/UNKNOWN sort-order states
- Dry-run results prohibited from becoming completed artifacts
- Lightweight BAM/CRAM output and optional index validation
- Conventional BAI, CSI, and CRAI discovery without header parsing
- Deterministic index strategy: CRAI for CRAM; BAI for ordinary BAM references;
  CSI when any reference contig exceeds the BAI 2^29 coordinate limit
- `AlignmentIndexRequest` refuses UNKNOWN sort order and silent overwrite
- Minimal dedicated `SamtoolsWrapper` limited to `samtools index`
- samtools executable/version/error/output handling through `CommandRunner`
- Dry-run indexing without a local samtools installation or filesystem writes
- Lightweight alignment QC for file size, source, declared sort order, and index
- Missing index and unknown sort represented as explicit QC warnings
- No BAM/CRAM record scan, header inference, sorting, coverage, or mapping metrics
- Unit and tiny post-processing integration coverage

## Phase 2.5: Phase 2 Integration — COMPLETE

- `Phase2Settings` loads the implemented pbmm2/BAM/indexing subset from YAML
- `run_phase2()` orchestrates ordered input QC, planning, alignment/reuse,
  output validation, indexing, alignment QC, and provenance
- FASTQ path: input QC → ALIGN → pbmm2 sorted BAM → samtools index → QC
- Existing BAM/CRAM path: input QC → REUSE original path → lightweight QC
- Existing BAM/CRAM never invokes pbmm2 and is never silently copied or realigned
- Existing alignment indexing is not rebuilt while sort order remains UNKNOWN
- Full dry-run plans pbmm2 and samtools commands without requiring either tool
- Frozen per-sample and run-level result models with deterministic QC aggregation
- Atomic UTF-8 JSON/YAML Phase 2 report with tool versions and command provenance
- Tiny mixed FASTQ/BAM/CRAM mock integration and all-FASTQ dry-run coverage
- Phase 0 Snakemake DAG deliberately remains infrastructure-only; no tool command
  is duplicated outside wrappers and `CommandRunner`

## Final Phase 2 behavior

| Primary input | Alignment action | Indexing action | Current checks |
| --- | --- | --- | --- |
| FASTQ | pbmm2 sorted BAM | samtools BAI/CSI | lightweight input/output/index QC |
| BAM | REUSE original BAM | retain discovered index; do not rebuild UNKNOWN sort | path, size, declared sort state, index presence |
| CRAM | REUSE original CRAM | retain discovered CRAI; do not rebuild UNKNOWN sort | path, size, declared sort state, index presence |

Lightweight means that HiFiVar checks paths, non-empty files, suffixes, planned
provenance, declared sorting, and conventional indexes. It does not parse
BAM/CRAM headers or records and does not calculate read counts, N50, QV, yield,
coverage, mapping quality, duplicates, or contig compatibility.

An unaligned BAM (uBAM) is still represented by the same `InputType.BAM` as a
coordinate alignment. Phase 2 therefore continues to REUSE BAM by explicit
contract rather than guessing from its suffix. Routing uBAM to pbmm2 remains
deferred as X-08 until a supported BAM-header parser or explicit uBAM input type
can distinguish it without risking accidental realignment of a valid BAM.

Real pbmm2 and samtools execution remains a Linux/HPC verification requirement.
Windows tests use dry-run and deterministic CommandRunner doubles. No variant
caller or Phase 3 module is implemented.

Phase 2.1 does not calculate read count, N50, QV, yield, GC, mapping statistics,
or coverage. Phase 2.2 does not parse BAM/CRAM headers, inspect or synthesize
read groups, create indexes, select tool-specific presets, construct concrete
pbmm2/minimap2 commands, execute external tools, register a CLI command, or add
a biological Snakemake rule. REUSE means only that alignment execution is
skipped; it does not claim that existing alignment headers or reference
compatibility have been validated.
