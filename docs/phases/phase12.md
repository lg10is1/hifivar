# Phase 12 — Cohort / Multi-sample Analysis

[中文阶段说明](../zh_CN/phases/phase12.md)

**Status: implemented; validation remains route- and candidate-specific.**

This page describes the cohort contracts. Its original implementation checkpoint
used Windows mocks; later identified Linux cases exercised GLnexus/bcftools,
sample-order handling and the independent SV/TR summaries. Those cases do not
certify every later package or a complete multi-sample whole-genome DAG.
See [current capabilities](../capabilities.md) and [project status](../project_status.md)
for the current development/release boundary. Benchmarking is documented
separately in [Phase 13](phase13.md); the old phase checkpoint is not its status.

## Scope and scientific invariants

Phase 12 has three independent tracks:

```text
DeepVariant gVCFs -> GLnexus -> cohort.small.bcf + cohort.small.vcf.gz
Phase 9 harmonized SVs -> source-native cohort SV sites + sample-state matrix
TRGT VCFs with one catalog -> cohort TR loci + sample-state/genotype matrix
```

Small variants, SVs, and TRs are never merged. Track state, outputs, commands,
versions, metrics, and provenance remain independent. `run_phase12()` isolates
expected HiFiVar failures per track.

States are `CALLED`, `NO_CALLS`, `NOT_RUN`, `FAILED`, `DISABLED`,
`MISSING_INPUT`, and record-level `NOT_OBSERVED`. A missing VCF record is not a
homozygous-reference genotype. `FAILED`, `DISABLED`, and `NOT_RUN` do not mean
`NO_CALLS`.

## Cohort contract and manifest

`CohortDefinition` retains validated cohort ID, explicit sample order, and full
`ReferenceGenome` identity. `CohortSampleInput` retains source path/index,
tool/version, state, build, and optional TR catalog identity. Duplicate/unsafe
IDs, reordered or missing rows, reference conflicts, and unknown samples are
rejected.

The long-form TSV has one row per sample for every enabled track:

```text
sample track state source_path index_path source_tool source_version reference_build catalog_id
```

Paths may be relative to the TSV. JSON/YAML manifests retain partial states;
no row is silently dropped.

## Small variants / GLnexus

`GLnexusWrapper` uses `CommandRunner` only and no shell pipeline. The official
GLnexus 1.4.1 CLI emits multi-sample BCF on stdout and supports `--dir`,
`--config`, `--threads`, and `--mem-gbytes`. HiFiVar captures binary stdout to
`{cohort}.small.bcf`, then runs discrete `bcftools index`, `bcftools view -Oz`,
and `bcftools index --tbi` commands.

The production target is GLnexus 1.4.1, bcftools 1.21, and `DeepVariantWGS`,
pinned in `workflow/envs/glnexus.yaml`. Input gVCFs follow declared cohort
order. Each must be indexed, contain exactly its declared sample, and use exact
reference contig names. Final VCF identity is validated by exact sample-set
equality: missing, extra, duplicate, and empty sample IDs are rejected. GLnexus
may emit physical VCF columns in a deterministic order different from the
manifest order. HiFiVar therefore retains both `declared_sample_order` and
`output_sample_order`, records `sample_set_match` and `sample_order_match`, and
does not rewrite the raw cohort VCF merely to make those orders agree.

An enabled small-variant cohort must set a positive
`cohort.small_variants.memory_gb`; the same value becomes the Snakemake memory
request and GLnexus `--mem-gbytes` cap. No universal default is inferred. A
three-sample WGS deployment peaked near 153 GB and used 192–200 GB planning,
which is evidence for that workload only.

Streaming QC reports sample/variant/multiallelic counts, FILTER distribution,
per-sample non-reference counts, missing rate, and call rate. It does not load
the complete VCF or perform GWAS, PCA, relatedness, ancestry, or clinical work.
All per-sample metrics are mapped by VCF header sample name rather than by the
manifest's column position.

## SV cohort boundary

Phase 12 prefers Phase 9 harmonized per-sample VCFs. It does not introduce an
unbenchmarked cross-sample clustering algorithm: every source record becomes a
lossless source-native site row and sample-state rows. Other callable samples
are `NOT_OBSERVED`, never `0/0`; unavailable samples retain their explicit run
state. Source sample/ID/VCF, coordinates, REF/ALT, SVTYPE, callable denominator,
support count, and `sample_support_fraction` are retained. AC/AN/AF are not
synthesized. This safely preserves BND, insertion, complex, and native evidence.

### Source-record keys (development policy v2)

`cohort_site_id` identifies one source record, not a cross-sample biological
locus. A unique legacy `sample:source_variant_id` key remains unchanged. If a
legacy key occurs more than once, **all** records in that collision group use
`@hifivar-sv-record-v2:<sample>:<ordinal>`, where ordinal is the 1-based VCF data
record number within that sample's exact source order; headers, comments and
blank lines are not records. The reserved `@` prefix cannot be a legal sample
ID, so generated keys cannot collide with any legacy `sample:ID` key. A disk
index enforces uniqueness before the two formal tables are published.

TSV columns, record/sample order, GT text and sample states are unchanged.
`source_variant_id` retains its historical meaning: native VCF ID when present,
otherwise `CHROM:POS:REF:ALT` as a descriptive fallback. A fallback is not a
caller-assigned ID or a guaranteed locus identity. Original VCFs are never
rewritten. Duplicate records and multi-ALT records remain separate source rows
and intact records respectively, not deduplicated or split alleles.

Track metrics record `site_key_policy`, `disambiguated_record_count`, and
`site_key_uniqueness_verified`. The policy is
`legacy_unique_else_reserved_record_ordinal_v2`. Unambiguous historical tables
remain byte-compatible; only collision-group keys change. An ordinal is stable
for the exact source file/order, not across sorting or replacement of that file.
Do not join different samples' keys as though they were clustered SV sites.

SV/TR single-sample source records must contain exactly ten VCF columns, matching
their one-sample headers. Extra unadvertised sample columns fail validation;
they are not silently ignored. The SQLite index is task-owned staging: success
cleans only unchanged owned files, and a failure retains diagnostics without
publishing a partial table pair. This is not physical two-file atomicity.

## TR cohort boundary

Callable TR inputs declare the same `catalog_id`. VCFs stream through a temporary
on-disk SQLite index; equal TRIDs must have identical contig/start/end/motif
representation. Locus and sample matrices retain genotype text and explicit
states. Absent loci are `NOT_OBSERVED`; no disease threshold or clinical
interpretation is inferred.

TR records require exactly one explicit `INFO/TRID` and `INFO/END`. `TRID`
must be a nonmissing single identifier; an unescaped comma, whitespace, duplicate
key or flag form is rejected. Percent-encoded identifier text is retained without
decoding. The native VCF ID column may be `.`. `POS` and `END` must be explicit
positive ASCII integer coordinates within the signed-64-bit storage range, with
`END >= POS`. Coordinates keep their native VCF/padding convention; missing or
malformed END is never replaced with POS or inferred from REF. This does not
introduce a new STRUC, allele normalization or repeat-genotyping policy.

Duplicate TRIDs within one sample fail. SQLite open, use and standalone close
failures are typed output-validation errors, allowing Phase 12 to retain other
successful tracks. A close failure during an existing failure adds diagnostics
without replacing that failure; process-control exceptions remain interrupts.
No formal table pair is published after these validation/storage failures.
Failed staging is preserved for inspection. Native GT text, source VCF/index
bytes, catalog identity and explicit missing-state semantics remain unchanged.

## Configuration and workflow

`cohort.enabled` defaults false. `cohort.small_variants.enabled`,
`cohort.sv.enabled`, and `cohort.tr.enabled` are independent. Deterministic
outputs live below `results/cohort/<cohort_id>/<small|sv|tr>/`.

`workflow/rules/cohort.smk` contains separate rules plus a downstream manifest.
None is an upstream dependency of per-sample calling. Use Snakemake
`--keep-going` when all enabled independent tracks should continue after one
failure.

## Linux/HPC verification

```bash
conda env create -f workflow/envs/glnexus.yaml
conda activate hifivar-glnexus-1.4.1
glnexus_cli --help
bcftools --version
python -m pytest -p no:cacheprovider tests/unit/test_cohort.py \
  tests/unit/test_glnexus.py tests/unit/test_cohort_tracks.py \
  tests/unit/test_phase12.py tests/integration/test_phase12_complete.py \
  tests/integration/test_snakemake_phase12.py
snakemake --snakefile workflow/Snakefile --configfile effective.phase12.yaml \
  --cores 8 --use-conda --keep-going --dry-run
snakemake --snakefile workflow/Snakefile --configfile effective.phase12.yaml \
  --cores 8 --use-conda --keep-going
bcftools query -l results/cohort/COHORT/small/COHORT.small.vcf.gz
```

`GLNEXUS_LINUX_REAL_TOOL_EXECUTION: PASS`

The historical three-sample full-WGS execution used direct external-tool
commands. It therefore does not prove end-to-end execution through the packaged
HiFiVar DAG. Later bounded DAG, regional cohort and installed-table cases retain
their own source/package identities; they do not turn that historical WGS run
into a packaged-DAG run or certify every current backend combination.

Set `HIFIVAR_GIT_COMMIT` to the validated checkout SHA before a production
Snakemake run; the cohort manifest records it explicitly (or `null` when the
deployment did not provide it). Reference SHA256 is likewise retained when the
upstream `ReferenceGenome` provenance computed/provided it, without forcing an
expensive checksum of every large input.

## Definition of Done

- Cohort/order/reference/state/manifest contracts: complete.
- GLnexus wrapper, dry-run, output validation, provenance: complete.
- Streaming small-variant QC: complete.
- Lossless native SV tables and explicit denominators: complete.
- Catalog-consistent TR locus/sample matrix: complete.
- Independent optional Snakemake tracks: complete.
- Platform-specific software checks and selected Linux workflow/native cases
  are recorded against their exact source/package identities. Production
  Snakemake and bioinformatics execution target Linux/HPC.
- Linux/HPC GLnexus execution and later scoped cohort cases exist; a successful
  component or regional case is not all-route whole-genome acceptance.
- Benchmark development and scientific evaluation are separate from this
  cohort implementation checklist; see [Phase 13](phase13.md).
