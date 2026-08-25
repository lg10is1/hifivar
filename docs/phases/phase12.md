# Phase 12 — Cohort / Multi-sample Analysis

**Status: COMPLETE (implementation and Windows mock validation)**

Phase 13 benchmark is **NOT STARTED**. GLnexus and bcftools real-tool execution
still requires Linux/HPC delta validation.

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
reference contig names. Final VCF sample order/set must match the cohort.

Streaming QC reports sample/variant/multiallelic counts, FILTER distribution,
per-sample non-reference counts, missing rate, and call rate. It does not load
the complete VCF or perform GWAS, PCA, relatedness, ancestry, or clinical work.

## SV cohort boundary

Phase 12 prefers Phase 9 harmonized per-sample VCFs. It does not introduce an
unbenchmarked cross-sample clustering algorithm: every source record becomes a
lossless source-native site row and sample-state rows. Other callable samples
are `NOT_OBSERVED`, never `0/0`; unavailable samples retain their explicit run
state. Source sample/ID/VCF, coordinates, REF/ALT, SVTYPE, callable denominator,
support count, and `sample_support_fraction` are retained. AC/AN/AF are not
synthesized. This safely preserves BND, insertion, complex, and native evidence.

## TR cohort boundary

Callable TR inputs declare the same `catalog_id`. VCFs stream through a temporary
on-disk SQLite index; equal TRIDs must have identical contig/start/end/motif
representation. Locus and sample matrices retain genotype text and explicit
states. Absent loci are `NOT_OBSERVED`; no disease threshold or clinical
interpretation is inferred.

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

`GLNEXUS_LINUX_REAL_VERIFICATION: NOT_RUN`

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
- Windows mock/unit/integration/Snakemake regression: complete.
- Linux/HPC real GLnexus execution: pending, not claimed.
- Phase 13: not started.
