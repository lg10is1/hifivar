# Benchmarking

[English](benchmarking.md) | [简体中文](benchmarking.zh-CN.md)

This document describes the verified HG002 benchmark results for HiFiVar
release candidate 0.1.0rc4. All metrics are derived from frozen,
audited artifacts and are specific to the exact tool versions, parameters,
reference build, and truth datasets documented below.

## Important scope statements

- Benchmark performance is **dataset-, version-, configuration- and
  region-specific**. Different configurations, samples, or technologies
  may yield different results.
- The CMRG read-SV benchmark covers **challenging medically relevant genes
  only** and is **not** a whole-genome SV truth benchmark.
- No unmatched-truth SCZ result is presented as an accuracy benchmark.
- Caller support count is descriptive evidence, not biological truth.
- Annotation or functional impact is not call confidence.
- No clinical interpretation is provided.

## 1. Whole-genome small-variant benchmark

### Configuration

| Item | Value |
|------|-------|
| Sample | GIAB HG002 |
| Data | PacBio HiFi CCS 11kb |
| Reference | GRCh38 (human_GRCh38_no_alt_analysis_set) |
| Evaluation scope | chr1–chr22 confident regions |
| Caller | DeepVariant 1.10.0 (PACBIO model) |
| Benchmark tool | hap.py 0.3.15 |
| hap.py runtime | Python 2.7.15 legacy validated environment |
| DeepVariant shards | 4 (reduced from 8 due to resource limits; no scientific impact) |

### Canonical metrics

| Variant class | Precision | Recall | F1 |
|---------------|-----------|--------|-----|
| SNP | 0.995444 | 0.977927 | 0.986608 |
| INDEL | 0.916605 | 0.893150 | 0.904726 |

The combined metric (precision=0.9843, recall=0.9665, F1=0.9753) is
**non-canonical** and must not be used as headline performance. It is
derived from SNP and INDEL rows for internal comparison only.

### Stratified small-variant benchmark

Results are stratified by GIAB/GA4GH v3.1 genome difficulty categories.
Only **all difficult regions** and **not in all difficult regions** may be
presented as complementary categories. Do not sum overlapping strata.

#### SNP by stratification

| Stratification | Precision | Recall | F1 |
|----------------|-----------|--------|-----|
| all_difficult_regions | 0.988451 | 0.971033 | 0.979665 |
| not_in_all_difficult_regions | 0.997079 | 0.979560 | 0.988242 |
| low_complexity | 0.983696 | 0.968712 | 0.976146 |
| homopolymer | 0.977721 | 0.963155 | 0.970384 |
| tandem_repeats | 0.991771 | 0.974935 | 0.983281 |
| segmental_duplications | 0.973572 | 0.956138 | 0.964776 |
| mappability | 0.929612 | 0.846068 | 0.885875 |
| GC_content | 0.995844 | 0.979658 | 0.987685 |
| functional_regions | 0.994127 | 0.973027 | 0.983464 |

#### INDEL by stratification

| Stratification | Precision | Recall | F1 |
|----------------|-----------|--------|-----|
| all_difficult_regions | 0.888783 | 0.860543 | 0.874435 |
| not_in_all_difficult_regions | 0.984558 | 0.970714 | 0.977587 |
| low_complexity | 0.881102 | 0.851154 | 0.865869 |
| homopolymer | 0.853250 | 0.823097 | 0.837902 |
| tandem_repeats | 0.958893 | 0.927645 | 0.943010 |
| segmental_duplications | 0.950951 | 0.938938 | 0.944906 |
| mappability | 0.866345 | 0.804769 | 0.834422 |
| GC_content | 0.923517 | 0.901508 | 0.912380 |
| functional_regions | 0.813596 | 0.937028 | 0.870961 |

Do not assert that a genomic property caused lower accuracy unless
supported by a separate analysis. Describe observed differences only.

## 2. Read-based SV benchmark (CMRG GRCh38)

### Configuration

| Item | Value |
|------|-------|
| Sample | HG002 |
| Truth | GIAB CMRG GRCh38 v1.00 StructuralVariant |
| Scope | Challenging medically relevant genes only |
| Benchmark tool | Truvari v5.4.0 |
| Minimum SV size | 50 bp |
| Original CMRG records | 250 |
| Eligible benchmark base count | 216 |
| Excluded | 34 records below 50 bp |

### Why 250 records but 216 base count?

The CMRG truth VCF contains 250 total SV records. Truvari applies
`--sizemin 50` which filters records by SV size. Only records with
SV size >= 50bp are eligible for benchmarking, reducing the benchmark
universe to 216 records. The 34 excluded records are small variants
(< 50bp) that fall outside the standard SV benchmark regime.

### Former NIST SV v0.6 benchmark

The previous benchmark using NIST SV Tier1 v0.6 is
**INVALID_REFERENCE_BUILD_MISMATCH**:

- Truth: GRCh37
- Query calls/reference: GRCh38
- Chromosome renaming did not perform coordinate liftover

All former per-caller Truvari metrics against NIST v0.6 are invalid
for GRCh38 performance interpretation and must not appear in
performance tables or figures.

### Caller metrics

| Caller | Version | TP | FP | FN | Precision | Recall | F1 | GT Concordance |
|--------|---------|-----|-----|-----|-----------|--------|------|----------------|
| Sawfish | 2.2.1 | 196 | 0 | 20 | 1.000 | 0.907 | 0.951 | 0.944 |
| Sniffles2 | 2.8.0 | 160 | 7 | 56 | 0.958 | 0.741 | 0.836 | 0.875 |
| pbsv | 2.11.0 | 142 | 32 | 74 | 0.816 | 0.657 | 0.728 | 0.873 |
| cuteSV | 2.1.4 | 34 | 3 | 182 | 0.919 | 0.157 | 0.269 | 0.882 |

Every public statement must say that these are **CMRG-region results**
and are **not whole-genome SV performance**. Do not claim one caller
is universally superior.

### Truvari parameters

```
--sizemin 50
--sizemax 100000
--refdist 500
--pctseq 0.5
--pctsize 0.5
--includebed HG002_GRCh38_CMRG_SV_v1.00.bed
```

### Caller input BAMs

| Caller | BAM | State |
|--------|-----|-------|
| Sawfish | HG002.nohardclip.sorted.bam | Hard-clips removed (tool-compatibility artifact) |
| Sniffles2 | HG002_ccs11kb.sorted.bam | Original (canonical raw alignment) |
| pbsv | HG002.nohardclip.rg.sorted.bam | Hard-clips removed + RG tags (tool-compatibility artifact) |
| cuteSV | HG002_ccs11kb.sorted.bam | Original (canonical raw alignment) |

## 3. Reproducibility

All benchmark artifacts are checksummed. The frozen evidence files are:

- `HG002_BENCHMARK_FROZEN_SUMMARY.tsv` — machine-readable summary
- `HG002_BENCHMARK_FROZEN_METADATA.yaml` — full metadata with checksums
- `HG002_BENCHMARK_FROZEN_CHECKSUMS.sha256` — all artifact checksums
- `HG002_BENCHMARK_INTERPRETATION.md` — interpretation guidelines

See [Reproducibility](reproducibility.md) for detailed reproducibility
documentation.

## 4. Limitations

See [Known Limitations](known_limitations.md) for the complete list of
limitations, including:

- Single HG002 dataset
- One PacBio CCS 11kb configuration
- CMRG restricted-region SV scope
- Absence of whole-genome GRCh38 SV truth benchmark
- Legacy hap.py environment
- No clinical interpretation

## 5. Citation

When citing benchmark results, cite the individual tools:

- **DeepVariant**: Poplin et al., "A universal SNP and small-indel variant
  caller using deep neural networks", Nature Biotechnology 36, 983–987
  (2018). doi:10.1038/nbt.4235
- **hap.py**: Li et al., "A synthetic-diploid benchmark for accurate
  variant-calling evaluation", Nature Methods 14, 587–590 (2017).
  doi:10.1038/nmeth.4326
- **Truvari**: Alonge et al., "Assessing and benchmarking structural
  variant calls with Truvari", Genome Biology 24, 168 (2023).
  doi:10.1186/s13059-023-02999-w
- **Sawfish**: Huang et al., "Accurate discovery of all types of
  structural variants from long-read sequencing", Genome Research 33,
  1–17 (2023). doi:10.1101/gr.277943.123
- **Sniffles2**: Sedlazeck et al., "Sniffles2: Accurate structural
  variant calling for single-sample and trio long-read sequencing",
  bioRxiv (2023). doi:10.1101/2023.12.27.573420
- **pbsv**: PacBio, "pbsv: PacBio structural variant calling", GitHub
  https://github.com/PacificBiosciences/devbio/wiki/pbsv
- **cuteSV**: Jiang et al., "LongReadSVPipeline: a structural variant
  calling pipeline for long reads", Nature Communications 11, 6251
  (2020). doi:10.1038/s41467-020-19985-w

## Citation verification

CITATION_VERIFICATION_REQUIRED: Verify all DOIs and URLs before
formal publication.
# Phase 15 low-depth experiments

The optional Phase 15 framework constructs deterministic replicated depth
matrices and evaluates user-declared thresholds using Phase 13 metric keys. It
supports an aligned-BAM engineering mode and a raw-FASTQ mode that assigns
PacBio molecules to exact nested subsets before independent realignment. It
does not reinterpret the frozen HG002 results below and does not yet establish
a generally applicable minimum HiFi depth. See
[`docs/phases/phase15.md`](phases/phase15.md) for the experiment contract and
Linux/HPC validation requirements.
