# Known Limitations

[简体中文](zh_CN/known_limitations.md)

This document describes the known limitations of HiFiVar release
candidate 0.1.0rc4 and its benchmark results. Historical statements below are
scoped to that payload, not the latest RC5 source.

## Current development limitations

The [capability/evidence matrix](capabilities.md) describes later staged raw-input,
assembly, annotation/cohort, benchmark and review routing. The latest main launcher
still defers CRAM, distributed Slurm execution and automatic all-tool deployment.
Scoped native VEP execution exists; separate failed deployments remain failed.
TR/phase/BND accuracy, independent biological replication, universal minimum-depth
and sequencing-cost reduction are not established by regional engineering tests.
Screenshots do not complete human review. Report generation remains explicit.

The RC4 SV table uses CMRG only. Later HG002 GRCh38 draft-truth evaluations are
separate fixed-scope studies, not retroactive changes to RC4 values or complete
all-type genome truth. Do not interpret historical NOT_AVAILABLE below as a
current inventory of all GIAB resources.

## 1. Benchmark scope limitations

### Single dataset

All benchmarks use a single HG002 PacBio HiFi CCS 11kb dataset.
Results may not generalize to other samples, ancestries, or
sequencing technologies.

### One configuration

Benchmarks use one specific configuration:
- PacBio HiFi CCS 11kb
- GRCh38 reference
- chr1–chr22 confident regions (small variants)
- CMRG challenging medically relevant genes (SVs)

Different configurations may yield different results.

### CMRG is not whole-genome

The CMRG read-SV benchmark covers **challenging medically relevant
genes only** (250 truth records, 216 eligible after size filtering).
It is **not** a whole-genome SV truth benchmark. Whole-genome GRCh38
read-SV benchmark: NOT_AVAILABLE.

### RC4 whole-genome SV truth boundary

The RC4 evaluation did not establish comprehensive whole-genome, all-SV-type
accuracy. Its CMRG resource covers a specific subset of medically relevant
regions. Later draft resources retain their own build, regions, variant classes
and limitations; availability alone is not a complete accuracy validation.

## 2. Tool and environment limitations

### Legacy hap.py environment

The hap.py benchmark tool (0.3.15) runs on Python 2.7.15, which
is a legacy environment. Future versions of hap.py may produce
different results.

### External tool dependencies

HiFiVar orchestrates external bioinformatics tools. Each tool has
its own:
- version-specific behavior
- container/Conda environment requirements
- reference and database dependencies
- resource requirements

Different tool versions may yield different results.

### Container and deployment

Some optional workflow branches require separate container or
Conda deployment. The core HiFiVar package does not bundle
external caller binaries, references, or databases.

## 3. Technical limitations

### DeepVariant shard count

The DeepVariant WGS run used 4 shards instead of 8 due to
process/thread resource limits (ulimit -u=4096). Shard count
affected parallel execution but not the scientific algorithm.

### Hard-clip removal

Sawfish and pbsv require hard-clip-free BAM inputs. The
hard-clip removal preprocessing is a tool-compatibility artifact
that may slightly alter the alignment compared to the original.

### No unified FASTQ-to-VCF pipeline

The current release candidate does not expose a single unified
`hifivar run` command or Snakemake alignment rule. FASTQ
alignment must be performed separately before entering the
variant calling pipeline.

## 4. Interpretation limitations

### Support metrics, not biological truth

TP/FP/FN counts and derived metrics (precision, recall, F1)
represent computational concordance between caller output and
truth dataset. They do not constitute biological truth or
clinical interpretation.

### No clinical claims

Benchmark metrics do not imply:
- clinical validity
- diagnostic utility
- biological significance
- pathogenicity
- medical actionability

### Configuration-specific

All results are specific to the exact tool versions, parameters,
reference build, and truth datasets used. They should not be
generalized to other configurations without re-benchmarking.

### Caller support is not truth

The number of callers supporting a variant is descriptive
evidence only. It does not establish biological truth,
confidence, or clinical significance.

### Annotation is not confidence

Functional annotation scores (e.g., from ANNOVAR, VEP, or
AlphaGenome) do not represent variant calling confidence or
truth status.

## 5. Future work recorded at the RC4 checkpoint

- Whole-genome GRCh38 SV truth benchmark (pending GIAB resources)
- Unified FASTQ-to-VCF pipeline
- Additional sample benchmarks
- Additional technology benchmarks
- Cohort-level benchmarks
- Additional stratification analyses
