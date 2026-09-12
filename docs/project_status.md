# HiFiVar Project Status

Status date: 2026-09-09  
Current release candidate: `0.1.0rc4`  
Production platform: Linux/HPC  
License: Apache-2.0

## Current boundary

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
and are not yet connected to a unified end-to-end CLI workflow. Neither the
release nor the current development line should therefore be described as a
one-command raw-read workflow.

## Validation summary

- Windows regression: 917 passed, 2 skipped, 0 failed on 2026-09-06; the two
  skips require Windows symbolic-link privileges.
- Independent Linux RC4 package, clean-install, workflow-resource, security,
  PAV SV-only, and release validation: PASS.
- Real Linux executions include pbmm2, DeepVariant, the four read-SV callers,
  TRGT, HiPhase, hifiasm, PAV, SVIM-asm, Jasmine, GLnexus, IGV, ANNOVAR,
  hap.py, and Truvari within their documented scopes.
- VEP remains environment-blocked for real execution; AlphaGenome is an
  explicit selected-candidate interface rather than a validated production
  backend.
- Phase 15 low-depth implementation has passed local and Linux engineering
  validation, including nested raw-FASTQ sampling and a tiny end-to-end run.
  The preregistered full low-depth benchmark matrix remains a separate
  completion gate. No minimum supported depth or sequencing-cost reduction
  claim has yet been established.

## Frozen public benchmark scopes

The HG002 small-variant benchmark uses PacBio HiFi CCS 11 kb data, GRCh38,
DeepVariant 1.10.0, hap.py 0.3.15, and chr1–chr22 confident regions. The read-SV
benchmark uses GIAB CMRG GRCh38 v1.00 and covers challenging medically relevant
genes only. CMRG values are not whole-genome SV performance. See
`benchmarking.md` and `known_limitations.md` before reusing any metric.

## Known product gaps

1. No unified `hifivar run` entry point or raw-read-to-results DAG.
2. No formal uBAM-to-workflow handoff.
3. The post-0.1.0rc4 development line provides a minimap2 wrapper and a
   standalone alignment rule, but unified CLI/DAG integration remains pending.
4. Cross-sample SV clustering is intentionally unsupported.
5. VEP and AlphaGenome lack real production-backend validation.
6. Phase 15 provides aligned-BAM engineering subsampling and raw-FASTQ
   molecule-level nested materialization. Its preregistered 4-depth ×
   3-replicate matrix (5×/10×/15×/20×) is tracked as a separate Linux/HPC
   validation activity. No minimum depth or optimized profile is scientifically
   validated or part of 0.1.0rc4. PCR/primer-design validation remains a future
   research extension.
7. Core Docker/Apptainer and Conda distribution claims remain separate from
   tool-specific real validation and must not be overstated.

## Documentation policy

Phase documents retain point-in-time design history. This page and the current
README describe the latest supported boundary. Historical audit failures remain
preserved in the private evidence repository and may be superseded only by an
explicit later validation record.

## Recommended next milestones

1. Publish the audited benchmark, limitations, reproducibility, manuscript, and
   figure materials without private evidence or local environments.
2. Start a separately versioned development line for unified execution and
   uBAM/FASTQ alignment integration.
3. Complete and independently review the preregistered Phase 15 matrix on
   Linux/HPC and independent GIAB data; do not publish a universal minimum
   coverage before the data support it.
4. Add an optional IGV-to-PCR evidence branch with at least three in-silico
   primer candidates, explicit specificity checks, predicted allele products,
   and independent wet-lab result records.
5. Extend validation to independent GIAB samples before making generalized
   performance claims.
