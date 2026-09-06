# HiFiVar Project Status

Status date: 2026-09-06  
Current release candidate: `0.1.0rc4`  
Production platform: Linux/HPC  
License: Apache-2.0

## Current boundary

HiFiVar implements the Phase 0–14 modular contracts for configuration,
provenance, QC, alignment, small variants, read- and assembly-based structural
variants, tandem repeats, phasing, harmonization, manual review, annotation,
cohort analysis, benchmarking, and reporting. External bioinformatics tools
remain independently provisioned and versioned.

The packaged 0.1.0rc4 Snakemake small/SV/TR branches consume an existing
indexed BAM or CRAM. The Phase 2 Python API supports pbmm2 alignment planning
and execution, but there is not yet a unified `hifivar run` command or packaged
FASTQ/uBAM-to-results alignment rule. This release must therefore not be
described as a one-command raw-read workflow.

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

## Frozen public benchmark scopes

The HG002 small-variant benchmark uses PacBio HiFi CCS 11 kb data, GRCh38,
DeepVariant 1.10.0, hap.py 0.3.15, and chr1–chr22 confident regions. The read-SV
benchmark uses GIAB CMRG GRCh38 v1.00 and covers challenging medically relevant
genes only. CMRG values are not whole-genome SV performance. See
`benchmarking.md` and `known_limitations.md` before reusing any metric.

## Known product gaps

1. No unified `hifivar run` entry point or packaged alignment rule.
2. No formal uBAM-to-workflow handoff.
3. No production minimap2 wrapper equivalent to the validated pbmm2 path.
4. Cross-sample SV clustering is intentionally unsupported.
5. VEP and AlphaGenome lack real production-backend validation.
6. Low-depth optimization and PCR/primer-design validation are future research
   extensions, not 0.1.0rc4 capabilities.
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
3. Add a preregistered, replicated low-depth HiFi benchmark; do not hard-code a
   universal minimum coverage before the data support it.
4. Add an optional IGV-to-PCR evidence branch with at least three in-silico
   primer candidates, explicit specificity checks, predicted allele products,
   and independent wet-lab result records.
5. Extend validation to independent GIAB samples before making generalized
   performance claims.
