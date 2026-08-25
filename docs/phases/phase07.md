# Phase 7  hifiasm single-sample haplotype assembly

## Status

**COMPLETE**

Phase 7 adds the independent HiFi FASTQ-to-assembly branch. It does not perform
assembly-based variant calling.

## 7.1 Input contract

- Input is one or more ordered PacBio HiFi FASTQ files from a single sample.
- BAM and CRAM are rejected; Phase 7 never extracts reads from alignments.
- Assembly is reference-independent. The AnalysisContext reference is retained
  only as run provenance and is not passed to hifiasm.
- Existing owned outputs are never silently overwritten.

## 7.2 hifiasm wrapper

`HifiasmWrapper` checks availability and `--version`, builds a deterministic
argument list, and delegates execution only to `CommandRunner`. The command is:

```text
hifiasm -o <work/sample.asm> -t <threads> <ordered FASTQ inputs>
```

Dry-run plans the command without requiring hifiasm or creating outputs.
The compatibility target used for official CLI review is hifiasm 0.25.0-r726;
the exact production build remains subject to Linux/HPC validation.

## 7.3 Raw and derived outputs

The wrapper requires and preserves hifiasm's HiFi-mode GFA outputs:

- `sample.asm.bp.p_ctg.gfa`
- `sample.asm.bp.hap1.p_ctg.gfa`
- `sample.asm.bp.hap2.p_ctg.gfa`

A streaming, atomic converter writes sequence-bearing `S` records to:

- `sample.primary.fa`
- `sample.hap1.fa`
- `sample.hap2.fa`

The converter rejects missing sequence-bearing segments. It does not load a
complete GFA into memory and never deletes or rewrites the raw GFA.

## 7.4 Artifacts and provenance

`AssemblyArtifact` keeps the sample, raw GFA paths, primary/haplotype FASTA
artifacts, source-GFA links, command, and hifiasm version. `Phase7RunReport`
records ordered sample results, settings, reference-independent status, tool
version, runtime, and failures in atomic JSON/YAML.

## 7.5 Workflow integration

`workflow/rules/assembly.smk` is config-driven and disabled by default.
When enabled it creates deterministic work, result, and log paths and calls the
Python bridge, which uses the same wrapper. It is an independent branch and does
not contaminate the alignment, small-variant, SV, TR, or phasing DAGs.

## 7.6 Verification

- IMPLEMENTATION_VERIFICATION: PASS
- MOCK_VERIFICATION: PASS
- OFFICIAL_CLI_VERIFICATION: PASS
- LINUX_REAL_VERIFICATION: NOT_RUN

Windows tests cover single/multiple FASTQ ordering, BAM rejection, version
handling, dry-run, fake execution, GFA preservation/conversion, invalid GFA,
overwrite protection, run provenance, and Snakemake dry-run/execution.

Linux/HPC real verification:

```bash
hifiasm --version
hifiasm -o work/assembly/TINY/TINY.asm -t 4 tiny.hifi.fastq.gz
python -m pytest -p no:cacheprovider tests/unit/test_hifiasm.py \
  tests/integration/test_phase7_complete.py \
  tests/integration/test_snakemake_phase7.py
```

## Intentionally out of scope

Phase 7 does not implement PAV, SVIM-asm, dipcall, assembly-to-reference
alignment, assembly-based VCF generation, cohort assembly, polishing, or Phase 8.
