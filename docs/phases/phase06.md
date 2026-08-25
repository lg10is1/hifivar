# Phase 6: HiPhase Read-based Phasing

**Phase 6 status: COMPLETE**

## Phase 6.1: Contract and data model -- COMPLETE

- `PhasingRequest` requires one validated, indexed BAM and one validated,
  indexed single-sample DeepVariant VCF artifact for the same sample/reference.
- Raw FASTQ and CRAM are not accepted by the current HiPhase boundary.
- The source BAM and source VCF are immutable handoffs and are never overwritten.
- Deterministic output is `{sample}.phased.vcf.gz` plus TBI.

## Phase 6.2: HiPhase wrapper -- COMPLETE

- `HiPhaseWrapper` uses `CommandRunner` exclusively and builds argument lists.
- The official single-small-VCF interface is `--bam`, `--vcf`,
  `--output-vcf`, `--reference`, `--threads`, and `--sample-name`.
- `--disable-global-realignment` is explicit because the Phase 6 minimum
  contract provides only DeepVariant small variants, matching official guidance.
- HiPhase version detection, executable checks, dry-run, logging, overwrite
  protection, nonzero failure propagation, and explicit tabix indexing are
  implemented.
- Official CLI review target: HiPhase 1.7.0. Linux runtime verification remains
  separate from Windows mock verification.

## Phase 6.3: Artifact validation -- COMPLETE

- Validation streams only the BGZF VCF header.
- It requires VCF format, FORMAT/PS, the exact sample column, reference-compatible
  contig declarations, and a readable TBI index with correct magic.
- No full-record scan, phase-accuracy benchmark, Mendelian benchmark, or
  biological quality inference is performed.

## Phase 6.4: Provenance and interoperability -- COMPLETE

- `PhasedVariantArtifact` records sample, reference build, source BAM, source
  small VCF, HiPhase version, shell-free command, backend, phased VCF, and index.
- `Phase6RunReport` records ordered `AnalysisContext` handoffs and supports
  atomic JSON/YAML serialization.
- Original DeepVariant products remain unchanged and available independently.

## Phase 6.5: Snakemake integration -- COMPLETE

- The modular `hiphase_phasing` rule is disabled by default, sample-driven,
  config-driven, resource-declared, logged, and deterministic.
- It depends on the existing DeepVariant rule and delegates execution to the
  Python wrapper bridge.
- No hifiasm or assembly rule is present at the Phase 6 gate.

## Phase 6.6: Integration and regression -- COMPLETE

The tested handoff is:

```text
AnalysisContext
  -> indexed BAM
  -> DeepVariant small VCF + TBI
  -> HiPhase
  -> sample.phased.vcf.gz + TBI
  -> Phase6RunReport
```

Windows tests use deterministic runner doubles and Snakemake dry-run.
`LINUX_REAL_VERIFICATION: NOT_RUN`.

## Explicitly out of scope

- haplotagged BAM output (supported by HiPhase but not enabled by this minimum path)
- joint small/SV/TR phasing
- cohort or pedigree phasing
- biological phasing accuracy benchmark
- hifiasm, PAV, SVIM-asm, dipcall, annotation, or Phase 8 behavior

## Official interface evidence

- https://github.com/PacificBiosciences/HiPhase/blob/main/docs/user_guide.md
- https://github.com/PacificBiosciences/HiPhase/releases/tag/v1.7.0
