# Phase 5: TRGT Tandem-Repeat Calling

**Phase 5 status: COMPLETE**

## Phase 5.1: Catalog and artifact boundary — COMPLETE

- `TandemRepeatCatalog` streams a UTF-8 TRGT BED without loading it into memory.
- Every record requires valid BED coordinates plus `ID`, `MOTIFS`, and `STRUC`.
- Duplicate IDs, unknown contigs, and an explicitly mismatched reference build
  are rejected; chromosome names are never rewritten.
- `TandemRepeatArtifact` keeps the TR result independent from small and SV VCFs.

## Phase 5.2: TRGT wrapper — COMPLETE

- `TrgtWrapper` delegates every external command to `CommandRunner` with
  argument lists and `shell=False` behavior.
- The deterministic genotype command contains reference FASTA, indexed aligned
  BAM, repeat catalog, output prefix, sample name, XX/XY karyotype, threads, and
  WGS/targeted preset.
- Real execution checks `trgt`, `bcftools`, and `samtools` availability and
  records parsed versions. Dry-run needs no installed executable and writes no
  output or work directory.
- Existing raw or final products are refused unless overwrite is explicit.

The official TRGT genotype interface accepts an aligned HiFi BAM and produces an
unsorted compressed VCF plus an unsorted spanning-read BAM. Phase 5 therefore
rejects raw FASTQ and CRAM rather than performing a hidden conversion.

## Phase 5.3: Sorting, indexing, and validation — COMPLETE

- `bcftools sort -Oz` followed by TBI indexing produces
  `{sample}.tr.vcf.gz` and `{sample}.tr.vcf.gz.tbi`.
  `bcftools sort` is intentionally single-threaded because bcftools 1.21 does
  not expose a threads option for that subcommand; `bcftools index --threads`
  remains enabled.
- `samtools sort` and `samtools index` produce
  `{sample}.tr.spanning.bam` and `{sample}.tr.spanning.bam.bai` for future TRGT
  visualization and review.
- A zero-exit TRGT process is still rejected when its captured or persisted
  stderr contains `[ERROR] - Locus processing:`; partial-catalog results are
  never promoted to validated artifacts.
- Validation streams only the VCF header and checks BGZF/TBI, one exact sample,
  INFO/TRID/MOTIFS/STRUC declarations, and reference-contig compatibility.
- BAM/BAI checks remain lightweight filesystem checks; Phase 5 does not compute
  coverage, allele-quality summaries, expansion classifications, or methylation
  interpretation.

## Phase 5.4: Integration and provenance — COMPLETE

- `run_phase5()` preserves `AnalysisContext` sample order and accepts completed
  Phase 2 alignment artifacts or existing indexed BAM inputs.
- `tr.karyotype: auto` maps only declared female/male metadata to XX/XY. Missing
  or unknown metadata fails clearly; HiFiVar does not infer biological sex.
- Strict YAML configuration controls catalog, tools, resources, preset,
  karyotype, and overwrite behavior; Phase 5 is disabled by default.
- `Phase5RunReport` atomically serializes context, catalog, per-sample karyotype,
  all commands, versions, runtime, and validated artifacts to JSON/YAML.
- The modular `tandem_repeat` Snakemake rule delegates to the Python wrapper and
  declares deterministic VCF/TBI/spanning-BAM/BAI outputs.

## Official interface evidence

- TRGT CLI: <https://github.com/PacificBiosciences/trgt/blob/main/docs/cli.md>
- Official tiny tutorial and required sorting/indexing:
  <https://github.com/PacificBiosciences/trgt/blob/main/docs/tutorial.md>
- TRGT VCF fields: <https://github.com/PacificBiosciences/trgt/blob/main/docs/vcf_files.md>

The wrapper was checked against the official interface available on 2026-08-21.
The initial Linux target is TRGT 5.1.0, the official release used for this CLI
review; bcftools/samtools remain on the already established 1.21 baseline. This
Phase does not alter the previously verified X-02 matrix or `workflow/envs/`.
Real TRGT 5.1.0 execution remains a Linux/HPC verification boundary; Windows
tests use deterministic runner doubles and Snakemake dry-run.

## Explicitly out of scope

- TRGT cohort merge
- pathogenic expansion thresholds or clinical interpretation
- TRGT plot/deepdive automation and IGV review
- HiPhase or any other phasing
- hifiasm, PAV, SVIM-asm, or other assembly workflows
- combining TR calls into small-variant or structural-variant VCFs
- annotation, benchmark, and Phase 6 behavior
