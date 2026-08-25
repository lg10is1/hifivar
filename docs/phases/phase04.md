# Phase 4: Read-based Structural Variant Calling

**Phase 4 status: COMPLETE**

## Phase 4.1: Sawfish Wrapper — COMPLETE

- Dedicated `SawfishWrapper` using `CommandRunner` exclusively
- Official two-step `discover` then single-sample `joint-call` execution
- Indexed BAM/CRAM and reference FASTA/FAI validation
- Deterministic resources, work paths, final `sample.sawfish.sv.vcf.gz`, log
  handoff, dry-run commands, version/runtime provenance, and overwrite refusal
- Native Sawfish `genotyped.sv.vcf.gz` plus TBI are atomically materialized under
  the HiFiVar output convention without changing records
- Windows unit tests use a deterministic runner double; real Sawfish remains a
  Linux/HPC verification boundary

Phase 4.1 does not merge callers, normalize SV records, call tandem repeats, or
create a general SV harmonization framework.

## Phase 4.2: Sniffles2 Wrapper — COMPLETE

- Independent `Sniffles2Wrapper` using one shell-free `CommandRunner` call
- Indexed coordinate-alignment contract for BAM or CRAM plus explicit reference
- Deterministic sample ID, threads, minimum support, and minimum SV length
- Native `.vcf.gz` output uses Sniffles2's documented BGZF/TBI behavior
- Version detection, dry-run, logging handoff, overwrite refusal, output checks,
  runtime provenance, Unicode paths, and fake-runner failure tests

Sniffles2 thresholds remain caller-specific settings and are not interpreted as
shared scientific truth or applied to other callers.

## Phase 4.3: pbsv Wrapper — COMPLETE

- Dedicated `PbsvWrapper` with separate `pbsv discover` and `pbsv call`
  `CommandRunner` invocations; no shell, pipe, or hidden subprocess
- Deterministic `.pbsv.svsig.gz` intermediate and `.pbsv.raw.vcf` native output
- HiFi/CCS calling mode, explicit call threads, indexed BAM/reference inputs,
  version/runtime/command provenance, dry-run, logging handoff, and overwrite
  refusal
- Discover-output validation gates the call step and failures from either step
  are preserved

The official pbsv interface is BAM-to-SVSIG-to-plain-VCF. CRAM is rejected
explicitly rather than silently converted. BGZF/TBI finalization belongs to the
Phase 4.5 common artifact boundary.

## Phase 4.4: cuteSV Wrapper — COMPLETE

- Independent `CuteSvWrapper` using the official positional BAM, reference,
  native VCF, and work-directory contract through `CommandRunner`
- Explicit threads, sample, support/size thresholds, and PacBio HiFi clustering
  parameters kept in the cuteSV request model
- Genotype generation is enabled by default with the official `--genotype`
  switch and can be disabled explicitly through `sv.cutesv.genotype: false`
- Deterministic `.cutesv.raw.vcf`, isolated work directory, version/runtime and
  command provenance, dry-run, log handoff, overwrite refusal, and output checks
- Explicit overwrite removes stale cuteSV work products only when a sibling
  ownership marker proves that HiFiVar created the directory; unmarked or
  symlinked directories are preserved and rejected
- BAM-only workflow contract is explicit; no hidden CRAM conversion is attempted

cuteSV writes a plain native VCF. Compression and indexing are deliberately
owned by the Phase 4.5 common artifact boundary.

## Phase 4.5: Common SV Artifact / Validation Boundary — COMPLETE

- `StructuralVariantArtifact` records caller, sample, reference FASTA/build,
  caller version, complete argument-list commands, BGZF VCF, and TBI
- Streaming header validation checks VCF fileformat, one exact sample, declared
  contigs as a subset of the exact reference names, and INFO/SVTYPE metadata
- BGZF and TBI magic are checked without loading the VCF records into memory
- `BgzipTabixWrapper` finalizes the native plain pbsv/cuteSV VCF with two
  explicit `CommandRunner` calls and records bgzip/tabix versions
- Sawfish and Sniffles2 retain their caller-native BGZF/TBI outputs

This boundary validates and packages raw caller outputs. It does not rewrite
coordinates, alleles, BND pairs, insertion sequences, SVTYPE, FILTER, genotype,
or caller IDs. A raw caller VCF is explicitly `harmonized: false`.

## Phase 4.6: Multi-caller Integration / Regression — COMPLETE

- `run_phase4()` consumes indexed aligned artifacts in deterministic
  `AnalysisContext` sample order and runs enabled callers independently
- Strict `sv.sawfish`, `sv.sniffles2`, `sv.pbsv`, `sv.cutesv`, and
  `sv.finalization` YAML sections reject unknown and invalid settings
- Every sample retains four separately validated outputs by default:
  `sample.sawfish.sv.vcf.gz`, `sample.sniffles2.sv.vcf.gz`,
  `sample.pbsv.sv.vcf.gz`, and `sample.cutesv.sv.vcf.gz`
- `Phase4RunReport` records alignment handoff, caller/finalizer commands,
  caller/bgzip/tabix versions, resources, artifacts, HiFiVar version, and UTC
  timestamp in atomic JSON or YAML
- The modular `read_based_sv` Snakemake rule is sample/caller-driven,
  config-driven, deterministic, disabled by default, and delegates all command
  construction/execution to Python wrappers
- Tiny Windows integration uses deterministic fake tools; Snakemake dry-run
  verifies four independent targets without creating biological outputs

## Final Phase 4 behavior

| Aligned input | Enabled caller behavior |
| --- | --- |
| Indexed BAM | Sawfish, Sniffles2, pbsv, and cuteSV may run independently |
| Indexed CRAM | Sawfish and Sniffles2 may run; pbsv/cuteSV are rejected explicitly |
| Raw FASTQ | Rejected; a completed Phase 2 `AlignmentArtifact` is required |
| Missing alignment index | Rejected before any caller executes |

Lightweight SV validation streams headers and checks BGZF, TBI magic, exact
sample identity, INFO/SVTYPE metadata, and contig-name compatibility. It does
not scan all records, validate every BND mate, calculate caller concordance,
normalize representation, merge calls, select a truth caller, or produce a
consensus/final SV VCF.

Real Sawfish, Sniffles2, pbsv, cuteSV, bgzip, and tabix execution remains a
Linux/HPC verification requirement. Phase 4 does not implement TRGT, phasing,
assembly calling, cohort SV merging, Jasmine, Truvari, annotation, benchmark,
review, or any Phase 5 feature.
