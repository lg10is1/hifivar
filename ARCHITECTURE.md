# HiFiVar Architecture

## Scope

HiFiVar will orchestrate established bioinformatics tools rather than develop new SNV, Indel, or structural-variant callers in v0.1. External tools will be isolated behind dedicated wrappers and executed through shared infrastructure added in later phases.

## Planned analysis flow

```text
Input
  ↓
QC / Alignment
  ↓
├── DeepVariant
├── Read-based SV
├── TRGT
├── hifiasm
│     ├── PAV
│     └── SVIM-asm
├── HiPhase
├── Cohort
├── Annotation
├── Functional prediction
├── IGV review
├── Benchmark
└── Report
```

Assembly-based variant calling will use PAV and SVIM-asm. Dipcall is not a default assembly caller for the current design.

Small variants (SNVs/Indels), structural variants, and tandem-repeat variants will remain in separate outputs. A future internal evidence model may relate these result types without forcing them into a single VCF.

## Current implementation boundary

Phase 0.1–0.6 provides the repository structure, Python packaging, base exception hierarchy, shared logging, layered YAML configuration, and safe external command execution. Future tool integration will use this fixed boundary:

```text
Tool wrapper
    ↓
CommandRunner
    ↓
subprocess (argument list, shell=False)
```

The analysis diagram above describes future modules, not currently implemented behavior.

The Phase 0.7 application boundary is:

```text
User
  ↓
CLI
  ↓
Packaged default + preset + user config + CLI overrides
  ↓
Logging
  ↓
Future workflow and tool layers
```

The CLI does not call `CommandRunner` until a future implemented action needs an
external program. Planned analysis subcommands are intentionally not registered
as empty placeholders.

The Phase 0.8 validation boundary is independent of configuration loading and
external command execution:

```text
Input/config value
    ↓
Validation
    ↓
Future workflow/tool wrapper
```

Validation means checking paths, lightweight text structure, index presence,
contig compatibility, and expected outputs. Validation is not normalization and
is not repair: it does not build indexes, sort alignments, rewrite VCF records,
or rename chromosomes. Query contigs must be a subset of reference contigs, so
`chr1` and `1` remain incompatible unless the user explicitly normalizes inputs.

FASTA, FASTQ, VCF, and BED checks use streaming reads. Checksums use fixed-size
chunks. BAM/CRAM validation is deliberately limited to path, suffix, non-empty,
and index-presence checks until a reliable parser or external-tool layer exists.

The Phase 0.9 workflow boundary keeps configuration ownership in Python and DAG
ownership in Snakemake:

```text
CLI / Config API
     ↓
effective_config.yaml
     ↓
Snakemake DAG
     ↓
Future tool wrappers
     ↓
CommandRunner
     ↓
External tools
```

Snakemake consumes an already merged and validated effective config. It checks
only the small set of keys needed to construct the current DAG and does not
reimplement preset merging or the Config schema. Relative defaults are `work/`
for intermediates and `results/` for retained outputs; both can be overridden by
the effective config. Workflow logs follow `logs/<module>/<rule>.log`.

The current DAG contains only deterministic Phase 0 preparation and smoke
markers with no wildcards. Rule names use `snake_case`, and resource fields use
`threads`, `mem_mb`, and `runtime_min`. Future biological rules will remain in
separate rule modules and use dedicated Conda environments where appropriate.
Container directives and HPC scheduler profiles are reserved for later phases.

## Phase 0 completed boundary

Phase 0.10 verifies the complete foundation through an isolated end-to-end test:

```text
                User
                  │
                  ▼
                 CLI
           ┌──────┼──────────┐
           ▼      ▼          ▼
        Config  Logging  Validation
           │
           ▼
   effective_config.yaml
           │
           ▼
       Snakemake
           │
           ▼
   Future tool wrapper
           │
           ▼
      CommandRunner
           │
           ▼
  Future external tool
```

Packaging, exceptions, logging, configuration, CommandRunner, CLI, validation,
and the infrastructure-only Snakemake smoke DAG are implemented. The workflow
launcher, tool wrappers, and external
bioinformatics tools shown as future components are intentionally absent.

## Phase 1.1 reference boundary

Phase 1.1 introduces the first shared bioinformatics data model without adding
a biological workflow rule:

```text
User/config reference path + explicit build
                    │
                    ▼
            ReferenceGenome
         ┌────────┼────────┐
         ▼        ▼        ▼
      FASTA      FAI contigs   provenance
   validation   names/lengths  build/SHA256
                    │
                    ▼
        Future analysis modules
```

`ReferenceGenome` is immutable and retains all FAI contigs in index order. It
requires a pre-existing FAI and uses the Phase 0 validation API; it never builds
an index, infers a build from contig names, or rewrites chromosome names. Thus
`chr1` and `1` remain incompatible. The primary workflow reference is currently
limited to uncompressed `.fa`, `.fasta`, or `.fna` files.

Ordinary construction validates only the FASTA prefix and streams the small FAI.
It does not scan the complete FASTA. SHA256 calculation is an explicit
provenance operation through `compute_checksum=True` or `with_checksum()`.
Reference summaries are JSON/YAML-friendly and omit full contig metadata unless
explicitly requested.

## Phase 1.2 sample-input boundary

Phase 1.2 separates stable sample identity from the files used as primary input:

```text
Sample
  ├── sample_id
  └── InputDataset
        ├── FASTQ (one or more long-read files)
        ├── BAM   (exactly one file)
        └── CRAM  (exactly one file)
```

Each sample has one primary input mode per analysis. FASTQ, BAM, and CRAM cannot
be mixed in one dataset, and derived alignments, VCFs, assemblies, and reports
do not become fields on `Sample`. Input paths retain caller order and remain
unresolved so Unicode paths, relative paths, and visible symlink spellings are
preserved.

The model reuses lightweight Phase 0 validation. FASTQ validation reads only the
first record of each input; BAM and CRAM validation establishes only readable,
non-empty paths, supported suffixes, and optional index presence. Binary
integrity, read groups, sample-name matching, QC metrics, and CRAM/reference
compatibility are not yet checked.

Phase 1.4 combines these objects explicitly:

```text
AnalysisContext
  ├── ReferenceGenome
  └── SampleRecord
        └── Sample
              └── InputDataset
```

Keeping `ReferenceGenome` outside `Sample` allows the same biological sample to
be analyzed against more than one reference without changing its identity.

## Phase 1.3 sample-sheet boundary

Phase 1.3 adds an ordered cohort-metadata layer while leaving the Phase 1.2
sample identity unchanged:

```text
SampleSheet
  └── SampleRecord
        ├── Sample
        │     └── InputDataset
        ├── declared sex
        ├── father / mother IDs
        ├── phenotype
        └── group

ReferenceGenome (independent analysis context)
```

The only input format is a strict UTF-8 or UTF-8-BOM TSV. Required columns are
`sample_id` and `input`; optional canonical columns are `input_type`, `sex`,
`father`, `mother`, `phenotype`, and `group`. Unknown or duplicate columns are
errors. Blank and comment lines are ignored, while sample order is retained.

Relative input paths are interpreted beside the sheet. Semicolons separate
multiple ordered HiFi FASTQ inputs, and Phase 1.2 remains responsible for input
type inference and file validation. Sample IDs and normalized input paths must
be unique across the sheet.

Pedigree metadata records declared relationships only. Parent IDs must exist in
the same sheet; self-parent, duplicate-parent, declared parent-sex conflicts,
and cycles are rejected. Partial and multi-generation pedigrees are supported,
and complete child/father/mother trios can be enumerated. HiFiVar does not infer
sex or perform Mendelian, trio-calling, or cohort-calling analysis here.

## Phase 1.4 analysis-context and provenance boundary

Phase 1 closes with a run-level integration object and a portable provenance
snapshot. Neither object plans or executes a biological workflow:

```text
ReferenceGenome ─┐
SampleSheet ──────┼──> AnalysisContext ──> RunManifest
EffectiveConfig ─┘                              │
                                               ▼
                                Future workflow launcher
                                      (not implemented)
```

`AnalysisContext` supports either one Python `Sample` or an ordered
`SampleSheet`. It requires a reference and at least one sample, rejects duplicate
sample IDs and reused primary input paths, and detects obvious conflicts between
the effective reference config and `ReferenceGenome`. Mixed FASTQ/BAM/CRAM modes
are allowed across different samples. A generic context does not require
BAM/CRAM indexes.

The current compatibility boundary is explicit: FASTQ is unaligned and makes no
reference-compatibility claim; BAM/CRAM status is `not_checked` because headers
are not parsed. Explicit query-contig collections can still be validated through
the exact-name `ReferenceGenome` API.

`RunManifest` records schema and HiFiVar versions, a UTC timestamp, reference,
samples, effective config, sample-sheet source, and per-input file metadata.
Runtime file paths are made absolute with `Path.absolute()` without resolving
symlinks or mutating input models. Input SHA256 calculation is off by default and
must be requested explicitly; reference SHA256 is recorded only when already
available. Exact secret keys such as `token`, `password`, `secret`, and `api_key`
are redacted recursively. JSON and YAML writers use a sibling temporary file,
atomic replacement, and refuse existing destinations unless overwrite is
explicit.

## Phase 2.1 lightweight input-QC boundary

Phase 2 begins with reusable QC result models and low-cost inspection of the
already validated Phase 1 inputs:

```text
AnalysisContext
      │
      ▼
Lightweight Input QC
      │
      ▼
 RunQCReport
```

Validation and QC have different control-flow semantics:

```text
Validation failure ──> exception ──> stop
QC finding         ──> PASS / WARN / FAIL / NOT_CHECKED
```

For example, an input file that disappeared after model construction is a
validation failure. A readable BAM or CRAM without an index remains a valid
generic Phase 1 input, but Phase 2.1 reports `ALIGNMENT_INDEX_MISSING` and WARN.
QC never wraps a missing input as a quality result.

`QCMetric`, `QCIssue`, and `QCResult` provide frozen, ordered, standard-type
sample results. `RunQCReport` aggregates status deterministically with
FAIL > WARN > PASS > NOT_CHECKED precedence and records the HiFiVar version,
UTC timestamp, reference build/contig/checksum-availability metadata, and status
counts. Reports use the same shared atomic UTF-8 JSON/YAML writers and overwrite
protection as `RunManifest`.

Current FASTQ QC records only input type, file count, aggregate/per-file size,
absolute path spellings, and suffix-level compression state. BAM/CRAM QC records
the same filesystem metadata plus conventional index presence. Existing
lightweight validation may inspect a FASTQ prefix, but QC does not scan all
reads, calculate checksums, parse alignment headers, or call an external tool.
Alignment and biological workflow integration are outside Phase 2.1.

## Phase 2.2 alignment-interface boundary

Phase 2.2 introduces a tool-neutral planning layer between the validated Phase 1
run context and future alignment wrappers:

```text
AnalysisContext (FASTQ samples)
            |
            v
 build_alignment_requests()
            |
            v
  AlignmentRequest per sample
            |
            v
future AlignmentBackend implementation
```

`AlignmentRequest` is a frozen description of one FASTQ sample, reference,
output path/format, selected alignment tool, explicit overwrite policy, and
scheduler-neutral resources. `AlignmentPlan` is the mixed-input decision layer:
FASTQ produces an ALIGN request, while an existing BAM or CRAM produces an
explicit REUSE plan that retains its original path. Request and plan order
follows the validated `AnalysisContext`, and deterministic new output names use
`{sample_id}.aligned.bam` or `{sample_id}.aligned.cram`.

`AlignmentTool` currently names the planned pbmm2 and minimap2 implementations,
while `AlignmentOutputFormat` limits this boundary to BAM and CRAM. The
`AlignmentBackend` protocol defines only shell-free command construction;
`AlignmentCommandPlan` adds a serializable, display-only preview, and
`AlignmentResult` distinguishes PLANNED, COMPLETED, and REUSED outcomes. Later
wrappers remain responsible for executable/version/input/output checks and must
delegate execution to `CommandRunner`.

Planning is side-effect free and refuses an existing generated output unless
overwrite is explicit. It does not create output directories, inspect complete
FASTQ files, build indexes, parse read groups, choose tool-specific presets,
construct a concrete tool command, or execute software. REUSE does not claim
header/reference compatibility. Those operations belong to Phase 2.3 and later
phases.

## Phase 2.3 pbmm2 wrapper boundary

`Pbmm2Wrapper` is the first concrete external bioinformatics wrapper. It accepts
only a pbmm2/BAM `AlignmentRequest`, revalidates the reference/FAI and each FASTQ,
checks and versions the executable for real execution, builds a `list[str]`, and
delegates exclusively to `CommandRunner`. Dry-run deliberately does not require
pbmm2 to be installed and does not create output directories or helper files.

Commands use the HiFi CCS/HIFI preset, on-the-fly coordinate sorting, a stable
sample read group, and explicit alignment threads. pbmm2's automatic BAM index
is disabled with `--bam-index NONE`; Phase 2.4 therefore has one explicit
indexing owner. One FASTQ is passed directly. Multiple ordered FASTQs use a
deterministic, atomically written `fastq.fofn`, supported by upstream pbmm2.

Real execution records the parsed pbmm2 version, command, return state, and
runtime in `AlignmentResult`. Non-zero commands remain `CommandExecutionError`;
missing/empty expected BAMs are `OutputValidationError`. The wrapper does not
parse BAM headers, calculate mapping statistics, run samtools, or call variants.

## Phase 2.4 post-processing and indexing boundary

`AlignmentArtifact` is the validated handoff after alignment or explicit reuse.
Generated pbmm2 BAMs carry a COORDINATE sort-order claim because the wrapper
always requests `--sort`. Existing BAM/CRAM inputs remain UNKNOWN: filenames and
index presence are never used to infer sorting. UNKNOWN artifacts cannot be
indexed automatically.

Indexing is an explicit `AlignmentIndexRequest` executed only by the minimal
`SamtoolsWrapper` through `CommandRunner`. CRAM uses CRAI. BAM normally uses BAI,
but the deterministic auto strategy selects CSI when a reference contig exceeds
the BAI 2^29 coordinate limit. Existing index paths are protected unless
overwrite is explicit. The wrapper validates the alignment and expected index,
records samtools version/runtime, and has a no-install/no-write dry-run path.

Alignment QC remains metadata-level: file size, container, provenance source,
declared sort order, and readable conventional index presence. UNKNOWN sorting
and a missing index are warnings. Phase 2.4 does not parse headers or records,
infer SM/RG/contigs, sort files, calculate coverage/mapping statistics, or invoke
pysam. Those deeper checks remain a future Linux/external-tool boundary.

## Phase 2.5 integrated boundary

The complete Python integration is:

```text
AnalysisContext
      |
      +--> input QC
      |
      +--> AlignmentPlan
              |
              +--> FASTQ: pbmm2 sorted BAM --> samtools index --> alignment QC
              |
              +--> BAM/CRAM: REUSE original path -------------> alignment QC
      |
      +--> Phase2RunReport (JSON/YAML provenance)
```

`run_phase2()` preserves context order and runs no caller. A real FASTQ run must
complete pbmm2 output validation and explicit samtools indexing before alignment
QC. A dry-run plans both shell-free commands without executable checks or output
creation. Existing BAM/CRAM inputs skip pbmm2; their discovered indexes are
retained, but missing indexes are not rebuilt while sorting remains UNKNOWN.

`Phase2RunReport` records effective Phase 2 settings, ordered input/alignment QC,
plans, commands, results, artifacts, external versions, runtimes, and deterministic
overall QC status. Writers are atomic and refuse replacement unless explicit.

The Phase 2 Python boundary remains independently tested and the original Phase
0 smoke target is preserved. Phase 3 adds its biological Snakemake rule without
copying wrapper command construction into the DAG; Linux/HPC executes the same
wrappers through `CommandRunner`.

## Phase 3 DeepVariant boundary

Phase 3 adds only single-sample SNV/Indel calling:

```text
AnalysisContext
      |
      +--> existing indexed BAM/CRAM
      |
      +--> completed Phase 2 AlignmentArtifact
                    |
                    v
            DeepVariantRequest
                    |
                    v
    DeepVariantRuntime (native/docker/apptainer)
                    |
                    v
          DeepVariantWrapper
                    |
                    v
             CommandRunner
                    |
                    v
 sample.small.vcf.gz + sample.g.vcf.gz + TBI
                    |
                    v
             Phase3RunReport
```

Raw FASTQ never enters DeepVariant. BAM/CRAM must have a readable conventional
index and carry the exact `AnalysisContext` reference metadata. The wrapper uses
the PACBIO model, separate VCF/gVCF paths, explicit resources/logs, shell-free
commands, version detection, output validation, and strict overwrite policy.
Runtime/container prefixes are owned only by `DeepVariantRuntime`.

Validation streams BGZF VCF/gVCF headers and checks sample, contigs, gVCF
`NON_REF`, and tabix magic. It does not scan records, normalize variants, or
claim BAM/CRAM header/reference compatibility. The modular Snakemake rule is
disabled by default and delegates to a Python bridge that calls the same wrapper.
Phase 3 introduces no SV/TR caller, joint genotyping, annotation, benchmark,
review, or Phase 4 behavior.

## Phase 4 read-based structural-variant boundary

Phase 4 fans one indexed alignment out to independent callers and never treats
caller count as biological truth:

```text
AnalysisContext + indexed BAM/CRAM
                 |
                 +--> Sawfish discover + joint-call --> caller BGZF/TBI
                 +--> Sniffles2 ----------------------> caller BGZF/TBI
                 +--> pbsv discover + call --> plain VCF --+
                 +--> cuteSV -------------> plain VCF -----+--> bgzip + tabix
                 |
                 v
       four StructuralVariantArtifact objects
                 |
                 v
          Phase4RunReport (no merge)
```

Every external command is an argument list owned by a dedicated wrapper and
executed by `CommandRunner`. Sawfish and pbsv preserve their official multi-step
interfaces as separate invocations. pbsv and cuteSV native plain VCFs pass
through the narrow `BgzipTabixWrapper`; this step compresses/indexes but does
not alter records.

`StructuralVariantArtifact` is a provenance and lightweight-validation model,
not a harmonized variant representation. Validation streams only headers and
checks BGZF/TBI, one exact sample, INFO/SVTYPE declaration, and exact contig-name
subset compatibility with the declared reference. No chromosome renaming,
breakpoint rewriting, BND repair, insertion rewriting, caller comparison,
merging, consensus selection, or truth inference occurs.

The output contract keeps callers separate as
`{sample}.{caller}.sv.vcf.gz`. Indexed BAM supports all four current callers.
Sawfish and Sniffles2 can consume indexed CRAM, while the current pbsv and
cuteSV workflow contracts require BAM and fail clearly for CRAM. FASTQ must be
aligned in Phase 2 before entering Phase 4.

## Phase 5 tandem-repeat boundary

Phase 5 keeps tandem repeats separate from SNV/Indel and structural-variant
outputs:

```text
AnalysisContext + indexed aligned BAM + reference-specific TRGT BED
                              |
                              v
                       TRGT genotype
                              |
             +----------------+----------------+
             |                                 |
       unsorted VCF.gz                 unsorted spanning BAM
             |                                 |
       bcftools sort/index              samtools sort/index
             |                                 |
             v                                 v
   sample.tr.vcf.gz + TBI       sample.tr.spanning.bam + BAI
             |
             v
        Phase5RunReport
```

The catalog validator streams BED rows, requires `ID`, `MOTIFS`, and `STRUC`,
rejects duplicate locus IDs, and enforces exact reference contig names and an
optional exact reference build. TRGT receives only indexed BAM because its
documented genotype interface specifies aligned HiFi BAM; no implicit CRAM or
FASTQ conversion occurs. `karyotype: auto` uses only declared sample sex and
fails when metadata is absent or unknown.

TRGT outputs are explicitly sorted/indexed because the caller documents both
native outputs as unsorted. Lightweight VCF validation reads only the header and
checks BGZF/TBI, exact sample, INFO/TRID/MOTIFS/STRUC declarations, and exact
contig-name compatibility. Phase 5 does not interpret pathogenic expansion
thresholds, merge cohorts, phase calls, generate TRGT plots, or combine TR calls
with small/SV VCFs.
## Phase 6 single-sample phasing boundary

Phase 6 consumes completed Phase 2/3 artifacts and keeps phasing separate from
calling:

```text
AnalysisContext + indexed BAM + indexed sample.small.vcf.gz
                              |
                              v
                     HiPhaseWrapper
                              |
                              v
              sample.phased.vcf.gz + TBI
                              |
                              v
                       Phase6RunReport
```

`PhasingRequest` requires one indexed BAM and a validated small-variant
artifact for the same sample and reference build. It does not consume FASTQ,
CRAM, SV VCF, or TR VCF. The dedicated wrapper owns the deterministic HiPhase
and tabix commands and delegates every process to `CommandRunner`. Lightweight
output validation checks BGZF/TBI, sample, contigs, and the phased-genotype
`PS` FORMAT declaration without scanning the complete VCF.

The modular rule is disabled by default. Phase 6 does not merge variant classes,
perform cohort phasing, infer pedigrees, call variants, or implement assembly.

## Phase 7 haplotype-assembly boundary

Phase 7 is a parallel, reference-independent branch from primary HiFi FASTQ:

```text
Sample FASTQ file(s), in declared order
                  |
                  v
           HifiasmWrapper
                  |
       +----------+-----------+
       |          |           |
   primary GFA  hap1 GFA    hap2 GFA
       |          |           |
       v          v           v
   primary FASTA hap1 FASTA hap2 FASTA
                  |
                  v
           Phase7RunReport
```

BAM/CRAM are rejected rather than implicitly converted. hifiasm's raw
`.bp.p_ctg.gfa`, `.bp.hap1.p_ctg.gfa`, and `.bp.hap2.p_ctg.gfa` outputs
remain immutable evidence. A streaming converter materializes explicit FASTA
artifacts from sequence-bearing GFA segments using atomic replacement.

The config-disabled Snakemake branch calls only the Python wrapper. Phase 7 does
not align assemblies to the reference or implement PAV, SVIM-asm, dipcall,
assembly-based VCFs, cohort assembly, or Phase 8.
## Phase 8 assembly-derived SV boundary

Phase 8 fans validated hifiasm haplotype FASTAs into two independent evidence
streams:

    hap1/hap2 FASTA
       +--> PAV workflow adapter --> PAV raw/final assembly SV artifact
       +--> minimap2/samtools --> SVIM-asm --> SVIM-asm raw/final artifact

PAV and SVIM-asm have separate wrappers, work directories, commands, raw VCFs,
indexes, intermediate files, versions, and provenance. Phase 8 performs no
clustering and neither artifact is labelled as harmonized.

## Phase 9 SV harmonization boundary

    Sawfish / Sniffles2 / pbsv / cuteSV
    PAV / SVIM-asm
                 |
                 v
          SVEvidenceCollection
                 |
                 v
              Jasmine
                 |
         harmonized VCF + evidence TSV
                 |
                 v
     Truvari per-source comparison summaries

The boundary preserves raw source VCFs and caller-native fields. Evidence classes
describe read/assembly origin only. Caller counts, clustering, and Truvari
concordance never create truth or confidence labels. VCF processing is streaming
and partial-run states are explicit.

## Phase 10 IGV/manual-review boundary

Phase 10 is an optional downstream branch. It consumes an explicit selection;
it never changes whether a caller or harmonization task succeeds:

```text
explicit selected variant IDs / loci
        + indexed BAM/CRAM
        + reference FASTA/FAI
        + immutable source VCF/evidence
                       |
                       v
                 ReviewTarget
                       |
                       v
        deterministic IGV batch + screenshots
                       |
                       v
       JSON/YAML/TSV manual-review manifest
```

SNV and insertion loci use the variant anchor, interval variants use their
declared span, and BND targets retain two independent breakend loci. The current
default is the configurable `review.flank_bp`; per-type values may be introduced
later. Screenshot names derive from the unique review ID and locus ordinal, so
variants and BND breakends cannot overwrite one another.

`IgvWrapper` owns IGV batch-mode invocation and delegates execution only to
`CommandRunner`. It does not automate GUI clicks. TRGT-specific plots remain
separate optional evidence metadata and do not block ordinary IGV review.
Review statuses describe one reviewer's visual assessment only: they do not
modify raw VCFs, establish truth, encode pathogenicity, or provide clinical
classification.

## Phase 11 annotation/functional-prioritization boundary

Phase 11 remains an optional downstream fan-out and keeps variant classes and
annotation sources independent:

```text
small VCF ----+--> ANNOVAR annotation --+
SV VCF -------+--> offline VEP ----------+--> Phase 11 provenance
TR VCF -------+--> explicit BED overlap -+
                                           |
                          explicit selected candidates only
                                           |
                                           v
                             functional backend interface
```

`AnnotationInput` retains sample, source variant IDs, caller, VCF, reference,
and build. ANNOVAR and VEP own separate commands, outputs, database provenance,
and errors, and both execute only through `CommandRunner`. ANNOVAR databases and
VEP caches are externally provisioned; core HiFiVar performs no download.

Region overlap is a conservative SV/TR boundary for versioned gene, exon,
regulatory, repeat, and segmental-duplication BED sources. It requires explicit
variant IDs and preserves the caller's original coordinates, breakpoint, and
SVTYPE. It does not normalize or rewrite a VCF.

The AlphaGenome-compatible `FunctionalBackend` accepts only a non-empty explicit
selection and named modalities. Backend results retain model and source-variant
provenance. A high functional-impact score never upgrades call confidence,
creates truth, or supplies clinical classification.

## Phase 12 cohort/multi-sample boundary

Phase 12 adds a downstream `cohort` layer without changing Phase 0–11 outputs.
`hifivar.cohort` owns ordered cohort/sample-state/manifests and streaming VCF
QC; `hifivar.glnexus` owns shell-free GLnexus/bcftools joint genotyping;
`hifivar.cohort_tracks` writes lossless source-native SV and catalog-consistent
TR matrices; and `hifivar.phase12` isolates track failures. Optional modular
Snakemake tracks are not dependencies of single-sample calling. Small variants,
SVs, and TRs remain separate artifact families.

## Phase 13 benchmark/truth-set boundary

Phase 13 is an optional fan-out from immutable variant artifacts:

```text
small VCF + pinned truth/confident BED --> hap.py --------+
read SV VCF + pinned truth/regions ----> Truvari --------+--> benchmark manifest
assembly SV VCF + pinned truth/regions -> Truvari -------+
TR VCF + identical truth catalog ------> exact TR compare+
```

`hifivar.benchmark` owns truth/region/status/metric/result/manifest contracts;
`hifivar.happy` owns the hap.py command and column-based summary parser; Phase 9
`hifivar.truvari` is reused for explicit benchmark policies and streaming
stratification. No benchmark engine modifies query artifacts or supplies
pathogenicity. Independent Snakemake targets are never upstream requirements.

## Phase 14 reporting/release boundary

Phase 14 consumes prior phase records and immutable artifact references:

```text
phase reports + config + provenance + selected artifacts
                         |
                         v
                 FinalRunReport
                  /           \
       JSON/YAML + MD/HTML   explicit bundle plan
                                  |
                    reports/manifests/configs/provenance
                    selected VCF/BCF/TSV + large-data pointers
```

Reporting never reclassifies variants or infers success from file existence.
`FinalStatus` keeps COMPLETE/PARTIAL/FAILED/NOT_RUN/DISABLED distinct. Bundle
selection is explicit; BAM/CRAM/FASTQ and assembly inputs are pointers unless a
caller opts in to copying. Machine reports redact semantic secret keys.

The wheel contains configs and modular workflow data, located through
`installed_workflow_root()` without assuming a checkout path. CI exercises only
Python, fake-tool, packaging, CLI and Snakemake dry-run layers. Real external
tools/WGS remain an independent Linux/HPC release gate. No version bump or
publication is automatic.

## Public distribution boundary

The public distribution follows a core-plus-tools model:

```text
HiFiVar core (wheel / Conda / Docker / Apptainer)
  ├── Python package and CLI
  ├── Snakemake 8/9
  ├── packaged configs/rules/scripts/env YAML
  └── no caller binaries, references, databases, caches or credentials

Enabled analysis branches
  ├── site-native or tool-specific Conda executable
  ├── independently pinned Apptainer image where validated
  └── explicit reference/database/catalog/truth resources
```

The Docker and Apptainer core images are not monolithic WGS environments and
must not launch nested PAV/DeepVariant containers. Linux/HPC should normally run
the core from Conda/venv and mount or invoke tool-specific deployments.

The `0.1.0rc2` packaged DAG consumes existing indexed BAM/CRAM for small/SV/TR
tracks. FASTQ-to-pbmm2 execution remains available through the Phase 2 Python
integration but is not a packaged alignment rule or unified CLI command. Public
documentation must preserve this distinction until that integration exists.
