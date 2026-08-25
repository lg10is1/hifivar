# Phase 3: DeepVariant Small-Variant Calling

**Phase 3 status: COMPLETE**

## Phase 3.1: DeepVariant Wrapper — COMPLETE

Phase 3 is limited to single-sample DeepVariant SNV/Indel calling. Small
variants remain separate from structural variants and tandem repeats.

The Phase 3.1 boundary introduces tool-neutral request/result models and a
dedicated `DeepVariantWrapper`. Native, Docker, and Apptainer launch strategies
are isolated in `DeepVariantRuntime`; business and workflow code do not compose
container strings. Every real command is a `list[str]` delegated to the shared
`CommandRunner` with `shell=False`.

The wrapper uses the DeepVariant `PACBIO` model, accepts only an indexed aligned
BAM/CRAM artifact, and plans separate deterministic VCF and gVCF outputs. A
dry-run validates inputs and serializes the complete reproducible command
without requiring DeepVariant or a container engine to be installed.

## Phase 3.2: Single-sample calling — COMPLETE

Single-sample execution accepts aligned BAM or CRAM only and requires a readable
adjacent BAI/CSI/CRAI. The declared `AlignmentArtifact.reference` is handed to
DeepVariant, but BAM/CRAM header/reference compatibility is not inferred on
Windows. Output policy is `{sample}.small.vcf.gz`, `{sample}.g.vcf.gz`, plus a
tabix index for each file.

Validation streams headers only and checks BGZF framing, VCF fileformat, exact
sample column, reference contig names, and TBI magic. DeepVariant 1.10 gVCFs do
not need a `NON_REF` ALT declaration or `<NON_REF>` records; gVCF identity is
validated from native `DeepVariant_version`, `RefCall`, and MIN_DP/MED_DP
reference-block header markers. Validation never loads variant records into
memory. Existing outputs and output races are refused unless overwrite is
explicit.

Before a real launch, Linux checks require an open-file soft limit of at least
4096 and recommend 65536. The error reports the observed soft/hard limits and
suggests raising `ulimit -n` or reducing `small.threads` when site policy limits
the hard value. This preflight is skipped for dry-run and non-POSIX development
platforms.

Writable container bind sources are created before Docker/Apptainer command
construction. When callers provide the workflow stderr log path, CommandRunner
also persists stdout to a sibling `*.stdout.log`; this retains DeepVariant
make_examples, call_variants, and postprocess diagnostics.

If external execution succeeds but output validation fails, any existing VCF,
gVCF, and tabix indexes are atomically moved under
`small/quarantine/{sample}.{UTC}/` with `VALIDATION_ERROR.txt`. The declared
outputs disappear so the run still fails, while Snakemake cannot delete the
quarantined evidence during failed-job cleanup.

## Phase 3.3: CLI / Snakemake integration — COMPLETE

The modular `workflow/rules/small.smk` rule is disabled by default and becomes
sample-driven only when `small.enabled` is true. It reads aligned BAM/CRAM paths
from the validated sample sheet, requires their indexes, and declares separate
deterministic VCF/gVCF/TBI outputs, threads, memory, runtime, work paths, and
logs. The rule calls a narrow Python bridge; only `DeepVariantWrapper` owns the
external command. The original Phase 0 smoke target remains active.

## Phase 3.4: Integration / provenance — COMPLETE

`run_phase3()` preserves `AnalysisContext` sample order and resolves indexed
existing BAM/CRAM inputs directly. FASTQ contexts must provide completed Phase 2
`AlignmentArtifact` objects; `collect_phase2_alignment_artifacts()` rejects
dry-run handoffs. Every artifact must match the context sample and exact
reference FASTA/build/contig metadata.

The versioned `Phase3RunReport` records effective resources, execution backend,
container image when applicable, input alignment provenance, reproducible
commands, DeepVariant version/runtime, validated VCF/gVCF artifacts, HiFiVar
version, and UTC timestamp. JSON/YAML writers are atomic and protect existing
reports unless overwrite is explicit.

## Final Phase 3 behavior

| Input at Phase 3 boundary | Behavior |
| --- | --- |
| Indexed BAM | call DeepVariant directly |
| Indexed CRAM | call DeepVariant with the declared reference |
| FASTQ plus completed Phase 2 artifact | call from the derived indexed alignment |
| Raw FASTQ without alignment artifact | fail clearly; DeepVariant never consumes FASTQ |
| BAM/CRAM without index | fail clearly before external execution |

Native execution requires `run_deepvariant`. Docker and Apptainer modes require
their launcher plus an explicit image; backend-specific prefixes and mounts stay
inside `DeepVariantRuntime`. Windows validation uses deterministic fakes and
Snakemake dry-runs. Real DeepVariant execution remains a Linux/HPC requirement.

Phase 3 produces only SNV/Indel VCF and gVCF outputs. It does not perform joint
genotyping, structural/tandem-repeat calling, phasing, assembly, annotation,
functional prediction, benchmark, or review.

No SV/TR caller, phasing, assembly, joint genotyping, annotation, functional
prediction, benchmark, or review feature belongs to Phase 3.
