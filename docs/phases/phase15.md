# Phase 15 — Low-depth HiFi benchmark and optimization

Status: **IMPLEMENTATION COMPLETE; LINUX/HPC SCIENTIFIC VALIDATION REQUIRED**  
Phase 16: **NOT STARTED**

## Scope

Phase 15 adds a reproducible experiment framework for estimating the lowest
assessed PacBio HiFi depth that satisfies an explicitly preregistered accuracy
policy. It does not define a universal minimum depth, change a caller algorithm,
or claim that one parameter profile is optimal.

Two explicit materialization modes are supported. The original
`aligned_bam` engineering mode accepts a coordinate-sorted aligned BAM.
`SamtoolsDownsampleWrapper` performs deterministic template/read-level
subsampling with an explicit seed, writes a new BAM, and creates its BAI through
`CommandRunner`. Before any explicitly authorized overwrite, the wrapper checks
that the source header declares `@HD SO:coordinate`; each result must also pass
`samtools quickcheck`. The source BAM is never edited.

The `raw_fastq` scientific-study mode accepts one or more HiFi FASTQ/FASTQ.GZ
files. `FastqNestedDownsampler` streams each source once per replicate and uses
a stable hash of the PacBio movie/ZMW molecule key. One replicate reuses the
same hash seed at every target depth, making lower fractions exact subsets of
higher fractions. All parameter profiles reuse those same subset files. This
mode is the correct starting point for independent pbmm2 realignment,
DeepVariant calling and hap.py comparison. Raw uBAM conversion remains an
explicit upstream preparation boundary; it is not silently interpreted as
FASTQ.

An additional disabled-by-default `low_depth.end_to_end` branch connects each
raw subset to the existing pbmm2, samtools, DeepVariant and hap.py wrappers.
Every depth/replicate receives an isolated BAM, DeepVariant TMPDIR, VCF/gVCF,
benchmark prefix, measured-coverage record and observation JSON. This changes
orchestration only; it does not change any caller algorithm or infer a depth.
The current end-to-end contract accepts exactly one profile named `baseline`;
candidate optimization profiles require a separately implemented and validated
parameter-application design.

`samtools coverage` output is parsed from captured stdout when no log path is
configured and from CommandRunner's persisted sibling `*.stdout.*` log when a
stderr log is configured. The hap.py executable may be a benchmark-local
executable shim (for example, a launcher for a validated legacy environment);
the same configured executable is used for both `--version` and the benchmark.
If that command returns no parseable version, an independently verified,
explicit numeric `happy_version` is required and recorded with
`version_source: config`; HiFiVar never fabricates a version.

## Experiment contract

`build_low_depth_plan()` expands this preregistered matrix:

```text
parameter profile × target coverage × replicate
                         |
                         v
deterministic seed + output BAM/BAI + immutable plan.json
                         |
          unchanged Phase 3/4 calling and Phase 13 truth benchmark
                         |
                         v
observed coverage + per-class metrics + replicated threshold assessment
```

The plan distinguishes target coverage from observed coverage. Every acceptance
threshold names a variant class and metric (for example `SNP.f1`). A depth is
supported only when the configured number and fraction of completed replicates
meet every threshold. The reported minimum is conservative: that depth and all
higher assessed depths must pass. It is specific to the sample, truth set,
reference, software versions, profiles and policy in the plan.

## Configuration

The branch is disabled by default. The aligned-BAM engineering mode uses:

```yaml
low_depth:
  enabled: true
  study_id: HG002-low-depth-v1
  sample_id: HG002
  input_mode: aligned_bam
  source_bam: /work/HG002.sorted.bam
  baseline_coverage: 30
  target_coverages: [5, 10, 15, 20, 30]
  replicates: 3
  base_seed: 1729
  profiles:
    - name: baseline
      settings: {}
    - name: candidate-a
      settings:
        small.example_parameter: explicit-value
  thresholds:
    - variant_class: SNP
      metric: f1
      minimum: 0.98
    - variant_class: INDEL
      metric: f1
      minimum: 0.90
  minimum_successful_replicates: 3
  minimum_pass_fraction: 1.0
  samtools_executable: samtools
  threads: 4
  memory_mb: 8000
  runtime_minutes: 120
  overwrite: false
```

For a scientifically independent raw-read matrix, replace the input fields
with:

```yaml
low_depth:
  enabled: true
  study_id: HG002-low-depth-raw-v1
  sample_id: HG002
  input_mode: raw_fastq
  source_bam: null
  source_fastqs:
    - /data/HG002/movie1.fastq.gz
    - /data/HG002/movie2.fastq.gz
  baseline_coverage: 30
  target_coverages: [5, 10, 15, 20, 30]
  replicates: 3
  base_seed: 1729
  profiles:
    - name: baseline
      settings: {}
  thresholds:
    - variant_class: SNP
      metric: f1
      minimum: 0.98
  minimum_successful_replicates: 3
  minimum_pass_fraction: 1.0
  overwrite: false
  end_to_end:
    enabled: true
    truth_vcf: /truth/HG002_GRCh38.vcf.gz
    truth_version: GIAB-v4.2.1
    truth_source: GIAB
    confident_bed: /truth/HG002_GRCh38_confident.bed
    confident_bed_version: GIAB-v4.2.1
    metrics_compression: gzip
    happy_version: 0.3.15
    alignment_threads: 8
    deepvariant_threads: 4
    happy_threads: 2
```

The numbers above illustrate syntax only. They are not HiFiVar recommendations.
Profile settings are recorded as provenance; the Phase 15 downsampling branch
does not inject arbitrary profile keys into caller commands.

## Outputs and integration

- `results/low_depth/<study_id>/low_depth_plan.json`
- deterministic BAM/BAI paths under `subsamples/<coverage>/`; the same
  depth/replicate subset is reused across profiles for paired comparison; one
  replicate keeps the same seed across depths so lower fractions are nested
  within higher fractions
- in `raw_fastq` mode, deterministic gzip FASTQs under
  `raw_subsamples/<coverage>/replicate-<n>.fastq.gz`; each replicate is one
  failure-isolated Snakemake job that fans a single streaming pass into all
  requested depths
- Python JSON/YAML study-result serialization after Phase 13 observations are
  supplied
- optional independent BAM/BAI, observed-coverage JSON, VCF/gVCF/TBI, hap.py
  result and observation JSON outputs for every raw depth/replicate

The base Snakemake branch materializes only experiment inputs. The optional
end-to-end branch delegates alignment, indexing, calling and truth comparison
to the unchanged Phase 2, Phase 3 and Phase 13 wrappers. Each
depth/replicate is an independent Snakemake job with its own log, so one failed
replicate can be diagnosed and retried without rerunning the matrix. Profiles
share that exact subset, preventing sampling noise from being confused with a
profile effect. Replicate seeds remain distinct, but each replicate keeps one
seed across target depths. With deterministic samtools fractional sampling,
this makes lower-depth material a nested subset of its higher-depth counterpart
instead of an unrelated draw. Raw FASTQ mode applies the same nesting contract
directly at the molecule level before alignment; aligned-BAM mode remains an
engineering sensitivity study and must not be used alone for a sequencing-cost
claim.

## Scientific safeguards

- A target depth is not measured coverage.
- Passing one replicate is insufficient unless the configured policy says so.
- Missing/failed runs never become zero-valued accuracy observations.
- Support at a lower depth is not reported when a higher assessed depth fails.
- Parameter profiles are explicit hypotheses, not automatically optimized
  algorithms.
- Accuracy does not establish clinical utility, pathogenicity, or biological
  truth.
- Any cost claim needs a separate, documented sequencing-cost model.

## Definition of Done

- deterministic plan, run, threshold, observation and result models: complete;
- CommandRunner-only samtools downsampling/index boundary: complete;
- streaming multi-file raw FASTQ molecule-level nested materialization:
  implementation complete; Linux/HPC validation required;
- config schema/defaults and disabled-by-default Snakemake branch: complete;
- serialization, overwrite, immutability, invalid-input and dry-run tests:
  complete;
- full Windows regression, compileall, build and Snakemake regression: complete;
- raw FASTQ -> pbmm2 -> DeepVariant -> hap.py replicated matrix on Linux/HPC:
  orchestration implemented; tiny Linux/HPC E2E and the full replicated matrix
  remain required
  before publishing any minimum-depth or optimization claim.

Phase 15 changes no Phase 0–14 caller output or scientific interpretation and
does not enter Phase 16.
