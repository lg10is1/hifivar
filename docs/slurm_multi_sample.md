# Multi-sample execution with Slurm

This guide describes a transparent Slurm deployment for HiFiVar `0.1.0rc2`.
It uses three PacBio HiFi uBAMs as a concrete example, but the same manifests
and array pattern scale to larger cohorts.

> HiFiVar `0.1.0rc2` does not include a packaged uBAM/FASTQ alignment rule.
> Coordinate-sorted, indexed BAM/CRAM is the input contract of the packaged
> small-variant, read-SV, TR, and cohort DAG. For PacBio uBAM, run a controlled
> pbmm2 preprocessing array first and record that boundary in provenance.

## Execution overview

```text
three immutable PacBio uBAMs
          |
          v
Slurm array: pbmm2 align + samtools index       one task per sample
          |
          v
three coordinate-sorted, indexed BAMs
          |
          v
one HiFiVar samples.tsv + one effective config
          |
          v
Snakemake Slurm executor
  |       |         |
  |       |         +-- TRGT jobs, if one compatible catalog is configured
  |       +------------ independent read-SV caller jobs
  +-------------------- DeepVariant VCF/gVCF jobs
          |
          v
explicit cohort input manifest
  |       |       |
  |       |       +-- TR cohort matrix
  |       +---------- SV cohort tables
  +------------------ GLnexus small-variant cohort
```

Alignment must finish and pass validation before calling is submitted. Cohort
jobs must consume completed per-sample outputs; missing results are explicit
states, never inferred as homozygous reference.

## Two supported scheduler patterns

### Pattern A: Slurm executor plugin (recommended for multiple nodes)

Snakemake 8 and newer use executor plugins for distributed execution. The
official Slurm plugin submits each ready Snakemake job with `sbatch`; `--jobs`
limits simultaneous submissions. See the
[official plugin documentation](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)
and [Snakemake executor documentation](https://snakemake.readthedocs.io/en/stable/executing/executors.html).

Run the Snakemake controller on the login/submission host if site policy permits.
Do not wrap a Slurm-executor controller inside another batch job unless cluster
administrators explicitly require and support that pattern.

### Pattern B: one `sbatch` allocation (fallback)

Submit one batch job that runs Snakemake with the local executor and
`--cores "$SLURM_CPUS_PER_TASK"`. This needs no executor plugin but all rules run
inside one node/allocation. It is useful for tiny validation, not efficient for
three WGS samples or many independent callers.

## Example inputs and immutable reference

Reference:

```text
/data/project/references/GRCh38.fa
```

Example sample identities and uBAMs:

```text
SAMPLE_A  /data/project/pacbio_smrtlink/data/example_data/batch_1/hifi_reads/sample_a.hifi_reads.bam
SAMPLE_B  /data/project/pacbio_smrtlink/data/example_data/batch_1/hifi_reads/sample_b.hifi_reads.bam
SAMPLE_C  /data/project/pacbio_smrtlink/data/example_data/batch_2/hifi_reads/sample_c.hifi_reads.bam
```

Confirm each header `@RG SM` before accepting these IDs. Do not infer sex,
phenotype, or relatedness from filenames.

## 1. Create an isolated run directory

Use a new directory for every run. Never place generated files beside the
original uBAMs.

```bash
export RUN_ROOT="/work/hifivar/realdata_validation/example_rc2_$(date -u +%Y%m%dT%H%M%SZ)"
export REFERENCE_FASTA="/data/project/references/GRCh38.fa"
export REFERENCE_BUILD="GRCh38"

mkdir -p "$RUN_ROOT"/{configs,slurm,logs/alignment,logs/workflow,alignments,work,results,validation}
printf 'RUN_ROOT=%s\n' "$RUN_ROOT"
```

Before scheduling WGS, inspect space and Slurm policy:

```bash
df -h "$RUN_ROOT"
sinfo -o '%P %a %l %c %m %G'
sacctmgr show assoc user="$USER" format=User,Account,Partition,QOS 2>/dev/null || true
```

Do not choose an account, partition, or QoS from this guide. Use values assigned
by cluster administrators.

## 2. Preflight the three uBAMs

```bash
for ubam in \
  /data/project/pacbio_smrtlink/data/example_data/batch_1/hifi_reads/sample_a.hifi_reads.bam \
  /data/project/pacbio_smrtlink/data/example_data/batch_1/hifi_reads/sample_b.hifi_reads.bam \
  /data/project/pacbio_smrtlink/data/example_data/batch_2/hifi_reads/sample_c.hifi_reads.bam
do
  test -r "$ubam"
  test -s "$ubam"
  samtools quickcheck -v "$ubam"
  samtools view -H "$ubam" | grep -E '^@HD|^@SQ|^@RG'
done
```

An unaligned uBAM normally lacks reference `@SQ` records and coordinate sort
state. If a file is already coordinate-sorted, audit it separately instead of
realigning it silently.

Validate the reference without rewriting it:

```bash
test -r "$REFERENCE_FASTA"
test -s "$REFERENCE_FASTA"
test -s "${REFERENCE_FASTA}.fai"
head "${REFERENCE_FASTA}.fai"
```

## 3. Build the alignment array manifest

Create `$RUN_ROOT/configs/ubam_samples.tsv` as a real tab-separated file:

```text
sample_id	ubam
SAMPLE_A	/data/project/pacbio_smrtlink/data/example_data/batch_1/hifi_reads/sample_a.hifi_reads.bam
SAMPLE_B	/data/project/pacbio_smrtlink/data/example_data/batch_1/hifi_reads/sample_b.hifi_reads.bam
SAMPLE_C	/data/project/pacbio_smrtlink/data/example_data/batch_2/hifi_reads/sample_c.hifi_reads.bam
```

Check it before submission:

```bash
awk -F '\t' 'NR == 1 {if ($1 != "sample_id" || $2 != "ubam") exit 1}
  NR > 1 {if (NF != 2 || $1 == "" || $2 == "") exit 2; print NR-1, $1, $2}' \
  "$RUN_ROOT/configs/ubam_samples.tsv"
```

Array task `1` maps to the first data row, task `2` to the second, and so on.

## 4. Create the pbmm2 alignment array script

Save the following as `$RUN_ROOT/slurm/align_ubam_array.sbatch`. Replace account,
partition, Conda initialization, and environment with site-validated values.

```bash
#!/usr/bin/env bash
#SBATCH --job-name=hifivar-align
#SBATCH --partition=<SLURM_PARTITION>
#SBATCH --account=<SLURM_ACCOUNT>
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/alignment/%A_%a.stdout.log
#SBATCH --error=logs/alignment/%A_%a.stderr.log

set -euo pipefail

: "${RUN_ROOT:?RUN_ROOT must be exported by sbatch}"
: "${REFERENCE_FASTA:?REFERENCE_FASTA must be exported by sbatch}"

source <CONDA_ROOT>/etc/profile.d/conda.sh
conda activate <VALIDATED_PBMM2_ENV>

command -v pbmm2 >/dev/null
command -v samtools >/dev/null
pbmm2 --version
samtools --version | head -1

manifest="$RUN_ROOT/configs/ubam_samples.tsv"
row="$(awk -F '\t' -v task="$SLURM_ARRAY_TASK_ID" 'NR == task + 1 {print; exit}' "$manifest")"
if [[ -z "$row" ]]; then
  echo "No manifest row for array task $SLURM_ARRAY_TASK_ID" >&2
  exit 10
fi

IFS=$'\t' read -r sample_id ubam <<<"$row"
if [[ -z "$sample_id" || -z "$ubam" ]]; then
  echo "Malformed manifest row: $row" >&2
  exit 11
fi
if [[ ! -r "$ubam" || ! -s "$ubam" ]]; then
  echo "Unreadable or empty uBAM: $ubam" >&2
  exit 12
fi

output="$RUN_ROOT/alignments/${sample_id}.aligned.bam"
index="${output}.bai"
if [[ -e "$output" || -e "$index" ]]; then
  echo "Refusing existing output: $output or $index" >&2
  exit 13
fi

samtools quickcheck -v "$ubam"
read_group=$'@RG\tID:'"$sample_id"$'\tSM:'"$sample_id"$'\tPL:PACBIO'

pbmm2 align \
  "$REFERENCE_FASTA" \
  "$ubam" \
  "$output" \
  --preset CCS \
  --sort \
  --bam-index NONE \
  --rg "$read_group" \
  -j "$SLURM_CPUS_PER_TASK" \
  --log-level INFO

samtools quickcheck -v "$output"
samtools index -@ 4 "$output"
test -s "$index"

header="$(samtools view -H "$output")"
grep -q '^@SQ' <<<"$header"
grep -q 'SO:coordinate' <<<"$header"
grep -q $'SM:'"$sample_id" <<<"$header"

mapped="$(samtools view -c -F 4 "$output")"
if [[ "$mapped" -le 0 ]]; then
  echo "No mapped reads for $sample_id" >&2
  exit 14
fi

printf '%s\t%s\t%s\t%s\n' \
  "$sample_id" "$output" "$index" "$mapped" \
  >"$RUN_ROOT/validation/${sample_id}.alignment.complete.tsv"
```

Compare all arguments with the installed `pbmm2 align --help` before submission.
The real CLI remains authoritative.

## 5. Submit and monitor the alignment array

Run from `$RUN_ROOT`, because the log paths are relative to the submission
directory:

```bash
cd "$RUN_ROOT"

ALIGN_JOB_ID="$(sbatch \
  --parsable \
  --array=1-3%2 \
  --export=ALL,RUN_ROOT="$RUN_ROOT",REFERENCE_FASTA="$REFERENCE_FASTA" \
  slurm/align_ubam_array.sbatch)"

printf 'ALIGN_JOB_ID=%s\n' "$ALIGN_JOB_ID" | tee validation/alignment_job_id.txt
```

`--array=1-3%2` creates three tasks but permits at most two to run at once.
Change `%2` only after considering filesystem bandwidth, account limits, memory,
and available nodes.

Monitor:

```bash
squeue -j "$ALIGN_JOB_ID" -o '%.18i %.9P %.28j %.8T %.10M %.6D %R'
sacct -j "$ALIGN_JOB_ID" \
  --format=JobID,JobName%28,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList
```

After every task is `COMPLETED`, require all completion records:

```bash
test -s "$RUN_ROOT/validation/SAMPLE_A.alignment.complete.tsv"
test -s "$RUN_ROOT/validation/SAMPLE_B.alignment.complete.tsv"
test -s "$RUN_ROOT/validation/SAMPLE_C.alignment.complete.tsv"
```

Do not submit calling if any array task failed, was cancelled, timed out, or ran
out of memory.

## 6. Create the aligned HiFiVar sample sheet

Create `$RUN_ROOT/configs/aligned_samples.tsv`:

```text
sample_id	input	input_type
SAMPLE_A	<RUN_ROOT>/alignments/SAMPLE_A.aligned.bam	bam
SAMPLE_B	<RUN_ROOT>/alignments/SAMPLE_B.aligned.bam	bam
SAMPLE_C	<RUN_ROOT>/alignments/SAMPLE_C.aligned.bam	bam
```

Replace `<RUN_ROOT>` with the actual absolute path. Revalidate each BAM and BAI
before generating the DAG.

## 7. Create the multi-sample calling config

Start with one enabled track. This example enables only DeepVariant; add read-SV
and TR only after executables/images/catalogs pass real preflight.

```yaml
project:
  name: example_3sample_rc2

reference:
  fasta: /data/project/references/GRCh38.fa
  build: GRCh38

samples:
  sheet: <RUN_ROOT>/configs/aligned_samples.tsv

runtime:
  threads: 16
  tmpdir: <RUN_ROOT>/work/tmp

paths:
  workdir: <RUN_ROOT>/work/calling
  outdir: <RUN_ROOT>/results/calling

logging:
  level: INFO
  file: <RUN_ROOT>/logs/workflow/hifivar.log

small:
  enabled: true
  execution_mode: apptainer
  deepvariant_image: <VALIDATED_DEEPVARIANT_1.10_IMAGE>
  model_type: PACBIO
  threads: 16
  memory_mb: 64000
  runtime_minutes: 2880
  overwrite: false

sv:
  enabled: false
tr:
  enabled: false
phasing:
  enabled: false
assembly:
  enabled: false
assembly_sv:
  enabled: false
review:
  enabled: false
annotation:
  enabled: false
cohort:
  enabled: false
benchmark:
  enabled: false
```

Replace every placeholder and validate:

```bash
conda activate hifivar-rc2-public-test
unset PYTHONPATH PYTHONHOME

hifivar --config "$RUN_ROOT/configs/calling.yaml" config validate
hifivar --config "$RUN_ROOT/configs/calling.yaml" config dump-effective \
  --output "$RUN_ROOT/configs/calling.effective.yaml"

WORKFLOW_ROOT="$(python -c 'from hifivar.package_resources import installed_workflow_root; print(installed_workflow_root())')"

snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile "$RUN_ROOT/configs/calling.effective.yaml" \
  --cores 1 \
  --dry-run \
  --printshellcmds
```

The three-row sheet expands `deepvariant_small` into three sample jobs. Enabling
four read-SV callers expands to twelve read-SV jobs. No per-sample shell loop is
required inside the HiFiVar calling DAG.

## 8. Install the Slurm executor without modifying the validated environment

The plugin requires Python 3.11 or newer. Clone the validated Python 3.12
environment instead of modifying it:

```bash
mamba create -n hifivar-rc2-slurm --clone hifivar-rc2-public-test
conda activate hifivar-rc2-slurm
python -m pip install 'snakemake-executor-plugin-slurm==2.7.1'

snakemake --version
python -m pip show snakemake-executor-plugin-slurm
snakemake --help | grep -A2 -B2 slurm
```

The `2.7.1` pin is a deployment candidate, not a caller dependency or a claim of
site validation. If HPC cannot reach PyPI, transfer the plugin wheel through the
same checksum-verified Windows/SCP process used for the HiFiVar release. Record
the installed plugin version in provenance.

## 9. Create a site-specific Slurm profile

Save this as `$RUN_ROOT/configs/slurm/profile.v9+.yaml` and replace account and
partition values:

```yaml
executor: slurm
jobs: 20
latency-wait: 120
rerun-incomplete: true
keep-going: false
printshellcmds: true
show-failed-logs: true

default-resources:
  slurm_account: "<SLURM_ACCOUNT>"
  slurm_partition: "<SLURM_PARTITION>"
  mem_mb: 4000
  runtime: 60

set-resources:
  deepvariant_small:
    runtime: 2880
  read_based_sv:
    runtime: 1440
  tandem_repeat:
    runtime: 720
  cohort_small_variants:
    runtime: 1440
  cohort_sv:
    runtime: 240
  cohort_tr:
    runtime: 240
```

HiFiVar rules already expose `threads` and `mem_mb`, which the plugin maps to
CPU and total `--mem`. HiFiVar also records `runtime_min`, but the official Slurm
executor maps the standard resource `runtime` to `sbatch --time`. This profile
therefore supplies explicit `runtime` values. Do not assume `runtime_min` is
translated automatically.

`deepvariant_slots` is a workflow-wide concurrency cap, not a per-job CPU
request. HiFiVar registers it from the safe default
`small.max_concurrent_samples: 1`. Raise that config value only after a tiny DAG
run demonstrates independent per-sample `TMPDIR` paths and the site can sustain
concurrent TensorFlow/JAX model loading.

Set rule-specific `slurm_partition` values in the profile if needed. Keep
cluster account, partition, and constraints out of the portable HiFiVar config.

## 10. Dry-run and execute with the Slurm plugin

```bash
conda activate hifivar-rc2-slurm
unset PYTHONPATH PYTHONHOME

snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile "$RUN_ROOT/configs/calling.effective.yaml" \
  --profile "$RUN_ROOT/configs/slurm" \
  --dry-run \
  --printshellcmds
```

Then submit from the login/submission host:

```bash
snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile "$RUN_ROOT/configs/calling.effective.yaml" \
  --profile "$RUN_ROOT/configs/slurm" \
  --printshellcmds \
  2>&1 | tee "$RUN_ROOT/logs/workflow/snakemake.controller.log"
```

Use a persistent terminal multiplexer if permitted. The controller must remain
alive to monitor jobs. Do not start a second controller on the same output tree.

Monitor:

```bash
squeue -u "$USER" -o '%.18i %.9P %.32j %.8T %.10M %.9l %.6D %R'
sacct -S "$(date +%F)" -u "$USER" \
  --format=JobID,JobName%36,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS
```

## 11. Single-node fallback

If the plugin cannot be provisioned, submit one controller allocation:

```bash
#!/usr/bin/env bash
#SBATCH --job-name=hifivar-calling
#SBATCH --partition=<SLURM_PARTITION>
#SBATCH --account=<SLURM_ACCOUNT>
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/workflow/%j.stdout.log
#SBATCH --error=logs/workflow/%j.stderr.log

set -euo pipefail
source <CONDA_ROOT>/etc/profile.d/conda.sh
conda activate hifivar-rc2-public-test
unset PYTHONPATH PYTHONHOME

ulimit -n 65536 || true
if [[ "$(ulimit -n)" -lt 4096 ]]; then
  echo 'DeepVariant requires at least 4096 file descriptors.' >&2
  exit 20
fi

WORKFLOW_ROOT="$(python -c 'from hifivar.package_resources import installed_workflow_root; print(installed_workflow_root())')"

snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile "$RUN_ROOT/configs/calling.effective.yaml" \
  --cores "$SLURM_CPUS_PER_TASK" \
  --resources mem_mb=120000 deepvariant_slots=1 \
  --printshellcmds \
  --rerun-incomplete
```

This runs on one node. Do not confuse `--cores 32` with 32 cluster jobs.

## 12. Resource starting points

These are RC2 configuration starting points, not universal guarantees. Revise
them only after reviewing `sacct`, tool logs, coverage, and site limits.

| Stage/rule | Jobs for three samples | CPUs/job | Memory/job | Initial walltime |
|---|---:|---:|---:|---:|
| pbmm2 alignment | 3 | 16 | 48 GB | 48 h |
| DeepVariant | 3 | 16 | 64 GB | 48 h |
| Sawfish | 3 | 16 | 32 GB | 24 h |
| Sniffles2 | 3 | 8 | 16 GB | 12 h |
| pbsv | 3 | 8 | 32 GB | 24 h |
| cuteSV | 3 | 8 | 16 GB | 12 h |
| TRGT | 3 | 8 | 16 GB | 12 h |
| GLnexus cohort | 1 | 8 | 192–200 GB for the validated three-sample WGS workload | 24 h |
| SV cohort tables | 1 | 1 | 8 GB | 4 h |
| TR cohort matrix | 1 | 1 | 8 GB | 4 h |

Do not enable a track merely to maximize parallelism. Every enabled tool needs
a validated executable/image and compatible scientific resources.

## 13. Key Slurm parameters

| Parameter | Meaning | Example |
|---|---|---|
| `--array=1-3%2` | Three tasks, no more than two running | alignment |
| `--cpus-per-task` | CPUs visible to one program | pbmm2 `-j` |
| `--mem` | Total memory for one job | `48G` |
| `--time` | Hard walltime | `2-00:00:00` |
| `--partition` | Site queue | administrator value |
| `--account` | Charging/fair-share account | administrator value |
| `--output/--error` | Persistent separate logs | `%A_%a` for arrays |
| `--export` | Values passed to a batch script | run/reference paths |
| `--dependency=afterok:<id>` | Start only after success | optional chaining |
| `--parsable` | Return a machine-readable job ID | submission scripts |

Use `%A` for array job ID, `%a` for array task ID, and `%j` for a normal job ID.

## 14. Cohort handoff

After all three per-sample outputs pass validation, create a long-form TSV with
this exact header:

```text
sample	track	state	source_path	index_path	source_tool	source_version	reference_build	catalog_id
```

Small-variant rows use each DeepVariant gVCF and TBI, state `CALLED`, actual
DeepVariant version, and build `GRCh38`. SV rows should prefer a per-sample
harmonized SV VCF. TR rows require the same catalog ID for all callable samples.

Create a separate cohort config:

```yaml
cohort:
  enabled: true
  cohort_id: EXAMPLE_3SAMPLE_RC2
  input_manifest: <RUN_ROOT>/configs/cohort_inputs.tsv
  overwrite: false
  small_variants:
    enabled: true
    glnexus_executable: glnexus_cli
    bcftools_executable: bcftools
    preset: DeepVariantWGS
    threads: 8
    # Explicit run-specific cap. The validated three-sample WGS run peaked
    # near 153 GB; 192 GB leaves headroom for that workload.
    memory_gb: 192
    runtime_minutes: 1440
  sv:
    enabled: false
  tr:
    enabled: false
```

Enable SV/TR cohort tracks only when upstream rows are valid. GLnexus must use
the validated `hifivar-glnexus-1.4.1` environment or an equivalent pinned
deployment visible to compute nodes.

The GLnexus value is not a universal cohort formula. Re-estimate it from sample
count, gVCF density, site limits, and recorded peak RSS before scaling.

## 15. Failure handling and resumability

1. Inspect the failed rule log and `sacct` state/exit code.
2. Correct only run-local config, resource, environment, or path problems.
3. Never delete original uBAM, alignment BAM, VCF, gVCF, or validated outputs.
4. Do not use `--keep-going` during the first single-sample acceptance run.
5. Use `--rerun-incomplete` only after understanding the previous failure.
6. Preserve quarantine and logs when validation fails after a tool succeeds.
7. A missing sample/caller is `FAILED`, `NOT_RUN`, or `MISSING_INPUT`, never a
   negative biological call.

Useful inspection:

```bash
snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile "$RUN_ROOT/configs/calling.effective.yaml" \
  --summary

find "$RUN_ROOT/logs" -type f -size +0c -print
```

## 16. Recommended rollout

1. 0.1% deterministic SAMPLE_A uBAM smoke.
2. Full SAMPLE_A alignment and DeepVariant.
3. Add one read-SV caller, then other validated callers.
4. Add TRGT only with an explicitly compatible GRCh38 catalog.
5. Align and call SAMPLE_B and SAMPLE_C.
6. Run GLnexus.
7. Add SV/TR cohort tracks when their upstream contracts are complete.
8. Review Slurm efficiency before scaling to a larger cohort.

This order validates deployment before consuming three-sample WGS resources.
