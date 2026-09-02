# HiFiVar 0.1.0rc4 complete quick start

[English](quickstart.md) | [简体中文](zh_CN/quickstart.md)

This guide covers installation, configuration, dry-run, single-sample and
multi-sample execution on Linux/HPC. Start with tiny data before chromosome- or
whole-genome-scale execution.

## 1. Current execution and FASTQ boundary

HiFiVar is a core workflow plus separately deployed scientific tools. The wheel
does not bundle callers, references, databases, catalogs, truth sets or WGS
data. The packaged small-variant, read-SV and TR workflows consume an existing
indexed BAM/CRAM. The Python API supports pbmm2 alignment planning/execution,
but `0.1.0rc4` has no unified `hifivar run` command and no packaged Snakemake
alignment rule. A raw FASTQ/uBAM therefore needs an explicit upstream alignment
step. FASTQ is also the input boundary for the separate hifiasm assembly branch.

## 2. Install and verify the release

Download all three assets from the
[`v0.1.0-rc4` release](https://github.com/lg10is1/hifivar/releases/tag/v0.1.0-rc4):

```bash
curl -fLO https://github.com/lg10is1/hifivar/releases/download/v0.1.0-rc4/hifivar-0.1.0rc4-py3-none-any.whl
curl -fLO https://github.com/lg10is1/hifivar/releases/download/v0.1.0-rc4/hifivar-0.1.0rc4.tar.gz
curl -fLO https://github.com/lg10is1/hifivar/releases/download/v0.1.0-rc4/SHA256SUMS
sha256sum -c SHA256SUMS

python3 -m venv hifivar-env
source hifivar-env/bin/activate
python -m pip install --upgrade pip
python -m pip install './hifivar-0.1.0rc4-py3-none-any.whl[workflow]'

hifivar --version
hifivar --help
hifivar config validate
hifivar doctor
```

The expected version is `hifivar 0.1.0rc4`. HiFiVar is not yet published to
PyPI, Bioconda or conda-forge; use the GitHub Release assets or a tagged source
checkout. See [Installation](installation.md) for source and Conda options.

## 3. Create a project directory

```text
analysis/
├── config.yaml
├── effective_config.yaml
├── samples.tsv
├── references/
├── containers/
├── databases/
├── work/
├── results/
└── logs/
```

Keep original reads and alignments outside disposable `work/` directories.

## 4. Prepare the reference

Use an uncompressed FASTA with an existing `.fai` index and an explicit build.
HiFiVar never silently converts `chr1` to `1` or vice versa.

```bash
samtools faidx /data/reference/GRCh38.fa
```

Confirm that the BAM/CRAM, reference, catalogs and downstream resources use the
same build and contig naming. Caller-specific dictionaries and indexes remain
owned by their deployments.

## 5. Prepare the sample sheet

Create a UTF-8 tab-separated `samples.tsv`:

```text
sample_id	input	input_type	sex
SAMPLE01	/data/alignments/SAMPLE01.bam	bam	unknown
SAMPLE02	/data/alignments/SAMPLE02.cram	cram	unknown
```

Each BAM must be non-empty, coordinate sorted where required and indexed. CRAM
requires its matching reference and CRAI. Input sample identity, read groups and
reference compatibility must be checked before expensive callers are enabled.

For raw HiFi data, align and index first. Do not label an unaligned PacBio uBAM
as a reusable aligned BAM. For assembly, derive or provide HiFi FASTQ explicitly.

## 6. Create and validate the configuration

From a source checkout, copy the minimal example. From an installed wheel, use
it as a schema reference and create equivalent local files.

```bash
cp examples/minimal/config.yaml config.yaml
cp examples/minimal/samples.tsv samples.tsv
```

Replace every example `/data/...` and `/work/...` path with a location valid on
your system. The minimal example enables no biological caller by default and is
safe for validation. Global CLI options precede the subcommand:

```bash
hifivar --config config.yaml config validate
hifivar --config config.yaml config dump-effective --output effective_config.yaml
```

Review `effective_config.yaml`; archive it with the run. Never put passwords,
tokens or credentials in the config.

## 7. Discover and dry-run the packaged workflow

```bash
WORKFLOW_ROOT="$(python -c 'from hifivar.package_resources import installed_workflow_root; print(installed_workflow_root())')"
test -f "$WORKFLOW_ROOT/Snakefile"

snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile effective_config.yaml \
  --cores 1 \
  --dry-run \
  --printshellcmds
```

A successful infrastructure dry-run does not prove that an external caller is
installed or scientifically configured.

## 8. Enable branches incrementally

Provision and verify only the tools needed by a branch, using the pinned
[deployment matrix](deployment.md). Do not replace a validated version with
`latest`. Enable one branch, regenerate the effective config, and dry-run again.

Typical artifact families remain separate:

- DeepVariant: `results/small/<sample>.small.vcf.gz` and gVCF;
- read-based SV: caller-specific VCFs plus harmonized evidence;
- TRGT: `results/tr/<sample>.tr.vcf.gz`;
- assembly: hifiasm outputs, then PAV/SVIM-asm artifacts;
- cohort: GLnexus small-variant output and separate SV/TR cohort tables;
- review, annotation, benchmark and report: optional downstream branches.

Caller support counts, harmonization, manual review and annotation impact are
evidence metadata, not truth, confidence or clinical pathogenicity.

## 9. Execute one sample with Bash

After a successful dry-run:

```bash
snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile effective_config.yaml \
  --cores 8 \
  --printshellcmds \
  --rerun-incomplete
```

Use a data-disk `runtime.tmpdir`. DeepVariant temporary files are isolated under
`<tmp-root>/deepvariant/<sample>/tmp`. A file-descriptor limit of at least 4096
is required; 65536 is recommended where site policy permits. See the
[Chinese single-sample Bash guide](zh_CN/single_sample_bash.md) for a detailed
operational template.

## 10. Execute multiple samples on Slurm

First prove the DAG locally, then use a site-approved Snakemake profile:

```bash
snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile effective_config.yaml \
  --profile /path/to/site-profile \
  --rerun-incomplete \
  --printshellcmds
```

The profile must map `threads`, `mem_mb` and `runtime_min`. DeepVariant sample
concurrency is controlled by `small.max_concurrent_samples`; increase it only
after tiny tests prove isolated temporary directories and adequate resources.
GLnexus memory must be explicitly sized for the cohort. See
[Multi-sample execution with Slurm](slurm_multi_sample.md) and the
[Chinese Slurm guide](zh_CN/slurm_multi_sample.md).

## 11. Monitor, restart and inspect outputs

- preserve complete stdout/stderr and the effective config;
- use `--rerun-incomplete` after diagnosing a failure;
- do not delete original FASTQ/BAM/CRAM/VCF files;
- do not use `--keep-going` until partial-track behavior is understood;
- validate BGZF/TBI, sample names, contigs and expected output contracts;
- archive tool versions, container digests, reference checksums and Git SHA.

See [Outputs](outputs.md), [Troubleshooting](troubleshooting.md) and
[Linux/HPC execution](linux_hpc.md).

## 12. What is and is not validated

`0.1.0rc4` package installation, packaged resources, Snakemake regressions and
selected real Linux/HPC tool paths have passed release validation. A new site,
reference, container, database or scheduler configuration still requires its
own tiny real-tool validation before WGS execution. The software is for research
use and does not provide clinical interpretation.
