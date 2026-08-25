# Quick start

This guide exercises the public `0.1.0rc1` contract on Linux. Use a tiny test
dataset before chromosome- or whole-genome-scale execution.

## 1. Install the core

From a release wheel:

```bash
python3 -m venv hifivar-env
source hifivar-env/bin/activate
python -m pip install --upgrade pip
python -m pip install './hifivar-0.1.0rc1-py3-none-any.whl[workflow]'
hifivar --version
```

External callers are installed separately. Enable only branches whose tools,
versions, references and databases have been provisioned.

## 2. Prepare the reference

The primary reference is an uncompressed FASTA with an existing `.fai` index.
The configured build is explicit; HiFiVar does not infer or rename contigs.

```bash
samtools faidx /data/reference/GRCh38.fa
```

Caller-specific dictionaries and indexes remain owned by their deployment.

## 3. Prepare an existing-alignment sample sheet

Create a UTF-8 tab-separated `samples.tsv`:

```text
sample_id	input	input_type	sex
HG002	/data/alignments/HG002.bam	bam	male
```

The BAM must be non-empty, coordinate sorted for callers that require it, and
indexed as `/data/alignments/HG002.bam.bai` or a supported alternative. CRAM
requires its matching reference and CRAI. Sample IDs should be ASCII-safe.

## 4. Create a user config

Copy `examples/minimal/config.yaml` and replace every example `/data/...` and
`/work/...` path with locations valid on your system. The example enables no
biological caller by default, so it is safe for config and DAG validation.

```bash
cp examples/minimal/config.yaml config.yaml
cp examples/minimal/samples.tsv samples.tsv
```

Global CLI options precede the subcommand:

```bash
hifivar --config config.yaml config validate
hifivar --config config.yaml config dump-effective --output effective_config.yaml
```

## 5. Inspect the DAG

For a source checkout:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile effective_config.yaml \
  --cores 1 \
  --dry-run \
  --printshellcmds
```

For a wheel/Conda installation:

```bash
WORKFLOW_ROOT="$(python -c 'from hifivar.package_resources import installed_workflow_root; print(installed_workflow_root())')"
snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile effective_config.yaml \
  --cores 1 \
  --dry-run \
  --printshellcmds
```

## 6. Enable one branch at a time

After installing and validating the corresponding external tool, change the
relevant config flag, regenerate `effective_config.yaml`, and repeat dry-run.

Examples:

```yaml
small:
  enabled: true
  execution_mode: apptainer
  deepvariant_image: /containers/deepvariant_1.10.0.sif
```

```yaml
sv:
  enabled: true
  sawfish:
    enabled: true
  sniffles2:
    enabled: false
  pbsv:
    enabled: false
  cutesv:
    enabled: false
```

Do not enable every branch merely because it exists. Each enabled branch must
have its own validated runtime, resource limits and scientific inputs.

## 7. Execute

Use an HPC profile or explicit local resources. A minimal local invocation is:

```bash
snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile effective_config.yaml \
  --cores 8 \
  --printshellcmds \
  --rerun-incomplete
```

Do not use `--keep-going` until partial-track behavior is understood. Never
delete original FASTQ/BAM/CRAM/VCF files to recover a failed workflow.

## FASTQ boundary

`0.1.0rc1` provides pbmm2 alignment through the Python Phase 2 API, but the
packaged DAG does not yet contain an alignment rule. For the public Snakemake
quick start, use an existing indexed BAM/CRAM. FASTQ remains valid for the
separate hifiasm assembly branch. A unified FASTQ-to-calling CLI is future work.
