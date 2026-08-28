# Linux/HPC execution

Linux/HPC is the production platform. The recommended deployment is a Conda or
virtual-environment HiFiVar core plus tool-specific Conda environments and
Apptainer images.

## Site prerequisites

- Python 3.10–3.12 and Mamba/Conda or a virtual environment;
- Snakemake 8 or 9;
- Apptainer where containerized callers are used;
- a shared filesystem visible on compute nodes;
- scheduler-specific CPU, memory, walltime and log configuration;
- versioned references, indexes, catalogs, truth resources and databases.

HiFiVar does not install scheduler profiles or modify module/Conda base state.

## Recommended layout

```text
project/
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

Keep source reads and alignments outside disposable `work/` directories.

## Environment verification

```bash
conda activate hifivar
hifivar --version
hifivar doctor
snakemake --version
apptainer --version
ulimit -n
```

DeepVariant sharded execution should have a file-descriptor limit of at least
4096; 65536 is recommended where site policy permits. HiFiVar reports a clear
preflight error for insufficient limits.

DeepVariant temporary storage is isolated per sample at
`<runtime.tmpdir-or-workdir>/deepvariant/<sample>/tmp`. Select a shared or local
data-disk root with enough capacity; do not point `runtime.tmpdir` at a small
system `/tmp`. Native execution receives this path as `TMPDIR`; Docker and
Apptainer also receive an explicit writable bind and container environment.

## Workflow discovery after wheel/Conda installation

```bash
WORKFLOW_ROOT="$(python -c 'from hifivar.package_resources import installed_workflow_root; print(installed_workflow_root())')"
test -f "$WORKFLOW_ROOT/Snakefile"
```

## Scheduler use

Scheduler flags and profiles are site-specific and are not hard-coded in the
repository. First prove the DAG locally:

```bash
snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile effective_config.yaml \
  --cores 1 \
  --dry-run \
  --printshellcmds
```

Then invoke the site-approved profile, for example:

```bash
snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile effective_config.yaml \
  --profile /path/to/site-profile \
  --rerun-incomplete \
  --printshellcmds
```

The profile must map `threads`, `mem_mb` and `runtime_min` without weakening
rule-level resource requests.

The workflow registers one global `deepvariant_slots` unit by default, so only
one DeepVariant sample runs at a time. `small.max_concurrent_samples` controls
that capacity. Increase it only after a tiny DAG run confirms independent
sample temporary directories and the site can sustain the multiplied resource
load.

## External tools

Review `docs/deployment.md` before enabling a branch. Never replace a validated
tool tag with `latest`. Record executable versions, container digests, effective
config, Git SHA, reference/database versions and complete stdout/stderr.

PAV uses its validated independent Apptainer deployment. The HiFiVar core
container is not intended to launch another container from inside itself.

## Validation sequence

1. core CLI/config smoke;
2. Snakemake dry-run;
3. dependency-free Phase 0 marker execution;
4. fake/mock regression;
5. tiny real-tool run for each enabled branch;
6. chromosome-scale test;
7. WGS only after all prior gates pass.

Results from a different reference, image, executable, cache or database do not
automatically validate a new deployment.
