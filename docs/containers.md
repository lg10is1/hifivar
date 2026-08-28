# Conda, Docker and Apptainer

HiFiVar uses a core-plus-tools deployment. Core distributions contain Python,
HiFiVar, Snakemake and packaged workflow resources. Heavy callers, references
and databases remain external.

## Conda/Mamba

From a source checkout:

```bash
mamba env create -f environment.yml
conda activate hifivar
hifivar --version
```

Build the local noarch recipe:

```bash
mamba create -n hifivar-build -c conda-forge conda-build conda-verify
conda activate hifivar-build
conda build conda-recipe
```

The recipe is a local release artifact. It has not been accepted into Bioconda
or conda-forge, so public channel install commands must not be advertised yet.

## Docker core image

Build from the repository root:

```bash
docker build --pull -t hifivar:0.1.0rc2 .
docker run --rm hifivar:0.1.0rc2 --version
docker run --rm hifivar:0.1.0rc2 config validate
```

The image runs as non-root UID/GID 10001. Mount a writable work directory:

```bash
docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v "$PWD:/work" \
  hifivar:0.1.0rc2 \
  --config /work/config.yaml config validate
```

The image entrypoint is `hifivar`. To invoke bundled Snakemake, override it:

```bash
docker run --rm \
  --entrypoint /bin/sh \
  -v "$PWD:/work" \
  hifivar:0.1.0rc2 \
  -c 'WORKFLOW_ROOT=$(python -c "from hifivar.package_resources import installed_workflow_root; print(installed_workflow_root())"); snakemake --snakefile "$WORKFLOW_ROOT/Snakefile" --configfile /work/effective_config.yaml --cores 1 --dry-run'
```

Docker support remains `NOT_RUN` until independently built on Linux or CI.

## Apptainer core image

Build from the repository root so `%files` sources resolve correctly:

```bash
apptainer build hifivar_0.1.0rc2.sif containers/hifivar.def
apptainer test hifivar_0.1.0rc2.sif
apptainer run hifivar_0.1.0rc2.sif --version
```

On sites where unprivileged build is disabled, build the OCI image elsewhere
and convert/pull it according to site policy. Never request administrator access
automatically.

Bind the project and data explicitly:

```bash
apptainer run \
  --bind "$PWD:/work" \
  --pwd /work \
  hifivar_0.1.0rc2.sif \
  --config /work/config.yaml config validate
```

The core SIF does not contain PAV, DeepVariant, ANNOVAR databases or other heavy
tools. Do not attempt nested PAV/DeepVariant containers from the core SIF.

## Image identity

For publication, record the Docker/OCI digest and SIF SHA256 in release notes.
An image tag alone is not immutable. Rebuild and revalidate after any base image,
Python dependency, workflow or tool change.
