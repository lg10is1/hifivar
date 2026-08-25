# HiFiVar

HiFiVar is a modular, reproducible workflow framework for PacBio HiFi
whole-genome variant analysis. It orchestrates established bioinformatics
tools; it does not replace their calling algorithms.

The current release candidate is **0.1.0rc1**. It is intended for research use
on Linux/HPC and is not a clinical diagnostic system.

## What is included

- safe Python wrappers and provenance for pbmm2, DeepVariant, four read-SV
  callers, TRGT, HiPhase, hifiasm, PAV, SVIM-asm, Jasmine, Truvari, GLnexus,
  ANNOVAR, VEP, IGV and hap.py boundaries;
- independent small-variant, SV and tandem-repeat artifact families;
- modular Snakemake rules for optional downstream branches;
- JSON/YAML/TSV manifests, QC, benchmark, review and final reporting models;
- Python wheel, source distribution, Conda recipe, Docker core image and
  Apptainer core definition;
- mock/lightweight CI plus separate Linux/HPC real-tool validation boundaries.

External executables, container images, references, databases, caches and WGS
data are **not** bundled in the Python package or core container.

## Current execution boundary

The packaged Snakemake small/SV/TR branches consume an existing indexed
BAM/CRAM. The Phase 2 Python API can plan and execute FASTQ alignment with
pbmm2, but this release candidate does not yet expose one unified `hifivar run`
command or a Snakemake alignment rule. Do not claim a one-command FASTQ-to-all-
variants workflow for `0.1.0rc1`.

## Requirements

- Linux for production execution;
- Python 3.10–3.12;
- Snakemake 8 or 9;
- an indexed reference FASTA;
- caller-specific executables/containers and databases for enabled branches;
- sufficient site-specific CPU, memory, storage and scheduler configuration.

Windows supports the Python, configuration, packaging, mock and Snakemake
dry-run layers, but not the formal external-tool production workflow.

## Install

### GitHub release wheel

```bash
python3 -m venv hifivar-env
source hifivar-env/bin/activate
python -m pip install --upgrade pip
python -m pip install './hifivar-0.1.0rc1-py3-none-any.whl[workflow]'
```

### Source checkout

```bash
git clone https://github.com/lg10is1/hifivar.git
cd hifivar
git checkout v0.1.0-rc1
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install '.[workflow]'
```

### Conda/Mamba source environment

```bash
git clone https://github.com/lg10is1/hifivar.git
cd hifivar
mamba env create -f environment.yml
conda activate hifivar
```

The canonical source repository is
<https://github.com/lg10is1/hifivar>. The project is not yet published to PyPI,
Bioconda or conda-forge.

## Verify the installation

```bash
hifivar --version
hifivar --help
hifivar config validate
hifivar doctor
```

Expected candidate version:

```text
hifivar 0.1.0rc1
```

## Five-minute smoke test

Locate the workflow bundled in the installation:

```bash
WORKFLOW_ROOT="$(python -c 'from hifivar.package_resources import installed_workflow_root; print(installed_workflow_root())')"
hifivar --preset standard config dump-effective --output effective_config.yaml
snakemake \
  --snakefile "$WORKFLOW_ROOT/Snakefile" \
  --configfile effective_config.yaml \
  --cores 1 \
  --dry-run
```

Remove `--dry-run` to execute only the dependency-free infrastructure marker.
This smoke test does not call a biological tool.

## Start a real analysis

1. Read [the quick start](docs/quickstart.md).
2. Copy and edit [the minimal example](examples/minimal/README.md).
3. Provision only the external tools needed by enabled branches using the
   [deployment matrix](docs/deployment.md).
4. Validate paths, indexes, reference build and contig naming.
5. Generate an effective config and inspect a Snakemake dry-run before running.

HiFiVar never silently converts `chr1` to `1`, overwrites raw caller artifacts,
or equates caller support, manual review, annotation impact or benchmark status
with biological truth or pathogenicity.

## Distribution options

- [Installation methods](docs/installation.md)
- [Conda, Docker and Apptainer](docs/containers.md)
- [Linux/HPC execution](docs/linux_hpc.md)
- [Configuration and inputs](docs/quickstart.md)
- [Outputs](docs/outputs.md)
- [Troubleshooting](docs/troubleshooting.md)
- [Validated deployment matrix](docs/deployment.md)
- [Public distribution validation](docs/distribution_validation.md)

## Development

```bash
python -m pip install -e '.[dev,workflow]'
python -m pytest -p no:cacheprovider
python -m compileall src
python -m build
```

Heavy external tools and real WGS data are intentionally excluded from ordinary
CI. Release candidates require an independent Linux/HPC validation pass.

## License

HiFiVar is licensed under the [Apache License 2.0](LICENSE). Third-party tools,
containers, references and databases retain their own licenses and are not
redistributed by HiFiVar.
