# Installation

HiFiVar targets Linux/HPC for production. Windows supports Python, config,
packaging, fake/mock tests and Snakemake dry-run, but not the formal real-tool
workflow.

## Supported core

- Python 3.10, 3.11 or 3.12;
- Snakemake 8 or 9 for workflow execution;
- PyYAML 6 or newer.

External tools and scientific databases are provisioned separately using
`docs/deployment.md`.

## GitHub release wheel

Download both release artifacts and published SHA256 checksums. Then:

```bash
python3 -m venv hifivar-env
source hifivar-env/bin/activate
python -m pip install --upgrade pip
python -m pip install './hifivar-0.1.0rc4-py3-none-any.whl[workflow]'
```

The expected artifacts are:

- `hifivar-0.1.0rc4-py3-none-any.whl`
- `hifivar-0.1.0rc4.tar.gz`

The project is not yet on PyPI. Do not use `pip install hifivar` until a
separate PyPI publication has been completed and verified.

## Source checkout

```bash
git clone https://github.com/lg10is1/hifivar.git
cd hifivar
git checkout v0.1.0-rc4
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install '.[workflow]'
```

Development installation:

```bash
python -m pip install -e '.[dev,workflow]'
python -m pytest -p no:cacheprovider
```

## Conda/Mamba source environment

```bash
git clone https://github.com/lg10is1/hifivar.git
cd hifivar
mamba env create -f environment.yml
conda activate hifivar
```

The local recipe depends on packages from both conda-forge and Bioconda. Build
it with:

```bash
conda build -c conda-forge -c bioconda conda-recipe
```

This validates the local recipe only; HiFiVar is not yet published to
Bioconda or conda-forge.

## Docker and Apptainer

See `docs/containers.md`. Both definitions package only the core framework and
workflow. Real Docker/Apptainer build results remain pending independent Linux
validation for this distribution delta.

## Installation verification

```bash
hifivar --help
hifivar --version
hifivar config validate
hifivar doctor
python -c 'from hifivar.package_resources import installed_workflow_root; print(installed_workflow_root())'
```

The final command must print a directory containing `Snakefile`, `rules/`,
`scripts/` and `envs/`.

## Version policy

`src/hifivar/__init__.py` is the canonical version source. Package metadata and
the CLI must both report the PEP 440 version `0.1.0rc4`.

## License

HiFiVar is Apache-2.0. External executables, containers, references, databases
and caches are not redistributed and retain their own terms.
