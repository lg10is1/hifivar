# Install HiFiVar 0.1.0rc5

[简体中文](zh_CN/installation.md)

Use Linux or an HPC compute node for native analysis. Python 3.12 with Snakemake 8 or 9 is recommended. Configure caller software, images, reference and catalog separately.

## GitHub wheel

Download the three assets from [the RC5 release](https://github.com/lg10is1/hifivar/releases/tag/v0.1.0-rc5) in a new directory.

~~~bash
base=https://github.com/lg10is1/hifivar/releases/download/v0.1.0-rc5
wheel=hifivar-0.1.0rc5-py3-none-any.whl
for asset in "$wheel" hifivar-0.1.0rc5.tar.gz SHA256SUMS; do
  curl -fL "$base/$asset" -o "$asset" || exit 1
done
sha256sum -c SHA256SUMS || exit 1
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install "./$wheel[workflow]"
hifivar --version
hifivar config validate
hifivar doctor
~~~

Expected version: `hifivar 0.1.0rc5`. Earlier development packages with the same version string have distinct checksums.

## Tagged source

~~~bash
git clone --branch v0.1.0-rc5 --single-branch https://github.com/lg10is1/hifivar.git
cd hifivar
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install '.[workflow]'
~~~

## Installed workflow check

~~~bash
WORKFLOW_ROOT="$(python -c 'from hifivar.package_resources import installed_workflow_root; print(installed_workflow_root())')"
hifivar --preset standard config dump-effective --output effective_config.yaml
snakemake --snakefile "$WORKFLOW_ROOT/Snakefile" --configfile effective_config.yaml --cores 1 --dry-run
~~~

This plan checks workflow infrastructure. Continue with a caller-configured sample in the [quick start](quickstart.md). [Conda and containers](containers.md) describes local builds. To build the recipe locally use `conda build -c conda-forge -c bioconda conda-recipe`. Optional extras `review`, `sv-benchmark` and `alphagenome` support the corresponding features.
