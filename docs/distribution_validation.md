# Public distribution validation

Candidate: `HiFiVar 0.1.0rc2`

This document separates locally completed package checks from distribution
paths that still require a clean Linux host. `NOT_RUN` is not a failure and must
not be rewritten as PASS without command/output evidence.

## Local public-copy results

| Check | Result |
|---|---|
| Public source/privacy audit | PASS |
| Production `src/hifivar` parity with validated baseline | PASS |
| Production `workflow` parity with validated baseline | PASS |
| Distribution targeted tests | 31 passed |
| Full pytest | 888 passed, 2 skipped |
| Snakemake regression | 29 passed |
| `compileall src` | PASS |
| wheel/sdist build | PASS |
| clean wheel install | PASS |
| installed CLI/version/config/resources | PASS |
| Apache-2.0 package metadata | PASS |
| Conda recipe render/build | NOT_RUN: local legacy Conda proxy/TLS failure; base was not modified |
| Docker core build | NOT_RUN: Docker unavailable locally |
| Apptainer core build | NOT_RUN: Apptainer unavailable locally |

Artifact hashes change whenever the source distribution is rebuilt. Generate a
separate `SHA256SUMS` release asset only from the final committed/tagged source;
do not embed self-referential archive hashes in packaged documentation.

## Linux revalidation required

Use a fresh clone or synchronized public tree with no inherited `PYTHONPATH`.

```bash
python3 -m venv /tmp/hifivar-public-venv
source /tmp/hifivar-public-venv/bin/activate
python -m pip install --upgrade pip
python -m pip install '.[dev,workflow]'
python -m pytest -p no:cacheprovider
python -m pytest tests/integration -k snakemake -p no:cacheprovider
python -m compileall src
python -m build
python scripts/audit_public_release.py . --artifacts-dir dist
```

Conda:

```bash
mamba env create -f environment.yml
conda run -n hifivar hifivar --version
conda run -n hifivar hifivar config validate
conda build -c conda-forge -c bioconda conda-recipe
```

Docker, when authorized:

```bash
docker build --pull -t hifivar:0.1.0rc2 .
docker run --rm hifivar:0.1.0rc2 --version
docker run --rm hifivar:0.1.0rc2 config validate
docker run --rm --entrypoint /usr/bin/id hifivar:0.1.0rc2 -u
```

Apptainer:

```bash
apptainer build hifivar_0.1.0rc2.sif containers/hifivar.def
apptainer test hifivar_0.1.0rc2.sif
apptainer run hifivar_0.1.0rc2.sif --version
apptainer run hifivar_0.1.0rc2.sif config validate
sha256sum hifivar_0.1.0rc2.sif
```

Documentation smoke must execute the README and quick-start commands rather
than merely inspect them.

## Publication blockers

- confirm the canonical repository URL remains
  `https://github.com/lg10is1/hifivar` (completed for this candidate);
- complete Linux Conda and available container delta validation;
- rebuild artifacts from the final commit and publish their new hashes;
- do not claim PyPI, Bioconda, conda-forge or a container registry until each
  upload and clean downstream install has been verified.
