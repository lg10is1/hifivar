# syntax=docker/dockerfile:1
FROM python:3.12-slim-bookworm

ARG HIFIVAR_VERSION=0.1.0rc2
LABEL org.opencontainers.image.title="HiFiVar" \
      org.opencontainers.image.version="${HIFIVAR_VERSION}" \
      org.opencontainers.image.licenses="Apache-2.0" \
      org.opencontainers.image.description="HiFiVar core and Snakemake workflow; external bioinformatics tools are not bundled"

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1

RUN groupadd --gid 10001 hifivar \
    && useradd --uid 10001 --gid hifivar --create-home hifivar

WORKDIR /opt/hifivar-src
COPY pyproject.toml README.md LICENSE MANIFEST.in ./
COPY src ./src
COPY workflow ./workflow
RUN python -m pip install --upgrade pip \
    && python -m pip install ".[workflow]" \
    && hifivar --version \
    && python -c "from hifivar.package_resources import installed_workflow_root; assert installed_workflow_root().joinpath('Snakefile').is_file()"

WORKDIR /work
USER 10001:10001
ENTRYPOINT ["hifivar"]
CMD ["--help"]
