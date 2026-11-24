# SynSlug container image
# Build:  docker build -t synslug .
# Run:    docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/data synslug synslug config.toml

FROM python:3.12-slim

LABEL org.opencontainers.image.title="SynSlug" \
      org.opencontainers.image.description="Synteny visualisation tool" \
      org.opencontainers.image.source="https://github.com/Nepomorpha/SynSlug" \
      org.opencontainers.image.licenses="MIT"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    MPLBACKEND=Agg \
    PYTHONPATH=/opt/synslug

# procps provides `ps`, which Nextflow and similar workflow managers expect
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/synslug

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY synslug.py functions.py loaders.py orthofinder.py ./
COPY scripts/ ./scripts/
COPY example_data/ ./example_data/
COPY examples/*.toml ./examples/

# `synslug` command on PATH
# Give matplotlib a per-user cache folder: under Singularity/Apptainer /tmp is
# the host's shared /tmp, and under `docker run -u` there is no writable $HOME.
RUN printf '%s\n' \
      '#!/bin/sh' \
      ': "${MPLCONFIGDIR:=${TMPDIR:-/tmp}/matplotlib-$(id -u)}"' \
      'export MPLCONFIGDIR' \
      'exec python /opt/synslug/synslug.py "$@"' \
      > /usr/local/bin/synslug \
 && chmod +x /usr/local/bin/synslug

# Configs are usually written with paths relative to the project folder,
# so mount that folder here.
WORKDIR /data

CMD ["synslug"]
