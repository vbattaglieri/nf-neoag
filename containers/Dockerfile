# Runtime image for the neo-from-SNV / neo-from-Indel steps.
# The analysis code is the pure-Python-3 `neoag` package (this repo); the only
# external runtime dependency is tcsh, needed to launch the netMHC-4.0 script.
# netMHC is NOT included here (license-restricted) — it is mounted at runtime
# and located via the $NETMHC environment variable.
FROM python:3.11-slim-bookworm

LABEL org.opencontainers.image.title="nf-neoag" \
      org.opencontainers.image.description="Neoantigen prediction (SNV + Indel) runtime — Python3 neoag package + netMHC-4 harness" \
      org.opencontainers.image.licenses="MIT"

# tcsh: required to run the netMHC-4.0 executable (its launcher is a csh script)
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        tcsh \
        ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# Install the neoag package (provides the neo-snv / neo-indel entrypoints).
# The build context is the repo root (see docker build / CI).
COPY pyproject.toml /opt/neoag/pyproject.toml
COPY src /opt/neoag/src
RUN pip install --no-cache-dir /opt/neoag

CMD ["/bin/bash"]
