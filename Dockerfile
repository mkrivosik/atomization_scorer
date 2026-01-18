# -------------------------
# Base image
# -------------------------
FROM python:3.11-slim

LABEL author="Matej Krivosik"
LABEL maintainer="krivosik7@uniba.sk"
LABEL version="1.0"
LABEL description="Atomization Scorer"

# -------------------------
# System dependencies
# -------------------------
RUN apt-get update && apt-get install -y \
    minimap2 \
    mash \
    && rm -rf /var/lib/apt/lists/*

# -------------------------
# Working directory
# -------------------------
WORKDIR /app

# -------------------------
# Install uv
# -------------------------
RUN pip install --no-cache-dir uv

# -------------------------
# Create venv (recommended in containers)
# -------------------------
ENV VIRTUAL_ENV=/opt/venv
RUN uv venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# -------------------------
# Copy dependency metadata first (better caching)
# -------------------------
COPY pyproject.toml uv.lock README.md ./

# -------------------------
# Copy source
# -------------------------
COPY src/ ./src/

# -------------------------
# Sync deps + install your project (from lock)
# --frozen ensures uv.lock is respected exactly
# --no-dev avoids installing dev extras into runtime image
# -------------------------
RUN uv sync --frozen --no-dev

# -------------------------
# Entry point
# -------------------------
# Use your console script name from [project.scripts]
# Example:
# [project.scripts]
# atomization-scorer = "atomization_scorer.cli:main"
ENTRYPOINT ["atomization-scorer"]

