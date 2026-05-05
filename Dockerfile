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
# Copy dependency metadata first (better caching)
# -------------------------
COPY pyproject.toml uv.lock README.md ./

# -------------------------
# Copy source
# -------------------------
COPY src/ ./src/

# -------------------------
# Sync deps + install your project (from lock)
# uv creates .venv/ in the project root (/app/.venv).
# --frozen ensures uv.lock is respected exactly.
# --no-dev avoids installing dev extras into the runtime image.
# -------------------------
RUN uv sync --frozen --no-dev

# Put the project venv on PATH so the console script is invokable directly.
ENV PATH="/app/.venv/bin:$PATH"

# -------------------------
# Entry point
# -------------------------
# Use your console script name from [project.scripts]
# Example:
# [project.scripts]
# atomization_scorer = "atomization_scorer.cli:main"
ENTRYPOINT ["atomization_scorer"]

