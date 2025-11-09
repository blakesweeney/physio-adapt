# Multi-stage build for temperature adaptation pipeline
FROM python:3.13-slim AS base

# Install system dependencies
# Note: infernal package includes both Infernal tools (cmscan, cmalign, cmfetch, etc.)
# and Easel utilities (esl-sfetch, esl-alimask, etc.)
# HMMER is also installed as it provides additional easel utilities
RUN apt-get update && apt-get install -y \
    curl \
    unzip \
    infernal \
    hmmer \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Verify easel utilities are available
RUN esl-sfetch -h > /dev/null 2>&1 || { echo "ERROR: esl-sfetch not found"; exit 1; }

# Install NCBI datasets tool
RUN curl -L https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets -o /usr/local/bin/datasets \
    && chmod +x /usr/local/bin/datasets

# Install uv for Python package management
RUN pip install --no-cache-dir uv

# Set working directory
WORKDIR /app

# Copy Python project files
COPY pyproject.toml uv.lock ./

# Install Python dependencies using uv
RUN uv pip install --system -r pyproject.toml

# Copy Python scripts and make them executable
COPY bin/*.py /usr/local/bin/
RUN chmod +x /usr/local/bin/*.py

# Set PATH to include our scripts
ENV PATH="/usr/local/bin:${PATH}"

# Default command
CMD ["/bin/bash"]
