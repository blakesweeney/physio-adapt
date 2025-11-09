# Multi-stage build for temperature adaptation pipeline
FROM python:3.13-slim AS base

# Install system dependencies and build tools
RUN apt-get update && apt-get install -y \
    curl \
    unzip \
    wget \
    build-essential \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Install Infernal from source (includes easel utilities)
# Infernal bundles a copy of the Easel library with esl-* tools
RUN cd /tmp && \
    wget http://eddylab.org/infernal/infernal-1.1.5.tar.gz && \
    tar xzf infernal-1.1.5.tar.gz && \
    cd infernal-1.1.5 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd easel && \
    make install && \
    cd /tmp && \
    rm -rf infernal-1.1.5*

# Verify all required tools are available
RUN cmscan -h > /dev/null 2>&1 && \
    cmalign -h > /dev/null 2>&1 && \
    cmfetch -h > /dev/null 2>&1 && \
    esl-sfetch -h > /dev/null 2>&1 && \
    echo "All Infernal and Easel tools successfully installed"

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
