FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive

# Base tools + deps needed to build STAR
RUN apt-get update && apt-get install -y \
    python3 python3-pip \
    r-base \
    git curl wget ca-certificates \
    build-essential \
    zlib1g-dev libbz2-dev liblzma-dev \
    tar gzip pigz \
    fastqc default-jre \
    vim-common \
    && rm -rf /var/lib/apt/lists/*

# STAR (build from source tag)
ARG STAR_VERSION=2.7.11b
RUN mkdir -p /opt && \
    curl -L --fail -o /tmp/STAR.tar.gz \
      https://github.com/alexdobin/STAR/archive/refs/tags/${STAR_VERSION}.tar.gz && \
    tar -xzf /tmp/STAR.tar.gz -C /opt && \
    make -C /opt/STAR-${STAR_VERSION}/source STAR && \
    ln -s /opt/STAR-${STAR_VERSION}/source/STAR /usr/local/bin/STAR && \
    rm -f /tmp/STAR.tar.gz

# Python tools
RUN pip3 install --no-cache-dir \
    snakemake==7.32.4 \
    pulp==2.7.0 \
    multiqc \
    cutadapt==4.9

WORKDIR /work
