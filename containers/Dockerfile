FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive

# Pin CRAN snapshot for deterministic restores (renv reads this)
ENV RENV_CONFIG_REPOS_OVERRIDE=https://packagemanager.posit.co/cran/2024-01-01

# Base tools + deps needed to build STAR + R packages
RUN apt-get update && apt-get install -y \
    python3 python3-pip \
    r-base \
    git curl wget ca-certificates \
    build-essential \
    zlib1g-dev libbz2-dev liblzma-dev \
    tar gzip pigz \
    fastqc default-jre \
    vim-common \
    libglpk40 \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
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

# Pin CRAN snapshot for deterministic restores
ENV RENV_CONFIG_REPOS_OVERRIDE=https://packagemanager.posit.co/cran/2024-01-01

WORKDIR /work

# renv restore inputs
COPY renv.lock /work/renv.lock
COPY renv/ /work/renv/
COPY renv/vendor/ /work/renv/vendor/

# Install renv (bootstrap) + restore project library from renv.lock
RUN R -q -e "options(repos=c(CRAN=Sys.getenv('RENV_CONFIG_REPOS_OVERRIDE'))); install.packages('renv')" \
 && R -q -e "renv::restore(prompt=FALSE)"
