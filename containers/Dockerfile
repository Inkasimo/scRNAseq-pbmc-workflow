FROM rocker/r-ver:4.3.3
ENV DEBIAN_FRONTEND=noninteractive

ENV RENV_CONFIG_REPOS_OVERRIDE=https://packagemanager.posit.co/cran/2024-01-01
ENV RENV_PATHS_LIBRARY=/opt/renv/library
ENV RENV_PATHS_CACHE=/opt/renv/cache
ENV RENV_CONFIG_CACHE_SYMLINKS=FALSE

RUN apt-get update && apt-get install -y \
    python3 python3-pip \
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

# Install renv globally
RUN R -q -e "install.packages('renv', repos=Sys.getenv('RENV_CONFIG_REPOS_OVERRIDE'))"

# Python tools
RUN pip3 install --no-cache-dir snakemake==7.32.4 pulp==2.7.0 multiqc cutadapt==4.9

WORKDIR /work

# Copy only what renv needs to restore
COPY renv.lock /work/renv.lock
COPY renv/vendor/ /work/renv/vendor/
COPY renv/activate.R /work/renv/activate.R
COPY .Rprofile /work/.Rprofile

# Restore into /opt/renv/library
RUN R -q -e "renv::restore(project='/work', lockfile='/work/renv.lock', prompt=FALSE, rebuild=TRUE)"
CMD ["bash"]
