FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    r-base \
    git \
    curl \
    wget \
    ca-certificates \
    build-essential \
    tar \
    gzip \
	fastqc \
	default-jre \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir \
    snakemake \
    pulp==2.7.0 \
	multiqc

WORKDIR /work
