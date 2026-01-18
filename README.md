# scRNA-seq PBMC workflow

Production-style scRNA-seq pipeline demo using Docker + Snakemake.
Focus: execution, structure, reproducibility.

Dataset: 10x PBMC 5k donors 1–4 (3').

## Wrapper script requirements (host-side only)

This repository includes an optional Python wrapper (run_analysis.py) that simplifies running selected parts of the Snakemake workflow inside Docker.

### Important:
These Python requirements are only needed for the wrapper script.
The core workflow itself runs entirely inside Docker and does not 
require any host-side Python dependencies beyond Docker.

## To use the wrapper:

pip install -r requirements-dev.txt


requirements-dev.txt contains only:

pyyaml — used by the wrapper to read config/config.yaml

If you do not use run_analysis.py, 
you can ignore this step completely and run Snakemake directly via Docker.

## scRNAseq PBMC Snakemake Workflow

This repository contains a Dockerized Snakemake workflow for processing 10x Genomics PBMC scRNA-seq data, including FASTQ validation/download, QC (FastQC + MultiQC), STARsolo alignment, and count matrix generation.

The workflow can be run directly with Snakemake or via a Python wrapper that provides section-based execution and cleaner CLI ergonomics.

Requirements
Docker (required)

All core tools are provided inside the Docker image:

Snakemake

STAR

FastQC

MultiQC

Cutadapt (reserved for future use)

Python (wrapper only)

## If you use the Python wrapper (run_analysis.py), install:

pip install -r requirements-dev.txt


This is only required for the wrapper.
The Snakemake workflow itself runs entirely inside Docker.

## Build Docker image

docker build -t scrnaseq-workflow -f containers/Dockerfile .

## Running without the wrapper (direct Snakemake)

### Dry run (inspect DAG only)

docker run --rm -it \
  -v "$(pwd)":/work \
  -w /work \
  scrnaseq-workflow \
  snakemake -n -p

### Run full pipeline (rule all)

docker run --rm -it \
  -v "$(pwd)":/work \
  -w /work \
  scrnaseq-workflow \
  snakemake --cores 8 --rerun-incomplete

### Run a specific target (example: align all 4 donors)


docker run --rm -it \
  -v "$(pwd)":/work \
  -w /work \
  scrnaseq-workflow \
  snakemake -s workflow/Snakefile \
    --cores 8 \
    --rerun-incomplete \
    results/alignment/starsolo/donor{1,2,3,4}/starsolo.done

### Running with the Python wrapper (recommended)

The wrapper enforces explicit execution choices and avoids running the full pipeline by accident.

### Available sections

download_data

download_data_and_qc

qc (FastQC + MultiQC)

ref (barcode whitelist + STAR index)

align

all (full workflow, explicit only)

### Dry run with the wrapper

python run_analysis.py all --dry-run

### Download data
python run_analysis.py download_data \
  --cpus 8 \
  --cores 8

### Download data and run qc (FastQC + MultiQC)

python run_analysis.py download_data_and_qc \
  --cpus 8 \
  --cores 8


### Run QC only (FastQC + MultiQC)

python run_analysis.py qc \
  --cpus 8 --cores 8
  
### Build STAR index

Run reference preparation (whitelist + STAR index)
python run_analysis.py ref \
  --cpus 8 --cores 8

### Run alignments


#### All donors

Run alignment for all donors
python run_analysis.py align \
  --donor all \
  --cpus 8 --cores 8 \
  --set-threads starsolo=8

#### Selected donors
Run alignment for selected donors
python run_analysis.py align \
  --donor donor1 \
  --donor donor3 \
  --cpus 8 --cores 8 \
  --set-threads starsolo=8

Run the full workflow (explicit)
python run_analysis.py all \
  --cpus 8 --cores 8



## Barcode whitelist

The 10x Genomics barcode whitelist (3M-3pgex-may-2023_TRU.txt) is provided under:

resources/barcodes/


This file is required by STARsolo and is referenced via config/config.yaml.

The workflow does not attempt to download barcode files automatically, ensuring reproducibility and avoiding reliance on unstable external URLs.

## Notes

Large input/output files are intentionally excluded via .gitignore

FASTQ downloading is controlled via io.download_fastqs in the config

The workflow is designed for section-based execution, not implicit end-to-end runs


## WSL2 memory configuration (required for STAR index)

STAR genome indexing for GRCh38 requires ~25–30 GB RAM.

If you are running on Windows with WSL2, you must increase WSL memory.
A template is provided at:

    docs/wslconfig.example

To enable it:
1. Copy the file to:
       C:\Users\<your-username>\.wslconfig
2. Adjust values if needed
3. Restart WSL:
       wsl --shutdown
4. Restart Docker Desktop

## Barcode whitelist

The 10x Genomics barcode whitelist (3M-3pgex-may-2023_TRU.txt) is included directly in this repository under resources/barcodes/.
This is intentional: upstream 10x download URLs for barcode whitelists are unstable and frequently return HTTP 403 errors, which breaks fully automated workflows.
Bundling the whitelist ensures the pipeline runs reproducibly without external dependencies.

### Running via wrapper
./run.sh build
./run.sh dry
./run.sh run -p -j 1
