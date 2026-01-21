#Structure

Purpose (what / what not)

Requirements (disk, RAM, WSL note)

Quick start (minimal)

Configuration (config.yaml)

Running specific steps (many examples)

Outputs (directory tree)

Troubleshooting (OOM, rerun-incomplete, etc.)

## 1. Purpose

One paragraph:

What the pipeline does

What it does not do

Example:

This pipeline processes 10x Genomics GEM-X 3′ v4 scRNA-seq data using STARsolo to produce gene–cell count matrices. It does not perform downstream analysis or visualization.

2. Requirements

Bullet list:

Docker

(Windows) WSL2 with increased memory (link to docs/wslconfig.example)

Disk space note (e.g. “~50–70 GB recommended”)

3. Quick start

Minimal, copy-pasteable:

docker build -t scrnaseq-workflow .

docker run --rm -it \
  -v "$(pwd)":/work \
  -w /work \
  scrnaseq-workflow \
  snakemake -j 1

4. Configuration

Explain only:

where config.yaml lives

what users are expected to edit:

FASTQ URLs / paths

number of threads (optional)

Do not document every config key.

5. Running specific steps (very useful)

Examples only:

# Build reference only
snakemake data/ref/star_index.done

# Run alignment for one donor
snakemake results/alignment/starsolo/donor1/starsolo.done


This saves users time.

6. Outputs

High-level only:

results/qc/ → FastQC / MultiQC

results/alignment/starsolo/<donor>/Solo.out/ → count matrices

logs/ → execution logs

## Windows (WSL2) memory configuration

STAR genome indexing for GRCh38 requires ~25–30 GB RAM.
The default WSL2 memory limit is often insufficient and may cause
out-of-memory errors during index construction.

A template configuration file is provided at:

docs/wslconfig.example

### Steps

1. Copy the file to:
   C:\Users\<your-username>\.wslconfig
2. Adjust values if needed
3. Restart WSL:
   wsl --shutdown
4. Restart Docker Desktop

### Example configuration

[wsl2]
memory=32GB
processors=8
swap=16GB


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

### Run the full workflow (explicit)
python run_analysis.py all \
  --cpus 8 --cores 8


### Run trimming only (all donors)
python3 run_analysis.py trim \
  --cpus 8 --cores 8

### Run trimming + trimmed QC (FastQC + MultiQC)
python3 run_analysis.py trim_and_qc \
  --cpus 8 --cores 8


Outputs:

data/trimmed/{donor}/trim.done

results/qc/fastqc/trimmed/{donor}/fastqc.done

results/qc/multiqc/trimmed/multiqc_report.html

### Run alignment on trimmed reads
All donors
python3 run_analysis.py align \
  --donor all \
  --trimmed \
  --cpus 8 --cores 8

### Selected donors
python3 run_analysis.py align \
  --donor donor1 \
  --donor donor3 \
  --trimmed \
  --cpus 8 --cores 8

###Run full workflow in trimmed mode

Runs trimming, trimmed QC, and trimmed alignment.

python3 run_analysis.py all \
  --trimmed \
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




## Barcode whitelist

The 10x Genomics barcode whitelist (3M-3pgex-may-2023_TRU.txt) is included directly in this repository under resources/barcodes/.
This is intentional: upstream 10x download URLs for barcode whitelists are unstable and frequently return HTTP 403 errors, which breaks fully automated workflows.
Bundling the whitelist ensures the pipeline runs reproducibly without external dependencies.

### Running via wrapper
./run.sh build
./run.sh dry
./run.sh run -p -j 1
