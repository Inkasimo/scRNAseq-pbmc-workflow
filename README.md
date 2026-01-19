# scRNA-seq PBMC Workflow

Production-style scRNA-seq pipeline using **Docker + Snakemake**.

Focus:
- reproducible execution
- explicit workflow structure
- containerized dependencies

Dataset:
- 10x Genomics PBMC 5k
- donors 1–4
- 3′ Gene Expression

This repository is intended as a **technical portfolio / learning project**.

## Requirements

### Required
- Docker

All core tools are provided inside the Docker image e.g.

-Snakemake

-STAR

-FastQC

-MultiQC

### Optional (wrapper only)
- Python ≥3.9

## Build Docker image

`docker build -t scrnaseq-workflow -f containers/Dockerfile .`

## How to run

You can run this workflow in **two ways**:

1. **Directly with Snakemake (inside Docker)**  
   – minimal, explicit, closest to “pure Snakemake”

2. **Via the Python wrapper (`run_analysis.py`)**  
   – recommended for interactive use  
   – enforces explicit execution choices  
   – avoids accidental full runs

### Option A: Direct Snakemake (no wrapper)

#### Dry run

```bash
docker run --rm -it \
  -v "$(pwd)":/work \
  -w /work \
  scrnaseq-workflow \
  snakemake -n -p
  ```


#### Run full workflow:

```bash
docker run --rm -it \
  -v "$(pwd)":/work \
  -w /work \
  scrnaseq-workflow \
  snakemake --cores 8 --rerun-incomplete
  ```


### Run alignment for one donor:

```bash
  docker run --rm -it \
  -v "$(pwd)":/work \
  -w /work \
  scrnaseq-workflow \
  snakemake results/alignment/starsolo/donor1/starsolo.done
  ```


### Option B: Python wrapper (recommended)

The wrapper provides:
- section-based execution
- explicit control (no accidental `rule all`)
- simpler CLI for interactive work


#### Wrapper requirements (host-side only)
```bash
pip install -r requirements-dev.txt
```

This installs:
- pyyaml (used only by the wrapper)

If you do not use the wrapper, you can ignore this step entirely.


List available sections:
```bash
python run_analysis.py --list-sections
```

List donors:

```bash
python run_analysis.py --list-donors
```

Download data:

```bash
python run_analysis.py download_data --cpus 8 --cores 8
```

QC only:

```bash
python run_analysis.py qc --cpus 8 --cores 8
```

Align all donors:

```bash
python run_analysis.py align --donor all --cpus 8 --cores 8 --set-threads starsolo=8
```

Dry run:

```bash
python run_analysis.py all --dry-run
```

### CPUs vs cores (Docker vs Snakemake)

The wrapper exposes two separate parameters:

- `--cpus`  
  Controls the **Docker container CPU limit** (`docker run --cpus`).  
  This is a hard upper bound on how many CPU cores *all processes inside the container combined* may use.

- `--cores` / `-j`  
  Controls **Snakemake scheduling** (how many rules/jobs may run in parallel).

- `-j / --jobs` (Snakemake):  
  Maximum number of rules (jobs) that may run concurrently.

In practice:
- `--cores` limits **CPU usage**
- `-j` limits **job-level parallelism**
- A single job may use multiple threads (e.g. STARsolo)

For STARsolo-heavy workloads, it is common to use:
- `-j 1`
- `--cores = --set-threads starsolo`

If the resources allow you can run eg alignment in parallel: 

```bash
python run_analysis.py align \
  --donor all \
  --cpus 16 \
  --cores 16 \
  -j 2 \
  --set-threads starsolo=8
```

## Barcode whitelist

The 10x Genomics barcode whitelist is bundled under:

resources/barcodes/3M-3pgex-may-2023_TRU.txt

This avoids reliance on unstable upstream URLs (frequent HTTP 403 errors)
and ensures fully reproducible execution.

## Repository structure

```text
containers/          # Dockerfile(s) for reproducible execution
workflow/            # Snakemake workflow (rules, DAG)
config/              # User-editable configuration (config.yaml)
resources/            # Static resources bundled with the workflow
  barcodes/           # 10x barcode whitelist(s)
data/                # Input data and references (not versioned)
  raw/                # FASTQ files (downloaded or user-provided)
  ref/                # Reference genome, GTF, STAR index
results/             # Outputs and logs (not versioned)
  qc/                 # FastQC / MultiQC reports
  alignment/           # STARsolo outputs
  logs/                # Execution logs
docs/                # Documentation (user manual, report, notes)
scripts/             # Helper scripts (currently empty / reserved)
run_analysis.py      # Optional Python wrapper for section-based execution
```

## Notes
- Large data files are excluded via `.gitignore`
- FASTQ downloading is controlled via `io.download_fastqs`
- STAR index requires ~25–30 GB RAM


