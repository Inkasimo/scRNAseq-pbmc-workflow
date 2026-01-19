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

docker build -t scrnaseq-workflow -f containers/Dockerfile .

## How to run

You can run this workflow in **two ways**:

1. **Directly with Snakemake (inside Docker)**  
   – minimal, explicit, closest to “pure Snakemake”

2. **Via the Python wrapper (`run_analysis.py`)**  
   – recommended for interactive use  
   – enforces explicit execution choices  
   – avoids accidental full runs

### Option A: Direct Snakemake (no wrapper)

Dry run:
docker run --rm -it -v "$(pwd)":/work -w /work scrnaseq-workflow \
  snakemake -n -p

Run full workflow:
docker run --rm -it -v "$(pwd)":/work -w /work scrnaseq-workflow \
  snakemake --cores 8 --rerun-incomplete

Run alignment for one donor:
docker run --rm -it -v "$(pwd)":/work -w /work scrnaseq-workflow \
  snakemake results/alignment/starsolo/donor1/starsolo.done

### Option B: Python wrapper (recommended)

The wrapper provides:
- section-based execution
- explicit control (no accidental `rule all`)
- simpler CLI for interactive work

#### Wrapper requirements (host-side only)

pip install -r requirements-dev.txt

This installs:
- pyyaml (used only by the wrapper)

If you do not use the wrapper, you can ignore this step entirely.


List available sections:
python run_analysis.py --list-sections

List donors:
python run_analysis.py --list-donors

Download data:
python run_analysis.py download_data --cpus 8 --cores 8

QC only:
python run_analysis.py qc --cpus 8 --cores 8

Align all donors:
python run_analysis.py align --donor all --cpus 8 --cores 8 --set-threads starsolo=8

Dry run:
python run_analysis.py all --dry-run

## Barcode whitelist

The 10x Genomics barcode whitelist is bundled under:

resources/barcodes/3M-3pgex-may-2023_TRU.txt

This avoids reliance on unstable upstream URLs (frequent HTTP 403 errors)
and ensures fully reproducible execution.

## Repository structure

├── containers/          # Dockerfile(s) for reproducible execution
├── workflow/            # Snakemake workflow (rules, DAG)
├── config/              # User-editable configuration (config.yaml)
├── resources/           # Static resources bundled with the workflow
│   └── barcodes/        # 10x barcode whitelist(s)
├── data/                # Input data and references (not versioned)
│   ├── raw/             # FASTQ files (downloaded or user-provided)
│   └── ref/             # Reference genome, GTF, STAR index
├── results/             # Outputs and logs (not versioned)
│   ├── qc/              # FastQC / MultiQC reports
│   ├── alignment/       # STARsolo outputs
│   └── logs/            # Execution logs
├── docs/                # Documentation (user manual, report, notes)
├── scripts/             # Helper scripts (currently empty / reserved)
└── run_analysis.py      # Optional Python wrapper for section-based execution


## Notes

- Large data files are excluded via `.gitignore`
- FASTQ downloading is controlled via `io.download_fastqs`
- STAR index requires ~25–30 GB RAM
- Windows users: see `docs/wslconfig.example`

