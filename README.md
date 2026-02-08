# scRNA-seq PBMC Workflow

Production-style scRNA-seq analysis pipeline
using **Docker + Snakemake** with a Python CLI wrapper.

This repository is intended as a **technical portfolio / learning project**.

## Status:
Draft — core pipeline implemented; downstream analyses in progress

## Focus
- Reproducible execution
- Explicit workflow structure
- Containerized dependencies
- Clear separation between engineering (upstream) and analysis (downstream)

## Dataset
- 10x Genomics PBMC 5k
- Donors 1–4
- 3′ Gene Expression

## What this pipeline does

### Upstream

- FASTQ acquisition (download or reuse of existing data)
- Raw and optional trimmed QC (FastQC + MultiQC)
- Optional read trimming (Cutadapt; non-default)
- Reference preparation (barcode whitelist, STAR index)
- STARsolo alignment and gene–cell count matrix generation

### Downstream (work in progress)

- Seurat objects
- Cell-level QC and annotation
- Differential expression (DESeq2) and equivalence testing (TOST)
- Enrichment analysis
- Co-expression network analysis (work in progress)

## Requirements

- Full run requires ~25–30 GB RAM for STAR index and alignment
  and several hours depending on cores.

### Required

- Docker

All core tools are provided inside the Docker image, including:

- Snakemake
- STAR/STARsolo
- FastQC
- MultiQC
- Cutadapt
- Seurat
- DESeq2
- igraph

### Optional (wrapper only)

- Python ≥3.9 (used only for the execution wrapper)

## Build Docker image

```bash
docker build -t scrnaseq-workflow -f containers/Dockerfile .
```

## How to run

This workflow can be run in **two ways**:

1. **Via the Python wrapper** (`run_analysis.py`) – recommended
   
- Section-based execution
- Explicit control over what runs
- Prevents accidental full runs

2. **Directly with Snakemake inside Docker**   

Most users should use the **wrapper**. 

### Option A: Python wrapper (recommended)


#### Wrapper requirements (host-side only)

```bash
pip install -r wrapper-requirements.txt
```
This installs:
- pyyaml (used only by the wrapper)

If you do not use the wrapper, you can ignore this step entirely.


#### Example runs:

****Inspect available sections:****

```bash
python3 run_analysis.py --list-sections
```

****Inspect donors:****

```bash
python3 run_analysis.py --list-donors
```

****Download data:****

```bash
python3 run_analysis.py download_data --cpus 8 --cores 8
```
****QC only:****

```bash
python3 run_analysis.py qc --cpus 8 --cores 8
```

****Align all donors:****

```bash
python3 run_analysis.py align \
  --donor all \
  --cpus 8 --cores 8 \
  -j 1 \
  --set-threads starsolo=8

```

****Dry run (no execution, sanity check):****

```bash
python3 run_analysis.py all --dry-run
```

#### Trimmed flag:

Use the `--trimmed` flag to enable read trimming and run all downstream steps
using trimmed reads instead of raw reads.

```bash
python3 run_analysis.py all --trimmed

```

### Option B: Direct Snakemake (no wrapper)

****Dry run****

```bash
docker run --rm -it \
  -v "$(pwd)":/work \
  -w /work \
  scrnaseq-workflow \
  snakemake -n -p
  ```


****Run a specific target (example: one donor alignment):****

```bash
docker run --rm -it \
  -v "$(pwd)":/work \
  -w /work \
  scrnaseq-workflow \
  snakemake results/alignment/starsolo/donor1/starsolo.done

  ```
  
## Resources and reproducibility

### Barcode whitelist

The 10x Genomics barcode whitelist is bundled directly in the repository:

`resources/barcodes/3M-3pgex-may-2023_TRU.txt`


This avoids reliance on unstable upstream URLs and ensures reproducible execution.

## Workflow

### Execution flow

```text

config/config.yaml
        |
        v
run_analysis.py
(section-based CLI)
        |
        v
Snakemake DAG
(workflow logic)
        |
        v
Docker container
(reproducible execution environment)

```

### Upstream (engineering)

```text

FASTQs
  |
  +-- download (optional)
  |
  +-- validate presence
  |
  v
QC (raw)
  - FastQC
  - MultiQC
  |
  +-- trim (optional)
  |     |
  |     v
  |   QC (trimmed)
  |     - FastQC
  |     - MultiQC
  |
  v
Reference
  - genome FASTA
  - GTF
  - barcode whitelist
  - STAR index
  |
  v
Alignment / Counting 
  - STARsolo
  - gene x cell matrix
  - alignment on raw OR trimmed reads depending on mode

```

### Downstream (analysis)

```text

Per-donor
  |
  +-- Seurat object
  +-- cell QC + filtering
  +-- normalization + HVGs
  +-- clustering + annotation
  |
  v
Cross-donor
  |
  +-- pseudobulk by donor & cell type
  +-- DESeq2 (DE)
  +-- TOST (equivalence test)
  +-- enrichment (GSEA / ORA)
  |
  v
(Work in progress)
  - co-expression network analysis

```

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
  trimmed/            # Trimmed FASTQ files (generated only if read trimming is enabled)
results/             # Outputs and logs (not versioned)
  qc/                 # FastQC / MultiQC reports
  alignment/           # STARsolo outputs
  logs/                # Execution logs
  downstream/         # Downstream analysis results
	/deg_and_tost     # DEG and TOST analysis results
	/seurat           # Seurat objects and related plots and tables
	/networks         # Network analysis results (Work in progress)
docs/                # Documentation (user manual, report, notes)
scripts/             # R-scripts and helpers
run_analysis.py      # Optional Python wrapper for section-based execution
```

## Representative Results Directory Layout (Example Run)

```text
results/
├── qc/
│   ├── fastqc/
│   │   └── raw/
│   │       ├── donor1/
│   │       ├── donor2/
│   │       ├── donor3/
│   │       └── donor4/
│   └── multiqc/
│       └── raw/
│           ├── multiqc_report.html
│           └── multiqc_report_data/
├── alignment/
│   └── starsolo/
│       └── raw/
│           ├── donor1/
│           ├── donor2/
│           ├── donor3/
│           └── donor4/
├── downstream/
│   ├── seurat/
│   │   └── untrimmed/
│   │       ├── donor1/
│   │       ├── donor2/
│   │       ├── donor3/
│   │       └── donor4/
│   ├── deg_and_tost/
│   │   └── untrimmed/
│   │       └── deg_and_tost/
│   └── networks/
│       └── untrimmed/
│           ├── consensus/
│           ├── per_donor/
│           ├── plots/
│           └── tables/
└── logs/
    ├── qc/
    ├── alignment/
    ├── download/
    └── ref/


```


## Notes
- Large data files are excluded via `.gitignore`
- FASTQ downloading is controlled via `io.download_fastqs` in `config/config.yaml`  
  (automatically set by the wrapper for relevant sections)
- STAR index requires ~25–30 GB RAM
- Biological analyses are included to validate 
  pipeline correctness and demonstrate statistically coherent downstream usage, 
  not to claim novel biological findings.
- R package versions inside the Docker image are managed with `renv`
  to ensure reproducible R environments.
  Users do not need to interact with `renv` directly.
- A representative execution artifact (multiqc_report.html) has been copied to 
 docs/example_outputs/ as a lightweight proof of pipeline execution (open locally in a browser after cloning);
 full outputs remain under results/ and are not version-controlled.
- Developer / maintainer notes (including environment maintenance,
STAR indexing details, and Docker CPU behavior) are kept in
`docs/DEVELOPER_NOTES.md` and are not required to run the workflow.

  
## Non-goals

- This pipeline is not intended to benchmark methods or claim novel biological findings.

## Data availability (planned)

Once the workflow is finalized, a stable snapshot of representative results will be archived on Zenodo with a DOI.

This archive will include:
- Representative output data produced by the full pipeline (alignment and downstream analysis).
- A small toy dataset (downsampled FASTQs) intended for demonstration and pipeline sanity-check runs only.

At the time of final release, the Docker image used to generate the archived results will also be published and referenced explicitly, 
ensuring reproducible execution of the archived snapshot.


## Citation

If you use this workflow in a publication, please cite the repository.
A Zenodo DOI for a stable release will be provided upon publication.
