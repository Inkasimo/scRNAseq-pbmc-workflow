# scRNA-seq PBMC Workflow

**Work in progress**
Production-style scRNA-seq analysis pipeline
using **Docker + Snakemake** with a Python CLI wrapper.

This repository is intended as a **technical portfolio / learning project**.

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

- FASTQ validation or download
- Raw and optional trimmed QC (FastQC + MultiQC)
- Optional read trimming (Cutadapt; non-default)
- Reference preparation (barcode whitelist, STAR index)
- STARsolo alignment and gene–cell count matrix generation

### Downstream (work in progress)

- Seurat objects
- Cell-level QC and annotation
- Differential expression and TOST analysis
- Enrichment analysis
- Co-expression network analysis (work in progress)

## Requirements

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

## Option A: Python wrapper (recommended)


### Wrapper requirements (host-side only)

```bash
pip install -r wrapper-requirements.txt
```
This installs:
- pyyaml (used only by the wrapper)

If you do not use the wrapper, you can ignore this step entirely.


### Inspect available sections:

```bash
python3 run_analysis.py --list-sections
```

### Inspect donors:

```bash
python3 run_analysis.py --list-donors
```
### Example runs

#### Download data:

```bash
python3 run_analysis.py download_data --cpus 8 --cores 8
```

#### QC only:

```bash
python3 run_analysis.py qc --cpus 8 --cores 8
```

#### Align all donors:

```bash
python3 run_analysis.py align \
  --donor all \
  --cpus 8 --cores 8 \
  -j 1 \
  --set-threads starsolo=8

```

#### Dry run (no execution):

```bash
python3 run_analysis.py all --dry-run
```

### Option B: Direct Snakemake (no wrapper)

#### Dry run

```bash
docker run --rm -it \
  -v "$(pwd)":/work \
  -w /work \
  scrnaseq-workflow \
  snakemake -n -p
  ```


#### Run a specific target (example: one donor alignment):

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

scRNA-seq PBMC Workflow (Docker + Snakemake + CLI wrapper)
=========================================================

                    +----------------------+
                    |   config/config.yaml |
                    +----------+-----------+
                               |
                               v
                    +----------------------+
                    |  run_analysis.py     |
                    |  (section runner)    |
                    +----------+-----------+
                               |
                               v
                    +----------------------+
                    |  Snakemake (DAG)     |
                    |  workflow/Snakefile  |
                    +----------+-----------+
                               |
   -------------------------------------------------------------------------
   UPSTREAM (engineering)                                               
   -------------------------------------------------------------------------

  [FASTQs]
    |-- download OR validate presence ------------------------------+
    |                                                               |
    v                                                               v
+------------------+                                       +------------------+
| FastQC (raw)     |                                       | Cutadapt (opt)   |
+--------+---------+                                       +--------+---------+
         |                                                          |
         v                                                          v
+------------------+                                       +------------------+
| MultiQC (raw)    |                                       | FastQC (trimmed) |
+--------+---------+                                       +--------+---------+
                                                           | 
                                                           v
                                                  +------------------+
                                                  | MultiQC (trimmed)|
                                                  +--------+---------+

  [REFERENCE]
    |
    |-- download genome+gtf (if needed)
    |-- build STAR index (or validate existing)
    |-- validate whitelist
    v
+------------------------------+
| Reference ready (STAR index) |
+---------------+--------------+
                |
                v

  [ALIGNMENT / COUNTING]
                |
                v
+---------------------------------------------+
| STARsolo (raw OR trimmed mode)              |
| -> gene x cell count matrix (MTX+features+barcodes)
+----------------------+----------------------+
                       |
                       v

   -------------------------------------------------------------------------
   DOWNSTREAM (analysis)  (per donor, then cross-donor)
   -------------------------------------------------------------------------

  Per-donor:
    v
+------------------------------+
| Build Seurat object + QC     |
| (CreateSeuratObject, QC plots|
|  qc_metrics.tsv)             |
+---------------+--------------+
                |
                v
+------------------------------+
| QC filter + normalization    |
| NormalizeData + HVGs + Scale |
+---------------+--------------+
                |
                v
+------------------------------+
| Clustering + annotation      |
| (markers / module scores)    |
+---------------+--------------+
                |
                v

  Cross-donor (all donors together):
                v
+----------------------------------------------+
| Pseudobulk by donor x cell-type-like group   |
| -> DESeq2 (donor blocking)                   |
| -> Markers (strict)                          |
| -> TOST equivalence (conserved)              |
+----------------------+-----------------------+
                       |
                       v
+----------------------------------------------+
| Enrichment                                    |
| - GSEA (fgsea) on ranked log2FC               |
| - ORA (clusterProfiler enricher) on gene sets |
+----------------------------------------------+

  (Next)
+----------------------------------------------+
| Co-expression network analysis (WIP)          |
+----------------------------------------------+


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
- FASTQ downloading is controlled via `io.download_fastqs` in `config/config.yaml`  
  (automatically set by the wrapper for relevant sections)
- STAR index requires ~25–30 GB RAM
- Biological analyses are included to validate 
  pipeline correctness and demonstrate statistically coherent downstream usage, 
  not to claim novel biological findings.


## Citation

If you use this workflow in a publication, please cite the repository.
A Zenodo DOI for a stable release will be provided upon publication.
