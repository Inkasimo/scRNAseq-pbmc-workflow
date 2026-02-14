[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18642101.svg)](https://doi.org/10.5281/zenodo.18642101)


# scRNA-seq PBMC Workflow

Reproducible, containerized single-cell RNA-seq workflow built with Snakemake + Docker, controlled via a Python CLI wrapper.

End-to-end execution:

FASTQ → QC → STARsolo → Seurat → DESeq2/TOST → enrichment → network analysis


## Quick Start (Toy Demonstration)

- Runs a chromosome 1 mini-reference with downsampled FASTQs.
- Execution time: ~5–10 minutes after image pull and toy data download.
- Download size: ~ 81.3 MB (toy bundle) 
- Toy runs upstream+Seurat object creation (No downstream from Seurat object creation)

### 1. Clone the repository

```bash
git clone https://github.com/inkasimo/scRNAseq-pbmc-workflow.git
```

### 2. Pull the versioned Docker image (in repository)

```bash
docker pull ghcr.io/inkasimo/scrnaseq-pbmc-workflow:v1.0.3
```

First pull may take several minutes (image ~1.5 GB).

### 3. Install wrapper dependency (host only)

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r wrapper-requirements.txt
```

Required only if using run_analysis.py.
Not needed if running Snakemake directly via Docker.

Activate the environment before using run_analysis.py.

### 4. Download toy bundle

```bash
python3 run_analysis.py download_toy
```
This downloads and extracts:

- `data/ref/toy/` (chr1 reference files)
- `data/toy/donor1/` (toy FASTQs)


### 5. Run toy workflow

**Dry run:**
```bash
python3 run_analysis.py toy --dry-run
```

**Raw mode:**
```bash
python3 run_analysis.py toy
```

**Trimmed mode:**

```bash
python3 run_analysis.py toy --trimmed
```

Outputs written to:

`results/`


Representative example outputs from real run are available under:

`docs/example_outputs/`


## What this pipeline does

### Upstream

- FASTQ acquisition (download or reuse of existing data)
- Raw and optional trimmed QC (FastQC + MultiQC)
- Optional read trimming (Cutadapt; non-default)
- Reference preparation (barcode whitelist, STAR index)
- STARsolo alignment and gene–cell count matrix generation

### Downstream 

- Seurat objects
- Cell-level QC and annotation
- Differential expression (DESeq2) and equivalence testing (TOST)
- Enrichment analysis
- Co-expression network analysis and network modules
- Module enrichment analysis

## Focus

- Reproducible execution
- Explicit workflow structure
- Containerized dependencies
- Clear separation between engineering (upstream) and analysis (downstream)


## Full Dataset Execution

### Dataset
- 10x Genomics PBMC 5k
- Donors 1–4
- 3′ Gene Expression

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


### 1. Pull docker image 

Use the same Docker image as shown in Quick Start.

### Execution with python wrapper (recommended)

#### Wrapper requirements (host-side only)

See Quick Start

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

### Execution directly with Snakemake (no wrapper)

This workflow can be run in **Directly with Snakemake inside Docker** 

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

## Repository structure

```text
containers/           # Dockerfile(s) for reproducible execution
workflow/             # Snakemake workflow (rules, DAG)
config/               # User-editable configuration (config.yaml)
resources/            # Static resources bundled with the workflow
  barcodes/           # 10x barcode whitelist(s)
  genesets/           # Hallmark and C7 .gmt files
data/                 # Input data and references (not versioned)
  raw/                # FASTQ files (downloaded or user-provided)
  ref/                # Reference genome, GTF, STAR index
  trimmed/            # Trimmed FASTQ files (generated only if read trimming is enabled)
results/              # Outputs and logs (not versioned)
  qc/                 # FastQC / MultiQC reports
  alignment/          # STARsolo outputs
  logs/               # Execution logs
  downstream/         # Downstream analysis results
	deg_and_tost/     # DEG and TOST analysis results
	seurat/           # Seurat objects and related plots and tables
	networks/         # Network analysis results (Work in progress)
docs/                 # Documentation (user manual, report, notes, results layout)
  example_outputs/    # A small set of representative execution artifacts
scripts/              # R-scripts and helpers
run_analysis.py       # Optional Python wrapper for section-based execution
```


 
## Workflow

### Execution flow

```text

config/config.yaml
        |
        v
run_analysis.py  (host-side CLI wrapper)
        |
        v
docker run -v $PWD:/work -w /work ...  (bind-mount repo)
        |
        v
+------------------------------------------------------+
|                  Docker Container                    |
|                                                      |
|  Snakemake (workflow/Snakefile) → rules → tools      |
|                                                      |
+------------------------------------------------------+


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
  +-- co-expression network analysis

```

## Resources and reproducibility

### Barcode whitelist

The 10x Genomics barcode whitelist is bundled directly in the repository:

`resources/barcodes/3M-3pgex-may-2023_TRU.txt`


This avoids reliance on unstable upstream URLs and ensures reproducible execution.

### Gene sets (MSigDB)

Gene sets used for enrichment analyses are stored locally under:

`resources/genesets/`


This includes:

MSigDB Hallmark gene sets:

`h.all.v2026.1.Hs.symbols.gmt`


MSigDB C7 (immunologic signatures):

`c7.all.v2026.1.Hs.symbols.gmt`


Hallmark (H) and Immunologic Signature (C7) gene sets were obtained from the Molecular Signatures Database 
(MSigDB, Broad Institute) and are included locally to ensure reproducible execution of the workflow. 
All enrichment steps (GSEA and ORA) explicitly reference these local files

Users should ensure compliance with MSigDB licensing terms when reusing these resources.

```text
Liberzon A et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection.
Cell Systems (2015).

Godec J et al. Compendium of Immune Signatures Identifies Conserved and Species-Specific
Biology in Response to Inflammation. Immunity (2016).
```

## Documentation

Project documentation is organized under the `docs/` directory:

- **Workflow summary** (`docs/workflow_summary.pdf`) — concise 6 page technical overview highlighting architecture, reproducibility strategy, statistical modeling, and representative results.  
  *(Recommended starting point.)*
- **Technical report (full)** (`docs/scRNAseq-pbmc-worflow-report_full.pdf`) — comprehensive documentation of architectural decisions, biological rationale, statistical methodology, and extended results.
- **User manual** (`docs/user_manual.md`) — execution instructions, Docker image usage, and wrapper configuration details.
- **Results layout** (`docs/results_layout.md`) — representative output directory structure for result discovery and navigation.

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
  
## Non-goals

- This pipeline is not intended to benchmark methods or claim novel biological findings.

## Data availability (planned)

A stable snapshot of representative results are archived on Zenodo with DOI. 

This archive includes:
- Representative output data produced by the full pipeline (alignment and downstream analysis) (Work in progress).
- A small toy dataset (downsampled FASTQs) intended for demonstration and pipeline sanity-check runs only.


## Citation

If you use this workflow in a publication, please cite the repository.