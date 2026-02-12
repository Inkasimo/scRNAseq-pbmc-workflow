# scRNA-seq PBMC Workflow

Production-style scRNA-seq analysis pipeline
using **Docker + Snakemake** with a Python CLI wrapper.

This repository is intended as a **technical portfolio / learning project**.

## Status:
- Draft — core pipeline implemented; downstream analyses in progress.
- A small set of representative execution artifacts 
 (MultiQC report and selected downstream plots) is included under `docs/example_outputs/`
 as lightweight proof of successful execution. Full outputs are written to `results/`
 and are intentionally not version-controlled.

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

## Docker image

### Versioned release (recommended)

Pull the versioned release

```bash
docker pull ghcr.io/inkasimo/scrnaseq-pbmc-workflow:v1.0.1
```

### Exact archive (bit-identical to published results)

Pull the archived Docker image used for the released workflow:

```bash
docker pull ghcr.io/inkasimo/scrnaseq-pbmc-workflow@sha256:80354b76e76405636c43e73902236e0399d26978a214227afbafa46fc0555bb8
```

This ensures you are running exactly the same environment used to generate the archived results.

### Local build (development only)

To build the image locally instead:

```bash
docker build -t scrnaseq-workflow -f containers/Dockerfile .
```

Local builds are not guaranteed to be bit-identical to the archived environment.

### Note

- If you are using exact archive or local build as an docker image instead of versioned release, 
with the wrapper you need to define the docker image with `--image` -flag

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
docs/                 # Documentation (user manual, report, notes)
  example_outputs/    # A small set of representative execution artifacts
scripts/              # R-scripts and helpers
run_analysis.py       # Optional Python wrapper for section-based execution
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
- Developer / maintainer notes (including environment maintenance,
  STAR indexing details, and Docker CPU behavior) are kept in
  `docs/DEVELOPER_NOTES.md` and are not required to run the workflow.
  
### Windows / WSL / Docker Desktop

This workflow runs inside Docker and bind-mounts the repository into the container.

On Windows, this requires Docker Desktop to have access to the repository path.
If Docker cannot bind-mount the path, the container will not be able to read
`workflow/Snakefile`.

If you encounter:

    Snakefile "workflow/Snakefile" not found

check Docker Desktop file-sharing settings, or run the workflow on a Linux or
macOS system where bind mounts are always available.


## Non-goals

- This pipeline is not intended to benchmark methods or claim novel biological findings.

## Data availability (planned)

Once the workflow is finalized, a stable snapshot of representative results will be archived on Zenodo with a DOI.

This archive will include:
- Representative output data produced by the full pipeline (alignment and downstream analysis).
- A small toy dataset (downsampled FASTQs) intended for demonstration and pipeline sanity-check runs only.


## Citation

If you use this workflow in a publication, please cite the repository.
A Zenodo DOI for a stable release will be provided upon publication.
