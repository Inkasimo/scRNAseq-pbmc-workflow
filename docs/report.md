# Reproducible scRNA-seq Analysis Pipeline (PBMC)

## 1. Overview
This project implements a reproducible, containerized scRNA-seq analysis
pipeline using Snakemake and Docker. The pipeline is demonstrated on
10x Genomics PBMC datasets (donors 1–4).

The primary purpose of this repository is to demonstrate the design,
implementation, and execution of a production-style bioinformatics
workflow, including reproducibility, dependency management, and
end-to-end data processing.

This project is intended as a technical portfolio example. The analysis
is not designed to generate novel biological findings or support
biological claims beyond demonstrating a complete and functional
pipeline.

\textit{Also add my computer specs on the report}

## 2. Execution environment

All computations were performed on a local workstation
running Windows 11 with WSL2 and Docker Desktop.

Due to the memory requirements of STAR genome indexing,
the default WSL2 configuration was modified to allow
up to 32 GB of RAM.

Hardware:
- CPU: …
- RAM: …
- Storage: …


## 3. Repository Organization

The repository is structured to separate execution logic, configuration,
containerization, and generated data, following common best practices for
reproducible computational workflows.

- `containers/`  
  Docker build context defining a fully self-contained execution environment,
  including all bioinformatics tools and their exact versions.

- `workflow/`  
  Snakemake workflow definition, including rules for data acquisition,
  quality control, reference preparation, and alignment.

- `config/`  
  User-editable configuration files controlling dataset paths, resource usage,
  and execution toggles (e.g. downloading FASTQs vs. validating local files).

- `resources/`  
  Static, versioned resources required by the workflow.  
  This includes the 10x Genomics barcode whitelist, which is bundled directly
  to avoid reliance on unstable external download URLs.

- `data/` *(not versioned)*  
  Large input data and reference files generated or downloaded at runtime.
  This includes raw FASTQs, reference genomes, annotations, and STAR indices.

- `results/` *(not versioned)*  
  All pipeline outputs, including QC reports, alignment results, and log files.
  These are excluded from version control by design.

- `docs/`  
  Supplementary documentation, including the user manual and this report.

- `scripts/`  
  Reserved for future downstream analysis scripts (not yet implemented).

- `run_analysis.py`  
  A lightweight Python wrapper that orchestrates Snakemake execution inside
  Docker, enabling section-based execution and preventing accidental full runs.
  
## R Environment and Reproducibility Strategy

### Overview

All R-based analyses in this workflow are executed inside a Docker container with a fully specified and reproducible R environment. 
Dependency management is handled using renv, while containerization ensures isolation from host-specific configurations. 
This design guarantees that analyses can be reproduced deterministically across machines and over time.

### R Version and Package Management

The R environment is defined using renv, which creates a project-local R library and records exact package versions in a lockfile 
`renv.lock`. The lockfile captures:

- The exact R version used (R 4.3.3)
- The CRAN snapshot date (https://packagemanager.posit.co/cran/2024-01-01)
- All R package versions required for the analysis

During Docker image build, the environment is restored only from this lockfile using renv::restore(). 
No package installation or snapshotting occurs at runtime, preventing accidental environment drift.

### Vendoring of Seurat

The core analysis depends on Seurat and SeuratObject, which are critical, 
fast-evolving packages that have historically introduced reproducibility issues due to 
upstream changes and repository availability.

To eliminate this risk, specific versions are vendored directly into the repository:

- `renv/vendor/seurat-v4.4.0.tar.gz`
- `renv/vendor/seurat-object-v4.1.4.tar.gz`

These tarballs are installed locally and recorded in the lockfile with Source = "Local". As a result:

Docker builds never rely on GitHub availability or branch state

Seurat versions are fully decoupled from CRAN/Bioconductor changes

The exact Seurat code used in the analysis is preserved with the project

## Docker Integration

The Docker image is built using a rocker-based R image, which provides a fixed R runtime 
while retaining a Debian/Ubuntu-compatible system environment.
System-level dependencies required for R packages and external tools 
(e.g. STAR, FastQC) are installed explicitly during the image build.

- The CRAN snapshot is pinned via an environment variable.
- The renv.lock file and renv/ directory are copied into the image.
- A bootstrap installation of renv is performed.
- renv::restore() installs all R packages into the container’s project library.

Importantly, the Docker build does not modify the lockfile. 
All dependency changes must be made explicitly on the host and committed to version control.

### Separation of Concerns

- Host system: used only to update dependencies and regenerate the lockfile when necessary.
- Docker image: consumes the lockfile in read-only fashion and provides a stable execution environment.
- IDE tooling: automatic IDE-driven package installation is disabled to prevent contamination of the project library.

#### Result

This setup ensures:

- Bitwise reproducibility of the R environment
- Long-term stability of critical dependencies (notably Seurat)
- Clear provenance of all software components
- Safe integration with Snakemake and container-based execution

In summary, the combination of renv, vendored critical packages, 
and Docker provides a robust and auditable foundation for reproducible single-cell RNA-seq analysis.

## 4. Execution Model
The pipeline is orchestrated using Snakemake and executed inside Docker
containers to ensure reproducibility and avoid local environment drift.
All analysis steps are executed via containerized tools.

## 5. Data
Raw FASTQ files are downloaded from 10x Genomics public PBMC datasets
(donors 1–4, 3' Gene Expression chemistry).

## 6. Quality Control (planned)
Raw FASTQ quality will be assessed using FastQC and summarized with MultiQC.
Based on these results, a decision will be made whether explicit read trimming
is required.

FastQC/MultiQC did not indicate meaningful adapter contamination in the cDNA read.
R2 “Adapter Content” is reported as PASS across the lanes/samples 
in the aggregated MultiQC status table in MultiQC report

The main flagged item is Overrepresented sequences in Donor4 R2 files,
where the top sequences match a SMARTer/SO-like oligonucleotide and each sits at ~0.10%
of reads (FastQC’s reporting threshold), consistent with expected
library-structure/short-insert artifacts rather than pervasive adapter read-through. 

For 10x 3′ scRNA-seq, trimming is typically unnecessary and can be counterproductive: 
standard pipelines (e.g., Cell Ranger / STARsolo) are designed to tolerate modest 3′ 
artifacts via soft-clipping while preserving the read structure for accurate alignment 
and UMI counting; therefore trimming was omitted unless strong adapter signal is observed.

## 7. Alignment and Quantification (planned)

### Choice of aligner: STARsolo vs Cell Ranger
<1–2 paragraphs explaining transparency, control, reproducibility>

### Reference genome and annotation
- GRCh38 primary assembly
- GENCODE v45
- STAR 2.7.11b
- Index built locally (not distributed)

### Read structure and chemistry assumptions
- 10x GEM-X 3′ v4
- R1 = CB+UMI, R2 = cDNA
- Explicit CB/UMI coordinates

### Counting and outputs
- STARsolo gene-level counting
- Raw and filtered gene–cell matrices
- No external counting (e.g. featureCounts)

### BAM output
- Disabled by default to reduce memory and disk usage
- STARsolo matrices are complete without BAMs
- Instructions provided to enable BAMs if needed



Reads will be aligned and quantified using STARsolo with UMI-aware counting.

Running on Windows with WSL2, I had to increase WSL memory.
A template is provided at docs/wslconfig.example.  
The default WSL settings caused an OOM error. Rename the file to 
".wslconfig" and copy in folder C:\Users\"Username"

[wsl2]
memory=32GB
processors=8
swap=16GB

\textit{Also add my computer specs on the report}

[RESULTS TO BE ADDED]

## Trimming (cutadapt) Not default and not applied on this data but a possibility


What trimming does in this workflow

Only Read 2 (R2 / cDNA) is trimmed

Read 1 (R1 / cell barcode + UMI) is left unchanged

Trimming is performed with cutadapt

If trimming is enabled, all downstream steps use trimmed FASTQs

Trimmed FastQC / MultiQC

STARsolo alignment into a separate trimmed/ results directory

Trimming parameters (from config/config.yaml)
trim:
  enabled: false        # default; overridden by wrapper flags
  adapter_r2: AGATCGGAAGAG
  q_r2: 20
  minlen_r2: 20


Meaning:

adapter_r2: adapter sequence removed from R2 only

q_r2: quality trimming applied to R2 only

minlen_r2: minimum length enforced on R2
(read pairs are dropped if R2 fails)

Output locations

Trimmed FASTQs:
data/trimmed/{donor}/*.fastq.gz

Trimming completion marker:
data/trimmed/{donor}/trim.done

Original filenames and lane structure are preserved.

## 8. Downstream Analysis (planned)
- Cell filtering and normalization
- Feature selection
- DEG analysis between T-cells and B-cells
- Cell-type–specific co-expression network construction
- Module detection
- Functional enrichment analysis

[RESULTS TO BE ADDED]

#Notes

Maybe a t-cell network made per donor...

then take the edges that appear > 3 of networks..


Think about it

## 9. Reproducibility
The full pipeline, including tool versions and execution steps, is available
in the accompanying GitHub repository.

# Downstream technical notes

1) How the R environment is created in Docker and “picked up” by the workflow

The container builds on rocker/r-ver:4.3.3 and uses renv to make the R package set reproducible. In the Dockerfile, you set RENV_CONFIG_REPOS_OVERRIDE (Posit Package Manager snapshot) plus RENV_PATHS_LIBRARY=/opt/renv/library and RENV_PATHS_CACHE=/opt/renv/cache, then install renv globally and run renv::restore() during the image build. That means the final image already contains a fully restored R library under /opt/renv/library, so runtime does not need to download/install packages.

Snakemake “reads” this environment implicitly because your workflow is executed inside that same Docker image (via docker run ... scrnaseq-workflow snakemake ...). When a rule calls Rscript ..., it uses the container’s R binary and the restored renv library paths. You also copy renv.lock, renv/activate.R, and .Rprofile into /work during build; with WORKDIR /work, Rscript runs in the project context where renv auto-activation is available (and the library path is already pinned via RENV_PATHS_LIBRARY). Net effect: Rscript calls in the Snakefile resolve packages consistently across hosts.

2) Downstream logic (Seurat object + QC)

Starsolo folders raw, filtered

vs upstream raw, trimmed

The downstream step is encapsulated in rule seurat_qc. It consumes STARsolo “filtered” count outputs (matrix.mtx*, barcodes.tsv*, and features.tsv*/genes.tsv*) and runs scripts/build_seurat_objects_qc.R, writing both an .rds Seurat object and a sentinel seurat_qc.done file. The sentinel is what Snakemake uses as the “this stage completed” marker; it prevents re-running unless inputs/metadata say otherwise.

Mode routing is handled by the trim_enabled() / is_trimmed() logic plus the LEGACY_BRANCH mapping. Conceptually:

If trim_enabled=false ⇒ “untrimmed” downstream uses alignment outputs in results/alignment/starsolo/raw/...

If trim_enabled=true ⇒ “trimmed” downstream uses results/alignment/starsolo/trimmed/...

That mapping is why downstream paths use trim_state (untrimmed vs trimmed) while STARsolo folders remain raw/trimmed. The top-level rule all also respects this: it only expands the downstream targets for the active mode (raw or trimmed), so you don’t accidentally build both branches unless you explicitly request both via targets/config.

Finally, your wrapper (run_analysis.py) implements “section → target set” selection. The downstream and build_seurat_object_qc sections simply emit the relevant results/downstream/seurat/{trim_state}/{donor}/seurat_qc.done targets (for all donors by default), and Snakemake computes the minimal DAG needed to produce them from whatever already exists upstream.

QC filtering and normalization rationale
Overview

This workflow performs donor-specific quality control and normalization of single-cell RNA-seq data prior to downstream analysis. The design prioritizes robustness, reproducibility, and interpretability, rather than maximizing sensitivity or applying complex statistical modeling. All decisions are motivated by the structure of the dataset (PBMCs from four donors) and the intended downstream analyses (cell-type–aware DEG, pseudobulk aggregation, and coexpression networks).

QC filtering strategy
Data-adaptive thresholds

Rather than applying fixed global cutoffs, QC thresholds are derived from donor-specific empirical distributions (qc_metrics.tsv). This avoids over-filtering donors with systematically higher or lower sequencing depth and ensures that filtering decisions are data-driven and reproducible.

The following rules are applied per donor:

Mitochondrial content
percent.mt ≤ min(mt_q90, 25)
This caps extreme mitochondrial outliers while preventing overly aggressive filtering in donors with modestly elevated MT fractions.

Gene complexity (nFeature_RNA)

Lower bound: nFeature_q05
Removes empty droplets and very low-quality cells.

Upper bound: min(nFeature_q99, 6000)
Acts as a conservative doublet/multiplet guardrail.

UMI counts (nCount_RNA)

Lower bound: nCount_q05

Upper bound: nCount_q99
These bounds are kept consistent with the gene complexity filters.

Hemoglobin content (PBMC-specific)
percent.hb ≤ 1
This removes residual erythrocyte contamination without affecting true PBMC populations.

Transparency and auditability

For each donor, the pipeline records:

exact QC thresholds used (qc_thresholds.tsv)

number and percentage of cells removed

per-rule and multi-rule failure counts (qc_filter_summary.tsv)

This makes QC decisions explicit, inspectable, and easy to compare across donors.

Normalization strategy
Choice of LogNormalize

After QC filtering, expression values are normalized using Seurat’s LogNormalize method (library-size scaling followed by log-transformation).

This choice is intentional:

Log-normalized expression preserves relative expression relationships and covariance structure, which is important for:

pseudobulk aggregation

coexpression network analysis

It avoids introducing model-based residuals that can complicate interpretation in small-n, donor-aware analyses.

It is simple, transparent, and widely understood, making the pipeline easier to reason about and explain.

More complex alternatives (e.g. SCTransform) were considered but deliberately not used in order to:

avoid implicit information sharing across donors,

preserve biologically interpretable expression values,

keep downstream aggregation and correlation analyses well-defined.

Variable feature selection

Highly variable genes (HVGs) are identified using the VST method (n = 2000).
This step is used to:

focus dimensionality reduction on informative genes,

reduce noise from low-variance features,

keep computation predictable.

HVG selection does not alter the normalized expression matrix used for downstream biological analyses.

Scaling and regression (visualization only)

Expression values are scaled and mitochondrial percentage is regressed out only for highly variable genes. This step is included exclusively to support dimensionality reduction and visualization (PCA / UMAP).

Key points:

Scaling is restricted to HVGs to avoid unnecessary transformation of low-information genes.

Regression is limited to percent.mt to reduce technical structure in embeddings.

Scaled data are not used for:

differential expression

pseudobulk analysis

coexpression networks

This separation ensures that visualization does not distort biologically meaningful signal used in inference.

Outputs and reproducibility

The final QC-filtered and normalized Seurat object contains:

raw counts

log-normalized expression (RNA@data)

HVGs

scaled data for visualization (RNA@scale.data)

All normalization parameters and summary statistics are written to disk (norm_metrics.tsv), and diagnostic plots are generated to confirm that normalization behaves as expected.

Summary

This QC and normalization strategy is intentionally conservative and explicit. It favors:

donor-aware filtering,

data-adaptive thresholds,

interpretable normalization,

clear separation between visualization and inference.

The result is a stable foundation for downstream analyses without overfitting assumptions or obscuring biological signal.

