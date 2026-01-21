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

