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

## 2. Repository Structure
- containers/: Docker images for isolated execution
- workflow/: Snakemake rules
- scripts/: R/Python analysis scripts
- config/: dataset and parameter configuration
- data/: raw and intermediate data (not versioned)
- results/: analysis outputs and logs (not versioned)
- docs/: documentation and reports

## 3. Execution Model
The pipeline is orchestrated using Snakemake and executed inside Docker
containers to ensure reproducibility and avoid local environment drift.
All analysis steps are executed via containerized tools.

## 4. Data
Raw FASTQ files are downloaded from 10x Genomics public PBMC datasets
(donors 1–4, 3' Gene Expression chemistry).

## 5. Quality Control (planned)
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

## 6. Alignment and Quantification (planned)
Reads will be aligned and quantified using STARsolo with UMI-aware counting.

The default docker settings caused an OOM error so needed to 
write a .wslconfig file in C:\Users\"Username"

[wsl2]
memory=32GB
processors=8
swap=16GB

\textit{Also add my computer specs on the report}

[RESULTS TO BE ADDED]

## 7. Downstream Analysis (planned)
- Cell filtering and normalization
- Feature selection
- DEG analysis between T-cells and B-cells
- Cell-type–specific co-expression network construction
- Module detection
- Functional enrichment analysis

[RESULTS TO BE ADDED]

## 8. Reproducibility
The full pipeline, including tool versions and execution steps, is available
in the accompanying GitHub repository.
