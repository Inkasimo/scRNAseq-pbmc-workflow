# User manual


## 1. Purpose and Scope

### Purpose

This pipeline implements a reproducible, containerized workflow for processing
10x Genomics 3′ scRNA-seq data, demonstrated on PBMC datasets.

The pipeline is designed to cover both upstream processing and downstream
biological analysis.  

### Upstream Processing 

The upstream component of the pipeline handles data ingestion, quality control,
reference preparation, and alignment/quantification, including:

- FASTQ validation or download
- Quality control (FastQC and MultiQC)
- Optional read trimming (non-default)
- Reference preparation (barcode whitelist and STAR index)
- STARsolo alignment and generation of gene–cell count matrices

### Downstream Analysis

- Cell-level filtering, normalization, and clustering
- Cell type annotation
- Differential expression analysis
- Equivalence test (TOST)
- Pathway or gene set enrichment analysis
- Co-expression and modularity analysis (Work in progress)
- Result summarization and visualization

### Project Goal

The primary goal of this project is to demonstrate disciplined workflow engineering:
clear execution structure, reproducibility via containerization, and restartable
execution using Snakemake.

### Intended Audience

This project is intended for technical users working in computational biology or
bioinformatics, including:

- Bioinformaticians and computational biologists
- Junior practitioners or trainees learning to build reproducible pipelines
- Wet-lab scientists seeking to develop practical computational workflow skills
- Software-oriented scientists interested in scientific workflow engineering

Familiarity with the command line and basic container concepts (Docker) is assumed,
but the project is designed to be readable and instructive rather than opaque.

The repository is structured as a technical portfolio and learning project, not
as a production or enterprise pipeline. Design decisions prioritize clarity,
reproducibility, and explicit control over maximal automation or fault tolerance.

## 2. Execution Model (How This Pipeline Runs)

### Containerized Execution

All workflow steps are executed inside a Docker container. The container provides
a fully defined software environment, including Snakemake and all bioinformatics
tools required by the pipeline.

**The host system is responsible only for:**
- Providing Docker
- Managing available CPU and memory resources

Input data are downloaded automatically by the workflow unless explicitly disabled.
A default configuration file is provided and does not require modification for
standard execution.

### Workflow Orchestration with Snakemake

Snakemake is used to define and orchestrate pipeline steps as a directed acyclic
graph (DAG) of rules. Each rule produces one or more output files and may depend
on outputs from upstream rules.

The workflow is designed for explicit, target-driven execution rather than
implicit end-to-end runs.

### Python Wrapper (run_analysis.py)

The repository includes an optional Python wrapper that invokes Snakemake inside
Docker and provides a simplified, section-based command-line interface.

**The wrapper:**

- Requires users to explicitly select which pipeline section to run
- Prevents accidental execution of the full workflow
- Sets configuration flags automatically based on the selected section
- Exposes common Snakemake options (cores, jobs, thread settings)

The wrapper does not replace Snakemake; it constrains and simplifies how Snakemake
is invoked.

The Python wrapper requires a local Python installation (≥ 3.9) and a small set
of host-side dependencies defined in `wrapper-requirements.txt`. These dependencies
are used only by the wrapper and are not required when running Snakemake directly
inside Docker.

### Success Indicators (.done Files)

Each major pipeline step creates a .done file only after successful completion.

.done files are used as the authoritative success indicators for:

- Data download or validation
- Quality control
- Trimming
- Reference preparation
- Alignment

Partial outputs may exist if a step is interrupted, but the absence of a .done
file indicates that the step did not complete successfully.

### Interruptions and Reruns

Workflow interruptions (e.g. Ctrl+C, system sleep, container termination) are
expected and supported.

If a run is interrupted:

- Partial outputs may remain on disk
- Corresponding .done files will not be created

Rerunning the same command with --rerun-incomplete enabled (the default behavior
when using the wrapper) will resume execution and rerun any incomplete steps.

Snakemake lock files may remain after abnormal termination and can be cleared
using the standard Snakemake --unlock mechanism.


## 3. System Requirements

### Docker

Docker is required to run the pipeline. All workflow steps are executed inside a
Docker container that provides the complete software environment.

No local installation of bioinformatics tools is required outside Docker.

### Disk Space

Significant disk space is required due to reference files, intermediate outputs,
and alignment results.

As a guideline:

- ~170 GB of free disk space is recommended for a full run on the PBMC dataset
(including reference files, FASTQs, and results).

Disk usage may increase if optional steps such as BAM output or trimming are enabled.

### Memory (RAM)

STAR genome indexing requires substantial memory.

- Building the GRCh38 STAR index typically requires ~25–30 GB RAM
- Alignment requires less memory once the index is built

Systems with insufficient memory may encounter out-of-memory (OOM) errors during
reference preparation.

### Windows / WSL2 Notes

On Windows systems using WSL2 and Docker Desktop, the default WSL2 memory limit is
often insufficient for STAR index construction.

A template WSL configuration file is provided at:

`docs/wslconfig.example`

Refer to the User Manual section on WSL2 configuration for details.


## 4. Quick Start (Minimal Path)

This section shows the simplest way to run the pipeline using the provided
defaults. No configuration changes are required.

### Docker images

This workflow can be executed using a versioned release, a digest-pinned archival image, or a local development build.

#### Versioned Release (recommended)

```bash
docker pull ghcr.io/Inkasimo/scrnaseq-pbmc-workflow:v1.0.7
```

#### Exact Archival Image (Digest-Pinned)

```bash
docker pull ghcr.io/Inkasimo/scrnaseq-pbmc-workflow@sha256:80354b76e76405636c43e73902236e0399d26978a214227afbafa46fc0555bb8
````

Use this image when strict reproducibility of archived results is required.

#### Local Development Build

```bash
docker build -t scrnaseq-workflow -f containers/Dockerfile .
```

Local builds are intended for development and testing. They are not guaranteed to be bit-identical to published container releases.

### 2. (Optional) Install wrapper dependencies

The Python wrapper requires a local Python installation (≥ 3.9) and one dependency.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r wrapper-requirements.txt
```

Required only if using run_analysis.py.
Not needed if running Snakemake directly via Docker.

Activate the environment every time before using run_analysis.py.

### 3. Run the pipeline (wrapper)

Run a dry run to inspect what would be executed:

```bash 
python3 run_analysis.py all --dry-run
```

Inspect the available sections: 

```bash 
python3 run_analysis.py --list-sections
```

Run the full upstream workflow:

```bash 
python3 run_analysis.py upstream --cpus 8 --cores 8
```

This will:
- Download the PBMC FASTQ data
- Run quality control
- Prepare references
- Perform STARsolo alignment
- Generate gene–cell count matrices

### 4. Check outputs

Pipeline outputs will appear under:

`results/`

Quality control reports and alignment outputs are organized by step and donor.

## 5. Configuration (config/config.yaml)

The pipeline is configured via a single YAML file located at:

`config/config.yaml`

A complete, working configuration is provided. No changes are required to run
the pipeline using the default PBMC dataset.

### What Users May Edit

Users may modify the configuration to adapt the pipeline to different datasets
or execution environments. Common use cases include:

1. Input data sources
- Enabling or disabling FASTQ downloading
- Pointing to locally available FASTQ files instead of downloading

2. Reference data paths
- Genome FASTA and annotation (GTF) locations
- STAR index reuse vs. rebuild

3. Resource-related parameters
- Thread counts for individual tools
- Memory-related options (where applicable)

These changes are optional and intended for users who understand the implications
for reproducibility and results.

### What Users Should Not Edit

Most configuration keys control internal workflow behavior and are not intended
to be modified casually.

In particular:
- Execution toggles set automatically by the wrapper
- Paths and flags that coordinate rule dependencies
- Parameters tightly coupled to workflow logic

Changing these values may require manual cleanup of intermediate outputs.

### Wrapper Overrides

When the Python wrapper is used, certain configuration values are set automatically
based on the selected execution section (e.g. data download, trimming, alignment).

These overrides are applied at runtime and do not permanently modify
config/config.yaml.


### Configuration and Reproducibility

The configuration file is part of the repository and should be version-controlled
when results need to be reproducible.

Changes to configuration values that affect results should be accompanied by
corresponding cleanup of generated outputs.

### Barcode Whitelist

The 10x Genomics barcode whitelist required by STARsolo is bundled directly with
the pipeline at:

resources/barcodes/3M-3pgex-may-2023_TRU.txt

This file is referenced via the configuration and is not downloaded at runtime.
Bundling the whitelist avoids reliance on unstable external URLs and ensures
reproducible execution.


## 6. Running the Workflow with the Python Wrapper

The Python wrapper `run_analysis.py` provides a simplified, section-based
interface for running the Snakemake workflow inside Docker.

It requires users to explicitly choose what to run, reducing the risk of
accidental full executions and simplifying interactive use.

### Wrapper Requirements (Host-side Only)

The wrapper requires:

Python ≥ 3.9

One dependency listed in `wrapper-requirements.txt`

```bash
pip install -r wrapper-requirements.txt
```

These requirements are only needed for the wrapper.
Running Snakemake directly inside Docker does not require any host-side Python
dependencies.

### Inspecting Available Sections

To see which pipeline sections can be executed:

```bash
python3 run_analysis.py --list-sections
```

Typical sections include:

- download_data
- download_data_and_qc
- qc
- trim
- trim_and_qc
- ref
- align
- upstream
- upstream_no_download
- downstream
- all
- all_no_download
- build_seurat_object_qc
- filter_and_normalize_seurat
- cluster_annotate_seurat
- deg_and_tost
- downstream
- network inference **Work in progress, not implemented yet**
- unlock



### Inspecting Available Donors

Donors are defined in config/config.yaml.

To list them:

``` bash
python3 run_analysis.py --list-donors
```

### Running Individual Sections

#### Download data

``` bash
python3 run_analysis.py download_data --cpus 8 --cores 8
```

Downloads FASTQ files (unless disabled) and marks completion with .done files.

#### Quality control only

``` bash
python3 run_analysis.py qc --cpus 8 --cores 8
```

Runs FastQC and MultiQC on raw FASTQs.

#### Reference preparation

``` bash
python3 run_analysis.py ref --cpus 8 --cores 8
```

Prepares static reference resources, including the barcode whitelist and STAR
genome index.

#### Alignment (all donors)

``` bash
python3 run_analysis.py align \
  --donor all \
  --cpus 8 --cores 8 \
  -j 1 \
  --set-threads starsolo=8
```

Runs STARsolo alignment and generates gene–cell count matrices.

#### Alignment (selected donors)

``` bash
python3 run_analysis.py align \
  --donor donor1 \
  --donor donor3 \
  --cpus 8 --cores 8 \
  -j 1 \
  --set-threads starsolo=8
```

#### Trimmed vs Untrimmed Mode

Read trimming is disabled by default.

To enable trimming for all applicable downstream steps:

``` bash
python3 run_analysis.py trim --cpus 8 --cores 8
```

To run alignment on trimmed reads:

``` bash
python3 run_analysis.py align \
  --donor all \
  --trimmed \
  --cpus 8 --cores 8
```

Trimmed outputs are written to separate directories and do not overwrite
untrimmed results.

#### Running upstream (download + QC + reference + alignment)

Runs the full upstream workflow for all donors defined in `config/config.yaml`.

```bash
python3 run_analysis.py upstream --cpus 8 --cores 8
```

To run upstream in trimmed mode: 

```bash
python3 run_analysis.py upstream --trimmed --cpus 8 --cores 8
```
 
#### Running upstream without downloading FASTQs

Runs the upstream workflow assuming FASTQs already exist locally (no downloads).

```bash
python3 run_analysis.py upstream_no_download --cpus 8 --cores 8
```

Trimmed mode is also supported:

```bash
python3 run_analysis.py upstream_no_download --trimmed --cpus 8 --cores 8
```
 

#### Running the Full Workflow Explicitly

To run the full workflow using the wrapper:

``` bash
python3 run_analysis.py all --cpus 8 --cores 8
```

A dry run can be performed without executing any steps:

``` bash
python3 run_analysis.py all --dry-run
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

**An example:** 
If the resources allow you can run alignment in parallel: 

```bash
python3 run_analysis.py align \
  --donor all \
  --cpus 16 \
  --cores 16 \
  -j 2 \
  --set-threads starsolo=8
```

## 7. Using a Custom Docker Image with the Wrapper

If you are using a digest-pinned image or a locally built image instead of the versioned release, 
you must explicitly define the Docker image when invoking the wrapper:

```bash
python3 run_analysis.py all --image <image_reference>

```

Example (digest-pinned image):

```bash
python3 run_analysis.py all --image ghcr.io/Inkasimo/scrnaseq-pbmc-workflow@sha256:80354b76e76405636c43e73902236e0399d26978a214227afbafa46fc0555bb8
```

Example (local build):

```bash
python3 run_analysis.py all --image scrnaseq-workflow
```


## 8. Running the Workflow Directly with Snakemake

The workflow can be executed directly using Snakemake inside the Docker container.
This mode provides maximal transparency and control but requires familiarity with
Snakemake’s command-line interface.

Direct Snakemake execution is most useful for:

- Inspecting the workflow DAG
- Running specific targets
- Advanced or ad hoc usage

Dry Run (Inspect the DAG)

``` bash
docker run --rm -it \
  -v "$(pwd)":/work \
  -w /work \
  scrnaseq-workflow \
  snakemake -n -p
```

This prints the planned execution steps without running any commands.

### Running a Specific Target

Snakemake targets correspond to output files produced by the workflow.

Example: run alignment for a single donor:

``` bash
docker run --rm -it \
  -v "$(pwd)":/work \
  -w /work \
  scrnaseq-workflow \
  snakemake results/alignment/starsolo/donor1/starsolo.done
```

Only the rules required to produce the specified target will be executed.

### Notes on Direct Execution

- Users are responsible for selecting appropriate Snakemake options
(e.g. --cores, -j, --rerun-incomplete)
- Direct execution bypasses the safeguards provided by the wrapper
- Partial outputs may remain after interruptions and should be handled carefully


## 9. Outputs and Directory Layout

The pipeline produces outputs under a small number of top-level directories.
Large data and generated results are intentionally excluded from version control.

### data/

The `data/` directory contains input data and reference files generated or
downloaded at runtime.

Typical contents include:
- `data/raw/`
Raw FASTQ files, either downloaded automatically or provided locally
- `data/trimmed/`
Trimmed FASTQ files (only present if trimming is enabled)
- `data/ref/`
Reference genome files, annotations, barcode whitelist, and STAR index

This directory is not versioned.

### results/

The `results/` directory contains all pipeline outputs and logs.

Key subdirectories include:

- `results/qc/`
FastQC outputs and MultiQC reports for raw and trimmed data

- `results/alignment/`
STARsolo alignment outputs organized by donor and mode (raw or trimmed)

- `results/logs/`
Execution logs for individual pipeline steps

This directory is not versioned.

## Completion Markers (.done Files)

Many directories contain `.done` files that indicate successful completion of
a pipeline step.

A .done file:

- Is created only after a step completes successfully
- Serves as the authoritative success indicator for that step
- Allows safe interruption and rerunning of the workflow

The presence or absence of a `.done` file determines whether a step will be
considered complete by the workflow.

## 10. Trimming (Optional)

Read trimming is optional and disabled by default in this pipeline.

When enabled, trimming is performed using Cutadapt and is applied in a controlled
and limited manner to preserve the structure required for scRNA-seq alignment.

### Trimming Behavior

- Only Read 2 (R2 / cDNA read) is trimmed
- Read 1 (R1 / cell barcode + UMI) is left unchanged
- Trimming is applied uniformly across all donors
- Read pairs are discarded if the trimmed R2 fails minimum length requirements

### Enabling Trimming

Trimming can be enabled explicitly via the Python wrapper:

```bash
python3 run_analysis.py trim --cpus 8 --cores 8
```

Once trimming is enabled, all downstream steps that depend on FASTQ inputs
(e.g. QC and alignment) will operate on the trimmed reads.

### Trimmed vs Untrimmed Outputs

Trimmed data are written to separate directories and do not overwrite untrimmed
outputs.

Typical locations include:

- `data/trimmed/<donor>/` for trimmed FASTQs
- `results/qc/fastqc/trimmed/` for trimmed QC reports
- `results/alignment/starsolo/trimmed/` for trimmed alignment results

This allows direct comparison between trimmed and untrimmed workflows if desired.

When running downstream steps explicitly, the `--trimmed` flag instructs the
wrapper to use trimmed FASTQs as inputs instead of raw reads. This flag must be
provided when running alignment or other steps that should operate on trimmed data.

### Trimming Parameters

Trimming parameters are defined in `config/config.yaml`, including:

- Adapter sequence removal for R2
- Quality trimming thresholds
- Minimum read length


## 11. Troubleshooting

### Snakemake Lock Errors

If execution is interrupted (e.g. Ctrl+C, system sleep, WSL shutdown), Snakemake may
leave a stale lock file that prevents subsequent runs.

To remove the lock:

```bash 
python3 run_analysis.py unlock
```

Alternatively, the lock can be cleared directly via Snakemake inside Docker.

### Incomplete Outputs After Interruption

If a run is interrupted, partial output files may remain on disk.

Incomplete steps are indicated by the absence of a corresponding .done file.
Rerunning the same command with --rerun-incomplete enabled will rerun any
incomplete steps.

The wrapper enables --rerun-incomplete by default.

### Out-of-Memory Errors During STAR Indexing

STAR genome index construction requires substantial memory.

Symptoms may include:
- Sudden termination during reference preparation
- Out-of-memory (OOM) errors reported by Docker or WSL

Ensure that at least 25–30 GB RAM is available to the Docker container.
On Windows systems using WSL2, the default memory limit may need to be increased
(see Section 3).

### Unexpected Re-execution of Steps

If configuration values that affect outputs are changed, Snakemake may rerun
downstream steps.

If outputs appear inconsistent with expectations:
- Verify the current configuration
- Remove affected output directories if necessary
- Rerun the workflow

### Stale or Unexpected Docker Containers

The pipeline runs Snakemake inside short-lived Docker containers. In normal
operation, containers exit automatically when a command finishes.

If a run is interrupted abruptly, a container may occasionally remain running.
This is not specific to the pipeline and reflects Docker runtime behavior.

If resource usage appears higher than expected after an interruption, verify that
no stale containers are running and stop them if necessary.

### Windows / WSL / Docker Desktop

This workflow runs inside Docker and bind-mounts the repository into the container.

On Windows, this requires Docker Desktop to have access to the repository path.
If Docker cannot bind-mount the path, the container will not be able to read
`workflow/Snakefile`.

If you encounter:

    Snakefile "workflow/Snakefile" not found

check Docker Desktop file-sharing settings, or run the workflow on a Linux or
macOS system where bind mounts are always available.

## 12. Limitations and Planned Extensions

### Scope Limitations

This pipeline is intentionally scoped to a specific class of scRNA-seq workflows.

The following limitations apply by design:
- Support is focused on 10x Genomics 3′ scRNA-seq data
- Reference configuration targets human datasets
- The pipeline is not designed for multi-project orchestration or large-scale
production deployment
- Interactive analysis environments (e.g. notebooks) are not embedded into the
execution workflow
- HPC execution is outside the scope of this project

These boundaries reflect deliberate design choices rather than missing
functionality.

#### Out of Scope

The pipeline does not attempt to:
- Serve as a general-purpose workflow framework
- Abstract away biological decision-making
- Automatically optimize parameters across heterogeneous datasets


### Running on HPC / Batch Systems (Notes)

This pipeline is designed to run in a containerized environment and can, in
principle, be executed on HPC systems that support Docker or compatible
runtimes (e.g. Singularity/Apptainer).

- Snakemake supports native integration with common schedulers (e.g. SLURM),
but scheduler-specific execution is not configured or tested in this project.
- The Python wrapper is designed for interactive and local execution and does
not provide scheduler submission logic.
- Users wishing to run this workflow on HPC systems may adapt the execution
model by:
1. Running Snakemake directly inside a scheduler job
2. Replacing Docker with a compatible container runtime if required
-Resource usage (CPU, memory) should be adapted to local cluster policies.

HPC execution is outside the scope of this project and intentionally left to
user customization.



## Final Notes

This repository is intended to serve both as a functional scRNA-seq workflow and
as a clear, defensible example of pipeline design and implementation.


