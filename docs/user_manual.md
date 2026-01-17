1. Purpose

One paragraph:

What the pipeline does

What it does not do

Example:

This pipeline processes 10x Genomics GEM-X 3′ v4 scRNA-seq data using STARsolo to produce gene–cell count matrices. It does not perform downstream analysis or visualization.

2. Requirements

Bullet list:

Docker

(Windows) WSL2 with increased memory (link to docs/wslconfig.example)

Disk space note (e.g. “~50–70 GB recommended”)

3. Quick start

Minimal, copy-pasteable:

docker build -t scrnaseq-workflow .

docker run --rm -it \
  -v "$(pwd)":/work \
  -w /work \
  scrnaseq-workflow \
  snakemake -j 1

4. Configuration

Explain only:

where config.yaml lives

what users are expected to edit:

FASTQ URLs / paths

number of threads (optional)

Do not document every config key.

5. Running specific steps (very useful)

Examples only:

# Build reference only
snakemake data/ref/star_index.done

# Run alignment for one donor
snakemake results/alignment/starsolo/donor1/starsolo.done


This saves users time.

6. Outputs

High-level only:

results/qc/ → FastQC / MultiQC

results/alignment/starsolo/<donor>/Solo.out/ → count matrices

logs/ → execution logs