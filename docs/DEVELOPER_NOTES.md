# Developer / maintainer notes (including environment maintenance)

STAR indexing details, and Docker CPU behavior) are kept in
`docs/DEVELOPER_NOTES.md` and are not required to run the workflow. 


NOTE:
STAR genome indexing is CPU-parallel but Snakemake defaults to 1 core
unless --cores / -j is specified at runtime.

When building the STAR index, always run Snakemake with multiple cores,
e.g.:

  snakemake -j 4 data/ref/star_index.done
  (or -j 8 if memory allows)

Otherwise STAR will run with --runThreadN 1 and indexing will be
unnecessarily slow.


#Run index

docker run --rm -it --cpus 4 -v "$(pwd)":/work scrnaseq-workflow \
  snakemake data/ref/whitelist.done data/ref/star_index.done -j 4 --rerun-incomplete
  
# SLURM

SEE SNAKEMAKE DOCUMENTATION ON CLUSTER PROFILES. DONT IMPLEMENT. 
YOU CAN WRAP THIS IN SLURM IF YOU CAN WORK IN SLURM
  
 

1) --cpus (Docker)
docker run --cpus 8 …


What it controls

Linux cgroup CPU quota

Hard upper limit on how many CPU cores the container can use

Scope

Entire container

All processes inside it combined

If you set

--cpus 8


Then:

STAR, FastQC, Java, Python, everything together cannot exceed ~8 cores

If you don’t set it → container can use all host CPUs

Think of this as:

“How much of my machine am I willing to give this container?”


Run the container with your repo mounted:

docker run --rm -it \
  -v /mnt/g/scRNAseq-pbmc-workflow:/work \
  -w /work \
  scrnaseq-workflow \
  bash

#Renv notes

A) Sanity: confirm Seurat is still vendor/local
R -q -e 'renv::activate("/work"); cat("Seurat:", as.character(packageVersion("Seurat")), "\n"); cat("SeuratObject:", as.character(packageVersion("SeuratObject")), "\n");'

B) Install CRAN msigdbr (safe)
R -q -e 'renv::activate("/work"); install.packages("msigdbr", repos=Sys.getenv("RENV_CONFIG_REPOS_OVERRIDE"))'

C) Install Bioconductor packages with updates blocked

This is the critical part—do NOT run the “BiocManager::install(version=...)” line now.

R -q -e 'renv::activate("/work"); BiocManager::install(version="3.18", ask=FALSE, update=FALSE)'
R -q -e 'renv::activate("/work"); BiocManager::install(c("DESeq2","fgsea","clusterProfiler"), ask=FALSE, update=FALSE)'

D) Verify packages are installed
R -q -e 'renv::activate("/work"); pkgs<-c("DESeq2","fgsea","clusterProfiler","msigdbr


# Docker bash

docker run --rm -it \
  -v "$(pwd)":/work \
  -w /work \
  scrnaseq-workflow \
  bash
