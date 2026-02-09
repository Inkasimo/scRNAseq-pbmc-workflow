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
  
  
  # Mock opt object for interactive testing (RStudio / docker shell)

opt <- list(
  seurat = paste(
    "results/downstream/seurat/untrimmed/donor1/seurat_cluster_annot/donor1_annotated_object.rds",
    "results/downstream/seurat/untrimmed/donor2/seurat_cluster_annot/donor2_annotated_object.rds",
    "results/downstream/seurat/untrimmed/donor3/seurat_cluster_annot/donor3_annotated_object.rds",
    "results/downstream/seurat/untrimmed/donor4/seurat_cluster_annot/donor4_annotated_object.rds",
    sep = ","
  ),

  outdir = "results/downstream/networks/test",

  celltype_sets = "scripts/celltype_sets.R",
  markers = "scripts/markers_pbmc.R",
  donor_names = "donor1,donor2,donor3,donor4",

  # ---- network defaults ----
  metacell_input = "linear_then_log",        # or "log"
  min_cells_per_donor_group = 200,
  metacell_size = 30,
  seed = 42,
  gene_detect_frac = 0.05,
  hvg_n = 5000,
  cor_method = "spearman",
  positive_only = TRUE,
  top_k = 100,
  min_abs_cor = 0.15,
  consensus_min_donors = 2,
  require_same_sign = TRUE,
  leiden_resolution = 1.0,
  deg_tables_dir ="results/downstream/deg_and_tost/untrimmed/deg_and_tost/tables"
)

str(opt)


## Maybe add MSigDB C7 (immunologic signatures) to enrichment at some point

ghcr.io/inkas/scrnaseq-pbmc-workflow:archive-2026-02-09


STEP 1 — Decide the GitHub image name (copy-paste this)

We will publish to GitHub Container Registry.

Use this exact name (lowercase is important):

ghcr.io/inkas/scrnaseq-pbmc-workflow:archive-2026-02-09

STEP 2 — Tag the existing image (NO rebuild)

This just gives your current image a GitHub name.

docker tag scrnaseq-workflow:latest \
  ghcr.io/inkas/scrnaseq-pbmc-workflow:archive-2026-02-09


Check:

docker images | grep ghcr.io


You should now see the same IMAGE ID (80354b76e764) with the new name.

STEP 3 — Login to GitHub Container Registry

You need a GitHub Personal Access Token with:

write:packages

read:packages

Then run:

echo YOUR_GITHUB_TOKEN | docker login ghcr.io -u inkas --password-stdin


If this succeeds, move on.

STEP 4 — Push the image to GitHub

This uploads the image. Takes a while (6 GB).

docker push ghcr.io/inkas/scrnaseq-pbmc-workflow:archive-2026-02-09


When this finishes: the image is published.

STEP 5 — Get the immutable digest (THIS IS IMPORTANT)

Run:

docker inspect ghcr.io/inkas/scrnaseq-pbmc-workflow:archive-2026-02-09 \
  --format='{{index .RepoDigests 0}}'


You’ll get something like:

ghcr.io/inkas/scrnaseq-pbmc-workflow@sha256:abc123deadbeef...


👉 Copy that exact string somewhere safe.

This is what guarantees reproducibility.

STEP 6 — How you (or anyone) runs it later

From any machine:

docker pull ghcr.io/inkas/scrnaseq-pbmc-workflow@sha256:abc123deadbeef...


Run your pipeline:

docker run --rm -it \
  -u "$(id -u):$(id -g)" \
  -v "$PWD:/work" -w /work \
  ghcr.io/inkas/scrnaseq-pbmc-workflow@sha256:abc123deadbeef... \
  snakemake -j 8