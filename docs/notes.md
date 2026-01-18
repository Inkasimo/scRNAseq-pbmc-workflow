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
yOU CAN WRAP THIS IN SLURM IF YOU CAN WORK IN SLURM
  
 
# Check process in container

docker ps
