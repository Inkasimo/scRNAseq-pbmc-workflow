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

# Make several config files where the booleans differ

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