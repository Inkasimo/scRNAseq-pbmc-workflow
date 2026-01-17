NOTE:
STAR genome indexing is CPU-parallel but Snakemake defaults to 1 core
unless --cores / -j is specified at runtime.

When building the STAR index, always run Snakemake with multiple cores,
e.g.:

  snakemake -j 4 data/ref/star_index.done
  (or -j 8 if memory allows)

Otherwise STAR will run with --runThreadN 1 and indexing will be
unnecessarily slow.
