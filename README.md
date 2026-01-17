# scRNA-seq PBMC workflow

Production-style scRNA-seq pipeline demo using Docker + Snakemake.
Focus: execution, structure, reproducibility.

Dataset: 10x PBMC 5k donors 1–4 (3').

## WSL2 memory configuration (required for STAR index)

STAR genome indexing for GRCh38 requires ~25–30 GB RAM.

If you are running on Windows with WSL2, you must increase WSL memory.
A template is provided at:

    docs/wslconfig.example

To enable it:
1. Copy the file to:
       C:\Users\<your-username>\.wslconfig
2. Adjust values if needed
3. Restart WSL:
       wsl --shutdown
4. Restart Docker Desktop
