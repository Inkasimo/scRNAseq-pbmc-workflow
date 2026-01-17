Below is a design note you can drop straight into your repo (e.g. docs/network_design.md or a commented block in the Snakefile). It is written as a forward-looking architectural note, not as something you are committing right now.

Design note: Config-driven, robust cell-type–specific coexpression networks
Scope and intent

This pipeline is primarily an alignment + counting workflow (STARsolo-based).
Network construction is downstream analysis, included to demonstrate extensibility and reproducibility, not to make statistical claims.

The goal of the network component is robustness, not scale or formal inference.

Specifically:

small number of donors (n = 4)

avoid pseudo-replication

avoid edge p-values

focus on reproducibility across donors

High-level idea

The network module is designed so that:

The cell type(s) used for network construction are specified in config.yaml

The pipeline loops over the requested cell types

A separate network is constructed per donor

Only edges reproducible across donors are retained

This shifts the definition of “significant edge” from statistical testing to cross-donor reproducibility, which is defensible for small n.

Configuration concept (not yet implemented)

The user will eventually specify something like:

network:
  enabled: true
  cell_types:
    - "T"
    - "B"
  top_k: 20
  consensus_min_donors: 3


Key idea:

The pipeline does not assume which cell type to analyze.

Cell types are treated as parameters, not hard-coded logic.

Adding/removing a cell type requires only a config change.

Data dependencies (future)

For this to work, the pipeline must eventually produce or consume:

STARsolo count matrices

Per donor

Typically Solo.out/Gene/filtered/

Per-cell metadata

Mapping of cell_barcode → cell_type

Could be produced by a downstream Scanpy/Seurat step

Must use consistent cell type labels across donors

The network module assumes that cell type labels already exist.
It does not attempt to infer cell types itself.

Network construction strategy (robust by design)
1. Cell-type subsetting per donor

For each donor:

Select only cells belonging to the requested cell type (e.g. T cells)

Use the same gene universe across donors (intersection)

This avoids:

donor-specific genes

edges driven by presence/absence artifacts

2. Avoiding cell-level pseudo-replication

Single cells are not treated as independent replicates.

Instead, one of the following strategies is used (implementation choice):

Meta-cells: randomly pool ~25 cells into one profile

Pseudo-bulk clusters: aggregate within donor-specific subclusters

This reduces noise and makes correlations more stable.

3. Donor-specific network construction

For each donor independently:

Compute gene–gene correlations

Sparsify aggressively using a top-k edges per gene rule

k = 10–30 typically

Store only the sparse edge list

This ensures:

bounded memory usage

comparable graph density across donors

4. Consensus network definition

The final network for a given cell type is defined as:

Edges that appear in at least N donors (e.g. ≥3 out of 4)

Optional constraints:

consistent sign across donors

average or median edge weight

This makes robustness explicit and easy to explain.

What this deliberately does NOT do

The design explicitly avoids:

edge p-values or permutation tests

pooled cell-level networks across donors

dense correlation matrices

claims of causality or regulation

This is intentional and aligned with the dataset size.

Why this fits a pipeline demo

This approach:

Demonstrates config-driven analysis

Keeps compute predictable and laptop-safe

Separates infrastructure (Snakemake/Docker) from analysis logic

Avoids overclaiming while still producing meaningful structure

It shows how a production pipeline can:

parameterize analysis goals

scale to additional cell types later

remain reproducible and transparent

Summary

The planned network module is:

optional

downstream

config-controlled

donor-aware

reproducibility-focused

It is not intended to compete with large-scale coexpression studies, but to illustrate how robust analytical components can be layered on top of a clean alignment/counting pipeline.