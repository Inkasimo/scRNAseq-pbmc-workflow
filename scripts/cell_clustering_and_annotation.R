#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
})

option_list <- list(
  make_option("--seurat",  type = "character", help = "Path to QC+normalized Seurat RDS (e.g., donorX_qcfilt_norm_object.rds)"),
  make_option("--markers", type = "character", help = "Path to markers script defining `markers_pbmc <- list(...)`"),
  make_option("--outdir",  type = "character", help = "Output directory (e.g., .../donorX/annot)"),
  make_option("--donor",   type = "character", default = "donor"),
  
  # Keep these fixed unless you have a reason to change:
  make_option("--dims",        type = "integer", default = 30, help = "Number of PCs to use for neighbors/UMAP"),
  make_option("--resolution",  type = "double",  default = 0.3, help = "Clustering resolution (PBMC: 0.2-0.6)"),
  make_option("--seed",        type = "integer", default = 1,   help = "Seed for deterministic UMAP")
)

opt <- parse_args(OptionParser(option_list = option_list))

# -------------------------
# Directories
# -------------------------
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
plots_dir <- file.path(opt$outdir, "plots_annot")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Load inputs
# -------------------------

obj <- readRDS(opt$seurat)

# Load marker sets
source(opt$markers)  # expected to define `markers_pbmc`
if (!exists("markers_pbmc")) stop("markers script must define object `markers_pbmc` (a named list of character vectors).")
if (!is.list(markers_pbmc) || is.null(names(markers_pbmc))) stop("`markers_pbmc` must be a *named* list.")

# Intersect markers with genes present
markers_pbmc <- lapply(markers_pbmc, function(g) intersect(g, rownames(obj)))
markers_pbmc <- markers_pbmc[sapply(markers_pbmc, length) > 0]
if (length(markers_pbmc) == 0) stop("None of the provided markers were found in the Seurat object feature names.")

# -------------------------
# Embeddings + clustering (for visualization + structure)
# -------------------------
#set.seed(42)
set.seed(opt$seed)

# Ensure required steps exist; rerun cheaply if missing
if (!"pca" %in% names(obj@reductions)) {
  obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)
}


obj <- FindNeighbors(obj, dims = 1:opt$dims, verbose = FALSE)
obj <- FindClusters(obj, resolution = opt$resolution, verbose = FALSE)

# Run UMAP if missing
if (!"umap" %in% names(obj@reductions)) {
  #obj <- RunUMAP(obj, dims = 1:dims, verbose = FALSE)
  set.seed(opt$seed)
  obj <- RunUMAP(obj, dims = 1:opt$dims, verbose = FALSE)
}

# -------------------------
# Marker-based scoring (module scores)
# -------------------------
# AddModuleScore creates columns: <name>1, <name>2, ...; so we use unique prefixes
score_prefixes <- paste0("score_", names(markers_pbmc))
names(score_prefixes) <- names(markers_pbmc)

# Add scores one marker-set at a time to keep naming predictable
for (ct in names(markers_pbmc)) {
  feats <- markers_pbmc[[ct]]
  if (length(feats) < 2) {
    warning(sprintf("Marker set '%s' has <2 genes after intersection; skipping module score.", ct))
    next
  }
  obj <- AddModuleScore(
    obj,
    features = list(feats),
    name = paste0("score_", ct, "_"),
    assay = DefaultAssay(obj),
    search = FALSE
  )
  # AddModuleScore creates score_<ct>_1
}

# Collect created score columns
score_cols <- grep("^score_.*_1$", colnames(obj@meta.data), value = TRUE)
if (length(score_cols) == 0) stop("No module score columns were created. Check marker sets and feature names.")

# -------------------------
# Per-cell label = argmax score
# -------------------------
score_mat <- obj@meta.data[, score_cols, drop = FALSE]
# Map score column back to cell type name
score_to_ct <- sub("^score_(.*)_1$", "\\1", score_cols)

max_idx <- apply(score_mat, 1, which.max)
cell_type_pred <- score_to_ct[max_idx]
obj$cell_type_pred <- factor(cell_type_pred, levels = names(markers_pbmc))

# -------------------------
# Per-cluster summary (majority vote)
# -------------------------
cl <- obj$seurat_clusters
cluster_ct_tbl <- table(cl, obj$cell_type_pred)
cluster_ct_df <- as.data.frame.matrix(cluster_ct_tbl)
cluster_ct_df$cluster <- rownames(cluster_ct_df)

# Majority label per cluster
majority_ct <- apply(cluster_ct_tbl, 1, function(x) {
  if (sum(x) == 0) return(NA_character_)
  names(which.max(x))
})
obj$cell_type_cluster_majority <- majority_ct[as.character(cl)]
obj$cell_type_cluster_majority <- factor(obj$cell_type_cluster_majority, levels = names(markers_pbmc))

# -------------------------
# Save tables
# -------------------------
# Per-cell counts
celltype_counts <- as.data.frame(table(obj$cell_type_pred), stringsAsFactors = FALSE)
colnames(celltype_counts) <- c("cell_type", "n_cells")
celltype_counts$donor <- opt$donor
#celltype_counts$donor <- "donor1"

write.table(
  celltype_counts,
  file = file.path(opt$outdir, "cell_type_counts.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Cluster ↔ predicted cell type table
cluster_ct_long <- as.data.frame(cluster_ct_tbl, stringsAsFactors = FALSE)
colnames(cluster_ct_long) <- c("cluster", "cell_type", "n_cells")
cluster_ct_long$donor <- opt$donor


write.table(
  cluster_ct_long,
  file = file.path(opt$outdir, "cluster_celltype_table.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Majority label per cluster
cluster_majority_df <- data.frame(
  donor = opt$donor,
  cluster = names(majority_ct),
  majority_cell_type = unname(majority_ct),
  stringsAsFactors = FALSE
)

write.table(
  cluster_majority_df,
  file = file.path(opt$outdir, "cluster_majority_cell_type.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# -------------------------
# Plots
# -------------------------
p_umap_cluster <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle(paste0(opt$donor, " - UMAP by cluster"))


ggsave(file.path(plots_dir, "umap_by_cluster.png"), p_umap_cluster, width = 6, height = 5, dpi = 300)

p_umap_ct <- DimPlot(obj, reduction = "umap", group.by = "cell_type_pred", label = TRUE) +
  ggtitle(paste0(opt$donor, " - UMAP by predicted cell type"))


ggsave(file.path(plots_dir, "umap_by_cell_type_pred.png"), p_umap_ct, width = 6, height = 5, dpi = 300)

# DotPlot of markers across clusters
# Flatten markers list in a stable order
marker_vec <- unique(unlist(markers_pbmc))
p_dot <- DotPlot(obj, features = marker_vec, group.by = "seurat_clusters") + RotatedAxis() +
  ggtitle(paste0(opt$donor, " - Marker DotPlot by cluster"))


ggsave(file.path(plots_dir, "dotplot_markers_by_cluster.png"), p_dot, width = 10, height = 6, dpi = 300)

# Optional: DotPlot by predicted cell type (often nicer)
p_dot_ct <- DotPlot(obj, features = marker_vec, group.by = "cell_type_pred") + RotatedAxis() +
  ggtitle(paste0(opt$donor, " - Marker DotPlot by predicted cell type"))



ggsave(file.path(plots_dir, "dotplot_markers_by_cell_type_pred.png"), p_dot_ct, width = 10, height = 6, dpi = 300)

# -------------------------
# Save annotated object
# -------------------------
saveRDS(obj, file.path(opt$outdir, paste0(opt$donor, "_annotated_object.rds")))

