#!/usr/bin/env Rscript

# ADD HARD CAPS TO CONFIG

#FILTER

#mt_quantile (0.90)
#mt_cap (25)
#nFeature_min_quantile (0.05)
#nFeature_max_quantile (0.99)
#nFeature_cap (6000)
#nCount_min_quantile (0.05)
#nCount_max_quantile (0.99)
#hb_max (1 or null/disable)

#normalization

#method ("LogNormalize")
#scale_factor (10000)
#n_hvg (2000)
#selection.method = "vst"

regress (["percent.mt"])
# QC-filter Seurat object 
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
})

option_list <- list(
  make_option("--seurat",     type = "character", help = "Path to base Seurat RDS (e.g., donorX_object.rds)"),
  make_option("--qc_metrics", type = "character", help = "Path to qc_metrics.tsv produced by base script"),
  make_option("--outdir",     type = "character", help = "Output directory (e.g., .../donorX/qcfilt_norm)"),
  make_option("--donor",      type = "character", default = "donor")
)

opt <- parse_args(OptionParser(option_list = option_list))

# -------------------------
# Directories
# -------------------------
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
plots_dir <- file.path(opt$outdir, "plots_qcfilt")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Load inputs
# -------------------------
#obj<-"G:/scRNAseq-pbmc-workflow-copy/scRNAseq-pbmc-workflow/results/downstream/seurat/untrimmed/donor1/donor1_object.rds"
#obj<-readRDS(obj)

obj <- readRDS(opt$seurat)

#metrics<-"G:/scRNAseq-pbmc-workflow-copy/scRNAseq-pbmc-workflow/results/downstream/seurat/untrimmed/donor1/qc_metrics.tsv"
#metrics <- read.table(metrics, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
#if (nrow(metrics) != 1) {
  #stop("qc_metrics.tsv must contain exactly one row for this donor.")
#}

metrics <- read.table(opt$qc_metrics, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(metrics) != 1) {
  stop("qc_metrics.tsv must contain exactly one row for this donor.")
}

# -------------------------
# Ensure QC fields exist (robust if object moved)
# -------------------------
gene_names <- rownames(obj)

if (!"percent.mt" %in% colnames(obj@meta.data) || all(is.na(obj@meta.data$percent.mt))) {
  if (any(grepl("^MT-", gene_names))) {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  } else if (any(grepl("^mt-", gene_names))) {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  } else {
    obj[["percent.mt"]] <- rep(NA_real_, ncol(obj))
    warning("No MT-/mt- genes detected; percent.mt set to NA.")
  }
}

if (!"percent.ribo" %in% colnames(obj@meta.data) || all(is.na(obj@meta.data$percent.ribo))) {
  if (any(grepl("^RPL|^RPS", gene_names))) {
    obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^(RPL|RPS)")
  } else if (any(grepl("^Rpl|^Rps", gene_names))) {
    obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^(Rpl|Rps)")
  } else {
    obj[["percent.ribo"]] <- rep(NA_real_, ncol(obj))
    warning("No ribosomal genes detected; percent.ribo set to NA.")
  }
}

if (!"percent.hb" %in% colnames(obj@meta.data) || all(is.na(obj@meta.data$percent.hb))) {
  if (any(grepl("^HB[ABDEGMQZ]", gene_names))) {
    obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^HB[ABDEGMQZ]")
  } else if (any(grepl("^Hba|^Hbb", gene_names))) {
    obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^(Hba|Hbb)")
  } else {
    obj[["percent.hb"]] <- rep(NA_real_, ncol(obj))
  }
}

# -------------------------
# Derive thresholds from qc_metrics.tsv
# Policy:
#  - mt: percent.mt <= min(mt_q90, 25)
#  - nFeature: >= nFeature_q05 and <= min(nFeature_q99, 6000)
#  - nCount: >= nCount_q05 and <= nCount_q99
#  - hb (PBMC): percent.hb <= 1
# -------------------------
get1 <- function(col) {
  if (!col %in% colnames(metrics)) stop(paste("Missing column in qc_metrics.tsv:", col))
  as.numeric(metrics[[col]][1])
}

thr <- data.frame(
  #donor="donor1",
  donor = opt$donor,
  mt_max = min(get1("mt_q90"), 25),
  nFeature_min = get1("nFeature_q05"),
  nFeature_max = min(get1("nFeature_q99"), 6000),
  nCount_min = get1("nCount_q05"),
  nCount_max = get1("nCount_q99"),
  hb_max = 1
)

write.table(thr, file = file.path(opt$outdir, "qc_thresholds.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# -------------------------
# Apply filtering + record reasons
# -------------------------
md <- obj@meta.data
n_before <- nrow(md)

fail_mt <- if (!all(is.na(md$percent.mt))) (md$percent.mt > thr$mt_max[1]) else rep(FALSE, n_before)
fail_nf_low <- md$nFeature_RNA < thr$nFeature_min[1]
fail_nf_high <- md$nFeature_RNA > thr$nFeature_max[1]
fail_nc_low <- md$nCount_RNA   < thr$nCount_min[1]
fail_nc_high <- md$nCount_RNA  > thr$nCount_max[1]
fail_hb <- if (!all(is.na(md$percent.hb))) (md$percent.hb > thr$hb_max[1]) else rep(FALSE, n_before)

fail_any <- fail_mt | fail_nf_low | fail_nf_high | fail_nc_low | fail_nc_high | fail_hb

summary_tbl <- data.frame(
  #donor="donor1",
  donor = opt$donor,
  n_before = n_before,
  n_after = sum(!fail_any),
  n_dropped = sum(fail_any),
  pct_dropped = (sum(fail_any) / n_before) * 100,
  n_fail_mt = sum(fail_mt),
  n_fail_nFeature_low = sum(fail_nf_low),
  n_fail_nFeature_high = sum(fail_nf_high),
  n_fail_nCount_low = sum(fail_nc_low),
  n_fail_nCount_high = sum(fail_nc_high),
  n_fail_hb = sum(fail_hb),
  n_fail_multiple = sum((fail_mt + fail_nf_low + fail_nf_high + fail_nc_low + fail_nc_high + fail_hb) > 1)
)

write.table(summary_tbl, file = file.path(opt$outdir, "qc_filter_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

obj_filt <- subset(obj, cells = rownames(md)[!fail_any])

# Save QC-filtered object
saveRDS(obj_filt, file.path(opt$outdir, paste0(opt$donor, "_qcfilt_object.rds")))

# -------------------------
# Plots AFTER filtering
# -------------------------
qc_feats <- c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo","percent.hb")
qc_feats <- qc_feats[qc_feats %in% colnames(obj_filt@meta.data)]
qc_feats <- qc_feats[!sapply(qc_feats, function(x) all(is.na(obj_filt@meta.data[[x]])))]

p_vln <- VlnPlot(obj_filt, features = qc_feats, pt.size = 0.1, ncol = length(qc_feats))
ggsave(file.path(plots_dir, "qc_violin.png"), plot = p_vln, width = 10, height = 4, dpi = 300)

p_sc1 <- FeatureScatter(obj_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(file.path(plots_dir, "scatter_nFeat_RNA.png"), plot = p_sc1, width = 5, height = 4, dpi = 300)

df <- obj_filt@meta.data
p_hist1 <- ggplot(df, aes(nFeature_RNA)) + geom_histogram(bins = 60) + scale_x_log10()
p_hist2 <- ggplot(df, aes(nCount_RNA))   + geom_histogram(bins = 60) + scale_x_log10()
ggsave(file.path(plots_dir, "hist_nFeat_RNA.png"), plot = p_hist1, width = 5, height = 4, dpi = 300)
ggsave(file.path(plots_dir, "hist_nCount_RNA.png"), plot = p_hist2, width = 5, height = 4, dpi = 300)

if ("percent.mt" %in% colnames(df) && !all(is.na(df$percent.mt))) {
  p_sc2 <- FeatureScatter(obj_filt, feature1 = "nCount_RNA", feature2 = "percent.mt")
  ggsave(file.path(plots_dir, "scatter_mt_RNA.png"), plot = p_sc2, width = 5, height = 4, dpi = 300)

  p_hist3 <- ggplot(df, aes(percent.mt)) + geom_histogram(bins = 60)
  ggsave(file.path(plots_dir, "hist_mt_RNA.png"), plot = p_hist3, width = 5, height = 4, dpi = 300)
}

# marker file
#file.create(file.path(opt$outdir, "seurat_qcfilt.done"))


####Normalization. Apply first lognormalize plus findvariable features vst. 
#### Think about SCtransfrom. It is an option. 

obj_filt <- NormalizeData(
  obj_filt,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)


obj_filt <- FindVariableFeatures(
  obj_filt,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

vars_to_regress <- c("percent.mt")  # optionally c("percent.mt", "nCount_RNA")

obj_filt <- ScaleData(
  obj_filt,
  features = VariableFeatures(obj_filt),
  vars.to.regress = vars_to_regress,
  verbose = FALSE
)


saveRDS(obj_filt, file.path(opt$outdir, paste0(opt$donor, "_qcfilt_norm_object.rds")))

#Normalization QC

p_hvg <- VariableFeaturePlot(obj_filt)
ggsave(file.path(plots_dir, "norm_variable_features.png"), p_hvg, width=6, height=4, dpi=300)

top10 <- head(VariableFeatures(obj_filt), 10)
p_hvg_lbl <- LabelPoints(plot = p_hvg, points = top10, repel = TRUE)
ggsave(file.path(plots_dir, "norm_variable_features_top10.png"), p_hvg_lbl, width=7, height=5, dpi=300)

qc_feats <- c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo","percent.hb")
qc_feats <- qc_feats[qc_feats %in% colnames(obj_filt@meta.data)]
qc_feats <- qc_feats[!sapply(qc_feats, function(x) all(is.na(obj_filt@meta.data[[x]])))]
p_vln_post <- VlnPlot(obj_filt, features = qc_feats, pt.size = 0.1, ncol = length(qc_feats))
ggsave(file.path(plots_dir, "qc_violin_postfilter.png"), p_vln_post, width = 10, height = 4, dpi = 300)

genes_check <- intersect(c("LST1","MS4A1","NKG7"), rownames(obj_filt))
if (length(genes_check) > 0) {
  p_feat <- VlnPlot(obj_filt, features = genes_check, pt.size = 0, ncol = length(genes_check))
  ggsave(file.path(plots_dir, "norm_marker_violin.png"), p_feat, width=10, height=4, dpi=300)
}


#normalization summary

meta <- obj_filt@meta.data

norm_metrics <- data.frame(
  donor = opt$donor,
  n_cells = ncol(obj_filt),
  n_genes = nrow(obj_filt),
  normalization_method = "LogNormalize",
  scale_factor = 10000,
  n_hvg = length(VariableFeatures(obj_filt)),
  vars_regressed = paste(vars_to_regress, collapse = ","),
  median_nCount_RNA = median(meta$nCount_RNA, na.rm=TRUE),
  median_percent_mt = if ("percent.mt" %in% colnames(meta)) median(meta$percent.mt, na.rm=TRUE) else NA_real_
)

write.table(
  norm_metrics,
  file = file.path(opt$outdir, "norm_metrics.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

