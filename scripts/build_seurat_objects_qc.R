#Build Seurat objects Scaffold
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

option_list <- list(
  make_option("--mtx", type = "character"),
  make_option("--features", type = "character"),
  make_option("--barcodes", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--donor", type = "character", default = "donor")
)

#base <- "/mnt/g/scRNAseq_pbmc_workflow/results/alignment/starsolo/raw/donor1/Solo.out/Gene/filtered"

#counts <- ReadMtx(
  #mtx = file.path(base, "matrix.mtx"),
  #features = file.path(base, "features.tsv"),
  #cells = file.path(base, "barcodes.tsv")
#)

opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
plots_dir <- file.path(opt$outdir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

counts <- ReadMtx(
  mtx = opt$mtx,
  features = opt$features,
  cells = opt$barcodes
)

obj <- CreateSeuratObject(counts)
saveRDS(obj, file.path(opt$outdir, paste0(opt$donor, "_object.rds")))


#####QC


gene_names <- rownames(obj)

# Mito
if (any(grepl("^MT-", gene_names))) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
} else if (any(grepl("^mt-", gene_names))) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
} else {
  obj[["percent.mt"]] <- rep(NA_real_, ncol(obj))
  warning("No MT-/mt- genes detected; percent.mt set to NA.")
}

# Ribo
if (any(grepl("^RPL|^RPS", gene_names))) {
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^(RPL|RPS)")
} else if (any(grepl("^Rpl|^Rps", gene_names))) {
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^(Rpl|Rps)")
} else {
  obj[["percent.ribo"]] <- rep(NA_real_, ncol(obj))
  warning("No ribosomal genes detected; percent.ribo set to NA.")
}

# Hb
if (any(grepl("^HB[ABDEGMQZ]", gene_names))) {
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^HB[ABDEGMQZ]")
} else if (any(grepl("^Hba|^Hbb", gene_names))) {
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^(Hba|Hbb)")
} else {
  obj[["percent.hb"]] <- rep(NA_real_, ncol(obj))
}

p_vln <- VlnPlot(
  obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"),
  pt.size = 0.1,
  ncol = 5
)
ggsave(file.path(plots_dir, "qc_violin.png"), plot = p_vln, width = 10, height = 4, dpi = 300)

p_sc1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p_sc2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
ggsave(file.path(plots_dir, "scatter_nFeat_RNA.png"), plot = p_sc1, width = 5, height = 4, dpi = 300)
ggsave(file.path(plots_dir, "scatter_mt_RNA.png"), plot = p_sc2, width = 5, height = 4, dpi = 300)

df <- obj@meta.data
p_hist1 <- ggplot(df, aes(nFeature_RNA)) + geom_histogram(bins = 60) + scale_x_log10()
p_hist2 <- ggplot(df, aes(nCount_RNA))   + geom_histogram(bins = 60) + scale_x_log10()
p_hist3 <- ggplot(df, aes(percent.mt))   + geom_histogram(bins = 60)

ggsave(file.path(plots_dir, "hist_nFeat_RNA.png"), plot = p_hist1, width = 5, height = 4, dpi = 300)
ggsave(file.path(plots_dir, "hist_nCount_RNA.png"), plot = p_hist2, width = 5, height = 4, dpi = 300)
ggsave(file.path(plots_dir, "hist_mt_RNA.png"), plot = p_hist3, width = 5, height = 4, dpi = 300)

# -------------------------
# QC TABLE 
# -------------------------
qc <- obj@meta.data

q_probs <- c(.01, .05, .10, .25, .50, .75, .90, .95, .99)

q_fmt <- function(x) {
  qs <- quantile(x, probs = q_probs, na.rm = TRUE, names = FALSE)
  names(qs) <- paste0("q", sprintf("%02d", q_probs * 100))
  qs
}

nf_q  <- q_fmt(qc$nFeature_RNA)
nc_q  <- q_fmt(qc$nCount_RNA)
mt_q  <- q_fmt(qc$percent.mt)

metrics <- data.frame(
  donor = opt$donor,
  n_cells = nrow(qc),

  median_nFeature_RNA = median(qc$nFeature_RNA),
  median_nCount_RNA   = median(qc$nCount_RNA),
  median_percent_mt   = median(qc$percent.mt, na.rm = TRUE),

  pct_cells_mt_gt_20  = mean(qc$percent.mt > 20, na.rm = TRUE) * 100,
  pct_cells_mt_gt_25  = mean(qc$percent.mt > 25, na.rm = TRUE) * 100,
  pct_cells_mt_gt_30  = mean(qc$percent.mt > 30, na.rm = TRUE) * 100,

  pct_cells_nFeature_lt_500  = mean(qc$nFeature_RNA < 500) * 100,
  pct_cells_nFeature_gt_6000 = mean(qc$nFeature_RNA > 6000) * 100,

  n_drop_mt20 = sum(qc$percent.mt > 20, na.rm = TRUE),
  n_drop_mt25 = sum(qc$percent.mt > 25, na.rm = TRUE),
  n_drop_mt30 = sum(qc$percent.mt > 30, na.rm = TRUE)
)

# bind quantiles with prefixes so it stays one row per donor
metrics <- cbind(
  metrics,
  setNames(as.data.frame(as.list(nf_q)), paste0("nFeature_", names(nf_q))),
  setNames(as.data.frame(as.list(nc_q)), paste0("nCount_", names(nc_q))),
  setNames(as.data.frame(as.list(mt_q)), paste0("mt_", names(mt_q)))
)


write.table(
  metrics,
  file = file.path(opt$outdir, "qc_metrics.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


#obj_filt <- subset(
  #obj,
  #subset =
    #nFeature_RNA >= 500 &
    #nFeature_RNA <= 6000 &
    #nCount_RNA   >= 1000 &
    #nCount_RNA   <= 50000 &
    #percent.mt   <= 20 &
    #percent.hb <= 1
#)


#
#
#
