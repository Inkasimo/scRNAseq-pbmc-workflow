library(Seurat)
library(Matrix)
library(ggplot2)

base <- "results/alignment/starsolo/raw/donor1/Solo.out/Gene/filtered"

counts <- ReadMtx(
  mtx = file.path(base, "matrix.mtx"),
  features = file.path(base, "features.tsv"),
  cells = file.path(base, "barcodes.tsv")
)

obj <- Seurat::CreateSeuratObject(counts)


# QC

# 1) Detect species-ish gene naming to pick sensible patterns
gene_names <- rownames(obj)

is_human_like <- any(grepl("^MT-", gene_names)) || any(grepl("^RPL|^RPS", gene_names))
is_mouse_like <- any(grepl("^mt-", gene_names)) || any(grepl("^Rp[sl]", gene_names))

# 2) Compute QC percentages
# Mito
if (any(grepl("^MT-", gene_names))) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")     # human
} else if (any(grepl("^mt-", gene_names))) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")     # mouse
} else {
  obj[["percent.mt"]] <- NA_real_
  warning("No MT-/mt- genes detected; percent.mt set to NA.")
}

# Ribosomal
# human: RPL/RPS ; mouse: Rpl/Rps
if (any(grepl("^RPL|^RPS", gene_names))) {
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^(RPL|RPS)")
} else if (any(grepl("^Rpl|^Rps", gene_names))) {
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^(Rpl|Rps)")
} else {
  obj[["percent.ribo"]] <- NA_real_
  warning("No ribosomal genes detected; percent.ribo set to NA.")
}

# Hemoglobin (useful for blood contamination)
# human: HBA/HBB... mouse: Hba/Hbb...
if (any(grepl("^HB[ABDEGMQZ]", gene_names))) {
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^HB[ABDEGMQZ]")
} else if (any(grepl("^Hba|^Hbb", gene_names))) {
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^(Hba|Hbb)")
} else {
  obj[["percent.hb"]] <- NA_real_
}

# 3) Core QC plots
# Violin plots (raw)
p_vln <- VlnPlot(
  obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"),
  pt.size = 0.1,
  ncol = 5
)
print(p_vln)


#obj_filt <- subset(
  #obj,
  #subset =
    #nFeature_RNA >= 500 &
    #nFeature_RNA <= 6000 &
    #nCount_RNA   >= 1000 &
    #nCount_RNA   <= 50000 &
    #percent.mt   <= 20 &
    #percent.hb <= 1
)

# Scatter: complexity and mito relationship
p_sc1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p_sc2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
print(p_sc1 + p_sc2)

# Distributions (useful for choosing cutoffs)
df <- obj@meta.data
p_hist1 <- ggplot(df, aes(nFeature_RNA)) + geom_histogram(bins = 60) + scale_x_log10()
p_hist2 <- ggplot(df, aes(nCount_RNA))   + geom_histogram(bins = 60) + scale_x_log10()
p_hist3 <- ggplot(df, aes(percent.mt))   + geom_histogram(bins = 60)
print(p_hist1)
print(p_hist2)
print(p_hist3)

# Re-plot after filtering (highly recommended)
print(VlnPlot(obj_qc, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo","percent.hb"), pt.size = 0.1, ncol = 5))
print(FeatureScatter(obj_qc, "nCount_RNA", "nFeature_RNA") + FeatureScatter(obj_qc, "nCount_RNA", "percent.mt"))

