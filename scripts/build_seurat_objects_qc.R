#Build Seurat objects Scaffold

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
})

option_list <- list(
  make_option("--mtx", type = "character"),
  make_option("--features", type = "character"),
  make_option("--barcodes", type = "character"),
  make_option("--outdir", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

counts <- ReadMtx(
  mtx = opt$mtx,
  features = opt$features,
  cells = opt$barcodes
)

obj <- CreateSeuratObject(counts)

saveRDS(obj, file.path(opt$outdir, "seurat.rds"))
