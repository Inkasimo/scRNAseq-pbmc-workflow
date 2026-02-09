#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(DESeq2)
  library(ggplot2)
  library(fgsea)
  library(clusterProfiler)
  library(stringr)
})

option_list <- list(
  make_option(
  c("--seurat"),
  type = "character",
  action = "store",
  default = NULL,
  help = "Comma-separated list of annotated Seurat RDS paths (no spaces)."
),
  make_option("--celltype_sets", type="character", default="scripts/celltype_sets.R",
              help="R script defining: deg_label_col, deg_sets, deg_contrasts"),
  make_option("--outdir", type="character", help="Output directory"),
  make_option("--assay", type="character", default="RNA", help="Assay to use (default: RNA)"),
  make_option("--min_baseMean", type="double", default=10,
              help="Filter low-expression genes (DESeq2 baseMean) for TOST/DE (default 10)"),
  make_option("--min_cells_per_group", type="integer", default=50,
              help="Minimum cells required per donor x group to create pseudobulk (default 50)"),
  make_option("--seed", type="integer", default=42),
  
  make_option("--marker_padj", type="double", default=1e-10,
              help="Very strict padj threshold for marker calls (default 1e-10)"),
  make_option("--marker_lfc", type="double", default=3,
              help="Absolute log2FC threshold for marker calls (default 3)"),
  
  make_option("--equiv_alpha", type="double", default=0.05,
              help="FDR threshold for equivalence calls (TOST) (default 0.05)"),
  make_option("--equiv_delta", type="double", default=0.75,
              help="Equivalence margin (abs log2FC < delta) (default 0.75)"),
  make_option("--hallmark_gmt", type="character", default=NULL,
            help="Local Hallmark GMT file path"),
  make_option("--c7_gmt", type="character", default=NULL,
            help="Local C7 (immunologic signatures) GMT file path")
  )



opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$seurat) || !nzchar(opt$seurat)) {
  stop("Provide >=2 donors via --seurat path1,path2,...")
}
opt$seurat <- strsplit(opt$seurat, ",", fixed = TRUE)[[1]]
opt$seurat <- opt$seurat[nzchar(opt$seurat)]
if (length(opt$seurat) < 2) stop("Provide >=2 donors via --seurat path1,path2,...")


if (is.null(opt$seurat) || length(opt$seurat) < 2) {
  stop("Provide >=2 donors via repeated --seurat flags.")
}

if (is.null(opt$hallmark_gmt) || !file.exists(opt$hallmark_gmt)) {
  stop("Missing/invalid --hallmark_gmt")
}
if (is.null(opt$c7_gmt) || !file.exists(opt$c7_gmt)) {
  stop("Missing/invalid --c7_gmt")
}


set.seed(opt$seed)
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(opt$outdir, "tables"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(opt$outdir, "plots"), recursive=TRUE, showWarnings=FALSE)

# -------------------------
# Load celltype sets (your file)
# -------------------------
source(opt$celltype_sets)
if (!exists("deg_label_col")) stop("celltype_sets.R must define deg_label_col")
if (!exists("deg_sets")) stop("celltype_sets.R must define deg_sets")
if (!exists("deg_contrasts")) stop("celltype_sets.R must define deg_contrasts")

# -------------------------
# Load donors, build per-cell group labels (T_like/B_like/Mono_like)
# -------------------------

strip_dup_suffix <- function(x) sub("\\.\\d+$", "", x)

objs <- list()
donors <- character()

seurat_paths <- opt$seurat
seurat_paths <- vapply(seurat_paths, function(p) {
  normalizePath(p, winslash = "/", mustWork = FALSE)
}, character(1))

for (p in seurat_paths) {
  donor <- basename(dirname(dirname(p)))
  p <- normalizePath(p, winslash = "/", mustWork = FALSE)
  
  if (!file.exists(p)) {
    stop(sprintf("Missing RDS for %s:\n  %s", donor, p))
  }
  
  obj <- readRDS(p)
  
  fine <- as.character(obj@meta.data[[deg_label_col]])
  coarse <- setNames(rep(NA_character_, length(fine)), rownames(obj@meta.data))
  
  for (grp in names(deg_sets)) {
    coarse[fine %in% deg_sets[[grp]]] <- grp
  }
  
  obj$deg_group <- factor(coarse, levels = names(deg_sets))
  obj$donor <- donor
  
  objs[[donor]] <- obj
  donors <- c(donors, donor)
}

# -------------------------
# Pseudobulk: sum raw counts per donor x deg_group
# -------------------------
message("Building pseudobulk counts (raw counts summed)...")

pb_counts_list <- list()
pb_coldata <- data.frame()

for (donor in names(objs)) {
  obj <- objs[[donor]]
  md <- obj@meta.data
  
  keep_cells <- rownames(md)[!is.na(md$deg_group)]
  obj2 <- subset(obj, cells=keep_cells)
  md2 <- obj2@meta.data
  
  tab <- table(md2$deg_group)
  valid_groups <- names(tab)[tab >= opt$min_cells_per_group]
  obj2 <- subset(obj2, cells=rownames(md2)[md2$deg_group %in% valid_groups])
  md2 <- obj2@meta.data
  
  if (ncol(obj2) == 0) {
    warning(sprintf("Donor %s: no cells after filtering by min_cells_per_group.", donor))
    next
  }
  
  counts <- GetAssayData(obj2, slot="counts")  # sparse dgCMatrix
  
  for (grp in levels(obj2$deg_group)) {
    cells_grp <- rownames(md2)[md2$deg_group == grp]
    if (length(cells_grp) < opt$min_cells_per_group) next
    
    mat <- counts[, cells_grp, drop=FALSE]
    pb <- Matrix::rowSums(mat)
    
    sample_id <- paste(donor, grp, sep="_")
    pb_counts_list[[sample_id]] <- pb
    
    pb_coldata <- rbind(pb_coldata, data.frame(
      sample = sample_id,
      donor = donor,
      group = grp,
      n_cells = length(cells_grp),
      stringsAsFactors = FALSE
    ))
  }
}

if (length(pb_counts_list) < 4) stop("Too few pseudobulk samples created; check min_cells_per_group and labels.")

pb_counts <- do.call(cbind, pb_counts_list)
colnames(pb_counts) <- names(pb_counts_list)

write.table(pb_coldata,
            file=file.path(opt$outdir, "tables", "pseudobulk_samples.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# -------------------------
# Helper: run DESeq2 with donor blocking for a 2-group contrast
# design: ~ donor + group
# -------------------------
run_deseq2_contrast <- function(pb_counts, coldata, g1, g2) {
  sub <- coldata[coldata$group %in% c(g1, g2), , drop=FALSE]
  sub$group <- droplevels(factor(sub$group, levels=c(g1, g2)))
  sub$donor <- factor(sub$donor)

  donors_ok <- names(which(tapply(sub$group, sub$donor, function(x) all(c(g1,g2) %in% x))))
  sub <- sub[sub$donor %in% donors_ok, , drop=FALSE]
  sub$donor <- droplevels(sub$donor)

  if (nrow(sub) < 4) stop(sprintf("Not enough samples for contrast %s vs %s after donor pairing.", g1, g2))

  m <- pb_counts[, sub$sample, drop=FALSE]

  keep <- rowSums(m) >= 5
  m <- m[keep, , drop = FALSE]

  m_mat <- as.matrix(m)
  if (anyDuplicated(rownames(m_mat)) > 0) {
    m_mat <- rowsum(m_mat, group = rownames(m_mat), reorder = FALSE)
  }

  keep_samples <- colSums(m_mat) > 0
  m_mat <- m_mat[, keep_samples, drop = FALSE]
  sub <- sub[keep_samples, , drop = FALSE]

  storage.mode(m_mat) <- "integer"

  dds <- DESeqDataSetFromMatrix(
    countData = m_mat,
    colData = sub,
    design = ~ donor + group
  )

  dds <- DESeq(dds, sfType = "poscounts", quiet = TRUE)

  res <- results(dds, contrast=c("group", g1, g2))
  res <- as.data.frame(res)
  res$gene <- rownames(res)
  res <- res[order(res$padj, res$pvalue), ]
  
  list(
  dds = dds,
  res = res,
  sub = sub,
  universe = rownames(m_mat)
)

}


# -------------------------
# TOST equivalence test on log2FC using DESeq2 SE
# -------------------------
tost_from_res <- function(df, delta, alpha, min_baseMean) {
  out <- df
  
  out$baseMean_ok <- !is.na(out$baseMean) & out$baseMean >= min_baseMean
  out$se_ok <- !is.na(out$lfcSE) & out$lfcSE > 0
  ok <- out$baseMean_ok & out$se_ok & !is.na(out$log2FoldChange)
  
  out$p_tost1 <- NA_real_
  out$p_tost2 <- NA_real_
  out$p_equiv <- NA_real_
  out$equiv <- FALSE
  
  z1 <- (out$log2FoldChange[ok] + delta) / out$lfcSE[ok]
  z2 <- (out$log2FoldChange[ok] - delta) / out$lfcSE[ok]
  
  p1 <- 1 - pnorm(z1)
  p2 <- pnorm(z2)
  
  p_equiv <- pmax(p1, p2)
  out$p_tost1[ok] <- p1
  out$p_tost2[ok] <- p2
  out$p_equiv[ok] <- p_equiv
  out$padj_equiv <- NA_real_
  out$padj_equiv[ok] <- p.adjust(out$p_equiv[ok], method="BH")
  
  out$equiv[ok] <- (out$padj_equiv[ok] < alpha) &
    (abs(out$log2FoldChange[ok]) < delta)
  
  out
}

# -------------------------
# GSEA/ORA utilities (explicit; error if deps missing)
# -------------------------
require_pkg <- function(x) {
  if (!requireNamespace(x, quietly = TRUE)) {
    stop(sprintf("Required package '%s' is not installed in this environment.", x))
  }
}

read_gmt_df <- function(gmt_path) {
  lines <- readLines(gmt_path, warn = FALSE)
  parts <- strsplit(lines, "\t", fixed = TRUE)

  # GMT format: NAME \t DESCRIPTION \t GENE1 \t GENE2 ...
  df_list <- lapply(parts, function(x) {
    if (length(x) < 3) return(NULL)
    gs <- x[[1]]
    genes <- unique(x[3:length(x)])
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) == 0) return(NULL)
    data.frame(gs_name = gs, gene_symbol = genes, stringsAsFactors = FALSE)
  })

  df <- do.call(rbind, df_list)
  if (is.null(df) || nrow(df) == 0) stop("GMT produced empty pathway table: ", gmt_path)
  df
}


run_gsea_fgsea <- function(ranks, pathways_df, out_tsv, out_png, label, top_n = 25) {
  require_pkg("fgsea")
  pw <- split(pathways_df$gene_symbol, pathways_df$gs_name)

  res <- fgsea::fgsea(pathways = pw, stats = ranks, eps = 0)
  res <- as.data.frame(res)
  res <- res[order(res$padj, res$pval), , drop = FALSE]

  is_list_col <- vapply(res, is.list, logical(1))
  if (any(is_list_col)) {
    for (nm in names(res)[is_list_col]) {
      res[[nm]] <- vapply(res[[nm]], function(x) if (is.null(x)) "" else paste(x, collapse=","), character(1))
    }
  }

  write.table(res, out_tsv, sep="\t", row.names=FALSE, quote=FALSE)

  df <- res[!is.na(res$padj), , drop = FALSE]
  df <- head(df, top_n)
  if (nrow(df) > 0) {

  # make long names readable
 df$pathway_plot <- make.unique(stringr::str_wrap(stringr::str_trunc(df$pathway, 80), width = 40))
 df$pathway_plot <- factor(df$pathway_plot, levels = rev(df$pathway_plot))

  p <- ggplot(df, aes(x = NES, y = pathway_plot, size = size, color = -log10(padj))) +
    geom_point(alpha = 0.9) +
    labs(title = paste0("GSEA (fgsea) ", label), x = "NES", y = "Pathway") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 8))

  h <- max(6, 0.25 * nrow(df))   # autoscale height
  ggsave(out_png, p, width = 10, height = h, dpi = 300)
  }

  invisible(res)
}



run_ora_enricher <- function(genes, universe, pathways_df, out_tsv, out_png, label, top_n = 25) {
  require_pkg("clusterProfiler")
  term2gene <- unique(pathways_df[, c("gs_name","gene_symbol")])
  colnames(term2gene) <- c("term","gene")

  genes <- unique(genes)
  genes <- genes[!is.na(genes) & nzchar(genes)]

  universe <- unique(universe)
  universe <- universe[!is.na(universe) & nzchar(universe)]

  genes <- intersect(genes, universe)

  if (length(genes) < 5) {
    write.table(data.frame(), out_tsv, sep="\t", row.names=FALSE, quote=FALSE)
    return(invisible(NULL))
  }

  enr <- clusterProfiler::enricher(genes, universe = universe, TERM2GENE = term2gene)
  df <- as.data.frame(enr)
  df <- df[order(df$p.adjust, df$pvalue), , drop = FALSE]
  write.table(df, out_tsv, sep="\t", row.names=FALSE, quote=FALSE)

  d2 <- df[!is.na(df$p.adjust), , drop = FALSE]
  d2 <- head(d2, top_n)
  if (nrow(d2) > 0) {
  d2$Description_plot <- stringr::str_wrap(stringr::str_trunc(d2$Description, 80), width = 40)
  d2$Description_plot <- make.unique(d2$Description_plot)
  d2$Description_plot <- factor(d2$Description_plot, levels = rev(d2$Description_plot))


p <- ggplot(d2, aes(x = Count, y = Description_plot, size = Count, color = -log10(p.adjust))) +
  geom_point(alpha = 0.9) +
  labs(title = paste0("ORA (enricher) ", label), x = "Gene count", y = "Pathway") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(size = 8))

 h <- max(6, 0.25 * nrow(d2))
 ggsave(out_png, p, width=10, height=h, dpi=300)
  }

  invisible(df)
}



pathways_h <- read_gmt_df(opt$hallmark_gmt)
pathways_c7 <- read_gmt_df(opt$c7_gmt)


# -------------------------
# Run contrasts
# -------------------------
all_results <- list()
marker_sets <- list()
conserved_sets <- list()
universe_sets <- list()

for (cc in deg_contrasts) {
  g1 <- cc[[1]]
  g2 <- cc[[2]]
  tag <- paste0(g1, "_vs_", g2)
  
  message(sprintf("DESeq2: %s", tag))
  fit <- run_deseq2_contrast(pb_counts, pb_coldata, g1, g2)
  
  res <- fit$res

  universe_sets[[tag]] <- strip_dup_suffix(fit$universe)

  # -------------------------
  # STRICT MARKER CALLS (use opt thresholds)
  # -------------------------
  res$is_marker <-
    !is.na(res$padj) &
    res$padj < opt$marker_padj &
    abs(res$log2FoldChange) > opt$marker_lfc
  
  # -------------------------
  # TOST equivalence
  # -------------------------
  res2 <- tost_from_res(
    res,
    delta = opt$equiv_delta,
    alpha = opt$equiv_alpha,
    min_baseMean = opt$min_baseMean
  )
  
  # carry marker flag into res2
  res2$is_marker <- res$is_marker[match(res2$gene, rownames(res))]
  
  # CONSERVED (use opt thresholds)
  res2$is_conserved <- res2$equiv & res2$baseMean >= opt$min_baseMean
  
  # -------------------------
  # WRITE TABLES
  # -------------------------
  write.table(
    res2,
    file = file.path(opt$outdir, "tables", paste0("deseq2_", tag, ".tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  
  markers <- subset(res2, is_marker)
  marker_sets[[tag]] <- unique(strip_dup_suffix(markers$gene))
  markers<-markers[order(as.numeric(markers$padj), decreasing=FALSE),]
  write.table(
    markers[, c("gene","log2FoldChange","lfcSE","stat","pvalue","padj","baseMean")],
    file = file.path(opt$outdir, "tables", paste0("markers_", tag, ".tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  
  conserved <- subset(res2, is_conserved)
  conserved_sets[[tag]] <- unique(strip_dup_suffix(conserved$gene))
  conserved<-conserved[order(as.numeric(conserved$padj_equiv), decreasing=FALSE),]
  write.table(
    conserved[, c("gene","log2FoldChange","lfcSE","p_equiv","padj_equiv","baseMean")],
    file = file.path(opt$outdir, "tables", paste0("conserved_", tag, ".tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  
  # -------------------------
  # GSEA (rank by log2FC) + ORA (gene lists)
  # -------------------------
  
  rnk_df <- res2[!is.na(res2$log2FoldChange) & !is.na(res2$gene), c("gene", "log2FoldChange")]
  # if rowsum collapsed duplicates upstream, this is already clean; keep robust anyway
  rnk_df$gene <- strip_dup_suffix(rnk_df$gene)
  rnk <- tapply(rnk_df$log2FoldChange, rnk_df$gene, function(x) x[1])
  rnk <- sort(rnk, decreasing = TRUE)
  
  # ranks already computed as `rnk`
# markers / conserved already computed

# ---- Hallmark ----
run_gsea_fgsea(
  ranks = rnk,
  pathways_df = pathways_h,
  out_tsv = file.path(opt$outdir, "tables", paste0("gsea_hallmark_", tag, ".tsv")),
  out_png = file.path(opt$outdir, "plots",  paste0("gsea_hallmark_", tag, "_bubble.png")),
  label="Hallmark"
)

run_ora_enricher(
  genes = strip_dup_suffix(markers$gene),
  universe = strip_dup_suffix(fit$universe),
  pathways_df = pathways_h,
  out_tsv = file.path(opt$outdir, "tables", paste0("ora_markers_hallmark_", tag, ".tsv")),
  out_png = file.path(opt$outdir, "plots",  paste0("ora_markers_hallmark_", tag, "_bubble.png")),
  label="Hallmark"
)

run_ora_enricher(
  genes = strip_dup_suffix(conserved$gene),
  universe = strip_dup_suffix(fit$universe),
  pathways_df = pathways_h,
  out_tsv = file.path(opt$outdir, "tables", paste0("ora_conserved_hallmark_", tag, ".tsv")),
  out_png = file.path(opt$outdir, "plots",  paste0("ora_conserved_hallmark_", tag, "_bubble.png")),
  label="Hallmark"
)

# ---- C7 (immunologic signatures) ----
run_gsea_fgsea(
  ranks = rnk,
  pathways_df = pathways_c7,
  out_tsv = file.path(opt$outdir, "tables", paste0("gsea_c7_", tag, ".tsv")),
  out_png = file.path(opt$outdir, "plots",  paste0("gsea_c7_", tag, "_bubble.png")),
  label="C7"
)

run_ora_enricher(
  genes = strip_dup_suffix(markers$gene),
  universe = strip_dup_suffix(fit$universe),
  pathways_df = pathways_c7,
  out_tsv = file.path(opt$outdir, "tables", paste0("ora_markers_c7_", tag, ".tsv")),
  out_png = file.path(opt$outdir, "plots",  paste0("ora_markers_c7_", tag, "_bubble.png")),
  label="C7"
)

run_ora_enricher(
  genes = strip_dup_suffix(conserved$gene),
  universe = strip_dup_suffix(fit$universe),
  pathways_df = pathways_c7,
  out_tsv = file.path(opt$outdir, "tables", paste0("ora_conserved_c7_", tag, ".tsv")),
  out_png = file.path(opt$outdir, "plots",  paste0("ora_conserved_c7_", tag, "_bubble.png")),
  label="C7"
)



  
  # -------------------------
  # VOLCANO
  # -------------------------
  dfp <- res2
  dfp$flag <- "other"
  dfp$flag[dfp$is_marker] <- "Marker"
  dfp$flag[dfp$is_conserved] <- "Conserved"
  
  p_vol <- ggplot(dfp, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(data=subset(dfp, flag=="other"),
               color="grey70", alpha=0.25, size=0.6) +
    geom_point(data=subset(dfp, flag=="Marker"),
               color="red", alpha=0.9, size=1.6) +
    geom_point(data=subset(dfp, flag=="Conserved"),
               color="blue", alpha=0.9, size=1.6) +
    geom_vline(xintercept=c(-opt$marker_lfc, opt$marker_lfc), linetype="dashed") +
    labs(title=paste0(tag, " volcano"),
         x="log2 fold change",
         y="-log10(padj)")
  
  ggsave(
    file.path(opt$outdir, "plots", paste0("volcano_", tag, ".png")),
    p_vol, width=7, height=4, dpi=300
  )
  
  # -------------------------
  # LFC histogram (ALL genes)
  # -------------------------
  p_lfc_all <- ggplot(dfp, aes(x = log2FoldChange)) +
    geom_histogram(bins = 80) +
    labs(title = paste0(tag, " log2FC (all genes)"), x = "log2 fold change", y = "genes")
  
  ggsave(
    file.path(opt$outdir, "plots", paste0("lfc_hist_", tag, ".png")),
    p_lfc_all, width = 7, height = 4, dpi = 300
  )
  
  # -------------------------
  # LFC histogram (Markers vs Conserved)
  # -------------------------
  df_lfc <- dfp[dfp$flag %in% c("Marker", "Conserved"), , drop = FALSE]
  
  if (nrow(df_lfc) > 0) {
    p_lfc_flag <- ggplot(df_lfc, aes(x = log2FoldChange)) +
      geom_histogram(bins = 60) +
      facet_wrap(~flag, scales = "free_y", ncol = 1) +
      labs(title = paste0(tag, " log2FC (Marker vs Conserved)"),
           x = "log2 fold change", y = "genes")
    
    ggsave(
      file.path(opt$outdir, "plots", paste0("lfc_hist_flags_", tag, ".png")),
      p_lfc_flag, width = 7, height = 6, dpi = 300
    )
  }
  
  # -------------------------
  # SUMMARY
  # -------------------------
  all_results[[tag]] <- data.frame(
    contrast = tag,
    n_samples = nrow(fit$sub),
    n_donors = length(unique(fit$sub$donor)),
    n_markers = nrow(markers),
    n_conserved = nrow(conserved),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, all_results)
write.table(
  summary_df,
  file = file.path(opt$outdir, "tables", "summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# -------------------------
# Cross-contrast intersections / uniqueness
# -------------------------
conserved_all <- Reduce(intersect, conserved_sets)
markers_any <- Reduce(union, marker_sets)
markers_all <- Reduce(intersect, marker_sets)
universe_all <- Reduce(intersect, lapply(universe_sets, unique))
universe_all <- universe_all[!is.na(universe_all) & nzchar(universe_all)]


# -------------------------
# Cross-contrast ANY universes (for union gene sets)
# -------------------------
universe_any <- Reduce(union, lapply(universe_sets, unique))
universe_any <- universe_any[!is.na(universe_any) & nzchar(universe_any)]

# optional: also do conserved_any if you want it
conserved_any <- Reduce(union, conserved_sets)
# conserved_any already stripped above; if not, do strip_dup_suffix(conserved_any)



markers_unique <- lapply(names(marker_sets), function(nm) {
  setdiff(marker_sets[[nm]], Reduce(union, marker_sets[names(marker_sets) != nm]))
})
names(markers_unique) <- names(marker_sets)

overlap_marker_conserved <- intersect(
  Reduce(union, conserved_sets),
  Reduce(union, marker_sets)
)

# -------------------------
# ORA on cross-contrast gene sets
# -------------------------

run_ora_enricher(
  genes = strip_dup_suffix(markers_all),
  universe = universe_all,
  pathways_df = pathways_h,
  out_tsv = file.path(opt$outdir, "tables", "ora_markers_all_contrasts.tsv"),
  out_png = file.path(opt$outdir, "plots", "ora_markers_all_contrasts_bubble.png"),
  label="Hallmark"
)

run_ora_enricher(
  genes = strip_dup_suffix(conserved_all),
  universe = universe_all,
  pathways_df = pathways_h,
  out_tsv = file.path(opt$outdir, "tables", "ora_conserved_all_contrasts.tsv"),
  out_png = file.path(opt$outdir, "plots", "ora_conserved_all_contrasts_bubble.png"),
  label = "Hallmark"
)


# C7 (new)
run_ora_enricher(
  genes = strip_dup_suffix(markers_all),
  universe = universe_all,
  pathways_df = pathways_c7,
  out_tsv = file.path(opt$outdir, "tables", "ora_markers_all_contrasts_c7.tsv"),
  out_png = file.path(opt$outdir, "plots",  "ora_markers_all_contrasts_c7_bubble.png"),
  label="C7"
)

run_ora_enricher(
  genes = strip_dup_suffix(conserved_all),
  universe = universe_all,
  pathways_df = pathways_c7,
  out_tsv = file.path(opt$outdir, "tables", "ora_conserved_all_contrasts_c7.tsv"),
  out_png = file.path(opt$outdir, "plots",  "ora_conserved_all_contrasts_c7_bubble.png"),
  label="C7"
)

# -------------------------
# ORA on cross-contrast ANY gene sets (union sets)
# -------------------------
run_ora_enricher(
  genes = strip_dup_suffix(markers_any),
  universe = universe_any,
  pathways_df = pathways_h,
  out_tsv = file.path(opt$outdir, "tables", "ora_markers_any_contrast.tsv"),
  out_png = file.path(opt$outdir, "plots", "ora_markers_any_contrast_bubble.png"),
  label="Hallmark"
)


run_ora_enricher(
  genes = strip_dup_suffix(conserved_any),
  universe = universe_any,
  pathways_df = pathways_h,
  out_tsv = file.path(opt$outdir, "tables", "ora_conserved_any_contrast.tsv"),
  out_png = file.path(opt$outdir, "plots", "ora_conserved_any_contrast_bubble.png"),
  label="Hallmark"
)


run_ora_enricher(
  genes = strip_dup_suffix(markers_any),
  universe = universe_any,
  pathways_df = pathways_c7,
  out_tsv = file.path(opt$outdir, "tables", "ora_markers_any_contrast_c7.tsv"),
  out_png = file.path(opt$outdir, "plots",  "ora_markers_any_contrast_c7_bubble.png"),
  label="C7"
)

run_ora_enricher(
  genes = strip_dup_suffix(conserved_any),
  universe = universe_any,
  pathways_df = pathways_c7,
  out_tsv = file.path(opt$outdir, "tables", "ora_conserved_any_contrast_c7.tsv"),
  out_png = file.path(opt$outdir, "plots",  "ora_conserved_any_contrast_c7_bubble.png"),
  label="C7"
)


# -------------------------
# Write outputs
# -------------------------
write.table(
  data.frame(gene = conserved_all),
  file = file.path(opt$outdir, "tables", "conserved_all_contrasts.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  data.frame(gene = conserved_any),
  file = file.path(opt$outdir, "tables", "conserved_any_contrasts.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  data.frame(gene = markers_any),
  file = file.path(opt$outdir, "tables", "markers_any_contrast.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  data.frame(gene = markers_all),
  file = file.path(opt$outdir, "tables", "markers_all_contrasts.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

for (nm in names(markers_unique)) {
  write.table(
    data.frame(gene = markers_unique[[nm]]),
    file = file.path(opt$outdir, "tables", paste0("markers_unique_", nm, ".tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

write.table(
  data.frame(gene = overlap_marker_conserved),
  file = file.path(opt$outdir, "tables", "overlap_marker_and_conserved.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

message("Done.")
