#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(igraph)
  library(clusterProfiler)
  library(ggplot2)
  library(readr)
  library(muumi)
  library(org.Hs.eg.db)
  })

# -------------------------
# CLI
# -------------------------
option_list <- list(
  make_option("--seurat", type="character", help="Comma-separated list of annotated Seurat RDS files (one per donor)."),
  make_option("--outdir", type="character", help="Output directory."),
  make_option("--celltype_sets", type="character", default="scripts/celltype_sets.R",
              help="Path to scripts/celltype_sets.R defining network_sets (or celltype_sets) (named list)."),
  make_option("--markers", type="character", default="",
              help="Optional: path to scripts/markers_pbmc.R defining markers_pbmc (named list)."),
  make_option("--donor_names", type="character", default="",
              help="Optional: comma-separated donor names matching --seurat order. If empty, uses filenames."),
  make_option("--metacell_input", type="character", default="log",
              help="Metacell aggregation scale: 'log' (RNA@data) or 'linear_then_log' (expm1 before pooling, log1p after)."),

  # Front-end
  make_option("--min_cells_per_donor_group", type="integer", default=200),
  make_option("--metacell_size", type="integer", default=20),
  make_option("--seed", type="integer", default=1),
  make_option("--gene_detect_frac", type="double", default=0.05),
  make_option("--hvg_n", type="integer", default=3000),

  # Muumi MINET params (passed into muumi::calculate_correlation_matrix via muumi::get_ranked_consensus_matrix)
  make_option("--muumi_iMethods", type="character", default="clr",
              help="Comma-separated MINET inference algorithms (e.g. clr,aracne,mrnet,mrnetb)."),
  make_option("--muumi_iEst", type="character", default="spearman",
              help="Comma-separated MINET estimators (e.g. spearman,pearson,kendall,mi.empirical,mi.shrink,mi.sg...)."),
  make_option("--muumi_iDisc", type="character", default="none",
              help="Comma-separated discretization methods (e.g. none,equalfreq,equalwidth,globalequalwidth)."),
  make_option("--muumi_ncores", type="integer", default=2),

  # Per-donor sparsification on Muumi weights (NO fixed global K)
  make_option("--muumi_top_k_gene", type="integer", default=100,
              help="Per-gene top-k neighbors on Muumi weights."),
  make_option("--muumi_quantile", type="double", default=NA,
              help="Optional global cutoff per donor: keep edges with weight >= quantile (e.g. 0.995). Applied before per-gene top-k."),

  # Cross-donor replicate-stable filter
  make_option("--consensus_min_donors", type="integer", default=2,
              help="Edge must appear in >= N donors (for N=4, use 2; also export >=3 overlay)."),

  # Modules (Muumi)
  make_option("--muumi_module_method", type="character", default="louvain",
              help="Muumi get_modules method: louvain|walktrap|spinglass|greedy"),

  # Existing options kept (annotations/enrichment)
  make_option("--deg_tables_dir", type="character", default="",
              help="Directory containing conserved_*.tsv and markers_*.tsv tables to annotate nodes."),
  make_option("--hallmark_gmt", type="character", default="",
              help="Local Hallmark GMT file path (required for module ORA)."),
  make_option("--c7_gmt", type="character", default="",
              help="Local C7 GMT file path (required for module ORA)."),
  make_option("--network_set_names", type="character", default="",
            help="Optional comma-separated subset of network set names to run (e.g. CD4,B_cells).")
)

opt <- parse_args(OptionParser(option_list = option_list))


# -------------------------
# Helpers
# -------------------------
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "per_donor"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "consensus"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "plots"), recursive = TRUE, showWarnings = FALSE)


read_rds_list <- function(x) {
  xs <- strsplit(x, ",", fixed=TRUE)[[1]]
  xs <- trimws(xs)
  xs <- xs[nzchar(xs)]
  xs
}

infer_donor_names <- function(paths) {
  nm <- basename(paths)
  nm <- sub("\\.rds$", "", nm, ignore.case=TRUE)
  nm
}

split_csv <- function(x) {
  xs <- trimws(strsplit(x, ",", fixed=TRUE)[[1]])
  xs[nzchar(xs)]
}

# Make random metacells: pool N cells -> average expression per gene
make_metacells <- function(mat_genes_x_cells, size=20, seed=1) {
  set.seed(seed)
  nc <- ncol(mat_genes_x_cells)
  idx <- sample(seq_len(nc))
  groups <- split(idx, ceiling(seq_along(idx)/size))
  m <- sapply(groups, function(cols) {
    if (length(cols) == 1) return(mat_genes_x_cells[, cols, drop=FALSE][,1])
    Matrix::rowMeans(mat_genes_x_cells[, cols, drop=FALSE])
  })
  if (!is.matrix(m)) m <- as.matrix(m)
  colnames(m) <- paste0("mc", seq_len(ncol(m)))
  m
}

filter_genes_detect <- function(expr_genes_x_mc, frac=0.05) {
  det <- rowMeans(expr_genes_x_mc > 0)
  expr_genes_x_mc[det >= frac, , drop=FALSE]
}

read_gmt_df <- function(gmt_path) {
  lines <- readLines(gmt_path, warn = FALSE)
  parts <- strsplit(lines, "\t", fixed = TRUE)

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

parse_term <- function(x) {
  x <- gsub("_", " ", x)
  x <- tolower(x)
  gsub("\\b([a-z])", "\\U\\1", x, perl = TRUE)
}

plot_hist_png <- function(x, file, main, xlab="weight", breaks=60, add_quantiles=TRUE) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(invisible(NULL))
  dir.create(dirname(file), recursive=TRUE, showWarnings=FALSE)
  png(file, width=1200, height=900, res=120)
  hist(x, breaks=breaks, main=main, xlab=xlab, col="grey80", border="grey40")
  abline(v=median(x), lwd=2)
  if (isTRUE(add_quantiles)) {
    qs <- quantile(x, c(0.05, 0.95))
    abline(v=qs, lty=2)
  }
  dev.off()
}

plot_support_png <- function(s, file, main="Edge support", xlab="support (# donors)", ylab="edge count") {
  s <- s[is.finite(s)]
  if (length(s) < 1) return(invisible(NULL))
  dir.create(dirname(file), recursive=TRUE, showWarnings=FALSE)
  png(file, width=1200, height=900, res=120)
  tab <- table(factor(s, levels=sort(unique(s))))
  barplot(tab, main=main, xlab=xlab, ylab=ylab, col="grey80", border="grey40")
  dev.off()
}

plot_module_sizes <- function(modules, out_png, title) {
  tab <- table(modules$module)
  df <- data.frame(module=as.integer(names(tab)), size=as.integer(tab))
  png(out_png, width=900, height=600)
  par(mar=c(5,5,4,2))
  barplot(df$size, names.arg=df$module, las=2, col="grey80", border="grey40",
          ylab="Number of genes", xlab="Module ID", main=title)
  abline(h=median(df$size), lwd=2)
  dev.off()
  df
}

plot_enrichment_dotplot <- function(enrich_tsv, out_png, top_n_per_module=5, padj_max=0.05,
                                    title="Enrichment by module (ORA)") {
  df <- readr::read_tsv(enrich_tsv, show_col_types = FALSE)
  term_col <- if ("Description" %in% names(df)) "Description" else "ID"

  df <- df %>%
    mutate(term_raw = .data[[term_col]],
           term = parse_term(term_raw),
           module = as.factor(module),
           neglog10_padj = -log10(p.adjust + 1e-300)) %>%
    filter(is.finite(neglog10_padj))

  if (!is.null(padj_max)) df <- df %>% filter(p.adjust <= padj_max)

  df <- df %>% group_by(module) %>% arrange(p.adjust, desc(Count)) %>% slice_head(n=top_n_per_module) %>% ungroup()
  if (nrow(df) == 0) return(invisible(NULL))

  term_order <- df %>% group_by(term) %>% summarise(best=max(neglog10_padj), .groups="drop") %>% arrange(best) %>% pull(term)
  df$term <- factor(df$term, levels=term_order)
  df$term_plot <- stringr::str_wrap(stringr::str_trunc(as.character(df$term), 80), width=40)
  df$term_plot <- factor(df$term_plot, levels=unique(df$term_plot))

  p <- ggplot(df, aes(x=module, y=term_plot)) +
    geom_point(aes(size=Count, color=neglog10_padj), alpha=0.9) +
    scale_color_continuous(name="-log10(adj p)") +
    scale_size_continuous(name="Overlap (Count)") +
    labs(x="Module", y=NULL, title=title) +
    theme_bw(base_size=12) +
    theme(panel.grid.major.y=element_blank(), axis.text.y=element_text(size=10),
          plot.title=element_text(hjust=0.5))

  ggsave(out_png, p, width=12, height=max(6, 0.30*nrow(df)), dpi=150)
  out_png
}

ora_gmt <- function(genes, universe, pathways_df) {
  term2gene <- pathways_df %>% dplyr::transmute(term=gs_name, gene=gene_symbol)
  gs_genes <- unique(term2gene$gene)
  genes2 <- intersect(unique(genes), gs_genes)
  universe2 <- intersect(unique(universe), gs_genes)
  if (length(genes2) < 10 || length(universe2) < 50) return(NULL)

  suppressMessages(suppressWarnings(
    enricher(gene=genes2, universe=universe2, TERM2GENE=term2gene,
             pAdjustMethod="BH", minGSSize=10, maxGSSize=500)
  ))
}

keep_largest_cc <- function(g) {
  comp <- igraph::components(g)
  giant <- which.max(comp$csize)
  igraph::induced_subgraph(g, vids = which(comp$membership == giant))
}

graph_stats <- function(g) {
  data.frame(
    nodes = igraph::vcount(g),
    edges = igraph::ecount(g),
    avg_degree = mean(igraph::degree(g)),
    density = igraph::edge_density(g),
    components = igraph::components(g)$no,
    stringsAsFactors = FALSE
  )
}

# Sparsify Muumi weight matrix into an edge list using:
# (1) optional donor-specific quantile cutoff, then (2) per-gene top-k cap.
sparsify_muumi <- function(W, k_gene = 25, quantile_cutoff = NA_real_) {
  stopifnot(nrow(W) == ncol(W))
  genes <- rownames(W)
  if (is.null(genes)) stop("Muumi matrix must have rownames/colnames.")
  W <- as.matrix(W)
  diag(W) <- NA_real_

  # Optional global cutoff
  if (is.finite(quantile_cutoff)) {
    off <- W[upper.tri(W)]
    off <- off[is.finite(off)]
    if (length(off) > 20) {
      thr <- as.numeric(stats::quantile(off, probs = quantile_cutoff, na.rm = TRUE))
      W[W < thr] <- NA_real_
    }
  }

  # Per-gene top-k
  from <- character(0); to <- character(0); weight <- numeric(0)

  for (i in seq_along(genes)) {
    v <- W[, i]
    v[i] <- NA_real_
    keep <- which(is.finite(v))
    if (!length(keep)) next
    ord <- keep[order(v[keep], decreasing = TRUE)]
    if (length(ord) > k_gene) ord <- ord[seq_len(k_gene)]

    from   <- c(from, rep(genes[i], length(ord)))
    to     <- c(to, genes[ord])
    weight <- c(weight, v[ord])
  }

  if (!length(from)) {
    return(data.frame(from=character(), to=character(), weight=numeric(), stringsAsFactors=FALSE))
  }

  df <- data.frame(from=from, to=to, weight=weight, stringsAsFactors=FALSE)

  # Undirected canonicalization
  a <- ifelse(df$from < df$to, df$from, df$to)
  b <- ifelse(df$from < df$to, df$to, df$from)
  df$a <- a
  df$b <- b

  # Deduplicate: keep max weight per (a,b)
  key <- paste(df$a, df$b, sep="__")
  o <- order(df$weight, decreasing=TRUE)
  df <- df[o, ]
  keep <- !duplicated(key[o])

  out <- df[keep, c("a","b","weight")]
  colnames(out) <- c("from","to","weight")
  rownames(out) <- NULL
  out
}

# -------------------------
# Load definitions
# -------------------------
if (!file.exists(opt$celltype_sets)) stop("Missing: ", opt$celltype_sets)
source(opt$celltype_sets)

if (exists("celltype_sets")) {
  ct_sets <- celltype_sets
} else if (exists("network_sets")) {
  ct_sets <- network_sets
} else {
  stop("celltype_sets file must define 'celltype_sets' or 'network_sets' (named list).")
}
if (!is.list(ct_sets) || is.null(names(ct_sets))) stop("Celltype set object must be a named list.")

if (nzchar(opt$network_set_names)) {
  keep_sets <- trimws(strsplit(opt$network_set_names, ",", fixed = TRUE)[[1]])
  keep_sets <- keep_sets[nzchar(keep_sets)]
  bad <- setdiff(keep_sets, names(ct_sets))
  if (length(bad) > 0) {
    stop("Unknown network set(s): ", paste(bad, collapse = ", "))
  }
  ct_sets <- ct_sets[keep_sets]
}

markers_pbmc <- NULL
if (nzchar(opt$markers)) {
  if (!file.exists(opt$markers)) stop("Missing: ", opt$markers)
  source(opt$markers)
  if (exists("markers_pbmc")) {
    if (!is.list(markers_pbmc) || is.null(names(markers_pbmc))) stop("'markers_pbmc' must be a named list.")
  }
}

deg_markers <- list()
deg_conserved <- list()
if (nzchar(opt$deg_tables_dir)) {
  if (!dir.exists(opt$deg_tables_dir)) stop("Missing deg_tables_dir: ", opt$deg_tables_dir)
  marker_files <- list.files(opt$deg_tables_dir, pattern="^markers_.*\\.tsv$", full.names=TRUE)
  conserved_files <- list.files(opt$deg_tables_dir, pattern="^conserved_.*\\.tsv$", full.names=TRUE)

  read_gene_col <- function(f) {
    x <- suppressWarnings(read.delim(f, sep="\t", header=TRUE, stringsAsFactors=FALSE))
    cand <- c("gene", "gene_symbol", "symbol", "Gene", "GeneSymbol")
    col <- cand[cand %in% colnames(x)][1]
    if (is.na(col)) stop("No gene column found in ", f)
    unique(as.character(x[[col]]))
  }

  for (f in marker_files) {
    nm <- tools::file_path_sans_ext(basename(f))
    deg_markers[[nm]] <- read_gene_col(f)
  }
  for (f in conserved_files) {
    nm <- tools::file_path_sans_ext(basename(f))
    deg_conserved[[nm]] <- read_gene_col(f)
  }
}

# -------------------------
# Load Seurat objects
# -------------------------
seurat_paths <- read_rds_list(opt$seurat)
if (length(seurat_paths) < 2) stop("--seurat must include >=2 donor RDS files.")
for (p in seurat_paths) if (!file.exists(p)) stop("Missing RDS: ", p)

objs <- lapply(seurat_paths, readRDS)

donors <- if (nzchar(opt$donor_names)) {
  dn <- trimws(strsplit(opt$donor_names, ",", fixed=TRUE)[[1]])
  dn <- dn[nzchar(dn)]
  if (length(dn) != length(objs)) stop("--donor_names must match number of --seurat files.")
  dn
} else {
  infer_donor_names(seurat_paths)
}
names(objs) <- donors

pick_label_col <- function(obj) {
  md <- obj@meta.data
  for (c in c("cell_type_cluster_majority", "cell_type_pred", "cell_type", "celltype", "seurat_clusters")) {
    if (c %in% colnames(md)) return(c)
  }
  stop("No usable cell type label column found in metadata.")
}

if (exists("deg_label_col") && !is.null(deg_label_col) && nzchar(deg_label_col)) {
  label_col <- deg_label_col
} else {
  label_col <- pick_label_col(objs[[1]])
}
if (!label_col %in% colnames(objs[[1]]@meta.data)) stop("Label column not found: ", label_col)
for (d in names(objs)) {
  if (!label_col %in% colnames(objs[[d]]@meta.data)) stop("Label column missing in donor ", d, ": ", label_col)
}
message("Using label column: ", label_col)

# Load gene sets for ORA (kept from your pipeline)
if (!nzchar(opt$hallmark_gmt) || !file.exists(opt$hallmark_gmt)) stop("Missing/invalid --hallmark_gmt")
if (!nzchar(opt$c7_gmt) || !file.exists(opt$c7_gmt)) stop("Missing/invalid --c7_gmt")
pathways_h  <- read_gmt_df(opt$hallmark_gmt)
pathways_c7 <- read_gmt_df(opt$c7_gmt)

# Muumi parameter parsing
muumi_iMethods <- split_csv(opt$muumi_iMethods)
muumi_iEst     <- split_csv(opt$muumi_iEst)
muumi_iDisc    <- split_csv(opt$muumi_iDisc)

# -------------------------
# Main
# -------------------------
for (set_name in names(ct_sets)) {
  message("=== Celltype set: ", set_name, " ===")

  out_set <- file.path(opt$outdir, "consensus", set_name)
  out_per <- file.path(opt$outdir, "per_donor", set_name)
  plots_set   <- file.path(opt$outdir, "plots", set_name)
  plots_donor <- file.path(plots_set, "per_donor")

  dir.create(out_set, recursive=TRUE, showWarnings=FALSE)
  dir.create(out_per, recursive=TRUE, showWarnings=FALSE)
  dir.create(plots_set, recursive=TRUE, showWarnings=FALSE)
  dir.create(plots_donor, recursive=TRUE, showWarnings=FALSE)

  donor_edges <- list()
  donor_genes_universe <- list()

  for (d in donors) {
    obj <- objs[[d]]
    md <- obj@meta.data

    keep_labels <- ct_sets[[set_name]]
    if (is.null(keep_labels) || length(keep_labels) == 0) next

    cells <- rownames(md)[md[[label_col]] %in% keep_labels]
    if (length(cells) < opt$min_cells_per_donor_group) {
      message("Skip donor ", d, " (", length(cells), " cells < ", opt$min_cells_per_donor_group, ")")
      next
    }

    sub <- subset(obj, cells=cells)

    expr_log <- GetAssayData(sub, slot="data")
    if (ncol(expr_log) < opt$metacell_size) {
      message("Skip donor ", d, " (cells < metacell_size).")
      next
    }

    expr_to_pool <- if (opt$metacell_input == "linear_then_log") expm1(expr_log) else expr_log
    mc_pool <- make_metacells(expr_to_pool, size=opt$metacell_size, seed=opt$seed)
    mc <- if (opt$metacell_input == "linear_then_log") log1p(mc_pool) else mc_pool

    mc_f <- filter_genes_detect(mc, frac=opt$gene_detect_frac)

    sub <- FindVariableFeatures(sub, selection.method="vst", nfeatures=max(opt$hvg_n, 2000), verbose=FALSE)
    hvgs <- intersect(VariableFeatures(sub), rownames(mc_f))
    if (length(hvgs) < 500) {
      message("Skip donor ", d, " (too few HVGs after filtering: ", length(hvgs), ")")
      next
    }
    hvgs <- hvgs[seq_len(min(opt$hvg_n, length(hvgs)))]
    mc_hvg <- mc_f[hvgs, , drop=FALSE]

    donor_genes_universe[[d]] <- hvgs

    # ============================================================
    # ✅ MUUMI NETWORK INFERENCE (per donor)
    # ============================================================
    # Muumi expects genes as rows; your mc_hvg already is genes x metacells.
    message("Muumi inference for donor ", d, " ...")
    rank_mat <- muumi::get_ranked_consensus_matrix(
      gx_table = mc_hvg,
      iMethods = muumi_iMethods,
      iEst     = muumi_iEst,
      iDisc    = muumi_iDisc,
      ncores   = opt$muumi_ncores,
      ensemble_strategy = "minet",
      mat_weights = "rank"
    )

    # Convert ranks to a weight matrix (higher = stronger) for sparsification
    # rank_mat: 1 is best; 0 means absent.
    # weight = 1/rank so top edges get largest weights.
    W <- rank_mat
    W[W <= 0] <- NA_real_
    W <- 1 / W
    diag(W) <- NA_real_
    rownames(W) <- rownames(rank_mat)
    colnames(W) <- colnames(rank_mat)

    # Sparsify per donor using Muumi weights
    ed <- sparsify_muumi(W, k_gene=opt$muumi_top_k_gene, quantile_cutoff=opt$muumi_quantile)
    if (nrow(ed) == 0) {
      message("No Muumi edges retained for donor ", d, " / set ", set_name)
      next
    }
    ed$donor <- d

    plot_hist_png(
      ed$weight,
      file.path(plots_donor, sprintf("muumi_edge_weight_hist_sparsified_%s.png", d)),
      main=sprintf("%s: donor %s Muumi weights (sparsified)", set_name, d),
      xlab="Muumi weight (1/rank)"
    )

  plot_hist_png(
    -log10(ed$weight + 1e-12),
    file.path(plots_donor, sprintf("muumi_edge_weight_hist_sparsified_log_%s.png", d)),
    main = sprintf("%s: donor %s Muumi weights (sparsified, -log10 scale)", set_name, d),
    xlab = "-log10(Muumi weight)"
  )
    # Build donor graph + LCC (kept)
    g_donor <- igraph::graph_from_data_frame(ed %>% dplyr::select(from, to, weight), directed=FALSE)
    g_donor <- igraph::simplify(g_donor, remove.multiple=TRUE, remove.loops=TRUE)
    g_donor <- keep_largest_cc(g_donor)

    ed2 <- igraph::as_data_frame(g_donor, what="edges")
    if ("weight" %in% names(ed2)) names(ed2)[names(ed2) == "weight"] <- "weight"
    ed2$donor <- d

    plot_hist_png(
      ed2$weight,
      file.path(plots_donor, sprintf("muumi_edge_weight_hist_postLCC_%s.png", d)),
      main=sprintf("%s: donor %s Muumi weights (post donor-LCC)", set_name, d),
      xlab="Muumi weight (1/rank)"
    )

    plot_hist_png(
      -log10(ed2$weight + 1e-12),
      file.path(plots_donor, sprintf("muumi_edge_weight_hist_postLCC_log_%s.png", d)),
      main = sprintf("%s: donor %s Muumi weights (post donor-LCC, -log10 scale)", set_name, d),
      xlab = "-log10(Muumi weight)"
    )

    donor_edges[[d]] <- ed2

    # Export per-donor
    write.table(
      ed2 %>% transmute(Source=from, Target=to, Weight=weight, donor=donor),
      file=file.path(out_per, paste0("edges_muumi_", d, ".tsv")),
      sep="\t", row.names=FALSE, quote=FALSE
    )
    saveRDS(g_donor, file=file.path(out_per, paste0("graph_muumi_", d, ".rds")))
    write.table(graph_stats(g_donor), file=file.path(out_per, paste0("stats_muumi_", d, ".tsv")),
                sep="\t", row.names=FALSE, quote=FALSE)
  }

  if (length(donor_edges) < length(donors)) {
  message("Not enough donor networks for set ", set_name,
          " (have ", length(donor_edges), ", require all ", length(donors), ")")
  next
  }

  # Universe: genes present in >= consensus_min_donors donors (same as your logic)
  gene_counts <- table(unlist(donor_genes_universe))
  universe <- names(gene_counts[gene_counts >= opt$consensus_min_donors])

  # Stack edges and compute support + median weight
  all_ed <- dplyr::bind_rows(donor_edges) %>%
  dplyr::filter(.data$from %in% universe, .data$to %in% universe) %>%
  dplyr::mutate(a = pmin(.data$from, .data$to),
                b = pmax(.data$from, .data$to)) %>%
  dplyr::group_by(.data$a, .data$b) %>%
  dplyr::summarise(
    support = dplyr::n(),
    median_weight = stats::median(.data$weight),
    .groups = "drop"
  ) %>%
  dplyr::transmute(from = .data$a, to = .data$b,
                   support = .data$support, median_weight = .data$median_weight)

  plot_hist_png(
    all_ed$median_weight,
    file.path(plots_set, "muumi_median_weight_ALL_edges.png"),
    main=sprintf("%s: ALL edges median Muumi weight (before support filter)", set_name),
    xlab="median Muumi weight"
  )
  plot_support_png(
    all_ed$support,
    file.path(plots_set, "muumi_support_ALL_edges.png"),
    main=sprintf("%s: ALL edges support (before filter)", set_name)
  )

  # Replicate-stable filter
  cons <- all_ed %>% filter(support >= opt$consensus_min_donors)
  if (nrow(cons) == 0) {
    message("Replicate-stable Muumi network empty for set ", set_name)
    next
  }

  # Build final consensus graph (Muumi replicate-stable network)
  g <- igraph::graph_from_data_frame(
  cons %>% dplyr::transmute(from=.data$from, to=.data$to, weight=.data$median_weight),
  directed=FALSE
)
  g <- igraph::simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb="max")
  g <- keep_largest_cc(g)

  cons2 <- igraph::as_data_frame(g, what="edges")
  if ("weight" %in% names(cons2)) names(cons2)[names(cons2) == "weight"] <- "weight"

  # Add support back
  cons2 <- cons2 %>%
  dplyr::mutate(a = pmin(.data$from, .data$to),
                b = pmax(.data$from, .data$to)) %>%
  dplyr::left_join(
    cons %>%
      dplyr::mutate(a = pmin(.data$from, .data$to),
                    b = pmax(.data$from, .data$to)) %>%
      dplyr::select(.data$a, .data$b, .data$support),
    by = c("a","b")
  ) %>%
  dplyr::transmute(from = .data$a, to = .data$b,
                   weight = .data$weight, support = .data$support)

  cons2 <- cons2 %>%
  mutate(
    support_frac = support / length(donor_edges),
    weight_consensus = weight * support_frac
  )

 # Gephi-friendly positive transformed weights: larger = stronger
wc_log <- -log10(cons2$weight_consensus + 1e-12)
  cons2 <- cons2 %>%
   mutate(
     weight_consensus_log = max(wc_log, na.rm = TRUE) - wc_log + 1
    )

  plot_hist_png(
    cons2$weight,
    file.path(plots_set, "muumi_weight_FINAL_postLCC.png"),
    main=sprintf("%s: FINAL Muumi weights (post consensus-LCC)", set_name),
    xlab="Muumi weight"
  )
  plot_support_png(
    cons2$support,
    file.path(plots_set, "muumi_support_FINAL_postLCC.png"),
    main=sprintf("%s: FINAL support (post consensus-LCC)", set_name)
  )

  # additional plot: log-transformed weights for readability
  plot_hist_png(
  -log10(cons2$weight + 1e-12),
  file.path(plots_set, "muumi_weight_FINAL_postLCC_log10.png"),
  main = sprintf("%s: FINAL Muumi weights (-log10)", set_name),
  xlab = "-log10(Muumi weight)"
)

plot_hist_png(
  cons2$weight_consensus_log,
  file.path(plots_set, "muumi_weight_consensus_log_FINAL_postLCC.png"),
  main = sprintf("%s: FINAL Muumi consensus weights (Gephi log-scale)", set_name),
  xlab = "Gephi-friendly consensus weight"
)

  saveRDS(g, file=file.path(out_set, "graph_muumi.rds"))
  write.table(graph_stats(g), file=file.path(out_set, "stats_muumi.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)

  # ============================================================
  # ✅ MUUMI MODULES (Louvain etc.)
  # ============================================================
  # Muumi get_modules expects an igraph
  message("Muumi modules (", opt$muumi_module_method, ") ...")
  muumi_comms <- muumi::get_modules(iGraph = g, method = opt$muumi_module_method)
  memb <- igraph::membership(muumi_comms)

  modules <- data.frame(
    gene = names(memb),
    module = as.integer(memb),
    stringsAsFactors = FALSE
  )

  mod_sizes <- plot_module_sizes(
    modules,
    out_png = file.path(plots_set, "module_size_distribution_muumi.png"),
    title = paste0(set_name, ": module size distribution (Muumi ", opt$muumi_module_method, ")")
  )
  write.table(mod_sizes, file=file.path(out_set, "module_sizes_muumi.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  write.csv(modules, file=file.path(out_set, "modules_muumi.csv"), row.names=FALSE, quote=TRUE)

  # ============================================================
  # Nodes table (same style as yours)
  # ============================================================
  deg <- degree(g)
  btw <- betweenness(g, normalized=TRUE)
  nodes <- data.frame(
    id = names(deg),
    label = names(deg),
    module = modules$module[match(names(deg), modules$gene)],
    degree = as.numeric(deg),
    betweenness = as.numeric(btw),
    stringsAsFactors = FALSE
  )

  if (length(deg_markers) > 0) {
    if ("markers_any_contrast" %in% names(deg_markers)) {
      nodes$is_deg_marker_any <- nodes$id %in% deg_markers[["markers_any_contrast"]]
    } else {
      nodes$is_deg_marker_any <- nodes$id %in% unique(unlist(deg_markers))
    }
    for (nm in names(deg_markers)) nodes[[paste0("is_", nm)]] <- nodes$id %in% deg_markers[[nm]]
  }

  if (length(deg_conserved) > 0) {
    if ("conserved_any_contrast" %in% names(deg_conserved)) {
      nodes$is_deg_conserved_any <- nodes$id %in% deg_conserved[["conserved_any_contrast"]]
    } else {
      nodes$is_deg_conserved_any <- nodes$id %in% unique(unlist(deg_conserved))
    }
    for (nm in names(deg_conserved)) nodes[[paste0("is_", nm)]] <- nodes$id %in% deg_conserved[[nm]]
  }

  if (!is.null(markers_pbmc)) {
    for (mk in names(markers_pbmc)) {
      coln <- paste0("is_marker_", mk)
      nodes[[coln]] <- nodes$id %in% unique(markers_pbmc[[mk]])
    }
  }

  write.csv(nodes, file=file.path(out_set, "nodes_muumi.csv"), row.names=FALSE, quote=TRUE)

  # Edges (Gephi)
  write.csv(
  cons2 %>% transmute(
    Source = from,
    Target = to,
    Weight = weight_consensus_log,
    weight_raw = weight,
    support,
    support_frac,
    weight_consensus
  ),
  file = file.path(out_set, "edges_muumi.csv"),
  row.names = FALSE,
  quote = TRUE
)

  # ============================================================
  # Enrichment (your GMT ORA retained)
  # ============================================================
  enr_list <- list()
  for (m in sort(unique(modules$module))) {
    genes_m <- modules$gene[modules$module == m]
    if (length(genes_m) < 10) next
    enr <- tryCatch(ora_gmt(genes=genes_m, universe=universe, pathways_df=pathways_h), error=function(e) NULL)
    if (is.null(enr) || nrow(as.data.frame(enr)) == 0) next
    df_enr <- as.data.frame(enr) %>% mutate(module=m, n_genes=length(genes_m)) %>% dplyr::select(module, n_genes, everything())
    enr_list[[as.character(m)]] <- df_enr
  }
  enr_tbl <- bind_rows(enr_list)
  if (nrow(enr_tbl) > 0) {
    write.table(enr_tbl, file=file.path(out_set, "enrichment_hallmark_ora_muumi.tsv"),
                sep="\t", row.names=FALSE, quote=FALSE)
    plot_enrichment_dotplot(
      enrich_tsv=file.path(out_set, "enrichment_hallmark_ora_muumi.tsv"),
      out_png=file.path(plots_set, "enrichment_hallmark_dotplot_muumi.png"),
      top_n_per_module=6, padj_max=0.10,
      title="Hallmark enrichment by module (ORA)"
    )
  }

  enr_list_c7 <- list()
  for (m in sort(unique(modules$module))) {
    genes_m <- modules$gene[modules$module == m]
    if (length(genes_m) < 10) next
    enr <- tryCatch(ora_gmt(genes=genes_m, universe=universe, pathways_df=pathways_c7), error=function(e) NULL)
    if (is.null(enr) || nrow(as.data.frame(enr)) == 0) next
    df_enr <- as.data.frame(enr) %>% mutate(module=m, n_genes=length(genes_m)) %>% dplyr::select(module, n_genes, everything())
    enr_list_c7[[as.character(m)]] <- df_enr
  }
  enr_tbl_c7 <- bind_rows(enr_list_c7)
  if (nrow(enr_tbl_c7) > 0) {
    write.table(enr_tbl_c7, file=file.path(out_set, "enrichment_c7_ora_muumi.tsv"),
                sep="\t", row.names=FALSE, quote=FALSE)
    plot_enrichment_dotplot(
      enrich_tsv=file.path(out_set, "enrichment_c7_ora_muumi.tsv"),
      out_png=file.path(plots_set, "enrichment_c7_dotplot_muumi.png"),
      top_n_per_module=6, padj_max=0.10,
      title="C7 enrichment by module (ORA)"
    )
  }

# ============================================================
# MUUMI REACTOME
# ============================================================
message("Muumi Reactome enrichment + bubbleplot ...")

# 1) Export Reactome results (CSV)
muumi::get_reactome_from_modules(
  modules = muumi_comms,
  geneID = "SYMBOL",
  pval_cutoff = 0.01,
  outPath = out_set,
  layout = "overall"
)

# 2) Bubble plot
bubble_png <- file.path(plots_set, "muumi_reactome_bubbleplot.png")
png(bubble_png, width = 1400, height = 1400, res = 140)

tryCatch({

  # run Muumi enrichment (returns compareCluster object)
  res <- muumi::get_bubbleplot_from_pathways(muumi_comms, geneID = "SYMBOL")

  # redraw cleaner plot
  p <- clusterProfiler::dotplot(
    res,
    showCategory = 3   # top 3 pathways per module
  ) +
    ggplot2::scale_y_discrete(
      labels = function(x) stringr::str_wrap(x, 40)
    )

  print(p)

}, error = function(e) {
  message("Muumi bubbleplot failed: ", conditionMessage(e))
})

dev.off()

  # Summary table
  summary_tbl <- data.frame(
    celltype_set = set_name,
    donors_used = length(donor_edges),
    genes_universe = length(universe),
    edges_final = igraph::ecount(g),
    nodes_final = igraph::vcount(g),
    modules = length(unique(modules$module)),
    stringsAsFactors = FALSE
  )
  write.table(summary_tbl,
              file=file.path(opt$outdir, "tables", paste0("network_summary_muumi_", set_name, ".tsv")),
              sep="\t", row.names=FALSE, quote=FALSE)
}