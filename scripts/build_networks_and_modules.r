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
})
# ADD THE MARKER GENES AND CONSERVED GENES TO THE GEPHI NODE OUTPUT!

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

  # Defaults you requested
  make_option("--min_cells_per_donor_group", type="integer", default=200,
              help="Minimum cells required in a donor×celltype-set to construct metacells."),
  make_option("--metacell_size", type="integer", default=20,
              help="Metacell size (random pooling)."),
  make_option("--seed", type="integer", default=1,
              help="Random seed for metacells."),
  make_option("--gene_detect_frac", type="double", default=0.05,
              help="Keep genes detected (>0) in at least this fraction of metacells."),
  make_option("--hvg_n", type="integer", default=3000,
              help="Number of HVGs used for correlation/network."),
  make_option("--cor_method", type="character", default="spearman",
              help="Correlation method: spearman|pearson"),
  make_option("--positive_only", type="logical", default=TRUE,
              help="Keep only positive correlations."),
  make_option("--top_k", type="integer", default=25,
              help="Top-k edges per gene (sparsification)."),
  make_option("--min_abs_cor", type="double", default=0.25,
              help="Minimum absolute correlation to keep an edge."),
  make_option("--consensus_min_donors", type="integer", default=3,
              help="Edge must appear in >= N donors to be retained."),
  make_option("--require_same_sign", type="logical", default=TRUE,
              help="Require consistent sign across donors."),
  make_option("--leiden_resolution", type="double", default=1.0,
              help="Leiden resolution parameter."),
  make_option("--deg_tables_dir", type="character", default="",
             help="Directory containing conserved_*.tsv and markers_*.tsv tables to annotate nodes."),
  make_option("--hallmark_gmt", type="character", default="",
            help="Local Hallmark GMT file path (required for module ORA)."),
  make_option("--c7_gmt", type="character", default="",
            help="Local C7 GMT file path (required for module ORA).")
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
out_plot <- file.path(opt$outdir, "plots")
  
read_rds_list <- function(x) {
  xs <- strsplit(x, ",", fixed=TRUE)[[1]]
  xs <- trimws(xs)
  xs <- xs[nzchar(xs)]
  xs
}

infer_donor_names <- function(paths) {
  # take basename without extension as donor label fallback
  nm <- basename(paths)
  nm <- sub("\\.rds$", "", nm, ignore.case=TRUE)
  nm
}

# Make random metacells: pool N cells -> average expression per gene
make_metacells <- function(mat_genes_x_cells, size=20, seed=1) {
  set.seed(seed)
  nc <- ncol(mat_genes_x_cells)
  idx <- sample(seq_len(nc))
  groups <- split(idx, ceiling(seq_along(idx)/size))
  # average per group
  m <- sapply(groups, function(cols) {
    if (length(cols) == 1) return(mat_genes_x_cells[, cols, drop=FALSE][,1])
    Matrix::rowMeans(mat_genes_x_cells[, cols, drop=FALSE])
  })
  # ensure matrix form: genes x metacells
  if (!is.matrix(m)) m <- as.matrix(m)
  colnames(m) <- paste0("mc", seq_len(ncol(m)))
  m
}

# Keep genes detected in >= frac of metacells
filter_genes_detect <- function(expr_genes_x_mc, frac=0.05) {
  det <- rowMeans(expr_genes_x_mc > 0)
  keep <- det >= frac
  expr_genes_x_mc[keep, , drop=FALSE]
}

# Sparsify correlation matrix into edge list using top_k per gene + threshold
sparsify_topk <- function(C, top_k=25, min_abs_cor=0.25, positive_only=TRUE) {
  genes <- rownames(C)
  stopifnot(length(genes) == nrow(C))
  edges <- vector("list", length(genes))
  names(edges) <- genes

  for (i in seq_along(genes)) {
    g <- genes[i]
    v <- C[, i]
    v[i] <- NA_real_

    if (positive_only) {
      v[v <= 0] <- NA_real_
      score <- v
    } else {
      score <- abs(v)
    }

    keep_idx <- which(!is.na(score) & abs(v) >= min_abs_cor)
    if (length(keep_idx) == 0) next

    # top-k by score
    ord <- keep_idx[order(score[keep_idx], decreasing=TRUE)]
    if (length(ord) > top_k) ord <- ord[seq_len(top_k)]

    edges[[i]] <- data.frame(
      from = g,
      to = genes[ord],
      cor = v[ord],
      stringsAsFactors = FALSE
    )
  }

  df <- bind_rows(edges)
  if (nrow(df) == 0) return(df)

  # undirected canonicalization
  df <- df %>%
    mutate(a = pmin(from, to), b = pmax(from, to)) %>%
    select(a, b, cor) %>%
    group_by(a, b) %>%
    summarise(cor = max(cor, na.rm=TRUE), .groups="drop")

  df %>% rename(from=a, to=b)
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


run_leiden <- function(g, resolution=1.0, seed=1) {
  set.seed(seed)
  if (!"cluster_leiden" %in% getNamespaceExports("igraph")) stop("igraph::cluster_leiden not available.")

  cl <- tryCatch(
    igraph::cluster_leiden(
      g,
      resolution_parameter = resolution,
      weights = E(g)$weight,
      objective_function = "modularity"
    ),
    error = function(e) {
      message("Leiden objective_function not supported; falling back to default.")
      igraph::cluster_leiden(g, resolution_parameter = resolution, weights = E(g)$weight)
    }
  )
  membership(cl)
}

parse_term <- function(x) {
  x <- gsub("_", " ", x)
  x <- tolower(x)
  x <- gsub("\\b([a-z])", "\\U\\1", x, perl = TRUE)
  x
}



plot_hist_png <- function(x, file, main, xlab="weight", breaks=60,
                          add_quantiles=TRUE) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(invisible(NULL))

  dir.create(dirname(file), recursive=TRUE, showWarnings=FALSE)
  png(file, width=1200, height=900, res=120)

  hist(
    x,
    breaks = breaks,
    main   = main,
    xlab   = xlab,
    col    = "grey80",
    border = "grey40"
  )

  abline(v = median(x), lwd = 2)

  if (isTRUE(add_quantiles)) {
    qs <- quantile(x, c(0.05, 0.95))
    abline(v = qs, lty = 2)
  }

  dev.off()
}



plot_support_png <- function(s, file,
                             main="Edge support",
                             xlab="support (# donors)",
                             ylab="edge count") {
  s <- s[is.finite(s)]
  if (length(s) < 1) return(invisible(NULL))

  dir.create(dirname(file), recursive=TRUE, showWarnings=FALSE)
  png(file, width=1200, height=900, res=120)

  tab <- table(factor(s, levels=sort(unique(s))))
  barplot(
    tab,
    main  = main,
    xlab  = xlab,
    ylab  = ylab,
    col   = "grey80",
    border= "grey40"
  )

  dev.off()
}

plot_module_sizes <- function(modules, out_png, title) {
  tab <- table(modules$module)
  df <- data.frame(
    module = as.integer(names(tab)),
    size = as.integer(tab)
  )

  png(out_png, width = 900, height = 600)
  par(mar = c(5, 5, 4, 2))
  barplot(
    df$size,
    names.arg = df$module,
    las = 2,
    col = "grey80",
    border = "grey40",
    ylab = "Number of genes",
    xlab = "Module ID",
    main = title
  )
  abline(h = median(df$size), lwd = 2)
  dev.off()

  return(df)
}

plot_enrichment_dotplot <- function(enrich_tsv, out_png,
                                    top_n_per_module = 5,
                                    padj_max = 0.05,
                                    title = "Enrichment by module (ORA)") {
  df <- readr::read_tsv(enrich_tsv, show_col_types = FALSE)

  # clusterProfiler uses: ID, Description, pvalue, p.adjust, Count, GeneRatio...
  # Your TERM names might be in Description or ID. Prefer Description.
  term_col <- if ("Description" %in% names(df)) "Description" else "ID"

  df <- df %>%
    mutate(
      term_raw = .data[[term_col]],
      term = parse_term(term_raw),
      module = as.factor(module),
      neglog10_padj = -log10(p.adjust + 1e-300)
    ) %>%
    filter(is.finite(neglog10_padj))

  # Optional filtering
  if (!is.null(padj_max)) df <- df %>% filter(p.adjust <= padj_max)

  # Keep top N terms per module (by p.adjust)
  df <- df %>%
    group_by(module) %>%
    arrange(p.adjust, desc(Count)) %>%
    slice_head(n = top_n_per_module) %>%
    ungroup()

  if (nrow(df) == 0) {
    message("No enrichment hits to plot for: ", enrich_tsv)
    return(invisible(NULL))
  }

  # Order y-axis by overall significance
  term_order <- df %>%
    group_by(term) %>%
    summarise(best = max(neglog10_padj), .groups = "drop") %>%
    arrange(best) %>%
    pull(term)

  df$term <- factor(df$term, levels = term_order)

  df$term_plot <- stringr::str_wrap(stringr::str_trunc(as.character(df$term), 80), width = 40)
  df$term_plot <- factor(df$term_plot, levels = unique(df$term_plot))


  p <- ggplot(df, aes(x = module, y = term_plot)) +
    geom_point(aes(size = Count, color = neglog10_padj), alpha = 0.9) +
    scale_color_continuous(name = "-log10(adj p)") +
    scale_size_continuous(name = "Overlap (Count)") +
    labs(x = "Module", y = NULL, title=title) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5)
    )

  ggsave(out_png, p, width = 12, height = max(6, 0.30 * nrow(df)), dpi = 150)
  out_png
}


ora_gmt <- function(genes, universe, pathways_df) {
  term2gene <- pathways_df %>% dplyr::transmute(term = gs_name, gene = gene_symbol)

  gs_genes <- unique(term2gene$gene)

  genes2 <- intersect(unique(genes), gs_genes)
  universe2 <- intersect(unique(universe), gs_genes)

  if (length(genes2) < 10 || length(universe2) < 50) return(NULL)

  suppressMessages(suppressWarnings(
    enricher(
      gene = genes2,
      universe = universe2,
      TERM2GENE = term2gene,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500
    )
  ))
}

# Largest connected component
keep_largest_cc <- function(g) {
  comp <- igraph::components(g)
  giant <- which.max(comp$csize)
  igraph::induced_subgraph(g, vids = which(comp$membership == giant))
}

#Graph stats
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


# -------------------------
# Load definitions
# -------------------------
if (!file.exists(opt$celltype_sets)) stop(paste("Missing:", opt$celltype_sets))
source(opt$celltype_sets)

# Accept either naming convention:
# - new: celltype_sets
# - your current: network_sets
if (exists("celltype_sets")) {
  ct_sets <- celltype_sets
} else if (exists("network_sets")) {
  ct_sets <- network_sets
} else {
  stop("celltype_sets file must define 'celltype_sets' or 'network_sets' (named list).")
}

if (!is.list(ct_sets) || is.null(names(ct_sets))) {
  stop("Celltype set object must be a named list (celltype_sets or network_sets).")
}


markers_pbmc <- NULL
if (nzchar(opt$markers)) {
  if (!file.exists(opt$markers)) stop(paste("Missing:", opt$markers))
  source(opt$markers) # should define markers_pbmc
  if (exists("markers_pbmc")) {
    if (!is.list(markers_pbmc) || is.null(names(markers_pbmc))) stop("'markers_pbmc' must be a named list.")
  }
}

deg_markers <- list()
deg_conserved <- list()

if (nzchar(opt$deg_tables_dir)) {
  if (!dir.exists(opt$deg_tables_dir)) stop(paste("Missing deg_tables_dir:", opt$deg_tables_dir))

  marker_files <- list.files(opt$deg_tables_dir, pattern="^markers_.*\\.tsv$", full.names=TRUE)
  conserved_files <- list.files(opt$deg_tables_dir, pattern="^conserved_.*\\.tsv$", full.names=TRUE)

  read_gene_col <- function(f) {
    x <- suppressWarnings(read.delim(f, sep="\t", header=TRUE, stringsAsFactors=FALSE))
    # pick a likely gene column
    cand <- c("gene", "gene_symbol", "symbol", "Gene", "GeneSymbol")
    col <- cand[cand %in% colnames(x)][1]
    if (is.na(col)) stop(paste("No gene column found in", f, "columns:", paste(colnames(x), collapse=", ")))
    unique(as.character(x[[col]]))
  }

  for (f in marker_files) {
    nm <- tools::file_path_sans_ext(basename(f))  # e.g. markers_T_like_vs_B_like
    deg_markers[[nm]] <- read_gene_col(f)
  }
  for (f in conserved_files) {
    nm <- tools::file_path_sans_ext(basename(f))
    deg_conserved[[nm]] <- read_gene_col(f)
  }
}


# -------------------------
# Load Seurat objects (annotated)
# -------------------------
seurat_paths <- read_rds_list(opt$seurat)
if (length(seurat_paths) < 2) stop("--seurat must include >=2 donor RDS files (comma-separated).")
for (p in seurat_paths) if (!file.exists(p)) stop(paste("Missing RDS:", p))

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

# Heuristic for label column: prefer cell_type_cluster_majority then cell_type_pred then seurat_clusters
pick_label_col <- function(obj) {
  md <- obj@meta.data
  for (c in c("cell_type_cluster_majority", "cell_type_pred", "cell_type", "celltype", "seurat_clusters")) {
    if (c %in% colnames(md)) return(c)
  }
  stop("No usable cell type label column found in metadata (expected one of: cell_type_cluster_majority, cell_type_pred, ...).")
}

# Decide label column
if (exists("deg_label_col") && !is.null(deg_label_col) && nzchar(deg_label_col)) {
  label_col <- deg_label_col
} else {
  label_col <- pick_label_col(objs[[1]])
}

# Validate
if (!label_col %in% colnames(objs[[1]]@meta.data)) {
  stop(paste0("Label column '", label_col, "' not found in Seurat metadata."))
}

message("Using label column: ", label_col)



for (d in names(objs)) {
  if (!label_col %in% colnames(objs[[d]]@meta.data)) {
    stop(paste0("Label column '", label_col, "' missing in donor '", d, "'."))
  }
}

#Load gene lists

if (!nzchar(opt$hallmark_gmt) || !file.exists(opt$hallmark_gmt)) stop("Missing/invalid --hallmark_gmt")
if (!nzchar(opt$c7_gmt) || !file.exists(opt$c7_gmt)) stop("Missing/invalid --c7_gmt")

pathways_h  <- read_gmt_df(opt$hallmark_gmt)  # columns: gs_name, gene_symbol
pathways_c7 <- read_gmt_df(opt$c7_gmt)


# -------------------------
# Main: for each celltype_set build per-donor + consensus
# -------------------------
all_outputs <- list()

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

  # Collect donor edge lists
  donor_edges <- list()
  donor_genes <- list()
  donor_genes_universe <- list()

  for (d in donors) {
    obj <- objs[[d]]
    md <- obj@meta.data

    # Which fine labels map into this set?
    keep_labels <- ct_sets[[set_name]]
    if (is.null(keep_labels) || length(keep_labels) == 0) next

    cells <- rownames(md)[md[[label_col]] %in% keep_labels]
    if (length(cells) < opt$min_cells_per_donor_group) {
      message("Skip donor ", d, " (", length(cells), " cells < min ", opt$min_cells_per_donor_group, ")")
      next
    }

    sub <- subset(obj, cells=cells)

  # Input expression
  expr_log <- GetAssayData(sub, slot="data")  # RNA@data (log1p normalized)
    if (ncol(expr_log) < opt$metacell_size) {
    message("Skip donor ", d, " (cells < metacell_size).")
    next
  }

  # Choose aggregation scale
  if (opt$metacell_input == "linear_then_log") {
    expr_to_pool <- expm1(expr_log)   # back to linear normalized
  } else if (opt$metacell_input == "log") {
    expr_to_pool <- expr_log                  # keep log scale (demo-fast)
  } else {
    stop("--metacell_input must be 'log' or 'linear_then_log'")
  }

  # Build metacells on chosen scale
  mc_pool <- make_metacells(expr_to_pool, size=opt$metacell_size, seed=opt$seed)

  # Ensure correlations use log scale
  mc <- if (opt$metacell_input == "linear_then_log") log1p(mc_pool) else mc_pool

  # Gene filter by detection fraction 
  mc_f <- filter_genes_detect(mc, frac=opt$gene_detect_frac)

    # HVGs on the subset (to pick ~3000)
    # (We run FindVariableFeatures on the subset to get a stable HVG set for this donor×set.)
    sub <- FindVariableFeatures(sub, selection.method="vst", nfeatures=max(opt$hvg_n, 2000), verbose=FALSE)
    hvgs <- VariableFeatures(sub)
    hvgs <- intersect(hvgs, rownames(mc_f))
    if (length(hvgs) < 500) {
      message("Skip donor ", d, " (too few HVGs after filtering: ", length(hvgs), ")")
      next
    }
    hvgs <- hvgs[seq_len(min(opt$hvg_n, length(hvgs)))]

    mc_hvg <- mc_f[hvgs, , drop=FALSE]

    # Correlation across genes (genes x metacells => cor(t()))
    # cor() expects observations in rows, variables in columns => transpose
    C <- suppressWarnings(cor(t(mc_hvg), method=opt$cor_method, use="pairwise.complete.obs"))
    diag(C) <- NA_real_
    rownames(C) <- hvgs
    colnames(C) <- hvgs

    # Sparsify
    ed <- sparsify_topk(C, top_k=opt$top_k, min_abs_cor=opt$min_abs_cor, positive_only=opt$positive_only)
    if (nrow(ed) == 0) {
      message("No edges retained for donor ", d, " / set ", set_name)
      next
    }

    ed$donor <- d

  plot_hist_png(
    ed$cor,
    file.path(plots_donor, sprintf("edge_weight_hist_sparsified_%s.png", d)),
    main=sprintf("%s / %s: donor %s edge weights (sparsified)", set_name, label_col, d),
    xlab="cor"
  )


    # genes present in this donor BEFORE LCC (for universe counting)
    donor_genes_universe[[d]] <- hvgs



    g_donor <- igraph::graph_from_data_frame(
      ed %>% dplyr::select(from, to, weight = cor),
      directed = FALSE
    )

    g_donor <- igraph::simplify(g_donor, remove.multiple = TRUE, remove.loops = TRUE)
    g_donor <- keep_largest_cc(g_donor)

    # Re-extract edges after LCC
    ed2 <- igraph::as_data_frame(g_donor, what = "edges")

    # Standardize name: weight -> cor
    if ("weight" %in% names(ed2)) names(ed2)[names(ed2) == "weight"] <- "cor"

    ed2$donor <- d

    plot_hist_png(
      ed2$cor,
      file.path(plots_donor, sprintf("edge_weight_hist_postLCC_%s.png", d)),
      main=sprintf("%s: donor %s edge weights (post donor-LCC)", set_name, d),
      xlab="cor"
    )


    # IMPORTANT: store the post-LCC edges for consensus
    donor_edges[[d]] <- ed2
    donor_genes[[d]] <- sort(unique(c(ed2$from, ed2$to)))

    # Write per-donor edges (Gephi-ready)
    write.table(
      ed2 %>% dplyr::transmute(
      Source = from,
      Target = to,
      Weight = cor,
      donor = donor
    ),
      file = file.path(out_per, paste0("edges_", d, ".tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE
    )

    saveRDS(g_donor, file = file.path(out_per, paste0("graph_", d, ".rds")))

    write.table(
      graph_stats(g_donor),
      file = file.path(out_per, paste0("stats_", d, ".tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE
    )


  }

  if (length(donor_edges) < opt$consensus_min_donors) {
    message("Not enough donors for consensus in set ", set_name, " (have ", length(donor_edges), ")")
    next
  }

  gene_counts <- table(unlist(donor_genes_universe))
  universe <- names(gene_counts[gene_counts >= opt$consensus_min_donors])



  # Stack edges and compute support + sign consistency + summary weight
  all_ed <- bind_rows(donor_edges) %>%
    filter(from %in% universe, to %in% universe)

  all_ed <- all_ed %>%
    mutate(sign = ifelse(cor >= 0, 1L, -1L)) %>%
    group_by(from, to) %>%
    summarise(
      support = n(),
      sign_sum = sum(sign),
      median_cor = median(cor),
      .groups="drop"
    )
  
  plot_hist_png(
    all_ed$median_cor,
    file.path(plots_set, "consensus_hist_median_cor_ALL_edges.png"),
    main=sprintf("%s: ALL edges median_cor (before support/sign filters)", set_name),
    xlab="median_cor"
  )

  plot_support_png(
    all_ed$support,
    file.path(plots_set, "consensus_bar_support_ALL_edges.png"),
    main=sprintf("%s: ALL edges support (before filters)", set_name)
  )


  cons <- all_ed %>%
    filter(support >= opt$consensus_min_donors)

  if (opt$require_same_sign) {
    cons <- cons %>% filter(abs(sign_sum) == support)
  }

  # Keep only desired sign if positive_only
  if (opt$positive_only) {
    cons <- cons %>% filter(median_cor > 0)
  }

  if (nrow(cons) == 0) {
    message("Consensus empty for set ", set_name)
    next
  }



  cons0 <- cons %>% mutate(weight = median_cor) %>% select(from, to, weight, support)
  cons <- cons0
  edges_consensus_preLCC <- nrow(cons0)
  



  g <- igraph::graph_from_data_frame(
    cons %>% select(from, to, weight),
    directed = FALSE
  ) 

  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "max")
  g <- keep_largest_cc(g)

  cons <- igraph::as_data_frame(g, what = "edges")
  if ("weight" %in% names(cons)) names(cons)[names(cons) == "weight"] <- "cor"

  # add support back (because igraph edges currently only have weight/cor)
  cons <- cons %>%
    mutate(a = pmin(from, to), b = pmax(from, to)) %>%
    select(a, b, cor) %>%
    left_join(
      cons0 %>%
        mutate(a = pmin(from, to), b = pmax(from, to)) %>%
        select(a, b, support),
      by = c("a", "b")
    ) %>%
    select(from = a, to = b, cor, support)

  
  cons <- cons %>%
  mutate(
    support_frac = support / length(donor_edges),
    weight_consensus = cor * support_frac
  )
  edges_consensus_postLCC <- nrow(cons)
  plot_hist_png(
    cons$cor,
    file.path(plots_set, "consensus_hist_cor_FINAL_postLCC.png"),
    main=sprintf("%s: FINAL consensus edge weights (post consensus-LCC)", set_name),
    xlab="cor"
  )

  plot_support_png(
    cons$support,
    file.path(plots_set, "consensus_bar_support_FINAL_postLCC.png"),
    main=sprintf("%s: FINAL consensus support (post consensus-LCC)", set_name)
  )


  saveRDS(g, file = file.path(out_set, "graph_consensus.rds"))

  write.table(
    graph_stats(g),
    file = file.path(out_set, "stats_consensus.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )


  memb <- run_leiden(g, resolution=opt$leiden_resolution, seed=opt$seed)
  modules <- data.frame(
    gene = names(memb),
    module = as.integer(memb),
    stringsAsFactors = FALSE
  )

  mod_sizes <- plot_module_sizes(
  modules,
  out_png = file.path(plots_set, "module_size_distribution.png"),
  title = paste0(set_name, ": module size distribution (consensus)")
)


  write.table(
    mod_sizes,
    file = file.path(out_set, "module_sizes.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  # Node metrics + annotations
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

  # Aggregate flags
if (length(deg_markers) > 0) {
  
  if ("markers_any_contrast" %in% names(deg_markers)) {
    nodes$is_deg_marker_any <- nodes$id %in% deg_markers[["markers_any_contrast"]]
  } else {
    nodes$is_deg_marker_any <- nodes$id %in% unique(unlist(deg_markers))
  }

  for (nm in names(deg_markers)) {
    nodes[[paste0("is_", nm)]] <- nodes$id %in% deg_markers[[nm]]
  }
}

if (length(deg_conserved) > 0) {

  if ("conserved_any_contrast" %in% names(deg_conserved)) {
  nodes$is_deg_conserved_any <- nodes$id %in% deg_conserved[["conserved_any_contrast"]]

} else {
  nodes$is_deg_conserved_any <- nodes$id %in% unique(unlist(deg_conserved))
}

  for (nm in names(deg_conserved)) {
    nodes[[paste0("is_", nm)]] <- nodes$id %in% deg_conserved[[nm]]
  }
}


  # Annotate marker membership (if markers_pbmc provided)
  if (!is.null(markers_pbmc)) {
    # For each marker set, add boolean column: is_marker_<name>
    for (mk in names(markers_pbmc)) {
      genes_mk <- unique(markers_pbmc[[mk]])
      coln <- paste0("is_marker_", mk)
      nodes[[coln]] <- nodes$id %in% genes_mk
    }
  }

  # Write Gephi tables
  # Nodes
  write.csv(
    nodes,
    file = file.path(out_set, "nodes.csv"),
    row.names = FALSE,
    quote = TRUE
  )

  # Edges
  # Edges (Gephi wants Source/Target)
write.csv(
  cons %>%
    transmute(
      Source = from,
      Target = to,
      Weight = cor,
      support,
      support_frac,
      weight_consensus
    ),
  file = file.path(out_set, "edges.csv"),
  row.names = FALSE,
  quote = TRUE
)


  # Modules
  write.csv(
    modules,
    file = file.path(out_set, "modules.csv"),
    row.names = FALSE,
    quote = TRUE
  )

  # Enrichment per module (Hallmark ORA)
  enr_list <- list()
  for (m in sort(unique(modules$module))) {
    genes_m <- modules$gene[modules$module == m]
    if (length(genes_m) < 10) next
    enr <- tryCatch(
      ora_gmt(genes = genes_m, universe = universe, pathways_df = pathways_h),
      error = function(e) NULL
    )
    if (is.null(enr) || nrow(as.data.frame(enr)) == 0) next
    df_enr <- as.data.frame(enr) %>%
      mutate(module = m, n_genes = length(genes_m)) %>%
      select(module, n_genes, everything())
    enr_list[[as.character(m)]] <- df_enr
  }

  # Hallmark ORA (independent)
enr_tbl <- bind_rows(enr_list)
if (nrow(enr_tbl) > 0) {
  write.table(enr_tbl, file=file.path(out_set, "enrichment_hallmark_ora.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)

  plot_enrichment_dotplot(
    enrich_tsv = file.path(out_set, "enrichment_hallmark_ora.tsv"),
    out_png    = file.path(plots_set, "enrichment_hallmark_dotplot.png"),
    top_n_per_module = 6,
    padj_max = 0.10,
    title = "Hallmark enrichment by module (ORA)"
  )
}

# C7 ORA (independent)
enr_list_c7 <- list()
for (m in sort(unique(modules$module))) {
  genes_m <- modules$gene[modules$module == m]
  if (length(genes_m) < 10) next
  enr <- tryCatch(
    ora_gmt(genes = genes_m, universe = universe, pathways_df = pathways_c7),
    error = function(e) NULL
  )
  if (is.null(enr) || nrow(as.data.frame(enr)) == 0) next
  df_enr <- as.data.frame(enr) %>%
    mutate(module = m, n_genes = length(genes_m)) %>%
    select(module, n_genes, everything())
  enr_list_c7[[as.character(m)]] <- df_enr
}

enr_tbl_c7 <- bind_rows(enr_list_c7)
if (nrow(enr_tbl_c7) > 0) {
  write.table(enr_tbl_c7, file=file.path(out_set, "enrichment_c7_ora.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)

  plot_enrichment_dotplot(
    enrich_tsv = file.path(out_set, "enrichment_c7_ora.tsv"),
    out_png    = file.path(plots_set, "enrichment_c7_dotplot.png"),
    top_n_per_module = 6,
    padj_max = 0.10,
    title = "C7 enrichment by module (ORA)"
  )
}

  


  # Optional: ORA with marker sets as "gene sets" (module -> enrichment in marker sets)
  if (!is.null(markers_pbmc)) {
    marker_sets <- stack(markers_pbmc) %>% transmute(term=ind, gene=values)
    ms_enr_list <- list()
    for (m in sort(unique(modules$module))) {
      genes_m <- modules$gene[modules$module == m]
      if (length(genes_m) < 10) next
      enr2 <- tryCatch(
        enricher(
          gene = unique(genes_m),
          universe = unique(universe),
          TERM2GENE = marker_sets,
          pAdjustMethod = "BH"
        ),
        error=function(e) NULL
      )
      if (is.null(enr2) || nrow(as.data.frame(enr2)) == 0) next
      df2 <- as.data.frame(enr2) %>%
        mutate(module=m, n_genes=length(genes_m)) %>%
        select(module, n_genes, everything())
      ms_enr_list[[as.character(m)]] <- df2
    }
    ms_tbl <- bind_rows(ms_enr_list)
    if (nrow(ms_tbl) > 0) {
      write.table(ms_tbl,
        file=file.path(out_set, "enrichment_marker_sets_ora.tsv"),
        sep="\t", row.names=FALSE, quote=FALSE
      )
      plot_enrichment_dotplot(
        enrich_tsv = file.path(out_set, "enrichment_marker_sets_ora.tsv"),
        out_png    = file.path(plots_set, "enrichment_marker_sets_dotplot.png"),
        top_n_per_module = 6,
        padj_max = 0.10,
        title = "Marker-set enrichment by module (ORA)"
     )

    }
  }

  # Summary
  summary_tbl <- data.frame(
  celltype_set = set_name,
  donors_used = length(donor_edges),
  genes_universe = length(universe),
  edges_consensus_preLCC = edges_consensus_preLCC,
  edges_consensus_postLCC = edges_consensus_postLCC,
  modules = length(unique(modules$module)),
  stringsAsFactors = FALSE
)

  write.table(summary_tbl,
    file=file.path(opt$outdir, "tables", paste0("network_summary_", set_name, ".tsv")),
    sep="\t", row.names=FALSE, quote=FALSE
  )
}

