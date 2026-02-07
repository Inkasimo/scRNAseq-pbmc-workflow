#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(igraph)
  library(msigdbr)
  library(clusterProfiler)
})

# -------------------------
# CLI
# -------------------------
option_list <- list(
  make_option("--seurat", type="character", help="Comma-separated list of annotated Seurat RDS files (one per donor)."),
  make_option("--outdir", type="character", help="Output directory."),
  make_option("--celltype_sets", type="character", default="scripts/celltype_sets.R",
              help="Path to scripts/celltype_sets.R defining celltype_sets (named list)."),
  make_option("--markers", type="character", default="",
              help="Optional: path to scripts/markers_pbmc.R defining markers_pbmc (named list)."),
  make_option("--donor_names", type="character", default="",
              help="Optional: comma-separated donor names matching --seurat order. If empty, uses filenames."),
  make_option("--metacell_input", type="character", default="log",
            help="Metacell aggregation scale: 'log' (RNA@data) or 'linear_then_log' (expm1 before pooling, log1p after).")

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
  make_option("--done", type="character", default="network.done",
              help="Sentinel filename written into outdir when complete.")
)

opt <- parse_args(OptionParser(option_list = option_list))

# -------------------------
# Helpers
# -------------------------
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "per_donor"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "consensus"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)

touch_done_atomic <- function(path) {
  d <- dirname(path)
  if (!dir.exists(d)) dir.create(d, recursive=TRUE, showWarnings=FALSE)
  tmp <- paste0(path, ".tmp")
  writeLines("ok", con=tmp)
  file.rename(tmp, path)
}

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

# Leiden clustering wrapper (igraph >= 1.3 has cluster_leiden)
run_leiden <- function(g, resolution=1.0, seed=1) {
  set.seed(seed)
  if (!"cluster_leiden" %in% getNamespaceExports("igraph")) {
    stop("igraph::cluster_leiden not available. Update igraph in container or use an alternative Leiden package.")
  }
  cl <- igraph::cluster_leiden(g, resolution_parameter=resolution, weights=E(g)$weight)
  membership(cl)
}

# ORA enrichment for one module gene set against MSigDB Hallmark
ora_hallmark <- function(genes, universe, hallmark_df) {
  # hallmark_df columns: gs_name, gene_symbol
  gs <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
  enricher(
    gene = unique(genes),
    universe = unique(universe),
    TERM2GENE = stack(gs) %>% transmute(term=ind, gene=values),
    pAdjustMethod = "BH"
  )
}

# -------------------------
# Load definitions
# -------------------------
if (!file.exists(opt$celltype_sets)) stop(paste("Missing:", opt$celltype_sets))
source(opt$celltype_sets)  # must define celltype_sets

if (!exists("celltype_sets")) stop("celltype_sets.R must define object 'celltype_sets' (named list).")
if (!is.list(celltype_sets) || is.null(names(celltype_sets))) stop("'celltype_sets' must be a named list.")

markers_pbmc <- NULL
if (nzchar(opt$markers)) {
  if (!file.exists(opt$markers)) stop(paste("Missing:", opt$markers))
  source(opt$markers) # should define markers_pbmc
  if (exists("markers_pbmc")) {
    if (!is.list(markers_pbmc) || is.null(names(markers_pbmc))) stop("'markers_pbmc' must be a named list.")
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

label_col <- pick_label_col(objs[[1]])

# -------------------------
# MSigDB Hallmark (H)
# -------------------------
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol) %>%
  distinct()

# -------------------------
# Main: for each celltype_set build per-donor + consensus
# -------------------------
all_outputs <- list()

for (set_name in names(celltype_sets)) {
  message("=== Celltype set: ", set_name, " ===")

  out_set <- file.path(opt$outdir, "consensus", set_name)
  out_per <- file.path(opt$outdir, "per_donor", set_name)
  dir.create(out_set, recursive=TRUE, showWarnings=FALSE)
  dir.create(out_per, recursive=TRUE, showWarnings=FALSE)

  # Collect donor edge lists
  donor_edges <- list()
  donor_genes <- list()

  for (d in donors) {
    obj <- objs[[d]]
    md <- obj@meta.data

    # Which fine labels map into this set?
    keep_labels <- celltype_sets[[set_name]]
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
    expr_to_pool <- Matrix::expm1(expr_log)   # back to linear normalized
  } else if (opt$metacell_input == "log") {
    expr_to_pool <- expr_log                  # keep log scale (demo-fast)
  } else {
    stop("--metacell_input must be 'log' or 'linear_then_log'")
  }

  # Build metacells on chosen scale
  mc_pool <- make_metacells(expr_to_pool, size=opt$metacell_size, seed=opt$seed)

  # Ensure correlations use log scale
  mc <- if (opt$metacell_input == "linear_then_log") Matrix::log1p(mc_pool) else mc_pool

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
    donor_edges[[d]] <- ed
    donor_genes[[d]] <- sort(unique(c(ed$from, ed$to)))

    # Write per-donor edges (Gephi-ready)
    write.table(
      ed %>% select(from, to, weight=cor, donor),
      file=file.path(out_per, paste0("edges_", d, ".tsv")),
      sep="\t", row.names=FALSE, quote=FALSE
    )
  }

  if (length(donor_edges) < opt$consensus_min_donors) {
    message("Not enough donors for consensus in set ", set_name, " (have ", length(donor_edges), ")")
    next
  }

  # Consensus gene universe: intersection across donors (you suggested this)
  universe <- Reduce(intersect, donor_genes)

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

  cons <- cons %>% mutate(weight = median_cor) %>% select(from, to, weight, support)

  # Build igraph + Leiden modules
  g <- graph_from_data_frame(cons %>% select(from, to, weight), directed=FALSE)
  g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb="max")

  memb <- run_leiden(g, resolution=opt$leiden_resolution, seed=opt$seed)
  modules <- data.frame(
    gene = names(memb),
    module = as.integer(memb),
    stringsAsFactors = FALSE
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

  # Optional: annotate marker membership (if markers_pbmc provided)
  if (!is.null(markers_pbmc)) {
    # For each marker set, add boolean column: is_marker_<name>
    for (mk in names(markers_pbmc)) {
      genes_mk <- unique(markers_pbmc[[mk]])
      coln <- paste0("is_marker_", mk)
      nodes[[coln]] <- nodes$id %in% genes_mk
    }
  }

  # Write Gephi tables
  write.table(nodes, file=file.path(out_set, "nodes.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
  write.table(cons,  file=file.path(out_set, "edges.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
  write.table(modules, file=file.path(out_set, "modules.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

  # Enrichment per module (Hallmark ORA)
  enr_list <- list()
  for (m in sort(unique(modules$module))) {
    genes_m <- modules$gene[modules$module == m]
    if (length(genes_m) < 10) next
    enr <- tryCatch(
      ora_hallmark(genes=genes_m, universe=universe, hallmark_df=hallmark),
      error=function(e) NULL
    )
    if (is.null(enr) || nrow(as.data.frame(enr)) == 0) next
    df_enr <- as.data.frame(enr) %>%
      mutate(module = m, n_genes = length(genes_m)) %>%
      select(module, n_genes, everything())
    enr_list[[as.character(m)]] <- df_enr
  }

  enr_tbl <- bind_rows(enr_list)
  if (nrow(enr_tbl) > 0) {
    write.table(enr_tbl,
      file=file.path(out_set, "enrichment_hallmark_ora.tsv"),
      sep="\t", row.names=FALSE, quote=FALSE
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
    }
  }

  # Summary
  summary_tbl <- data.frame(
    celltype_set = set_name,
    donors_used = length(donor_edges),
    genes_universe = length(universe),
    edges_consensus = nrow(cons),
    modules = length(unique(modules$module)),
    stringsAsFactors = FALSE
  )
  write.table(summary_tbl,
    file=file.path(opt$outdir, "tables", paste0("network_summary_", set_name, ".tsv")),
    sep="\t", row.names=FALSE, quote=FALSE
  )
}

# Sentinel
touch_done_atomic(file.path(opt$outdir, opt$done))
