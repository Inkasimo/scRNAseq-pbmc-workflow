# Canonical PBMC marker genes curated from standard Seurat / 10x Genomics
# PBMC workflows and immunology literature; used for manual annotation
# Satija et al., Nat Biotechnol 2015 (Seurat)
# Stuart et al., Cell 2019 (Seurat v3 integration / PBMC examples)
# 10x Genomics PBMC 3k / 10k datasets documentation


markers_pbmc <- list(
  T_cells = c("TRAC", "CD3D", "CD3E", "LTB"),
  CD4_T   = c("IL7R", "CCR7", "LTB", "MALAT1"),   # MALAT1 optional; can be noisy
  CD8_T   = c("CD8A", "CD8B", "NKG7", "GZMK"),    # avoids NK-only marker set
  NK      = c("NKG7", "GNLY", "PRF1", "FCGR3A", "XCL1"),

  B_cells = c("MS4A1", "CD79A", "CD74", "HLA-DRA"),
  Plasma  = c("MZB1", "XBP1", "JCHAIN"),          # optional, if present

  CD14_Mono = c("LYZ", "S100A8", "S100A9", "LGALS3", "CTSS"),
  FCGR3A_Mono = c("LYZ", "FCGR3A", "MS4A7", "LGALS3"),

  DC      = c("FCER1A", "CST3", "CLEC10A"),
  Platelets = c("PPBP", "PF4", "NRGN")            # NRGN sometimes helps, optional
)
