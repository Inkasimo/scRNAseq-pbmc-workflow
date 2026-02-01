# Canonical PBMC marker genes curated from standard Seurat / 10x Genomics
# PBMC workflows and immunology literature; used for manual annotation
# Satija et al., Nat Biotechnol 2015 (Seurat)
# Stuart et al., Cell 2019 (Seurat v3 integration / PBMC examples)
# 10x Genomics PBMC 3k / 10k datasets documentation


markers_pbmc <- list(
  T_cells = c("CD3D", "CD3E", "IL7R"),
  CD4_T = c("CCR7", "IL7R", "LTB"),
  CD8_T = c("NKG7", "GNLY"),
  B_cells = c("MS4A1", "CD79A", "CD74"),
  NK = c("NKG7", "GNLY", "XCL1"),
  Monocytes = c("LYZ", "S100A8", "S100A9"),
  DC = c("FCER1A", "CST3"),
  Platelets = c("PPBP", "PF4")
)
