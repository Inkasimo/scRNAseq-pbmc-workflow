deg_label_col <- "cell_type_cluster_majority"

deg_sets <- list(
  T_like    = c("T_cells","CD4_T","CD8_T","NK"),
  B_like    = c("B_cells","Plasma"),
  Mono_like = c("CD14_Mono","FCGR3A_Mono")
)

deg_contrasts <- list(
  c("T_like","B_like"),
  c("T_like","Mono_like"),
  c("B_like","Mono_like")
)

network_sets <- list(
  CD4    = c("CD4_T"),
  B_cells    = c("B_cells"),
  CD14_mono = c("CD14_Mono")
)
