# =============================
# Malignancy validation using CopyKAT
# DLBCL tumor-like states
# =============================

library(Seurat)
library(dplyr)
library(ggplot2)

# Install once if needed:
# install.packages("remotes")
# remotes::install_github("navinlabcode/copykat")

library(copykat)

# =============================
# Load scored object
# =============================

obj <- readRDS("data/processed/GSE182434_ecosystem_scored.rds")

obj$cell_annotation <- as.character(Idents(obj))

# =============================
# Subset B-cell / tumor-like compartments
# =============================

b_obj <- subset(
  obj,
  subset = ecosystem_compartment %in% c(
    "B_cells",
    "Cycling_B",
    "Plasma_like"
  )
)

# =============================
# Extract raw counts
# =============================

counts <- GetAssayData(
  b_obj,
  assay = "RNA",
  layer = "counts"
)
counts_mat <- as.matrix(counts)
# =============================
# Run CopyKAT
# =============================
# This can take time.
# ngene.chr = 5 is commonly used.
# KS.cut can be adjusted if needed.

copykat_result <- copykat(
  rawmat = counts_mat,
  id.type = "S",
  ngene.chr = 5,
  win.size = 25,
  KS.cut = 0.1,
  sam.name = "GSE182434_Bcell",
  distance = "euclidean",
  norm.cell.names = "",
  output.seg = "FALSE",
  plot.genes = "TRUE",
  genome = "hg20"
)

# =============================
# Extract CopyKAT prediction
# =============================

copykat_pred <- copykat_result$prediction

head(copykat_pred)

# CopyKAT usually returns columns:
# cell.names
# copykat.pred

colnames(copykat_pred)

# =============================
# Add predictions to Seurat object
# =============================

pred_df <- copykat_pred

rownames(pred_df) <- pred_df$cell.names

b_obj$copykat_prediction <- pred_df[
  colnames(b_obj),
  "copykat.pred"
]
table(
  b_obj$cell_annotation,
  b_obj$copykat_prediction
)
# =============================
# Visualize CopyKAT calls
# =============================

DimPlot(
  b_obj,
  reduction = "umap",
  group.by = "copykat_prediction",
  label = TRUE
)

DimPlot(
  b_obj,
  reduction = "umap",
  group.by = "cell_annotation",
  label = TRUE
)

# =============================
# Cross-tabulate annotations vs CopyKAT
# =============================

copykat_table <- table(
  b_obj$cell_annotation,
  b_obj$copykat_prediction
)

copykat_table

write.csv(
  as.data.frame(copykat_table),
  "results/markers/CopyKAT_annotation_table.csv",
  row.names = FALSE
)

# =============================
# Focus on tumor-like states
# =============================

tumor_like <- subset(
  b_obj,
  subset = cell_annotation %in% c(
    "Cycling_B",
    "Highly_Cycling_B",
    "B_cells_STAT1",
    "Plasma_like"
  )
)

DimPlot(
  tumor_like,
  reduction = "umap",
  group.by = "copykat_prediction",
  label = TRUE
)

table(
  tumor_like$cell_annotation,
  tumor_like$copykat_prediction
)

# =============================
# Save object
# =============================

saveRDS(
  b_obj,
  "data/processed/GSE182434_Bcell_CopyKAT.rds"
)

saveRDS(
  tumor_like,
  "data/processed/GSE182434_tumor_like_CopyKAT.rds"
)

