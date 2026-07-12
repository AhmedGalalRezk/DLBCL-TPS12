# =========================================
# Single-cell TPS integration
# Which cells express the therapy-persistence program?
# =========================================

library(Seurat)
library(dplyr)
library(ggplot2)

# =========================================
# Load discovery single-cell object
# =========================================

obj <- readRDS("data/processed/GSE182434_ecosystem_scored.rds")

obj$cell_annotation <- as.character(Idents(obj))

# =========================================
# Define locked TPS panel
# =========================================

TPS_genes_12 <- c(
  "MKI67",
  "RRM2",
  "CHEK1",
  "FOXP1",
  "TCL1A",
  "MCM2",
  "MCM3",
  "MCM4",
  "MCM5",
  "MCM6",
  "PCNA",
  "TYMS"
)

keep_present <- function(genes, obj) {
  genes[genes %in% rownames(obj)]
}

TPS_present <- keep_present(TPS_genes_12, obj)
TPS_present

# =========================================
# Score TPS in all single cells
# =========================================

obj <- AddModuleScore(
  obj,
  features = list(TPS_present),
  name = "TPS12"
)

# =========================================
# Visualize TPS across whole ecosystem
# =========================================

FeaturePlot(
  obj,
  features = "TPS121",
  reduction = "umap",
  cols = c("lightgrey", "red")
)

VlnPlot(
  obj,
  features = "TPS121",
  group.by = "cell_annotation",
  pt.size = 0
)

# =========================================
# Summarize TPS by annotated cell state
# =========================================

TPS_cell_summary <- obj@meta.data %>%
  group_by(cell_annotation) %>%
  summarise(
    n_cells = n(),
    mean_TPS12 = mean(TPS121, na.rm = TRUE),
    median_TPS12 = median(TPS121, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_TPS12))

TPS_cell_summary

dir.create("results/singlecell_TPS", recursive = TRUE, showWarnings = FALSE)

write.csv(
  TPS_cell_summary,
  "results/singlecell_TPS/TPS12_by_cell_state.csv",
  row.names = FALSE
)

# =========================================
# Compare malignant/anueploid CopyKAT states if available
# =========================================

copykat_path <- "data/processed/GSE182434_Bcell_CopyKAT.rds"

if (file.exists(copykat_path)) {
  
  b_copykat <- readRDS(copykat_path)
  
  b_copykat$cell_annotation <- as.character(Idents(b_copykat))
  
  TPS_present_b <- keep_present(TPS_genes_12, b_copykat)
  
  b_copykat <- AddModuleScore(
    b_copykat,
    features = list(TPS_present_b),
    name = "TPS12"
  )
  
  VlnPlot(
    b_copykat,
    features = "TPS121",
    group.by = "copykat_prediction",
    pt.size = 0
  )
  
  TPS_copykat_summary <- b_copykat@meta.data %>%
    group_by(copykat_prediction) %>%
    summarise(
      n_cells = n(),
      mean_TPS12 = mean(TPS121, na.rm = TRUE),
      median_TPS12 = median(TPS121, na.rm = TRUE),
      .groups = "drop"
    )
  
  TPS_copykat_summary
  
  write.csv(
    TPS_copykat_summary,
    "results/singlecell_TPS/TPS12_by_CopyKAT_prediction.csv",
    row.names = FALSE
  )
}

# =========================================
# Tumor-like compartment TPS analysis
# =========================================

tumor_like <- subset(
  obj,
  subset = cell_annotation %in% c(
    "Cycling_B",
    "Highly_Cycling_B",
    "B_cells_STAT1",
    "Plasma_like"
  )
)

FeaturePlot(
  tumor_like,
  features = "TPS121",
  reduction = "umap",
  cols = c("lightgrey", "red")
)

VlnPlot(
  tumor_like,
  features = "TPS121",
  group.by = "cell_annotation",
  pt.size = 0
)

# =========================================
# Statistics across tumor-like states
# =========================================

kruskal.test(
  TPS121 ~ cell_annotation,
  data = tumor_like@meta.data
)

pairwise.wilcox.test(
  tumor_like$TPS121,
  tumor_like$cell_annotation,
  p.adjust.method = "BH"
)

# =========================================
# Save single-cell TPS object
# =========================================

saveRDS(
  obj,
  "data/processed/GSE182434_singlecell_TPS12_scored.rds"
)

VlnPlot(
  b_copykat,
  features = "TPS121",
  group.by = "copykat_prediction",
  pt.size = 0
)

wilcox.test(
  TPS121 ~ copykat_prediction,
  data = b_copykat@meta.data
)

kruskal.test(
  TPS121 ~ copykat_prediction,
  data = b_copykat@meta.data
)

copykat_clean <- subset(
  b_copykat@meta.data,
  copykat_prediction %in% c("aneuploid", "diploid")
)

wilcox.test(
  TPS121 ~ copykat_prediction,
  data = copykat_clean
)