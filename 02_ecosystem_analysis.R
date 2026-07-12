

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/module_scores", recursive = TRUE, showWarnings = FALSE)
dir.create("results/markers", recursive = TRUE, showWarnings = FALSE)


# =============================
# Load libraries and object
# =============================

library(Seurat)
library(dplyr)
library(ggplot2)

obj <- readRDS("data/processed/GSE182434_annotated.rds")


# =============================
# Visualize initial annotations
# =============================

DimPlot(
  obj,
  reduction = "umap",
  label = TRUE
)

table(Idents(obj))


# =============================
# Define ecosystem compartments cycling=proliferative and highly cycling = stress proliferative
# =============================

obj$ecosystem_compartment <- case_when(
  Idents(obj) %in% c("CD4_T_IFN","CD4_T","T_cells","T_cells_2",
                     "CD4_T_IFN_2","T_cells_IFN","T_IFN","T_IFN_2") ~ "T_cells",
  Idents(obj) %in% c("T_NK","NK_cytotoxic","NK_cells_IFN","Activated_T_NK") ~ "NK_cytotoxic",
  Idents(obj) %in% c("Inflammatory_Myeloid") ~ "Myeloid",
  Idents(obj) %in% c("B_cells","B_cells_2","B_cells_IFN","B_cells_STAT1") ~ "B_cells",
  Idents(obj) %in% c("Plasmablast","Plasmablast_2","Plasma_like","Plasma_CTSS","Plasma_IFN") ~ "Plasma_like",
  Idents(obj) %in% c("Cycling_B","Cycling_B_2","Highly_Cycling_B") ~ "Cycling_B",
  TRUE ~ "Other"
)

DimPlot(
  obj,
  reduction = "umap",
  group.by = "ecosystem_compartment",
  label = TRUE
)


# =============================
# Define biologic gene programs
# =============================

prolif_genes <- c("MKI67","TOP2A","STMN1","HMGB2","PCNA","TYMS","UBE2C","BIRC5")

antigen_genes <- c(
  "B2M","HLA-A","HLA-B","HLA-C",
  "HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1",
  "TAP1","TAP2"
)

myeloid_genes <- c(
  "LYZ","FCN1","S100A8","S100A9","CTSS",
  "C1QA","C1QB","C1QC","APOE"
)

plasma_genes <- c(
  "JCHAIN","MZB1","XBP1","PRDM1"
)


# =============================
# Keep genes present in dataset
# =============================

keep_present <- function(genes) {
  genes[genes %in% rownames(obj)]
}


# =============================
# Calculate module scores
# =============================

obj <- AddModuleScore(obj, list(keep_present(prolif_genes)), name = "ProliferationScore")
obj <- AddModuleScore(obj, list(keep_present(antigen_genes)), name = "AntigenPresentationScore")
obj <- AddModuleScore(obj, list(keep_present(myeloid_genes)), name = "MyeloidInflammationScore")
obj <- AddModuleScore(obj, list(keep_present(plasma_genes)), name = "PlasmaLikeScore")


# =============================
# Build proliferation(persistence)-associated score
# =============================

obj$PersistenceProgramScore <- with(
  obj@meta.data,
  ProliferationScore1 +
    MyeloidInflammationScore1 +
    PlasmaLikeScore1 -
    AntigenPresentationScore1
)


# =============================
# Visualize ecosystem scores
# =============================

FeaturePlot(
  obj,
  features = c(
    "ProliferationScore1",
    "AntigenPresentationScore1",
    "MyeloidInflammationScore1",
    "PlasmaLikeScore1",
    "PersistenceProgramScore"
  ),
  ncol = 2
)

VlnPlot(
  obj,
  features = "PersistenceProgramScore",
  group.by = "ecosystem_compartment",
  pt.size = 0
)

DotPlot(
  obj,
  features = c(
    "MKI67","TOP2A",
    "B2M","HLA-DRA","TAP1",
    "LYZ","FCN1","S100A8",
    "JCHAIN","MZB1"
  ),
  group.by = "ecosystem_compartment"
) + RotatedAxis()


# =============================
# Save scored Seurat object
# =============================

saveRDS(
  obj,
  file = "data/processed/GSE182434_ecosystem_scored.rds"
)


# =============================
# Create cell annotation metadata
# =============================

obj$cell_annotation <- as.character(Idents(obj))


# =============================
# Create B-cell focused subset
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
# Compare persistence scores
# =============================

VlnPlot(
  b_obj,
  features = "PersistenceProgramScore",
  group.by = "cell_annotation",
  pt.size = 0
)


# =============================
# Export figures and metadata
# =============================

ggsave(
  "results/figures/Bcell_PersistenceProgramScore_violin.pdf",
  width = 8,
  height = 5
)

write.csv(
  b_obj@meta.data,
  "results/module_scores/Bcell_subset_persistence_scores.csv",
  row.names = TRUE
)


# =============================
# Compare module scores
# =============================

VlnPlot(
  b_obj,
  features = c(
    "ProliferationScore1",
    "AntigenPresentationScore1",
    "PlasmaLikeScore1",
    "PersistenceProgramScore"
  ),
  group.by = "cell_annotation",
  pt.size = 0,
  ncol = 2
)


# =============================
# Visualize lineage markers
# =============================

FeaturePlot(
  b_obj,
  features = c(
    "MS4A1",
    "CD79A",
    "CD74",
    "MKI67",
    "B2M",
    "JCHAIN"
  ),
  ncol = 3
)


# =============================
# Define tumor-like populations
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
  label = TRUE
)

Idents(tumor_like) <- tumor_like$cell_annotation


# =============================
# Identify Highly_Cycling_B markers
# =============================

highcycling_markers <- FindMarkers(
  tumor_like,
  ident.1 = "Highly_Cycling_B",
  ident.2 = c("Cycling_B","B_cells_STAT1","Plasma_like"),
  logfc.threshold = 0.25,
  min.pct = 0.1
)

head(highcycling_markers, 20)

write.csv(
  highcycling_markers,
  "results/markers/Highly_Cycling_B_markers.csv"
)


# =============================
# Visualize Highly_Cycling_B genes
# =============================

FeaturePlot(
  tumor_like,
  features = c(
    "MKI67",
    "MACROD2",
    "SNHG6",
    "PSMA2",
    "NIBAN3",
    "IGHV3-74"
  ),
  ncol = 3
)


# =============================
# Compare Highly_Cycling_B vs B_cells_STAT1
# =============================

hc_vs_stat1 <- FindMarkers(
  tumor_like,
  ident.1 = "Highly_Cycling_B",
  ident.2 = "B_cells_STAT1",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

head(hc_vs_stat1, 20)


# =============================
# Visualize proliferation-associated genes
# =============================

FeaturePlot(
  tumor_like,
  features = c(
    "BIRC5",
    "RRM2",
    "UBE2C",
    "TYMS",
    "BCL2L12",
    "STAT1"
  ),
  ncol = 3
)


# =============================
# Visualize immunoglobulin-associated genes
# =============================

FeaturePlot(
  tumor_like,
  features = c(
    "IGHV3-74",
    "IGLV4-69",
    "IGHV3-49",
    "MKI67",
    "BIRC5"
  ),
  ncol = 3
)


# =============================
# Visualize IFN / antigen-presentation genes
# =============================

FeaturePlot(
  tumor_like,
  features = c(
    "STAT1",
    "IFI6",
    "ISG15",
    "HLA-DRA",
    "B2M",
    "TAP1"
  ),
  ncol = 3
)


# =============================
# Compare proliferation vs antigen presentation
# =============================

FeatureScatter(
  tumor_like,
  feature1 = "ProliferationScore1",
  feature2 = "AntigenPresentationScore1"
)


# =============================
# Compare proliferation vs persistence
# =============================

FeatureScatter(
  tumor_like,
  feature1 = "ProliferationScore1",
  feature2 = "PersistenceProgramScore"
)

saveRDS(obj, "data/processed/GSE182434_ecosystem_scored.rds")
saveRDS(b_obj, "data/processed/GSE182434_Bcell_subset.rds")
saveRDS(tumor_like, "data/processed/GSE182434_tumor_like_subset.rds")

write.csv(hc_vs_stat1, "results/markers/HighlyCyclingB_vs_BcellsSTAT1_markers.csv")

persistence_signature <- c(
  "MKI67",
  "TOP2A",
  "NUSAP1",
  "CDK1",
  "UBE2S",
  "CCNA2",
  "CENPF",
  "TPX2",
  "UBE2C",
  "RRM2",
  "CDCA3",
  "BIRC5",
  "TYMS"
)