library(Seurat)
library(data.table)
library(Matrix)
library(dplyr)
library(ggplot2)
library(harmony)

# -----------------------------
# Load counts
# -----------------------------

counts <- fread(
  "/Users/ahmedrezk/Desktop/DLBCL_singlecell_RT_project/data/raw/GSE182434_raw_count_matrix.txt.gz",
  data.table = FALSE
)

# -----------------------------
# Set rownames
# -----------------------------

rownames(counts) <- counts[,1]
counts <- counts[,-1]

# -----------------------------
# Convert to sparse matrix
# -----------------------------

counts_matrix <- as(as.matrix(counts), "dgCMatrix")

rm(counts)
gc()

# -----------------------------
# Create Seurat object
# -----------------------------

obj <- CreateSeuratObject(
  counts = counts_matrix,
  project = "GSE182434",
  min.cells = 3,
  min.features = 200
)

# -----------------------------
# QC metrics
# -----------------------------

obj[["percent.mt"]] <- PercentageFeatureSet(
  obj,
  pattern = "^MT-"
)

# -----------------------------
# QC filtering
# -----------------------------

obj <- subset(
  obj,
  subset =
    nFeature_RNA > 300 &
    nFeature_RNA < 6000 &
    percent.mt < 15
)

# -----------------------------
# Normalize
# -----------------------------
options(future.globals.maxSize = 8 * 1024^3)

library(future)
plan("sequential")


obj <- SCTransform(
  obj,
  vars.to.regress = "percent.mt",
  verbose = FALSE
)

# -----------------------------
# PCA
# -----------------------------

obj <- RunPCA(obj)
ElbowPlot(obj)

# -----------------------------
# Neighbors
# -----------------------------

obj <- FindNeighbors(
  obj,
  dims = 1:15
)

# -----------------------------
# Clustering
# -----------------------------

obj <- FindClusters(
  obj,
  resolution = 0.5
)
# -----------------------------
# UMAP
# -----------------------------

obj <- RunUMAP(
  obj,
  dims = 1:15
)

DimPlot(
  obj,
  reduction = "umap",
  label = TRUE
)
# -----------------------------
# Save UMAP
# -----------------------------

pdf("results/umap/GSE182434_umap.pdf")

DimPlot(
  obj,
  label = TRUE
)

dev.off()

# Find cluster markers
markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

markers %>%
  filter(
    gene %in% c(
      "CD3D","CD3E","TRBC1","IL7R","LTB",
      "MS4A1","CD79A","CD74","HLA-DRA",
      "NKG7","GNLY","KLRD1",
      "LYZ","FCN1","S100A8","CTSS",
      "MKI67","TOP2A",
      "IFI6","ISG15","IFITM3","STAT1",
      "JCHAIN","MZB1","SDC1"
    )
  ) %>%
  arrange(cluster, desc(avg_log2FC))

new.cluster.ids <- c(
  "CD4_T_IFN",
  "CD4_T",
  "Activated_T_NK",
  "B_cells",
  "T_cells",
  "NK_cells_IFN",
  "Plasmablast",
  "T_cells_2",
  "B_cells_IFN",
  "B_cells_2",
  "CD4_T_IFN_2",
  "T_cells_IFN",
  "Plasmablast_2",
  "T_IFN",
  "T_NK",
  "T_IFN_2",
  "Inflammatory_Myeloid",
  "Cycling_B",
  "Cycling_B_2",
  "Plasma_like",
  "B_cells_STAT1",
  "Plasma_CTSS",
  "NK_cytotoxic",
  "Plasma_IFN",
  "Highly_Cycling_B"
)

names(new.cluster.ids) <- levels(obj)

obj <- RenameIdents(obj, new.cluster.ids)


FeaturePlot(obj, features = c("STAT1","ISG15","IFI6","IFITM3"))
FeaturePlot(obj, features = c("CD3D","MS4A1","NKG7","LYZ"))

DimPlot(
  obj,
  reduction = "umap",
  label = TRUE
)

# -----------------------------
# Save object
# -----------------------------

saveRDS(
  obj,
  file = "data/processed/GSE182434_processed.rds"
)

saveRDS(
  obj,
  file = "data/processed/GSE182434_annotated.rds"
)
