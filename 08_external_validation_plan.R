# =====================================
# 08_external_validation_plan.R
# External validation: LBCL atlas
# =====================================

# Run installs only once if needed:
# install.packages("remotes")
# remotes::install_github("mojaveazure/seurat-disk")
# BiocManager::install("zellkonverter")

library(zellkonverter)
library(Seurat)
library(SeuratObject)
library(SingleCellExperiment)
library(dplyr)

# =====================================
# Load h5ad
# =====================================

sce_ext <- readH5AD(
  "data/external_validation/LBCL_atlas_validation.h5ad"
)

# =====================================
# Convert to Seurat
# =====================================

ext_obj <- as.Seurat(
  sce_ext,
  counts = "X",
  data = "X"
)

# =====================================
# Convert Ensembl IDs to gene symbols
# =====================================

gene_symbols <- as.character(rowData(sce_ext)$feature_name)
names(gene_symbols) <- rownames(sce_ext)

expr_ext <- GetAssayData(
  ext_obj,
  assay = "originalexp",
  layer = "counts"
)

rownames(expr_ext) <- make.unique(
  gene_symbols[rownames(expr_ext)]
)

ext_obj_symbol <- CreateSeuratObject(
  counts = expr_ext,
  meta.data = ext_obj@meta.data
)

DefaultAssay(ext_obj_symbol) <- "RNA"

# =====================================
# Inspect metadata
# =====================================

colnames(ext_obj_symbol@meta.data)

table(ext_obj_symbol$sample_type)
table(ext_obj_symbol$histology)
table(ext_obj_symbol$histology_immune)
table(ext_obj_symbol$cell_type)
unique(ext_obj_symbol$sample_id)

grep(
  "relapse|refractory|treat|response|diagnosis|sample|patient",
  colnames(ext_obj_symbol@meta.data),
  value = TRUE,
  ignore.case = TRUE
)

# =====================================
# Subset malignant cells only
# =====================================

ext_malignant <- subset(
  ext_obj_symbol,
  subset = cell_type == "malignant cell"
)

ext_malignant

head(rownames(ext_malignant), 20)

grep(
  "BCL2|CHEK1|RRM2|HLA-A|MKI67",
  rownames(ext_malignant),
  value = TRUE
)

# =====================================
# Normalize before module scoring
# =====================================

DefaultAssay(ext_malignant) <- "RNA"

ext_malignant <- NormalizeData(
  ext_malignant,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# =====================================
# Therapy-persistence programs
# =====================================

dna_damage_genes <- c(
  "ATM","ATR","CHEK1","CHEK2",
  "TP53","BRCA1","BRCA2",
  "RAD51","XRCC5","XRCC6","PRKDC"
)

replication_stress_genes <- c(
  "MCM2","MCM3","MCM4","MCM5",
  "MCM6","PCNA","TYMS","RRM2"
)

proteostasis_genes <- c(
  "HSPA1A","HSP90AA1","PSMB5",
  "PSMD1","XBP1"
)

immune_visibility_genes <- c(
  "B2M","HLA-A","HLA-B","HLA-C"
)

keep_present <- function(genes, obj) {
  present <- genes[genes %in% rownames(obj)]
  print(present)
  return(present)
}

# =====================================
# Score programs
# =====================================

ext_malignant <- AddModuleScore(
  ext_malignant,
  features = list(keep_present(dna_damage_genes, ext_malignant)),
  name = "DNADamageCheckpoint"
)

ext_malignant <- AddModuleScore(
  ext_malignant,
  features = list(keep_present(replication_stress_genes, ext_malignant)),
  name = "ReplicationStress"
)

ext_malignant <- AddModuleScore(
  ext_malignant,
  features = list(keep_present(proteostasis_genes, ext_malignant)),
  name = "ProteostasisStress"
)

ext_malignant <- AddModuleScore(
  ext_malignant,
  features = list(keep_present(immune_visibility_genes, ext_malignant)),
  name = "ImmuneVisibility"
)

grep(
  "DNA|Replication|Proteostasis|Immune",
  colnames(ext_malignant@meta.data),
  value = TRUE
)

# =====================================
# Composite external validation score
# Reduced score: apoptosis excluded
# =====================================

ext_malignant$TherapyPersistenceScore <-
  as.numeric(scale(ext_malignant$DNADamageCheckpoint1)) +
  as.numeric(scale(ext_malignant$ReplicationStress1)) +
  as.numeric(scale(ext_malignant$ProteostasisStress1)) -
  as.numeric(scale(ext_malignant$ImmuneVisibility1))

# =====================================
# Visualize external validation
# =====================================

# Check available reductions
Reductions(ext_malignant)



VlnPlot(
  ext_malignant,
  features = "TherapyPersistenceScore",
  group.by = "sample_id",
  pt.size = 0
)

VlnPlot(
  ext_malignant,
  features = "TherapyPersistenceScore",
  group.by = "sample_type",
  pt.size = 0
)

# =====================================
# Quantify by sample
# =====================================

external_summary <- ext_malignant@meta.data %>%
  group_by(sample_id, sample_type) %>%
  summarise(
    n_cells = n(),
    mean_DNA_damage_checkpoint = mean(DNADamageCheckpoint1, na.rm = TRUE),
    mean_replication_stress = mean(ReplicationStress1, na.rm = TRUE),
    mean_proteostasis_stress = mean(ProteostasisStress1, na.rm = TRUE),
    mean_immune_visibility = mean(ImmuneVisibility1, na.rm = TRUE),
    mean_therapy_persistence = mean(TherapyPersistenceScore, na.rm = TRUE),
    .groups = "drop"
  )

external_summary

dir.create(
  "results/external_validation",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  external_summary,
  "results/external_validation/LBCL_external_therapy_persistence_summary.csv",
  row.names = FALSE
)

saveRDS(
  ext_malignant,
  "data/processed/LBCL_external_malignant_therapy_persistence.rds"
)

ext_malignant <- FindVariableFeatures(ext_malignant)
ext_malignant <- ScaleData(ext_malignant)
ext_malignant <- RunPCA(ext_malignant)
ext_malignant <- RunUMAP(ext_malignant, dims = 1:20)
ext_malignant <- FindNeighbors(ext_malignant, dims = 1:20)
ext_malignant <- FindClusters(ext_malignant, resolution = 0.4)

FeaturePlot(
  ext_malignant,
  features = "TherapyPersistenceScore",
  cols = c("lightgrey", "red")
)

VlnPlot(
  ext_malignant,
  features = "TherapyPersistenceScore",
  group.by = "seurat_clusters",
  pt.size = 0
)

# Find markers of high-persistence cluster

Idents(ext_malignant) <- "seurat_clusters"

high_cluster_markers <- FindMarkers(
  ext_malignant,
  ident.1 = "1",
  only.pos = TRUE,
  logfc.threshold = 0.25
)

head(high_cluster_markers, 20)

FeaturePlot(
  ext_malignant,
  features = c(
    "MKI67",
    "RRM2",
    "CHEK1",
    "FOXP1",
    "TCL1A",
    "HLA-A"
  ),
  reduction = "umap",
  cols = c("lightgrey", "red"),
  ncol = 3
)
