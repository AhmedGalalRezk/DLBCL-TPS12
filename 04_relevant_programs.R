# =============================
#  mechanistic program overlays
# =============================

library(Seurat)
library(dplyr)
library(ggplot2)

obj <- readRDS("data/processed/GSE182434_ecosystem_scored.rds")

dir.create("results/RT_programs", recursive = TRUE, showWarnings = FALSE)

# =============================
# Recreate tumor-like subset
# =============================

obj$cell_annotation <- as.character(Idents(obj))

b_obj <- subset(
  obj,
  subset = ecosystem_compartment %in% c(
    "B_cells",
    "Cycling_B",
    "Plasma_like"
  )
)

tumor_like <- subset(
  b_obj,
  subset = cell_annotation %in% c(
    "Cycling_B",
    "Highly_Cycling_B",
    "B_cells_STAT1",
    "Plasma_like"
  )
)

# =============================
# Define relevant gene programs
# =============================

dna_damage_genes <- c(
  "ATM","ATR","CHEK1","CHEK2","TP53","BRCA1","BRCA2",
  "RAD51","XRCC5","XRCC6","PRKDC","PARP1","MRE11","NBN"
)

hypoxia_genes <- c(
  "HIF1A","VEGFA","CA9","SLC2A1","LDHA","ENO1","PGK1"
)

apoptosis_genes <- c(
  "BAX","BAK1","CASP3","CASP7","BCL2","BCL2L1","MCL1","BIRC5"
)

keep_present <- function(genes) {
  genes[genes %in% rownames(obj)]
}

# =============================
# Score relevant programs
# =============================

obj <- AddModuleScore(obj, list(keep_present(dna_damage_genes)), name = "DNADamageRepairScore")
obj <- AddModuleScore(obj, list(keep_present(hypoxia_genes)), name = "HypoxiaScore")
obj <- AddModuleScore(obj, list(keep_present(apoptosis_genes)), name = "ApoptosisResistanceScore")

# =============================
# Transfer scores to tumor-like subset
# =============================

tumor_like <- subset(
  obj,
  subset = cell_annotation %in% c(
    "Cycling_B",
    "Highly_Cycling_B",
    "B_cells_STAT1",
    "Plasma_like"
  )
)

# =============================
# Visualize relevant programs
# =============================

VlnPlot(
  tumor_like,
  features = c(
    "DNADamageRepairScore1",
    "HypoxiaScore1",
    "ApoptosisResistanceScore1"
  ),
  group.by = "cell_annotation",
  pt.size = 0,
  ncol = 3
)

FeaturePlot(
  tumor_like,
  features = c(
    "DNADamageRepairScore1",
    "HypoxiaScore1",
    "ApoptosisResistanceScore1"
  ),
  ncol = 3
)


# =============================
# Compare relevant programs
# across tumor states
# =============================

VlnPlot(
  tumor_like,
  features = c(
    "DNADamageRepairScore1",
    "HypoxiaScore1",
    "ApoptosisResistanceScore1"
  ),
  group.by = "cell_annotation",
  pt.size = 0,
  ncol = 3
)

FeatureScatter(
  tumor_like,
  feature1 = "ProliferationScore1",
  feature2 = "DNADamageRepairScore1"
)

FeatureScatter(
  tumor_like,
  feature1 = "ProliferationScore1",
  feature2 = "HypoxiaScore1"
)

FeatureScatter(
  tumor_like,
  feature1 = "ProliferationScore1",
  feature2 = "ApoptosisResistanceScore1"
)

# ====================================
# Stress-adapted proliferation analysis
# Highly_Cycling_B vs Cycling_B
# ====================================

hc_vs_cycling <- FindMarkers(
  tumor_like,
  ident.1 = "Highly_Cycling_B",
  ident.2 = "Cycling_B",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

head(hc_vs_cycling, 30)

FeaturePlot(
  tumor_like,
  features = c(
    "BIRC5",
    "MCL1",
    "BCL2",
    "HIF1A",
    "VEGFA",
    "CHEK1",
    "CHEK2",
    "PARP1",
    "XRCC5",
    "ATF4",
    "DDIT3"
  ),
  ncol = 4
)

 VlnPlot(
  tumor_like,
  features = c(
    "CHEK1",
    "PARP1",
    "XRCC5",
    "BIRC5"
  ),
  group.by = "cell_annotation",
  pt.size = 0,
  ncol = 2
)
