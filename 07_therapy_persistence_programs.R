# =============================
# Therapy-persistent malignant ecosystem programs
# =============================

library(Seurat)
library(dplyr)
library(ggplot2)

# =============================
# Load tumor-like object
# =============================

tumor_like <- readRDS(
  "data/processed/GSE182434_tumor_like_CopyKAT.rds"
)

# =============================
# Define therapy-persistence programs
# =============================

apoptosis_resistance_genes <- c(
  "BCL2","BCL2L1","MCL1","BIRC5","BCL2L12",
  "XIAP","BAX","BAK1","CASP3","CASP7"
)

dna_damage_checkpoint_genes <- c(
  "CHEK1","CHEK2","ATR","ATM","PARP1",
  "RAD51","XRCC5","XRCC6","PRKDC","MRE11","NBN"
)

replication_stress_genes <- c(
  "RRM2","TYMS","MCM2","MCM4","MCM6",
  "PCNA","CDC6","CDK1","CCNA2","UBE2C","TOP2A"
)

proteostasis_stress_genes <- c(
  "PSMA2","PSMB5","PSMD1","HSPA1A","HSPA8",
  "HSP90AA1","XBP1","ATF4","DDIT3"
)

immune_escape_genes <- c(
  "B2M","HLA-A","HLA-B","HLA-C",
  "HLA-DRA","HLA-DRB1","TAP1","TAP2",
  "CD274","PDCD1LG2"
)

keep_present <- function(genes) {
  genes[genes %in% rownames(tumor_like)]
}

# =============================
# Score therapy-persistence programs
# =============================

tumor_like <- AddModuleScore(
  tumor_like,
  list(keep_present(apoptosis_resistance_genes)),
  name = "ApoptosisResistance"
)

tumor_like <- AddModuleScore(
  tumor_like,
  list(keep_present(dna_damage_checkpoint_genes)),
  name = "DNADamageCheckpoint"
)

tumor_like <- AddModuleScore(
  tumor_like,
  list(keep_present(replication_stress_genes)),
  name = "ReplicationStress"
)

tumor_like <- AddModuleScore(
  tumor_like,
  list(keep_present(proteostasis_stress_genes)),
  name = "ProteostasisStress"
)

tumor_like <- AddModuleScore(
  tumor_like,
  list(keep_present(immune_escape_genes)),
  name = "ImmuneVisibility"
)

# =============================
# Composite therapy-persistence score
# =============================

tumor_like$TherapyPersistenceScore <- with(
  tumor_like@meta.data,
  ApoptosisResistance1 +
    DNADamageCheckpoint1 +
    ReplicationStress1 +
    ProteostasisStress1 -
    ImmuneVisibility1
)

# =============================
# Visualize therapy-persistence programs
# =============================

VlnPlot(
  tumor_like,
  features = c(
    "ApoptosisResistance1",
    "DNADamageCheckpoint1",
    "ReplicationStress1",
    "ProteostasisStress1",
    "ImmuneVisibility1",
    "TherapyPersistenceScore"
  ),
  group.by = "cell_annotation",
  pt.size = 0,
  ncol = 2
)

FeaturePlot(
  tumor_like,
  features = c(
    "ApoptosisResistance1",
    "DNADamageCheckpoint1",
    "ReplicationStress1",
    "ProteostasisStress1",
    "ImmuneVisibility1",
    "TherapyPersistenceScore"
  ),
  ncol = 2
)

VlnPlot(
  tumor_like,
  features = "TherapyPersistenceScore",
  group.by = "cell_annotation",
  pt.size = 0
)

kruskal.test(
  TherapyPersistenceScore ~ cell_annotation,
  data = tumor_like@meta.data
)

pairwise.wilcox.test(
  tumor_like$TherapyPersistenceScore,
  tumor_like$cell_annotation,
  p.adjust.method = "BH"
)

FeaturePlot(
  tumor_like,
  features = "TherapyPersistenceScore",
  cols = c("lightgrey", "red")
)
# =============================
# Compare aneuploid vs diploid cells
# =============================

VlnPlot(
  tumor_like,
  features = "TherapyPersistenceScore",
  group.by = "copykat_prediction",
  pt.size = 0
)

# =============================
# Quantify by malignant state
# =============================

therapy_summary <- tumor_like@meta.data %>%
  group_by(cell_annotation, copykat_prediction) %>%
  summarise(
    n_cells = n(),
    mean_apoptosis_resistance = mean(ApoptosisResistance1, na.rm = TRUE),
    mean_DNA_damage_checkpoint = mean(DNADamageCheckpoint1, na.rm = TRUE),
    mean_replication_stress = mean(ReplicationStress1, na.rm = TRUE),
    mean_proteostasis_stress = mean(ProteostasisStress1, na.rm = TRUE),
    mean_immune_visibility = mean(ImmuneVisibility1, na.rm = TRUE),
    mean_therapy_persistence = mean(TherapyPersistenceScore, na.rm = TRUE),
    .groups = "drop"
  )

therapy_summary

dir.create("results/therapy_persistence", recursive = TRUE, showWarnings = FALSE)

write.csv(
  therapy_summary,
  "results/therapy_persistence/tumor_state_therapy_persistence_summary.csv",
  row.names = FALSE
)

# =============================
# Pseudotime association
# =============================

if ("pseudotime_1" %in% colnames(tumor_like@meta.data)) {
  
  pt_df <- tumor_like@meta.data
  pt_df <- pt_df[!is.na(pt_df$pseudotime_1), ]
  
  cor.test(
    pt_df$pseudotime_1,
    pt_df$TherapyPersistenceScore,
    method = "spearman"
  )
  
  ggplot(
    pt_df,
    aes(
      x = pseudotime_1,
      y = TherapyPersistenceScore,
      color = cell_annotation
    )
  ) +
    geom_point(size = 0.5, alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE) +
    theme_classic()
}

# =============================
# Save final object
# =============================

saveRDS(
  tumor_like,
  "data/processed/GSE182434_tumor_like_therapy_persistence.rds"
)