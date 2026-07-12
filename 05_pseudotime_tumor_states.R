# =============================
# Pseudotime / trajectory analysis
# Tumor-like DLBCL states
# =============================

library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(ggplot2)

obj <- readRDS("data/processed/GSE182434_ecosystem_scored.rds")

obj$cell_annotation <- as.character(Idents(obj))

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
# Convert to SingleCellExperiment
# =============================

sce <- as.SingleCellExperiment(tumor_like)

# Add UMAP coordinates
reducedDims(sce)$UMAP <- Embeddings(tumor_like, "umap")

# =============================
# Run Slingshot
# =============================

sce <- slingshot(
  sce,
  clusterLabels = tumor_like$cell_annotation,
  reducedDim = "UMAP",
  start.clus = "Cycling_B"
)

# =============================
# Plot trajectory
# =============================

plot(
  reducedDims(sce)$UMAP,
  col = as.numeric(factor(tumor_like$cell_annotation)),
  pch = 16,
  asp = 1,
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

lines(SlingshotDataSet(sce), lwd = 2)

legend(
  "topright",
  legend = levels(factor(tumor_like$cell_annotation)),
  col = seq_along(levels(factor(tumor_like$cell_annotation))),
  pch = 16
)

# =============================
# Extract pseudotime
# =============================

pt <- slingPseudotime(sce)

head(pt)

tumor_like$pseudotime_1 <- pt[, 1]


pt <- slingPseudotime(sce)

tumor_like$pseudotime_1 <- pt[, 1]
tumor_like$pseudotime_2 <- pt[, 2]

head(tumor_like@meta.data[, c("pseudotime_1", "pseudotime_2")])
# =============================
# Visualize pseudotime
# =============================

FeaturePlot(
  tumor_like,
  features = c("pseudotime_1", "pseudotime_2"),
  reduction = "umap",
  slot = "data"
)
VlnPlot(
  tumor_like,
  features = "pseudotime_1",
  group.by = "cell_annotation",
  pt.size = 0
)

#========================
# trouble shoot
#====================
# =============================
# Recreate RT-relevant scores
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
  genes[genes %in% rownames(tumor_like)]
}

tumor_like <- AddModuleScore(tumor_like, list(keep_present(dna_damage_genes)), name = "DNADamageRepairScore")
tumor_like <- AddModuleScore(tumor_like, list(keep_present(hypoxia_genes)), name = "HypoxiaScore")
tumor_like <- AddModuleScore(tumor_like, list(keep_present(apoptosis_genes)), name = "ApoptosisResistanceScore")

pt <- slingPseudotime(sce)

tumor_like$pseudotime_1 <- pt[, 1]
tumor_like$pseudotime_2 <- pt[, 2]

table(is.na(tumor_like$pseudotime_1))
table(is.na(tumor_like$pseudotime_2))

pt_df <- tumor_like@meta.data
pt_df <- pt_df[!is.na(pt_df$pseudotime_1), ]

dim(pt_df)

cor.test(pt_df$pseudotime_1, pt_df$DNADamageRepairScore1, method = "spearman")
cor.test(pt_df$pseudotime_1, pt_df$ApoptosisResistanceScore1, method = "spearman")
cor.test(pt_df$pseudotime_1, pt_df$HypoxiaScore1, method = "spearman")
# =============================
# Save object
# =============================

saveRDS(
  tumor_like,
  "data/processed/GSE182434_tumor_like_pseudotime.rds"
)


# =============================
# Test RT-relevant programs along pseudotime
# =============================

library(ggplot2)

pt_df <- tumor_like@meta.data

pt_df <- pt_df[!is.na(pt_df$pseudotime_1), ]

# DNA repair vs pseudotime
ggplot(
  pt_df,
  aes(
    x = pseudotime_1,
    y = DNADamageRepairScore1,
    color = cell_annotation
  )
) +
  geom_point(size = 0.5, alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_classic()

# Apoptosis resistance vs pseudotime
ggplot(pt_df, aes(x = pseudotime_1, y = ApoptosisResistanceScore1, color = cell_annotation)) +
  geom_point(size = 0.5, alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_classic()

# Hypoxia/stress vs pseudotime
ggplot(pt_df, aes(x = pseudotime_1, y = HypoxiaScore1, color = cell_annotation)) +
  geom_point(size = 0.5, alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_classic()


cor.test(
  pt_df$pseudotime_1,
  pt_df$DNADamageRepairScore1,
  method = "spearman"
)

cor.test(
  pt_df$pseudotime_1,
  pt_df$ApoptosisResistanceScore1,
  method = "spearman"
)

cor.test(
  pt_df$pseudotime_1,
  pt_df$HypoxiaScore1,
  method = "spearman"
)

lm_dna <- lm(DNADamageRepairScore1 ~ pseudotime_1, data = pt_df)
summary(lm_dna)

lm_apop <- lm(ApoptosisResistanceScore1 ~ pseudotime_1, data = pt_df)
summary(lm_apop)

lm_hypoxia <- lm(HypoxiaScore1 ~ pseudotime_1, data = pt_df)
summary(lm_hypoxia)

lm_dna <- lm(DNADamageRepairScore1 ~ pseudotime_1, data = pt_df)
summary(lm_dna)

lm_apop <- lm(ApoptosisResistanceScore1 ~ pseudotime_1, data = pt_df)
summary(lm_apop)

lm_dna <- lm(DNADamageRepairScore1 ~ pseudotime_1, data = pt_df)
summary(lm_dna)