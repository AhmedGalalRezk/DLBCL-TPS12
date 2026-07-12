# =============================
# Bulk validation: GSE31312
# Therapy-proliferation signature
# =============================

library(GEOquery)
library(limma)
library(survival)
library(survminer)

gse <- getGEO("GSE31312", GSEMatrix = TRUE)
gse <- gse[[1]]

expr <- exprs(gse)
pheno <- pData(gse)
features <- fData(gse)

colnames(features)

symbol_col <- grep(
  "symbol|gene.symbol|gene symbol",
  colnames(features),
  ignore.case = TRUE,
  value = TRUE
)[1]

features$gene_symbol <- features[[symbol_col]]

valid <- !is.na(features$gene_symbol) & features$gene_symbol != ""

expr_valid <- expr[valid, ]
features_valid <- features[valid, ]

features_valid$gene_symbol <- sapply(
  strsplit(as.character(features_valid$gene_symbol), "///"),
  function(x) trimws(x[1])
)

expr_gene <- avereps(
  expr_valid,
  ID = features_valid$gene_symbol
)

# =============================
# Therapy-persistence signature
# =============================

therapy_genes <- c(
  "MKI67","RRM2","CHEK1","FOXP1","TCL1A","TCL1B",
  "HMGB2","FBL","TUBB","PGAM1","NIBAN3","SPIB",
  "MCM2","MCM3","MCM4","MCM5","MCM6","PCNA","TYMS"
)

sig_present <- therapy_genes[therapy_genes %in% rownames(expr_gene)]
sig_missing <- setdiff(therapy_genes, sig_present)

sig_present
sig_missing

expr_sig <- expr_gene[sig_present, , drop = FALSE]
expr_sig_z <- t(scale(t(expr_sig)))

pheno$TherapyPersistenceScore <- colMeans(expr_sig_z, na.rm = TRUE)

pheno$TherapyPersistenceGroup <- ifelse(
  pheno$TherapyPersistenceScore >= median(pheno$TherapyPersistenceScore, na.rm = TRUE),
  "High",
  "Low"
)

table(pheno$TherapyPersistenceGroup)

# =============================
# Inspect clinical columns
# =============================

colnames(pheno)

grep(
  "surv|os|pfs|progress|death|status|follow|time|event|response|relapse",
  colnames(pheno),
  ignore.case = TRUE,
  value = TRUE
)

write.csv(
  pheno,
  "results/external_validation/GSE31312_metadata_with_therapy_persistence_score.csv",
  row.names = TRUE
)

# =====================================
# Extract COO subtype
# =====================================

pheno$COO <- gsub(
  "immunohistochemistry subgroup: ",
  "",
  pheno$characteristics_ch1.1
)

table(pheno$COO)

# =====================================
# Compare therapy persistence by COO
# =====================================

boxplot(
  pheno$TherapyPersistenceScore ~ pheno$COO,
  main = "TherapyPersistenceScore by COO subtype",
  ylab = "TherapyPersistenceScore",
  xlab = "COO subtype"
)

# Statistical test
wilcox.test(
  TherapyPersistenceScore ~ COO,
  data = subset(pheno, COO %in% c("ABC", "GCB"))
)


# =====================================
# Molecular COO subtype
# =====================================

pheno$molecular_COO <- gsub(
  "gene expression profiling subgroup: ",
  "",
  pheno$characteristics_ch1
)

table(pheno$molecular_COO)

# =====================================
# Compare persistence score
# =====================================

boxplot(
  TherapyPersistenceScore ~ molecular_COO,
  data = pheno,
  main = "TherapyPersistenceScore by molecular COO",
  ylab = "TherapyPersistenceScore",
  xlab = "Molecular COO"
)

wilcox.test(
  TherapyPersistenceScore ~ molecular_COO,
  data = subset(pheno, molecular_COO %in% c("ABC", "GCB"))
)

# =========================================
# Publication-quality COO comparison
# =========================================

library(ggplot2)

ggplot(
  pheno,
  aes(
    x = molecular_COO,
    y = TherapyPersistenceScore,
    fill = molecular_COO
  )
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  theme_classic(base_size = 14) +
  labs(
    title = "Therapy Persistence Score by Molecular COO",
    x = "Molecular COO subtype",
    y = "Therapy Persistence Score"
  )


aggregate(
  TherapyPersistenceScore ~ molecular_COO,
  data = pheno,
  FUN = mean
)

aggregate(
  TherapyPersistenceScore ~ molecular_COO,
  data = pheno,
  FUN = median
)

kruskal.test(
  TherapyPersistenceScore ~ molecular_COO,
  data = pheno
)

pairwise.wilcox.test(
  pheno$TherapyPersistenceScore,
  pheno$molecular_COO,
  p.adjust.method = "BH"
)
