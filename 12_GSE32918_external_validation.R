# =========================================
# GSE32918 external survival validation
# TPS / reduced panel
# =========================================

library(GEOquery)
library(limma)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

dir.create("results/GSE32918", recursive = TRUE, showWarnings = FALSE)

# Download
gse <- getGEO("GSE32918", GSEMatrix = TRUE)
gse <- gse[[1]]

expr <- exprs(gse)
pheno <- pData(gse)
features <- fData(gse)

# Inspect survival / annotation columns
colnames(pheno)
grep(
  "surv|death|event|status|follow|time|os|pfs|efs|dfs|progress|months|years",
  colnames(pheno),
  ignore.case = TRUE,
  value = TRUE
)

# Find gene symbol column
colnames(features)

symbol_col <- grep(
  "symbol|gene.symbol|gene symbol",
  colnames(features),
  ignore.case = TRUE,
  value = TRUE
)[1]

symbol_col

features$gene_symbol <- features[[symbol_col]]

valid <- !is.na(features$gene_symbol) & features$gene_symbol != ""

expr_valid <- expr[valid, ]
features_valid <- features[valid, ]

features_valid$gene_symbol <- sapply(
  strsplit(as.character(features_valid$gene_symbol), " /// "),
  function(x) trimws(x[1])
)

expr_gene <- avereps(
  expr_valid,
  ID = features_valid$gene_symbol
)

# =========================================
# Locked 12-gene TPS panel
# =========================================

TPS12 <- c(
  "MKI67","RRM2","CHEK1","FOXP1","TCL1A",
  "MCM2","MCM3","MCM4","MCM5","MCM6",
  "PCNA","TYMS"
)

present <- intersect(TPS12, rownames(expr_gene))
missing <- setdiff(TPS12, present)

present
missing

expr_panel <- expr_gene[present, , drop = FALSE]
expr_panel_z <- t(scale(t(expr_panel)))

pheno$TPS12_Score <- colMeans(expr_panel_z, na.rm = TRUE)

pheno$TPS12_Group <- ifelse(
  pheno$TPS12_Score >= median(pheno$TPS12_Score, na.rm = TRUE),
  "High",
  "Low"
)

table(pheno$TPS12_Group)

write.csv(
  pheno,
  "results/GSE32918/GSE32918_metadata_with_TPS12.csv",
  row.names = TRUE
)