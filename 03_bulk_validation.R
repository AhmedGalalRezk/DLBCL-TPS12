# =============================
# Bulk validation of single-cell-derived persistence signature
# Dataset: GSE10846
# =============================

# Install once if needed:
# install.packages("BiocManager")
# BiocManager::install(c("GEOquery", "limma"))
# install.packages(c("survival", "survminer", "dplyr", "ggplot2"))

library(GEOquery)
library(limma)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)

# =============================
# Create folders
# =============================

dir.create("data/bulk", recursive = TRUE, showWarnings = FALSE)
dir.create("results/bulk_validation", recursive = TRUE, showWarnings = FALSE)

# =============================
# Define single-cell-derived signature
# =============================

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

# =============================
# Download GSE10846
# =============================

gse <- getGEO(
  "GSE10846",
  GSEMatrix = TRUE,
  getGPL = TRUE
)

gse <- gse[[1]]

expr <- exprs(gse)
pheno <- pData(gse)
features <- fData(gse)

# =============================
# Inspect metadata
# =============================

dim(expr)
dim(pheno)
head(pheno)
colnames(pheno)
head(features)

write.csv(
  pheno,
  "results/bulk_validation/GSE10846_clinical_metadata_raw.csv",
  row.names = TRUE
)

write.csv(
  features,
  "results/bulk_validation/GSE10846_feature_annotation_raw.csv",
  row.names = TRUE
)

# =============================
# Map probes to gene symbols
# =============================

colnames(features)

# Common symbol columns may include:
# "Gene Symbol", "GENE_SYMBOL", "Symbol", or similar.
# Check the feature table above and edit this if needed.

symbol_col <- grep(
  "symbol|gene.symbol|gene symbol",
  colnames(features),
  ignore.case = TRUE,
  value = TRUE
)[1]

symbol_col

features$gene_symbol <- features[[symbol_col]]

# Remove missing symbols
valid <- !is.na(features$gene_symbol) & features$gene_symbol != ""

expr_valid <- expr[valid, ]
features_valid <- features[valid, ]

# Some probes map to multiple symbols separated by ///.
# Keep first symbol for now.
features_valid$gene_symbol <- sapply(
  strsplit(as.character(features_valid$gene_symbol), "///"),
  function(x) trimws(x[1])
)

# Collapse probes to gene level by mean expression
expr_gene <- avereps(
  expr_valid,
  ID = features_valid$gene_symbol
)

# =============================
# Score persistence signature
# =============================

sig_present <- persistence_signature[persistence_signature %in% rownames(expr_gene)]
sig_missing <- setdiff(persistence_signature, sig_present)

sig_present
sig_missing

expr_sig <- expr_gene[sig_present, , drop = FALSE]

# z-score each gene across patients
expr_sig_z <- t(scale(t(expr_sig)))

# patient-level score = mean z-score of signature genes
signature_score <- colMeans(expr_sig_z, na.rm = TRUE)

pheno$PersistenceSignatureScore <- signature_score[rownames(pheno)]

# Define high vs low by median
pheno$PersistenceSignatureGroup <- ifelse(
  pheno$PersistenceSignatureScore >= median(pheno$PersistenceSignatureScore, na.rm = TRUE),
  "High",
  "Low"
)

table(pheno$PersistenceSignatureGroup)

# =============================
# Save scored metadata
# =============================

write.csv(
  pheno,
  "results/bulk_validation/GSE10846_metadata_with_persistence_score.csv",
  row.names = TRUE
)

# =============================
# Parse survival metadata
# =============================

extract_field <- function(x, field) {
  sub(
    paste0(".*", field, ": ([^;]+).*"),
    "\\1",
    x
  )
}

pheno$follow_up_status <- extract_field(
  pheno$`Clinical info:ch1`,
  "Follow up status"
)

pheno$follow_up_years <- as.numeric(
  extract_field(
    pheno$`Clinical info:ch1`,
    "Follow up years"
  )
)

pheno$survival_event <- ifelse(
  pheno$follow_up_status == "DEAD",
  1,
  0
)

table(pheno$follow_up_status)
summary(pheno$follow_up_years)
# =============================
# Kaplan-Meier survival analysis
# =============================

fit <- survfit(
  Surv(follow_up_years, survival_event) ~ PersistenceSignatureGroup,
  data = pheno
)

ggsurvplot(
  fit,
  data = pheno,
  pval = TRUE,
  risk.table = TRUE
)

# =============================
# Cox proportional hazards model
# =============================

cox_model <- coxph(
  Surv(follow_up_years, survival_event) ~ PersistenceSignatureGroup,
  data = pheno
)

summary(cox_model)