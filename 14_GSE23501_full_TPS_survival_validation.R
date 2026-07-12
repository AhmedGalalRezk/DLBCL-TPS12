# =========================================================
# GSE23501 external survival validation
# Full TPS panel vs reduced TPS12 panel
# =========================================================

library(GEOquery)
library(limma)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(AnnotationDbi)
library(hgu133plus2.db)

dir.create("results/GSE23501_full_TPS", recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Load GSE23501
# =========================================================

gse <- getGEO("GSE23501", GSEMatrix = TRUE)
gse <- gse[[1]]

expr <- exprs(gse)
pheno <- pData(gse)

# =========================================================
# Probe-to-symbol conversion
# =========================================================

probe_ids <- rownames(expr)

probe_to_symbol <- mapIds(
  hgu133plus2.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

valid <- !is.na(probe_to_symbol) & probe_to_symbol != ""

expr_valid <- expr[valid, , drop = FALSE]
symbols_valid <- probe_to_symbol[valid]

expr_gene <- avereps(
  expr_valid,
  ID = symbols_valid
)

# =========================================================
# Align expression and phenotype samples
# =========================================================

common_samples <- intersect(colnames(expr_gene), rownames(pheno))

expr_gene <- expr_gene[, common_samples, drop = FALSE]
pheno <- pheno[common_samples, , drop = FALSE]

dim(expr_gene)
dim(pheno)

# =========================================================
# Define panels
# =========================================================

TPS12 <- c(
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

TPS_full <- c(
  "MKI67",
  "RRM2",
  "CHEK1",
  "FOXP1",
  "TCL1A",
  "TCL1B",
  "HMGB2",
  "FBL",
  "TUBB",
  "PGAM1",
  "NIBAN3",
  "SPIB",
  "MCM2",
  "MCM3",
  "MCM4",
  "MCM5",
  "MCM6",
  "PCNA",
  "TYMS"
)

# =========================================================
# Helper function to score a gene panel
# =========================================================

score_panel <- function(expr_gene, genes, score_name) {
  
  present <- intersect(genes, rownames(expr_gene))
  missing <- setdiff(genes, present)
  
  message(score_name, " present genes: ", paste(present, collapse = ", "))
  message(score_name, " missing genes: ", paste(missing, collapse = ", "))
  
  expr_panel <- expr_gene[present, , drop = FALSE]
  expr_panel_z <- t(scale(t(expr_panel)))
  
  score <- colMeans(expr_panel_z, na.rm = TRUE)
  
  return(list(
    score = score,
    present = present,
    missing = missing
  ))
}

# =========================================================
# Score TPS12 and full TPS
# =========================================================

TPS12_result <- score_panel(
  expr_gene = expr_gene,
  genes = TPS12,
  score_name = "TPS12"
)

TPS_full_result <- score_panel(
  expr_gene = expr_gene,
  genes = TPS_full,
  score_name = "TPS_full"
)

pheno$TPS12_Score <- TPS12_result$score
pheno$TPS_full_Score <- TPS_full_result$score

pheno$TPS12_Group <- ifelse(
  pheno$TPS12_Score >= median(pheno$TPS12_Score, na.rm = TRUE),
  "High",
  "Low"
)

pheno$TPS_full_Group <- ifelse(
  pheno$TPS_full_Score >= median(pheno$TPS_full_Score, na.rm = TRUE),
  "High",
  "Low"
)

pheno$TPS12_Group <- factor(pheno$TPS12_Group, levels = c("Low", "High"))
pheno$TPS_full_Group <- factor(pheno$TPS_full_Group, levels = c("Low", "High"))

table(pheno$TPS12_Group)
table(pheno$TPS_full_Group)

# =========================================================
# Survival variables
# =========================================================

pheno$OS_time <- as.numeric(
  pheno$`overall_survival (years):ch1`
)

pheno$OS_event <- ifelse(
  tolower(pheno$`code_os:ch1`) == "dead",
  1,
  ifelse(tolower(pheno$`code_os:ch1`) == "alive", 0, NA)
)

pheno$PFS_time <- as.numeric(
  pheno$`progression-free_survival (years):ch1`
)

pheno$PFS_event <- ifelse(
  tolower(pheno$`code_pfs:ch1`) == "dead",
  1,
  ifelse(tolower(pheno$`code_pfs:ch1`) == "alive", 0, NA)
)

table(pheno$OS_event, useNA = "ifany")
table(pheno$PFS_event, useNA = "ifany")

# =========================================================
# OS complete cases
# =========================================================

os_df <- pheno %>%
  filter(
    !is.na(OS_time),
    !is.na(OS_event),
    !is.na(TPS12_Score),
    !is.na(TPS_full_Score)
  )

dim(os_df)

# =========================================================
# OS Kaplan-Meier: TPS12
# =========================================================

fit_os_tps12 <- survfit(
  Surv(OS_time, OS_event) ~ TPS12_Group,
  data = os_df
)

p_os_tps12 <- ggsurvplot(
  fit_os_tps12,
  data = os_df,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("steelblue", "firebrick"),
  title = "GSE23501 OS by TPS12",
  xlab = "Overall survival time (years)",
  ylab = "Overall survival probability"
)

print(p_os_tps12)

# =========================================================
# OS Kaplan-Meier: full TPS
# =========================================================

fit_os_full <- survfit(
  Surv(OS_time, OS_event) ~ TPS_full_Group,
  data = os_df
)

p_os_full <- ggsurvplot(
  fit_os_full,
  data = os_df,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("steelblue", "firebrick"),
  title = "GSE23501 OS by full TPS panel",
  xlab = "Overall survival time (years)",
  ylab = "Overall survival probability"
)

print(p_os_full)

pdf(
  "results/GSE23501_full_TPS/GSE23501_OS_TPS12_KM.pdf",
  width = 7,
  height = 6
)
print(p_os_tps12)
dev.off()

pdf(
  "results/GSE23501_full_TPS/GSE23501_OS_full_TPS_KM.pdf",
  width = 7,
  height = 6
)
print(p_os_full)
dev.off()

# =========================================================
# OS Cox models
# =========================================================

cox_os_tps12 <- coxph(
  Surv(OS_time, OS_event) ~ TPS12_Score,
  data = os_df
)

cox_os_full <- coxph(
  Surv(OS_time, OS_event) ~ TPS_full_Score,
  data = os_df
)

summary(cox_os_tps12)
summary(cox_os_full)

# =========================================================
# PFS complete cases
# =========================================================

pfs_df <- pheno %>%
  filter(
    !is.na(PFS_time),
    !is.na(PFS_event),
    !is.na(TPS12_Score),
    !is.na(TPS_full_Score)
  )

dim(pfs_df)

# =========================================================
# PFS Kaplan-Meier: TPS12
# =========================================================

fit_pfs_tps12 <- survfit(
  Surv(PFS_time, PFS_event) ~ TPS12_Group,
  data = pfs_df
)

p_pfs_tps12 <- ggsurvplot(
  fit_pfs_tps12,
  data = pfs_df,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("steelblue", "firebrick"),
  title = "GSE23501 PFS by TPS12",
  xlab = "Progression-free survival time (years)",
  ylab = "Progression-free survival probability"
)

print(p_pfs_tps12)

# =========================================================
# PFS Kaplan-Meier: full TPS
# =========================================================

fit_pfs_full <- survfit(
  Surv(PFS_time, PFS_event) ~ TPS_full_Group,
  data = pfs_df
)

p_pfs_full <- ggsurvplot(
  fit_pfs_full,
  data = pfs_df,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("steelblue", "firebrick"),
  title = "GSE23501 PFS by full TPS panel",
  xlab = "Progression-free survival time (years)",
  ylab = "Progression-free survival probability"
)

print(p_pfs_full)

pdf(
  "results/GSE23501_full_TPS/GSE23501_PFS_TPS12_KM.pdf",
  width = 7,
  height = 6
)
print(p_pfs_tps12)
dev.off()

pdf(
  "results/GSE23501_full_TPS/GSE23501_PFS_full_TPS_KM.pdf",
  width = 7,
  height = 6
)
print(p_pfs_full)
dev.off()

# =========================================================
# PFS Cox models
# =========================================================

cox_pfs_tps12 <- coxph(
  Surv(PFS_time, PFS_event) ~ TPS12_Score,
  data = pfs_df
)

cox_pfs_full <- coxph(
  Surv(PFS_time, PFS_event) ~ TPS_full_Score,
  data = pfs_df
)

summary(cox_pfs_tps12)
summary(cox_pfs_full)

# =========================================================
# Save Cox summary table
# =========================================================

cox_table <- data.frame(
  Dataset = "GSE23501",
  Endpoint = c("OS", "OS", "PFS", "PFS"),
  Panel = c("TPS12", "TPS_full", "TPS12", "TPS_full"),
  HR = c(
    exp(coef(cox_os_tps12)),
    exp(coef(cox_os_full)),
    exp(coef(cox_pfs_tps12)),
    exp(coef(cox_pfs_full))
  ),
  P_value = c(
    summary(cox_os_tps12)$coefficients[,"Pr(>|z|)"],
    summary(cox_os_full)$coefficients[,"Pr(>|z|)"],
    summary(cox_pfs_tps12)$coefficients[,"Pr(>|z|)"],
    summary(cox_pfs_full)$coefficients[,"Pr(>|z|)"]
  )
)

cox_table

write.csv(
  cox_table,
  "results/GSE23501_full_TPS/GSE23501_TPS12_vs_fullTPS_cox_results.csv",
  row.names = FALSE
)

write.csv(
  pheno,
  "results/GSE23501_full_TPS/GSE23501_metadata_TPS12_fullTPS.csv",
  row.names = TRUE
)

saveRDS(
  pheno,
  "results/GSE23501_full_TPS/GSE23501_metadata_TPS12_fullTPS.rds"
)