# ============================================================
# GSE193566 — TPS12 vs Blood Advances 30-gene relapse signature
# ============================================================

library(tidyverse)
library(ggpubr)
library(pheatmap)

set.seed(123)

outdir <- "results/GSE193566_TPS12_vs_relapse_signature"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. Load expression matrix
# -----------------------------
expr_file <- "data/raw/GSE193566/GSE193566_all_batchCorrectedCounts.txt.gz"

expr_raw <- read.delim(
  "data/raw/GSE193566/GSE193566_all_batchCorrectedCounts.txt.gz",
  check.names = FALSE
)

# Fix blank first-column name
colnames(expr_raw)[1] <- "gene_symbol"

expr_raw$gene_symbol <- as.character(expr_raw$gene_symbol)

expr <- expr_raw %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  column_to_rownames("gene_symbol") %>%
  as.matrix()

mode(expr) <- "numeric"

expr_log <- expr

cat("Expression matrix:", nrow(expr_log), "genes x", ncol(expr_log), "samples\n")
# -----------------------------
# 2. Define signatures
# -----------------------------
TPS12 <- c(
  "MKI67","RRM2","CHEK1","FOXP1","TCL1A",
  "MCM2","MCM3","MCM4","MCM5","MCM6",
  "PCNA","TYMS"
)

relapse_weights <- c(
  ANKRD6   =  0.534454,
  CASP7    =  0.788655,
  CCDC141  = -0.027820,
  CNIH4    =  0.236360,
  DHRS3    = -0.002580,
  EGLN1    =  0.524393,
  FRK      =  0.147857,
  KCNQ5    =  0.309994,
  KEL      = -0.117840,
  KLRB1    =  0.156239,
  MYC      =  0.445559,
  MYEF2    = -0.156940,
  PDE4DIP  =  0.094356,
  PEAK1    =  0.259342,
  PGRMC2   =  0.270161,
  PSD3     = -0.008670,
  PTK2     =  0.163725,
  PYGL     =  0.044382,
  SCP2     = -0.258140,
  SEC31A   = -0.003390,
  SIRPA    = -0.301470,
  SLAIN2   =  0.023591,
  SNORA71A = -0.224750,
  SNORA71C =  0.256894,
  SSBP2    =  0.207043,
  TCN2     =  0.107918,
  TNFRSF9  = -0.163550,
  TPM1     =  0.036197,
  UBB      =  0.068660,
  WIPI1    =  0.408337
)

RELAPSE30 <- names(relapse_weights)

# -----------------------------
# 3. Scoring functions
# -----------------------------
score_mean_z <- function(mat, genes, score_name) {
  genes_use <- intersect(genes, rownames(mat))
  message(score_name, ": ", length(genes_use), "/", length(genes), " genes found")
  
  if (length(genes_use) < 3) {
    stop(score_name, " has too few genes detected.")
  }
  
  z <- t(scale(t(mat[genes_use, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

score_weighted <- function(mat, weights, score_name) {
  genes_use <- intersect(names(weights), rownames(mat))
  message(score_name, ": ", length(genes_use), "/", length(weights), " genes found")
  
  if (length(genes_use) < 10) {
    stop(score_name, " has too few genes detected.")
  }
  
  z <- t(scale(t(mat[genes_use, , drop = FALSE])))
  weights_use <- weights[genes_use]
  
  as.numeric(crossprod(weights_use, z))
}

scores <- tibble(
  sample_id = colnames(expr_log),
  TPS12 = score_mean_z(expr_log, TPS12, "TPS12"),
  Relapse30_unweighted = score_mean_z(expr_log, RELAPSE30, "Relapse30 unweighted"),
  Relapse30_weighted = score_weighted(expr_log, relapse_weights, "Relapse30 weighted")
)

# -----------------------------
# 4. Infer diagnosis/relapse pairing
# -----------------------------
scores <- scores %>%
  mutate(
    timepoint = case_when(
      str_detect(sample_id, "D$") ~ "Diagnosis",
      str_detect(sample_id, "R$") ~ "Relapse",
      TRUE ~ NA_character_
    ),
    patient_id = case_when(
      str_detect(sample_id, "D$|R$") ~ str_remove(sample_id, "[DR]$"),
      TRUE ~ sample_id
    )
  )

print(table(scores$timepoint, useNA = "ifany"))

scores %>%
  count(patient_id, timepoint) %>%
  filter(n > 1 | is.na(timepoint)) %>%
  print(n = 100)

# -----------------------------
# 5. Paired diagnosis-relapse table
# -----------------------------
paired <- scores %>%
  filter(!is.na(timepoint)) %>%
  select(patient_id, timepoint, TPS12, Relapse30_unweighted, Relapse30_weighted) %>%
  pivot_wider(
    names_from = timepoint,
    values_from = c(TPS12, Relapse30_unweighted, Relapse30_weighted)
  ) %>%
  drop_na(
    TPS12_Diagnosis, TPS12_Relapse,
    Relapse30_unweighted_Diagnosis, Relapse30_unweighted_Relapse,
    Relapse30_weighted_Diagnosis, Relapse30_weighted_Relapse
  ) %>%
  mutate(
    delta_TPS12 = TPS12_Relapse - TPS12_Diagnosis,
    delta_Relapse30_unweighted = Relapse30_unweighted_Relapse - Relapse30_unweighted_Diagnosis,
    delta_Relapse30_weighted = Relapse30_weighted_Relapse - Relapse30_weighted_Diagnosis
  )

write_csv(paired, file.path(outdir, "GSE193566_TPS12_Relapse30_paired.csv"))

# -----------------------------
# 6. Statistics
# -----------------------------
safe_wilcox <- function(x, y) {
  wilcox.test(x, y, paired = TRUE, exact = FALSE)
}

safe_cor <- function(x, y) {
  cor.test(x, y, method = "spearman", exact = FALSE)
}

w_tps <- safe_wilcox(paired$TPS12_Relapse, paired$TPS12_Diagnosis)
w_rel_unw <- safe_wilcox(paired$Relapse30_unweighted_Relapse, paired$Relapse30_unweighted_Diagnosis)
w_rel_w <- safe_wilcox(paired$Relapse30_weighted_Relapse, paired$Relapse30_weighted_Diagnosis)

c_all_unw <- safe_cor(scores$TPS12, scores$Relapse30_unweighted)
c_all_w <- safe_cor(scores$TPS12, scores$Relapse30_weighted)

c_delta_unw <- safe_cor(paired$delta_TPS12, paired$delta_Relapse30_unweighted)
c_delta_w <- safe_cor(paired$delta_TPS12, paired$delta_Relapse30_weighted)

results <- tibble(
  analysis = c(
    "TPS12 relapse vs diagnosis",
    "Relapse30 unweighted relapse vs diagnosis",
    "Relapse30 weighted relapse vs diagnosis",
    "TPS12 vs Relapse30 unweighted, all samples",
    "TPS12 vs Relapse30 weighted, all samples",
    "Delta TPS12 vs Delta Relapse30 unweighted",
    "Delta TPS12 vs Delta Relapse30 weighted"
  ),
  statistic = c(
    unname(w_tps$statistic),
    unname(w_rel_unw$statistic),
    unname(w_rel_w$statistic),
    unname(c_all_unw$estimate),
    unname(c_all_w$estimate),
    unname(c_delta_unw$estimate),
    unname(c_delta_w$estimate)
  ),
  p_value = c(
    w_tps$p.value,
    w_rel_unw$p.value,
    w_rel_w$p.value,
    c_all_unw$p.value,
    c_all_w$p.value,
    c_delta_unw$p.value,
    c_delta_w$p.value
  )
)

write_csv(results, file.path(outdir, "GSE193566_TPS12_Relapse30_summary_statistics.csv"))
print(results)

# -----------------------------
# 7. Plots
# -----------------------------
plot_df <- scores %>%
  filter(!is.na(timepoint)) %>%
  mutate(timepoint = factor(timepoint, levels = c("Diagnosis", "Relapse")))

p1 <- ggplot(plot_df, aes(timepoint, Relapse30_weighted)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.7) +
  theme_classic(base_size = 12) +
  labs(
    title = "Weighted relapse-discriminator score in GSE193566",
    x = NULL,
    y = "Weighted Relapse30 score"
  )

ggsave(file.path(outdir, "Relapse30_weighted_diagnosis_vs_relapse.png"),
       p1, width = 4.5, height = 4, dpi = 300)

p2 <- ggplot(scores, aes(TPS12, Relapse30_weighted)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 12) +
  labs(
    title = "TPS12 vs weighted relapse-discriminator score",
    x = "TPS12 score",
    y = "Weighted Relapse30 score"
  )

ggsave(file.path(outdir, "TPS12_vs_Relapse30_weighted.png"),
       p2, width = 4.5, height = 4, dpi = 300)

p3 <- ggplot(paired, aes(delta_TPS12, delta_Relapse30_weighted)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 12) +
  labs(
    title = "Paired change: TPS12 vs relapse-discriminator",
    x = "Relapse - diagnosis TPS12",
    y = "Relapse - diagnosis weighted Relapse30"
  )

ggsave(file.path(outdir, "Delta_TPS12_vs_Delta_Relapse30_weighted.png"),
       p3, width = 4.5, height = 4, dpi = 300)

# -----------------------------
# 8. Session info
# -----------------------------
sink(file.path(outdir, "sessionInfo.txt"))
sessionInfo()
sink()