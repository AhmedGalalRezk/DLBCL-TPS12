# ============================================================
# GSE193566 — TPS12 dissection and relapse-evolution analysis
# Focus: Is TPS12 more than generic proliferation?
# ============================================================

library(tidyverse)
library(ggpubr)
library(patchwork)

set.seed(123)

outdir <- "results/GSE193566_TPS12_dissection"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Load expression matrix
# ============================================================

expr_file <- "data/raw/GSE193566/GSE193566_all_batchCorrectedCounts.txt.gz"

expr_raw <- read.delim(expr_file, check.names = FALSE)

# First column contains gene symbols but has blank name
colnames(expr_raw)[1] <- "gene_symbol"
expr_raw$gene_symbol <- as.character(expr_raw$gene_symbol)

expr <- expr_raw %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  column_to_rownames("gene_symbol") %>%
  as.matrix()

mode(expr) <- "numeric"

# IMPORTANT:
# This file is batch-corrected and contains negative values.
# Do NOT log2-transform.
expr_mat <- expr

cat("Expression matrix:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n")

# ============================================================
# 2. Define gene programs
# ============================================================

TPS12 <- c(
  "MKI67", "RRM2", "CHEK1", "FOXP1", "TCL1A",
  "MCM2", "MCM3", "MCM4", "MCM5", "MCM6",
  "PCNA", "TYMS"
)

TPS12_proliferation <- c(
  "MKI67", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6"
)

TPS12_replication_stress <- c(
  "CHEK1", "RRM2", "TYMS"
)

TPS12_bcell_regulatory <- c(
  "FOXP1", "TCL1A"
)

# Published Blood Advances weighted 30-gene ABC-DLBCL relapse discriminator
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

# ============================================================
# 3. Scoring functions
# ============================================================

score_mean_z <- function(mat, genes, score_name) {
  genes_use <- intersect(genes, rownames(mat))
  message(score_name, ": ", length(genes_use), "/", length(genes), " genes found")
  
  if (length(genes_use) < 2) {
    stop(score_name, " has too few genes detected.")
  }
  
  z <- t(scale(t(mat[genes_use, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

score_weighted_z <- function(mat, weights, score_name) {
  genes_use <- intersect(names(weights), rownames(mat))
  message(score_name, ": ", length(genes_use), "/", length(weights), " genes found")
  
  if (length(genes_use) < 10) {
    stop(score_name, " has too few genes detected.")
  }
  
  z <- t(scale(t(mat[genes_use, , drop = FALSE])))
  weights_use <- weights[genes_use]
  
  as.numeric(crossprod(weights_use, z))
}

safe_cor <- function(x, y) {
  ok <- complete.cases(x, y)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  tibble(
    n = sum(ok),
    rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

safe_wilcox <- function(x, y) {
  ok <- complete.cases(x, y)
  wt <- suppressWarnings(wilcox.test(x[ok], y[ok], paired = TRUE, exact = FALSE))
  tibble(
    n = sum(ok),
    statistic = unname(wt$statistic),
    p_value = wt$p.value
  )
}

# ============================================================
# 4. Score samples
# ============================================================

scores <- tibble(
  sample_id = colnames(expr_mat),
  TPS12 = score_mean_z(expr_mat, TPS12, "TPS12"),
  Proliferation_component = score_mean_z(expr_mat, TPS12_proliferation, "TPS12 proliferation component"),
  Replication_stress_component = score_mean_z(expr_mat, TPS12_replication_stress, "TPS12 replication-stress component"),
  Bcell_regulatory_component = score_mean_z(expr_mat, TPS12_bcell_regulatory, "TPS12 B-cell regulatory component"),
  Relapse30_unweighted = score_mean_z(expr_mat, RELAPSE30, "Relapse30 unweighted"),
  Relapse30_weighted = score_weighted_z(expr_mat, relapse_weights, "Relapse30 weighted")
)

# ============================================================
# 5. Assign diagnosis/relapse labels
# ============================================================

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

cat("\nTimepoint table:\n")
print(table(scores$timepoint, useNA = "ifany"))

cat("\nUnclassified samples:\n")
scores %>%
  filter(is.na(timepoint)) %>%
  select(sample_id) %>%
  print(n = 50)

write_csv(scores, file.path(outdir, "GSE193566_TPS12_dissection_sample_scores.csv"))

# ============================================================
# 6. Create paired table
# ============================================================

paired <- scores %>%
  filter(!is.na(timepoint)) %>%
  select(
    patient_id, timepoint,
    TPS12,
    Proliferation_component,
    Replication_stress_component,
    Bcell_regulatory_component,
    Relapse30_unweighted,
    Relapse30_weighted
  ) %>%
  group_by(patient_id, timepoint) %>%
  summarise(
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = timepoint,
    values_from = c(
      TPS12,
      Proliferation_component,
      Replication_stress_component,
      Bcell_regulatory_component,
      Relapse30_unweighted,
      Relapse30_weighted
    )
  ) %>%
  drop_na(TPS12_Diagnosis, TPS12_Relapse) %>%
  mutate(
    delta_TPS12 = TPS12_Relapse - TPS12_Diagnosis,
    delta_Proliferation = Proliferation_component_Relapse - Proliferation_component_Diagnosis,
    delta_Replication_stress = Replication_stress_component_Relapse - Replication_stress_component_Diagnosis,
    delta_Bcell_regulatory = Bcell_regulatory_component_Relapse - Bcell_regulatory_component_Diagnosis,
    delta_Relapse30_weighted = Relapse30_weighted_Relapse - Relapse30_weighted_Diagnosis
  )

cat("\nNumber of paired cases used:", nrow(paired), "\n")

write_csv(paired, file.path(outdir, "GSE193566_TPS12_dissection_paired_scores.csv"))

# ============================================================
# 7. Paired diagnosis-relapse tests
# ============================================================

paired_tests <- bind_rows(
  safe_wilcox(paired$TPS12_Relapse, paired$TPS12_Diagnosis) %>%
    mutate(analysis = "TPS12 relapse vs diagnosis"),
  
  safe_wilcox(paired$Proliferation_component_Relapse, paired$Proliferation_component_Diagnosis) %>%
    mutate(analysis = "Proliferation component relapse vs diagnosis"),
  
  safe_wilcox(paired$Replication_stress_component_Relapse, paired$Replication_stress_component_Diagnosis) %>%
    mutate(analysis = "Replication-stress component relapse vs diagnosis"),
  
  safe_wilcox(paired$Bcell_regulatory_component_Relapse, paired$Bcell_regulatory_component_Diagnosis) %>%
    mutate(analysis = "B-cell regulatory component relapse vs diagnosis"),
  
  safe_wilcox(paired$Relapse30_weighted_Relapse, paired$Relapse30_weighted_Diagnosis) %>%
    mutate(analysis = "Published weighted Relapse30 relapse vs diagnosis")
) %>%
  select(analysis, n, statistic, p_value)

write_csv(paired_tests, file.path(outdir, "GSE193566_TPS12_dissection_paired_tests.csv"))

print(paired_tests)

# ============================================================
# 8. Correlation analyses
# ============================================================

cor_results <- bind_rows(
  safe_cor(scores$TPS12, scores$Proliferation_component) %>%
    mutate(analysis = "TPS12 vs proliferation component"),
  
  safe_cor(scores$TPS12, scores$Replication_stress_component) %>%
    mutate(analysis = "TPS12 vs replication-stress component"),
  
  safe_cor(scores$TPS12, scores$Bcell_regulatory_component) %>%
    mutate(analysis = "TPS12 vs B-cell regulatory component"),
  
  safe_cor(scores$Proliferation_component, scores$Replication_stress_component) %>%
    mutate(analysis = "Proliferation vs replication-stress component"),
  
  safe_cor(scores$TPS12, scores$Relapse30_weighted) %>%
    mutate(analysis = "TPS12 vs published weighted Relapse30"),
  
  safe_cor(scores$Replication_stress_component, scores$Relapse30_weighted) %>%
    mutate(analysis = "Replication-stress component vs published weighted Relapse30"),
  
  safe_cor(scores$Proliferation_component, scores$Relapse30_weighted) %>%
    mutate(analysis = "Proliferation component vs published weighted Relapse30"),
  
  safe_cor(paired$delta_TPS12, paired$delta_Relapse30_weighted) %>%
    mutate(analysis = "Delta TPS12 vs delta weighted Relapse30"),
  
  safe_cor(paired$delta_Replication_stress, paired$delta_Relapse30_weighted) %>%
    mutate(analysis = "Delta replication-stress vs delta weighted Relapse30"),
  
  safe_cor(paired$delta_Proliferation, paired$delta_Relapse30_weighted) %>%
    mutate(analysis = "Delta proliferation vs delta weighted Relapse30")
) %>%
  select(analysis, n, rho, p_value) %>%
  mutate(FDR = p.adjust(p_value, method = "BH"))

write_csv(cor_results, file.path(outdir, "GSE193566_TPS12_dissection_correlations.csv"))

print(cor_results)

# ============================================================
# 9. Main plots
# ============================================================

plot_df <- scores %>%
  filter(!is.na(timepoint)) %>%
  mutate(timepoint = factor(timepoint, levels = c("Diagnosis", "Relapse")))

p_tps <- ggplot(plot_df, aes(timepoint, TPS12)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.75) +
  theme_classic(base_size = 12) +
  labs(
    title = "TPS12 is not enriched at relapse",
    x = NULL,
    y = "TPS12 score"
  )

p_prolif <- ggplot(plot_df, aes(timepoint, Proliferation_component)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.75) +
  theme_classic(base_size = 12) +
  labs(
    title = "Proliferation component",
    x = NULL,
    y = "Score"
  )

p_repl <- ggplot(plot_df, aes(timepoint, Replication_stress_component)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.75) +
  theme_classic(base_size = 12) +
  labs(
    title = "Replication-stress component",
    x = NULL,
    y = "Score"
  )

p_rel30 <- ggplot(plot_df, aes(timepoint, Relapse30_weighted)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.75) +
  theme_classic(base_size = 12) +
  labs(
    title = "Published weighted Relapse30",
    x = NULL,
    y = "Score"
  )

combined_box <- (p_tps | p_prolif) / (p_repl | p_rel30)

ggsave(
  file.path(outdir, "GSE193566_TPS12_components_diagnosis_vs_relapse.png"),
  combined_box,
  width = 8.5,
  height = 7,
  dpi = 300
)

# ============================================================
# 10. Correlation plots
# ============================================================

p_corr1 <- ggplot(scores, aes(Proliferation_component, Replication_stress_component)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 12) +
  labs(
    title = "Replication-stress component vs proliferation component",
    x = "Proliferation component",
    y = "Replication-stress component"
  )

p_corr2 <- ggplot(scores, aes(TPS12, Relapse30_weighted)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 12) +
  labs(
    title = "TPS12 vs published relapse discriminator",
    x = "TPS12 score",
    y = "Weighted Relapse30 score"
  )

p_corr3 <- ggplot(paired, aes(delta_TPS12, delta_Relapse30_weighted)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 12) +
  labs(
    title = "Paired change: TPS12 vs relapse discriminator",
    x = "Relapse - diagnosis TPS12",
    y = "Relapse - diagnosis Relapse30"
  )

ggsave(
  file.path(outdir, "GSE193566_Proliferation_vs_ReplicationStress.png"),
  p_corr1,
  width = 5,
  height = 4.5,
  dpi = 300
)

ggsave(
  file.path(outdir, "GSE193566_TPS12_vs_Relapse30_weighted.png"),
  p_corr2,
  width = 5,
  height = 4.5,
  dpi = 300
)

ggsave(
  file.path(outdir, "GSE193566_Delta_TPS12_vs_Delta_Relapse30_weighted.png"),
  p_corr3,
  width = 5,
  height = 4.5,
  dpi = 300
)

# ============================================================
# 11. Export manuscript-ready summary
# ============================================================

summary_df <- bind_rows(
  paired_tests %>%
    mutate(type = "paired diagnosis-relapse test") %>%
    rename(metric = statistic) %>%
    select(type, analysis, n, metric, p_value),
  
  cor_results %>%
    mutate(type = "correlation") %>%
    rename(metric = rho) %>%
    select(type, analysis, n, metric, p_value)
)

write_csv(summary_df, file.path(outdir, "GSE193566_TPS12_dissection_summary.csv"))

sink(file.path(outdir, "sessionInfo.txt"))
sessionInfo()
sink()

cat("\nDone. Results saved to:", outdir, "\n")