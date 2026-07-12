# ============================================================
# GSE193566 — TPS12 vs relapse-associated biology
# Diagnosis-relapse DLBCL paired cohort
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(tibble)
library(patchwork)

set.seed(123)

outdir <- "results/GSE193566_TPS12_relapse_biology"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

theme_blood <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size + 1),
      axis.text = element_text(size = base_size),
      plot.title = element_text(size = base_size + 2, face = "bold"),
      plot.subtitle = element_text(size = base_size),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

safe_z <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  if (sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

score_module <- function(expr_mat, genes) {
  genes_present <- intersect(genes, rownames(expr_mat))
  genes_missing <- setdiff(genes, rownames(expr_mat))
  
  message("Genes present: ", paste(genes_present, collapse = ", "))
  message("Genes missing: ", paste(genes_missing, collapse = ", "))
  
  if (length(genes_present) < 2) {
    warning("Too few genes present for this module.")
    return(rep(NA_real_, ncol(expr_mat)))
  }
  
  zmat <- t(scale(t(expr_mat[genes_present, , drop = FALSE])))
  colMeans(zmat, na.rm = TRUE)
}

# ============================================================
# Load expression matrix
# ============================================================

expr_file <- "data/raw/GSE193566/GSE193566_all_batchCorrectedCounts.txt"

expr_df <- read.delim(
  expr_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(expr_df)[1] <- "GeneSymbol"

expr_gene <- expr_df %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
  group_by(GeneSymbol) %>%
  summarise(
    across(
      where(is.numeric),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  column_to_rownames("GeneSymbol") %>%
  as.matrix()

# ============================================================
# Infer sample metadata from sample IDs
# ============================================================

sample_meta <- data.frame(
  sample_id = colnames(expr_gene),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Timepoint = case_when(
      grepl("(_d$|D$)", sample_id) ~ "Diagnosis",
      grepl("(R$|R[0-9]+$)", sample_id) ~ "Relapse",
      TRUE ~ NA_character_
    ),
    Patient_ID = sample_id %>%
      str_replace("(_d$|D$)", "") %>%
      str_replace("(R[0-9]+$|R$)", "")
  )

write.csv(
  sample_meta,
  file.path(outdir, "GSE193566_inferred_sample_metadata.csv"),
  row.names = FALSE
)

print(table(sample_meta$Timepoint, useNA = "ifany"))

# ============================================================
# Define gene modules
# ============================================================

TPS12 <- c(
  "MKI67", "RRM2", "CHEK1", "FOXP1", "TCL1A",
  "MCM2", "MCM3", "MCM4", "MCM5", "MCM6",
  "PCNA", "TYMS"
)

DNA_REPAIR_LIKE <- c(
  "BRCA1", "BRCA2", "RAD51", "CHEK1", "CHEK2",
  "ATR", "ATM", "TOPBP1", "RPA1", "RPA2",
  "FANCD2", "MRE11", "NBN"
)

G2M_LIKE <- c(
  "MKI67", "TOP2A", "UBE2C", "BIRC5", "CENPF",
  "CDK1", "CCNB1", "CCNA2", "AURKA", "TPX2"
)

MYC_LIKE <- c(
  "MYC", "NCL", "NPM1", "HSPD1", "HSPE1",
  "LDHA", "ENO1", "GAPDH", "PKM", "TPI1"
)

# IMPORTANT:
# Replace this placeholder with the exact 30-gene ABC relapse discriminator
# from the Mareschal / Blood Advances supplementary table.
# Do NOT use this placeholder for final manuscript claims.
MARESCHAL_30GENE <- c(
  "MYC"
  # paste the remaining 29 genes here once extracted
)

modules <- list(
  "TPS12" = TPS12,
  "DNA repair-like" = DNA_REPAIR_LIKE,
  "G2M-like" = G2M_LIKE,
  "MYC-like" = MYC_LIKE,
  "Mareschal 30-gene relapse score" = MARESCHAL_30GENE
)

# ============================================================
# Score modules
# ============================================================

score_df <- data.frame(
  sample_id = colnames(expr_gene),
  stringsAsFactors = FALSE
)

for (nm in names(modules)) {
  clean_nm <- make.names(nm)
  message("\nScoring: ", nm)
  score_df[[clean_nm]] <- score_module(expr_gene, modules[[nm]])
}

analysis_df <- score_df %>%
  left_join(sample_meta, by = "sample_id") %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("Diagnosis", "Relapse"))
  )

write.csv(
  analysis_df,
  file.path(outdir, "GSE193566_TPS12_relapse_biology_scores.csv"),
  row.names = FALSE
)

# ============================================================
# Paired diagnosis-relapse changes
# ============================================================

score_cols <- names(score_df)[names(score_df) != "sample_id"]

paired_long <- analysis_df %>%
  filter(!is.na(Patient_ID), !is.na(Timepoint)) %>%
  select(Patient_ID, Timepoint, all_of(score_cols)) %>%
  pivot_longer(
    cols = all_of(score_cols),
    names_to = "Program",
    values_to = "Score"
  ) %>%
  group_by(Patient_ID, Timepoint, Program) %>%
  summarise(
    Score = mean(Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Patient_ID, Program) %>%
  filter(all(c("Diagnosis", "Relapse") %in% as.character(Timepoint))) %>%
  ungroup()

paired_wide <- paired_long %>%
  pivot_wider(
    names_from = Timepoint,
    values_from = Score
  ) %>%
  mutate(
    Delta_Relapse_minus_Diagnosis = Relapse - Diagnosis
  )

paired_tests <- paired_wide %>%
  group_by(Program) %>%
  summarise(
    n_pairs = n(),
    median_diagnosis = median(Diagnosis, na.rm = TRUE),
    median_relapse = median(Relapse, na.rm = TRUE),
    median_delta = median(Delta_Relapse_minus_Diagnosis, na.rm = TRUE),
    wilcoxon_p = wilcox.test(
      Relapse,
      Diagnosis,
      paired = TRUE,
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(wilcoxon_p, method = "BH")
  )

print(paired_tests)

write.csv(
  paired_wide,
  file.path(outdir, "GSE193566_paired_program_changes_wide.csv"),
  row.names = FALSE
)

write.csv(
  paired_tests,
  file.path(outdir, "GSE193566_paired_program_change_tests.csv"),
  row.names = FALSE
)

# ============================================================
# Correlations with TPS12
# ============================================================

cor_df <- analysis_df %>%
  select(sample_id, Patient_ID, Timepoint, all_of(score_cols))

tps_col <- "TPS12"

cor_tests <- setdiff(score_cols, tps_col) %>%
  lapply(function(prog) {
    tmp <- cor_df %>%
      filter(!is.na(.data[[tps_col]]), !is.na(.data[[prog]]))
    
    ct <- cor.test(
      tmp[[tps_col]],
      tmp[[prog]],
      method = "spearman",
      exact = FALSE
    )
    
    data.frame(
      Program = prog,
      n = nrow(tmp),
      spearman_rho = as.numeric(ct$estimate),
      p_value = ct$p.value,
      stringsAsFactors = FALSE
    )
  }) %>%
  bind_rows() %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH")
  )

print(cor_tests)

write.csv(
  cor_tests,
  file.path(outdir, "GSE193566_TPS12_program_correlations.csv"),
  row.names = FALSE
)

# ============================================================
# Plot A: paired TPS12
# ============================================================

p1 <- paired_long %>%
  filter(Program == "TPS12") %>%
  ggplot(
    aes(
      x = Timepoint,
      y = Score,
      group = Patient_ID
    )
  ) +
  geom_line(color = "grey55", linewidth = 0.35, alpha = 0.8) +
  geom_point(
    aes(fill = Timepoint),
    shape = 21,
    size = 2.3,
    color = "grey25",
    stroke = 0.25
  ) +
  scale_fill_manual(
    values = c("Diagnosis" = "grey70", "Relapse" = "#B2182B")
  ) +
  labs(
    title = "TPS12 in paired diagnosis-relapse samples",
    subtitle = paste0(
      "Paired Wilcoxon P = ",
      format.pval(
        paired_tests$wilcoxon_p[paired_tests$Program == "TPS12"],
        digits = 3,
        eps = 1e-4
      )
    ),
    x = NULL,
    y = "TPS12 score"
  ) +
  theme_blood(base_size = 9) +
  theme(legend.position = "none")

# ============================================================
# Plot B: program delta summary
# ============================================================

delta_plot_df <- paired_wide %>%
  mutate(
    Program = factor(
      Program,
      levels = paired_tests$Program[order(paired_tests$median_delta)]
    )
  )

p2 <- ggplot(
  delta_plot_df,
  aes(
    x = Program,
    y = Delta_Relapse_minus_Diagnosis
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_boxplot(outlier.shape = NA, width = 0.55, fill = "grey80") +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.6) +
  coord_flip() +
  labs(
    title = "Relapse-minus-diagnosis program shifts",
    x = NULL,
    y = "Δ score"
  ) +
  theme_blood(base_size = 9)

# ============================================================
# Plot C: TPS12 correlation with DNA repair-like score
# ============================================================

p3 <- ggplot(
  analysis_df,
  aes(
    x = DNA.repair.like,
    y = TPS12
  )
) +
  geom_point(
    aes(fill = Timepoint),
    shape = 21,
    size = 2,
    color = "grey25",
    stroke = 0.25,
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.5,
    color = "#B2182B"
  ) +
  scale_fill_manual(
    values = c("Diagnosis" = "grey70", "Relapse" = "#B2182B")
  ) +
  labs(
    title = "TPS12 tracks DNA repair-like relapse biology",
    x = "DNA repair-like score",
    y = "TPS12 score",
    fill = NULL
  ) +
  theme_blood(base_size = 9) +
  theme(legend.position = "bottom")

fig <- (p1 | p3) / p2 +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

ggsave(
  file.path(outdir, "GSE193566_TPS12_relapse_biology_summary.pdf"),
  fig,
  width = 9,
  height = 7,
  units = "in",
  limitsize = FALSE
)

ggsave(
  file.path(outdir, "GSE193566_TPS12_relapse_biology_summary.tiff"),
  fig,
  width = 9,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw",
  limitsize = FALSE
)

message("GSE193566 TPS12 relapse-biology analysis complete.")