# ============================================================
# Supplementary Figure — Gene-level longitudinal analysis
# of the TPS12 B-cell regulatory module
# GSE193566 paired diagnosis–relapse cohort
#
# A: FOXP1 diagnosis vs relapse
# B: TCL1A diagnosis vs relapse
# C: Longitudinal changes in FOXP1 and TCL1A
# D: Delta FOXP1 vs delta B-cell regulatory module
# ============================================================

library(tidyverse)
library(patchwork)

set.seed(123)

# ============================================================
# 1. Files and output directory
# ============================================================

inputdir <- "results/GSE193566_FOXP1_TCL1A_remodeling"
outdir  <- "results/TPS12_additional_figures"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

paired_file <- file.path(
  inputdir,
  "GSE193566_paired_FOXP1_TCL1A.csv"
)

tests_file <- file.path(
  inputdir,
  "GSE193566_paired_tests_FOXP1_TCL1A.csv"
)

if (!file.exists(paired_file)) {
  stop("Paired file not found: ", paired_file)
}

if (!file.exists(tests_file)) {
  stop("Paired-test file not found: ", tests_file)
}

paired_gene <- read_csv(
  paired_file,
  show_col_types = FALSE
)

paired_tests_original <- read_csv(
  tests_file,
  show_col_types = FALSE
)

cat("Paired cases loaded:", nrow(paired_gene), "\n")

# ============================================================
# 2. Confirm required columns
# ============================================================

required_columns <- c(
  "patient_id",
  "FOXP1_Diagnosis",
  "FOXP1_Relapse",
  "TCL1A_Diagnosis",
  "TCL1A_Relapse",
  "Bcell_regulatory_Diagnosis",
  "Bcell_regulatory_Relapse",
  "delta_FOXP1",
  "delta_TCL1A",
  "delta_Bcell_regulatory"
)

missing_columns <- setdiff(
  required_columns,
  colnames(paired_gene)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns:\n",
    paste(missing_columns, collapse = "\n")
  )
}

# ============================================================
# 3. Color palette
# ============================================================

timepoint_colors <- c(
  "Diagnosis" = "#0072B2",
  "Relapse"   = "#D55E00"
)

gene_colors <- c(
  "FOXP1" = "#CC79A7",
  "TCL1A" = "#009E73"
)

paired_line_color  <- "#8C8C8C"
correlation_color  <- "#7B3294"
reference_line_col <- "#777777"

# ============================================================
# 4. Helper functions
# ============================================================

fmt_p <- function(p) {
  case_when(
    is.na(p)  ~ "P = NA",
    p < 0.001 ~ "P < 0.001",
    TRUE ~ paste0(
      "P = ",
      formatC(
        p,
        digits = 2,
        format = "fg"
      )
    )
  )
}

fmt_median <- function(x) {
  paste0(
    "Median \u0394 = ",
    formatC(
      x,
      digits = 2,
      format = "fg"
    )
  )
}

paired_wilcox <- function(diagnosis, relapse) {
  
  keep <- is.finite(diagnosis) & is.finite(relapse)
  
  if (sum(keep) < 3) {
    return(
      list(
        n = sum(keep),
        statistic = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  result <- wilcox.test(
    relapse[keep],
    diagnosis[keep],
    paired = TRUE,
    exact = FALSE
  )
  
  list(
    n = sum(keep),
    statistic = unname(result$statistic),
    p_value = result$p.value
  )
}

# ============================================================
# 5. Recalculate paired tests
# ============================================================

test_foxp1 <- paired_wilcox(
  paired_gene$FOXP1_Diagnosis,
  paired_gene$FOXP1_Relapse
)

test_tcl1a <- paired_wilcox(
  paired_gene$TCL1A_Diagnosis,
  paired_gene$TCL1A_Relapse
)

test_regulatory <- paired_wilcox(
  paired_gene$Bcell_regulatory_Diagnosis,
  paired_gene$Bcell_regulatory_Relapse
)

gene_statistics <- tibble(
  analysis = c(
    "FOXP1 diagnosis vs relapse",
    "TCL1A diagnosis vs relapse",
    "B-cell regulatory module diagnosis vs relapse"
  ),
  n = c(
    test_foxp1$n,
    test_tcl1a$n,
    test_regulatory$n
  ),
  statistic = c(
    test_foxp1$statistic,
    test_tcl1a$statistic,
    test_regulatory$statistic
  ),
  p_value = c(
    test_foxp1$p_value,
    test_tcl1a$p_value,
    test_regulatory$p_value
  )
) %>%
  mutate(
    FDR_within_3_tests = p.adjust(
      p_value,
      method = "BH"
    )
  )

print(gene_statistics)

# ============================================================
# 6. Long-format paired expression data
# ============================================================

gene_long <- bind_rows(
  
  paired_gene %>%
    transmute(
      patient_id,
      Diagnosis = FOXP1_Diagnosis,
      Relapse = FOXP1_Relapse,
      gene = "FOXP1"
    ),
  
  paired_gene %>%
    transmute(
      patient_id,
      Diagnosis = TCL1A_Diagnosis,
      Relapse = TCL1A_Relapse,
      gene = "TCL1A"
    )
  
) %>%
  pivot_longer(
    cols = c(Diagnosis, Relapse),
    names_to = "timepoint",
    values_to = "expression"
  ) %>%
  mutate(
    gene = factor(
      gene,
      levels = c("FOXP1", "TCL1A")
    ),
    timepoint = factor(
      timepoint,
      levels = c("Diagnosis", "Relapse")
    )
  ) %>%
  filter(
    !is.na(patient_id),
    is.finite(expression)
  )

# ============================================================
# 7. Longitudinal delta data
# ============================================================

delta_long <- paired_gene %>%
  transmute(
    patient_id,
    FOXP1 = delta_FOXP1,
    TCL1A = delta_TCL1A
  ) %>%
  pivot_longer(
    cols = c(FOXP1, TCL1A),
    names_to = "gene",
    values_to = "delta_expression"
  ) %>%
  mutate(
    gene = factor(
      gene,
      levels = c("FOXP1", "TCL1A")
    )
  ) %>%
  filter(is.finite(delta_expression))

delta_regulatory_df <- paired_gene %>%
  transmute(
    patient_id,
    delta_FOXP1,
    delta_Bcell_regulatory
  ) %>%
  filter(
    is.finite(delta_FOXP1),
    is.finite(delta_Bcell_regulatory)
  )

# ============================================================
# 8. Correlation analysis
# ============================================================

cor_foxp1_regulatory <- cor.test(
  delta_regulatory_df$delta_FOXP1,
  delta_regulatory_df$delta_Bcell_regulatory,
  method = "spearman",
  exact = FALSE
)

cat(
  "\nDelta FOXP1 vs delta B-cell regulatory module:\n",
  "Spearman rho =",
  unname(cor_foxp1_regulatory$estimate),
  "\nP =",
  cor_foxp1_regulatory$p.value,
  "\n"
)

# ============================================================
# 9. Panel C annotation table
# ============================================================

delta_test_labels <- tibble(
  gene = factor(
    c("FOXP1", "TCL1A"),
    levels = c("FOXP1", "TCL1A")
  ),
  p_value = c(
    test_foxp1$p_value,
    test_tcl1a$p_value
  )
)

delta_summary <- delta_long %>%
  group_by(gene) %>%
  summarise(
    median_delta = median(
      delta_expression,
      na.rm = TRUE
    ),
    ymin = min(
      delta_expression,
      na.rm = TRUE
    ),
    ymax = max(
      delta_expression,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  left_join(
    delta_test_labels,
    by = "gene"
  ) %>%
  mutate(
    span = if_else(
      ymax > ymin,
      ymax - ymin,
      1
    ),
    y = ymax + 0.10 * span,
    median_label = fmt_median(median_delta),
    p_label = fmt_p(p_value),
    label = paste(
      median_label,
      p_label,
      sep = "\n"
    )
  )

# ============================================================
# 10. Figure theme
# ============================================================

theme_supp <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust = 0
    ),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    plot.tag = element_text(
      face = "bold",
      size = 16
    ),
    legend.position = "none",
    plot.margin = margin(8, 8, 8, 8)
  )

# ============================================================
# 11. Paired gene-plot function
# ============================================================

make_gene_plot <- function(
    gene_name,
    p_value
) {
  
  plot_data <- gene_long %>%
    filter(gene == gene_name)
  
  y_range <- range(
    plot_data$expression,
    na.rm = TRUE
  )
  
  y_span <- diff(y_range)
  
  if (!is.finite(y_span) || y_span == 0) {
    y_span <- 1
  }
  
  y_annotation <- y_range[2] + 0.10 * y_span
  
  ggplot(
    plot_data,
    aes(
      x = timepoint,
      y = expression,
      group = patient_id
    )
  ) +
    geom_line(
      color = paired_line_color,
      alpha = 0.60,
      linewidth = 0.38
    ) +
    geom_point(
      aes(color = timepoint),
      size = 2.1,
      alpha = 0.92
    ) +
    geom_boxplot(
      aes(
        group = timepoint,
        fill = timepoint,
        color = timepoint
      ),
      width = 0.42,
      outlier.shape = NA,
      alpha = 0.18,
      linewidth = 0.65
    ) +
    annotate(
      "text",
      x = 1.5,
      y = y_annotation,
      label = fmt_p(p_value),
      size = 3.8
    ) +
    coord_cartesian(
      ylim = c(
        y_range[1] - 0.05 * y_span,
        y_range[2] + 0.18 * y_span
      )
    ) +
    scale_color_manual(
      values = timepoint_colors
    ) +
    scale_fill_manual(
      values = timepoint_colors
    ) +
    theme_supp +
    labs(
      title = gene_name,
      x = NULL,
      y = "Gene-expression score"
    )
}

# ============================================================
# 12. Panel A — FOXP1
# ============================================================

pA <- make_gene_plot(
  gene_name = "FOXP1",
  p_value = test_foxp1$p_value
)

# ============================================================
# 13. Panel B — TCL1A
# ============================================================

pB <- make_gene_plot(
  gene_name = "TCL1A",
  p_value = test_tcl1a$p_value
)

# ============================================================
# 14. Panel C — Longitudinal gene changes
# ============================================================

pC <- ggplot(
  delta_long,
  aes(
    x = gene,
    y = delta_expression,
    fill = gene,
    color = gene
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = reference_line_col
  ) +
  geom_boxplot(
    width = 0.46,
    outlier.shape = NA,
    alpha = 0.24,
    linewidth = 0.65
  ) +
  geom_jitter(
    width = 0.10,
    height = 0,
    size = 2.0,
    alpha = 0.85
  ) +
  geom_text(
    data = delta_summary,
    aes(
      x = gene,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    color = "black",
    size = 3.4,
    lineheight = 0.95
  ) +
  scale_color_manual(
    values = gene_colors
  ) +
  scale_fill_manual(
    values = gene_colors
  ) +
  theme_supp +
  labs(
    title = "Longitudinal regulatory-gene changes",
    x = NULL,
    y = expression(
      Delta * " gene-expression score"
    )
  )

# ============================================================
# 15. Panel D — Delta FOXP1 vs regulatory module
# ============================================================

pD <- ggplot(
  delta_regulatory_df,
  aes(
    x = delta_FOXP1,
    y = delta_Bcell_regulatory
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.40,
    color = reference_line_col
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.40,
    color = reference_line_col
  ) +
  geom_point(
    size = 2.4,
    alpha = 0.88,
    color = correlation_color
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.95,
    color = correlation_color
  ) +
  annotate(
    "text",
    x = min(
      delta_regulatory_df$delta_FOXP1,
      na.rm = TRUE
    ),
    y = max(
      delta_regulatory_df$delta_Bcell_regulatory,
      na.rm = TRUE
    ),
    hjust = 0,
    vjust = 1,
    label = paste0(
      "Spearman \u03c1 = ",
      round(
        unname(cor_foxp1_regulatory$estimate),
        2
      ),
      "\n",
      fmt_p(
        cor_foxp1_regulatory$p.value
      )
    ),
    size = 3.8
  ) +
  theme_supp +
  labs(
    title = "FOXP1 contribution to regulatory-module change",
    x = expression(
      Delta * "FOXP1"
    ),
    y = expression(
      Delta * "B-cell regulatory module"
    )
  )

# ============================================================
# 16. Combine figure
# ============================================================

Supplementary_FOXP1_TCL1A <- (
  pA | pB
) / (
  pC | pD
) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 16
    )
  )

# ============================================================
# 17. Save figure
# ============================================================

tiff_file <- file.path(
  outdir,
  "Supplementary_Figure_Regulatory_Module_Gene_Analysis.tiff"
)

png_file <- file.path(
  outdir,
  "Supplementary_Figure_Regulatory_Module_Gene_Analysis_preview.png"
)

ggsave(
  filename = tiff_file,
  plot = Supplementary_FOXP1_TCL1A,
  device = "tiff",
  width = 10,
  height = 9,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = png_file,
  plot = Supplementary_FOXP1_TCL1A,
  width = 10,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

# ============================================================
# 18. Export figure statistics
# ============================================================

correlation_statistics <- tibble(
  analysis = "Delta FOXP1 vs delta B-cell regulatory module",
  n = nrow(delta_regulatory_df),
  statistic = unname(
    cor_foxp1_regulatory$estimate
  ),
  p_value = cor_foxp1_regulatory$p.value,
  FDR_within_3_tests = NA_real_
)

figure_statistics <- bind_rows(
  gene_statistics,
  correlation_statistics
)

statistics_file <- file.path(
  outdir,
  "Supplementary_Figure_Regulatory_Module_Gene_Analysis_statistics.csv"
)

write_csv(
  figure_statistics,
  statistics_file
)

# ============================================================
# 19. Export plotting values
# ============================================================

plot_values_file <- file.path(
  outdir,
  "Supplementary_Figure_Regulatory_Module_Gene_Analysis_values.csv"
)

write_csv(
  paired_gene %>%
    select(
      patient_id,
      FOXP1_Diagnosis,
      FOXP1_Relapse,
      TCL1A_Diagnosis,
      TCL1A_Relapse,
      Bcell_regulatory_Diagnosis,
      Bcell_regulatory_Relapse,
      delta_FOXP1,
      delta_TCL1A,
      delta_Bcell_regulatory
    ),
  plot_values_file
)

# ============================================================
# 20. Final output messages
# ============================================================

cat("\nSupplementary TIFF saved to:\n")
cat(tiff_file, "\n")

cat("\nPNG preview saved to:\n")
cat(png_file, "\n")

cat("\nStatistics saved to:\n")
cat(statistics_file, "\n")

cat("\nPlotting values saved to:\n")
cat(plot_values_file, "\n")