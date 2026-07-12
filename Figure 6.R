# ============================================================
# FIGURE 6
# Therapeutic vulnerability programs associated with TPS-high DLBCL
#
# Input:
#   results/GSE10846/patient_level_continuous_vulnerability_axes.csv
#
# Outputs:
#   results/main_figures/Figure_6_TPS_therapeutic_vulnerability_overlay.pdf
#   results/main_figures/Figure_6_TPS_therapeutic_vulnerability_overlay.tiff
#   results/main_figures/Table_Figure6_TPS_vulnerability_wilcoxon.csv
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(stringr)

set.seed(123)

dir.create(
  "results/main_figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Theme
# ============================================================

theme_blood <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size + 1),
      axis.text = element_text(size = base_size),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# ============================================================
# Load patient-level vulnerability scores
# ============================================================

vuln_path <- "results/GSE10846/patient_level_continuous_vulnerability_axes.csv"

if (!file.exists(vuln_path)) {
  stop("Missing file: ", vuln_path)
}

vuln <- read.csv(
  vuln_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c(
  "geo_accession",
  "TherapyPersistenceScore",
  "TherapyPersistenceGroup",
  "RT_DNA_Damage",
  "Chemo_Resistance",
  "CAR_T_Bispecific_Visibility",
  "Polivy_ADC_Target",
  "Targeted_ABC_Biology",
  "Dominant_axis"
)

missing_cols <- setdiff(required_cols, colnames(vuln))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

vuln$TherapyPersistenceGroup <- factor(
  vuln$TherapyPersistenceGroup,
  levels = c("Low", "High")
)

# ============================================================
# Display labels
# ============================================================

axis_labels <- c(
  "RT_DNA_Damage" = "RT/DNA damage program",
  "Chemo_Resistance" = "Chemoresistance program",
  "CAR_T_Bispecific_Visibility" = "Immune target visibility program",
  "Polivy_ADC_Target" = "ADC target-expression program",
  "Targeted_ABC_Biology" = "Targetable ABC biology"
)

axis_order <- unname(axis_labels)

dominant_labels <- c(
  "RT_DNA_Damage" = "RT/DNA damage",
  "Chemo_Resistance" = "Chemoresistance",
  "CAR_T_Bispecific_Visibility" = "Immune target visibility",
  "Polivy_ADC_Target" = "ADC target expression",
  "Targeted_ABC_Biology" = "Targetable ABC biology"
)

dominant_order <- c(
  "RT/DNA damage",
  "Chemoresistance",
  "Immune target visibility",
  "ADC target expression",
  "Targetable ABC biology"
)

# ============================================================
# Long table
# ============================================================

vuln_long <- vuln %>%
  dplyr::select(
    geo_accession,
    TherapyPersistenceScore,
    TherapyPersistenceGroup,
    all_of(names(axis_labels))
  ) %>%
  tidyr::pivot_longer(
    cols = all_of(names(axis_labels)),
    names_to = "Axis_raw",
    values_to = "Score"
  ) %>%
  dplyr::mutate(
    Axis = unname(axis_labels[Axis_raw]),
    Axis = factor(Axis, levels = axis_order)
  )

# ============================================================
# Statistics: TPS-high vs TPS-low
# ============================================================

stats_df <- vuln_long %>%
  dplyr::group_by(Axis, Axis_raw) %>%
  dplyr::summarise(
    n_low = sum(TherapyPersistenceGroup == "Low" & !is.na(Score)),
    n_high = sum(TherapyPersistenceGroup == "High" & !is.na(Score)),
    median_low = median(Score[TherapyPersistenceGroup == "Low"], na.rm = TRUE),
    median_high = median(Score[TherapyPersistenceGroup == "High"], na.rm = TRUE),
    mean_low = mean(Score[TherapyPersistenceGroup == "Low"], na.rm = TRUE),
    mean_high = mean(Score[TherapyPersistenceGroup == "High"], na.rm = TRUE),
    delta_median_high_minus_low = median_high - median_low,
    delta_mean_high_minus_low = mean_high - mean_low,
    p_value = tryCatch(
      wilcox.test(Score ~ TherapyPersistenceGroup)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    FDR = p.adjust(p_value, method = "BH"),
    FDR_label = dplyr::case_when(
      is.na(FDR) ~ "FDR = NA",
      FDR < 0.001 ~ "FDR < 0.001",
      TRUE ~ paste0("FDR = ", signif(FDR, 2))
    ),
    Significant = ifelse(FDR < 0.05, "FDR < 0.05", "NS")
  )

write.csv(
  stats_df,
  "results/main_figures/Table_Figure6_TPS_vulnerability_wilcoxon.csv",
  row.names = FALSE
)

# ============================================================
# Panel A: TPS-high vs TPS-low distributions
# ============================================================

stats_labels <- stats_df %>%
  dplyr::mutate(
    label_y = sapply(
      as.character(Axis),
      function(a) {
        max(vuln_long$Score[vuln_long$Axis == a], na.rm = TRUE) + 0.08
      }
    )
  )

pA <- ggplot(
  vuln_long,
  aes(
    x = TherapyPersistenceGroup,
    y = Score,
    fill = TherapyPersistenceGroup
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    color = "grey30",
    linewidth = 0.25
  ) +
  geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.25
  ) +
  geom_text(
    data = stats_labels,
    aes(
      x = 1.5,
      y = label_y,
      label = FDR_label
    ),
    inherit.aes = FALSE,
    size = 2.4
  ) +
  facet_wrap(
    ~ Axis,
    ncol = 3,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c("Low" = "grey70", "High" = "#B2182B")
  ) +
  labs(
    title = "Therapeutic vulnerability programs by TPS group",
    x = "TPS group",
    y = "Module score"
  ) +
  theme_blood(base_size = 8) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 7, face = "bold")
  )

# ============================================================
# Panel B: Effect-size summary
# ============================================================

stats_plot <- stats_df %>%
  dplyr::mutate(
    Axis = factor(Axis, levels = rev(axis_order))
  )

pB <- ggplot(
  stats_plot,
  aes(
    x = delta_median_high_minus_low,
    y = Axis,
    fill = Significant
  )
) +
  geom_col(width = 0.7) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.35
  ) +
  geom_text(
    aes(label = FDR_label),
    hjust = ifelse(stats_plot$delta_median_high_minus_low >= 0, -0.05, 1.05),
    size = 2.5
  ) +
  scale_fill_manual(
    values = c("FDR < 0.05" = "#B2182B", "NS" = "grey70")
  ) +
  coord_cartesian(
    xlim = c(
      min(stats_plot$delta_median_high_minus_low, na.rm = TRUE) - 0.10,
      max(stats_plot$delta_median_high_minus_low, na.rm = TRUE) + 0.45
    )
  ) +
  labs(
    title = "TPS-high enrichment effect size",
    x = "Median difference: TPS-high minus TPS-low",
    y = NULL,
    fill = NULL
  ) +
  theme_blood(base_size = 8) +
  theme(
    legend.position = "right"
  )

# ============================================================
# Panel C: Dominant vulnerability axis among TPS-high tumors
# ============================================================

dominant_df <- vuln %>%
  dplyr::filter(TherapyPersistenceGroup == "High") %>%
  dplyr::mutate(
    Dominant_axis_label = unname(dominant_labels[Dominant_axis]),
    Dominant_axis_label = ifelse(
      is.na(Dominant_axis_label),
      Dominant_axis,
      Dominant_axis_label
    )
  ) %>%
  dplyr::group_by(Dominant_axis_label) %>%
  dplyr::summarise(
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Dominant_axis_label = factor(
      Dominant_axis_label,
      levels = rev(dominant_order)
    ),
    percent = 100 * n / sum(n),
    label = paste0(n, " (", round(percent), "%)")
  )

pC <- ggplot(
  dominant_df,
  aes(
    x = n,
    y = Dominant_axis_label
  )
) +
  geom_col(
    fill = "#B2182B",
    width = 0.7
  ) +
  geom_text(
    aes(label = label),
    hjust = -0.05,
    size = 2.6
  ) +
  coord_cartesian(
    xlim = c(0, max(dominant_df$n, na.rm = TRUE) * 1.35),
    clip = "off"
  ) +
  labs(
    title = "Dominant vulnerability axis among TPS-high tumors",
    x = "TPS-high patients (%)",
    y = NULL
  ) +
  theme_blood(base_size = 8) +
  theme(
    plot.margin = margin(5.5, 35, 5.5, 5.5)
  )

# ============================================================
# Assemble Figure 6
# ============================================================

fig6 <- pA /
  (pB | pC) +
  plot_layout(
    heights = c(1.35, 1),
    widths = c(1, 1.25)
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 10, face = "bold")
  )

# ============================================================
# Save Figure 6
# ============================================================

ggsave(
  filename = "results/main_figures/Figure_6_TPS_therapeutic_vulnerability_overlay.pdf",
  plot = fig6,
  width = 12,
  height = 8.5,
  units = "in",
  limitsize = FALSE
)

ggsave(
  filename = "results/main_figures/Figure_6_TPS_therapeutic_vulnerability_overlay.tiff",
  plot = fig6,
  width = 12,
  height = 8.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  limitsize = FALSE
)

message("Figure 6 saved successfully.")
message("Statistics saved to results/main_figures/Table_Figure6_TPS_vulnerability_wilcoxon.csv")