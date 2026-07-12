# ==========================================================
# FIGURE 3
# External single-cell validation of TPS-associated biology
#
# Updates:
#   - Shorter Panel A title to avoid clipping
#   - Keeps DLBCL_1–DLBCL_6 labels for cross-panel traceability
#   - Orders heatmap rows from component programs to composite TPS
#   - Gives Panel A more vertical space
#   - Saves both PDF and TIFF
# ==========================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(scales)

set.seed(123)

dir.create("results/main_figures", recursive = TRUE, showWarnings = FALSE)

# ==========================================================
# Theme
# ==========================================================

theme_blood <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size + 1),
      axis.text = element_text(size = base_size),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# ==========================================================
# Load external validation data
# ==========================================================

ext_obj_path <- "data/processed/LBCL_external_malignant_therapy_persistence.rds"
summary_path <- "results/external_validation/LBCL_external_therapy_persistence_summary.csv"

if (!file.exists(ext_obj_path)) {
  stop("Missing file: ", ext_obj_path)
}

if (!file.exists(summary_path)) {
  stop("Missing file: ", summary_path)
}

ext_obj <- readRDS(ext_obj_path)

external_summary <- read.csv(
  summary_path,
  stringsAsFactors = FALSE
)

sample_levels <- c(
  "DLBCL_1",
  "DLBCL_2",
  "DLBCL_3",
  "DLBCL_4",
  "DLBCL_5",
  "DLBCL_6"
)

patient_map <- c(
  DLBCL_1 = "Patient 1",
  DLBCL_2 = "Patient 2",
  DLBCL_3 = "Patient 3",
  DLBCL_4 = "Patient 4",
  DLBCL_5 = "Patient 5",
  DLBCL_6 = "Patient 6"
)

external_summary$patient_id <- patient_map[
  as.character(external_summary$sample_id)
]

external_summary$patient_id <- factor(
  external_summary$patient_id,
  levels = paste("Patient", 1:6)
)
# ==========================================================
# Panel A: Program heatmap
# ==========================================================

heat_df <- external_summary %>%
  select(
    sample_id,
    mean_replication_stress,
    mean_DNA_damage_checkpoint,
    mean_proteostasis_stress,
    mean_immune_visibility,
    mean_therapy_persistence
  ) %>%
  pivot_longer(
    cols = -sample_id,
    names_to = "Program",
    values_to = "Score"
  )

heat_df$Program <- factor(
  heat_df$Program,
  levels = c(
    "mean_replication_stress",
    "mean_DNA_damage_checkpoint",
    "mean_proteostasis_stress",
    "mean_immune_visibility",
    "mean_therapy_persistence"
  ),
  labels = c(
    "Replication stress",
    "DNA damage checkpoint",
    "Proteostasis stress",
    "Immune visibility",
    "Composite TPS"
  )
)

heat_df <- heat_df %>%
  group_by(Program) %>%
  mutate(
    z = as.numeric(scale(Score))
  ) %>%
  ungroup()

pA <- ggplot(
  heat_df,
  aes(x = sample_id, y = Program, fill = z)
) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    name = "Mean score"
  ) +
  labs(
    title = "External single-cell validation of TPS programs",
    x = "Sample",
    y = NULL
  ) +
  theme_blood(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(size = 12, face = "bold", hjust = 0)
  )

# ==========================================================
# Panel B: Composite TPS by sample
# ==========================================================

pB <- ggplot(
  external_summary,
  aes(
    x = reorder(sample_id, n_cells),
    y = mean_therapy_persistence
  )
) +
  geom_col(fill = "#B2182B", width = 0.75) +
  coord_flip() +
  labs(
    title = "Sample-level composite TPS",
    x = NULL,
    y = "Mean composite TPS"
  ) +
  theme_blood(base_size = 9)

# ==========================================================
# Panel C: Malignant cell counts by sample
# ==========================================================

pC <- ggplot(
  external_summary,
  aes(
    x = reorder(patient_id, n_cells),
    y = n_cells
  )
) +
  geom_col(fill = "grey55", width = 0.75) +
  coord_flip() +
  labs(
    title = "Malignant cells per sample",
    x = NULL,
    y = "Number of malignant cells"
  ) +
  theme_blood(base_size = 9)

# ==========================================================
# Assemble figure
# ==========================================================

fig3 <- pA /
  (pB | pC) +
  plot_layout(
    heights = c(1.25, 1)
  ) +
  plot_annotation(
    tag_levels = "A"
  )

# ==========================================================
# Save PDF and TIFF
# ==========================================================

ggsave(
  filename = "results/main_figures/Figure_3_external_single_cell_validation_REVISED.pdf",
  plot = fig3,
  width = 9.5,
  height = 6,
  units = "in"
)

ggsave(
  filename = "results/main_figures/Figure_3_external_single_cell_validation_REVISED.tiff",
  plot = fig3,
  width = 9.5,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

message("Figure 3 saved successfully as PDF and TIFF in results/main_figures/")