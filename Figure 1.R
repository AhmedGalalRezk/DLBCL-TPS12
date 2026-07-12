# ============================================================
# FIGURE 1 ONLY — TIFF VERSION
# Option A: adds Panel D showing refinement to TPS12
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(forcats)
library(ggforce)

set.seed(123)

dir.create("results/main_figures", recursive = TRUE, showWarnings = FALSE)

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

# ----------------------------
# Load object
# ----------------------------

obj <- readRDS("data/processed/GSE182434_ecosystem_scored.rds")
obj$cell_annotation <- as.character(Idents(obj))

# ----------------------------
# Rename states for manuscript
# ----------------------------

obj$figure_cell_state <- obj$cell_annotation

obj$figure_cell_state[obj$cell_annotation %in% c("Cycling_B", "Cycling_B_2")] <-
  "Proliferative B"

obj$figure_cell_state[obj$cell_annotation == "Highly_Cycling_B"] <-
  "Stress-Adaptive\nProliferative B"

obj$figure_cell_state[obj$cell_annotation %in% c("B_cells_IFN", "B_cells_STAT1")] <-
  "IFN-Responsive B"

obj$figure_cell_state[obj$cell_annotation %in% c("B_cells", "B_cells_2")] <-
  "B cells"

obj$figure_cell_state[obj$cell_annotation == "Plasma_like"] <-
  "Plasma-like"

obj$figure_cell_state[obj$cell_annotation == "Plasma_IFN"] <-
  "IFN-Responsive\nPlasma-like"

obj$figure_cell_state[obj$cell_annotation == "Plasma_CTSS"] <-
  "Antigen-Presenting\nPlasma-like"

obj$figure_cell_state[obj$cell_annotation %in% c("Plasmablast", "Plasmablast_2")] <-
  "Plasmablasts"

obj$figure_cell_state[obj$cell_annotation %in% c(
  "T_cells", "T_cells_2", "CD4_T", "CD4_T_IFN", "CD4_T_IFN_2",
  "T_IFN", "T_IFN_2", "T_cells_IFN"
)] <- "T cells"

obj$figure_cell_state[obj$cell_annotation %in% c("Activated_T_NK", "T_NK")] <-
  "Activated T/NK"

obj$figure_cell_state[obj$cell_annotation == "NK_cytotoxic"] <-
  "Cytotoxic NK"

obj$figure_cell_state[obj$cell_annotation == "NK_cells_IFN"] <-
  "IFN-Responsive NK"

obj$figure_cell_state[obj$cell_annotation == "Inflammatory_Myeloid"] <-
  "Inflammatory Myeloid"

# ----------------------------
# TPS12 score
# ----------------------------

TPS12_genes <- c(
  "MKI67", "RRM2", "CHEK1", "FOXP1", "TCL1A",
  "MCM2", "MCM3", "MCM4", "MCM5", "MCM6",
  "PCNA", "TYMS"
)

TPS12_present <- intersect(TPS12_genes, rownames(obj))

if (!"TPS121" %in% colnames(obj@meta.data)) {
  obj <- AddModuleScore(
    obj,
    features = list(TPS12_present),
    name = "TPS12"
  )
}

obj$TPS12_score <- obj$TPS121

# ----------------------------
# Panel A: annotated UMAP
# ----------------------------

p1a <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "figure_cell_state",
  label = TRUE,
  repel = TRUE,
  label.size = 2.1,
  raster = TRUE
) +
  labs(
    title = "Annotated DLBCL single-cell ecosystem",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Cell state"
  ) +
  theme_blood() +
  theme(legend.position = "none")

# ----------------------------
# Panel B: annotation marker DotPlot
# ----------------------------

annotation_markers <- c(
  "MKI67", "TOP2A", "CDK1",
  "RRM2", "TYMS", "CHEK1", "BIRC5",
  "MCM2", "MCM3", "MCM4", "MCM5", "MCM6",
  "STAT1", "ISG15", "IFI6", "IFITM3",
  "JCHAIN", "MZB1", "SDC1", "CTSS", "HLA-DRA", "CD74"
)

annotation_markers <- intersect(annotation_markers, rownames(obj))

dot_states <- c(
  "Stress-Adaptive\nProliferative B",
  "Proliferative B",
  "IFN-Responsive B",
  "B cells",
  "Plasma-like",
  "IFN-Responsive\nPlasma-like",
  "Antigen-Presenting\nPlasma-like",
  "Plasmablasts",
  "T cells",
  "Activated T/NK",
  "Cytotoxic NK",
  "IFN-Responsive NK",
  "Inflammatory Myeloid"
)

dot_states <- intersect(dot_states, unique(obj$figure_cell_state))

obj_dot <- subset(obj, subset = figure_cell_state %in% dot_states)
obj_dot$figure_cell_state <- factor(obj_dot$figure_cell_state, levels = rev(dot_states))
Idents(obj_dot) <- obj_dot$figure_cell_state

p1b <- DotPlot(
  obj_dot,
  features = annotation_markers,
  group.by = "figure_cell_state",
  dot.scale = 4.3
) +
  RotatedAxis() +
  scale_color_gradient2(
    low = "#4E79A7",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  labs(
    title = "Marker programs defining B-cell states",
    x = "Annotation markers",
    y = NULL,
    color = "Average\nexpression",
    size = "Percent\nexpressed"
  ) +
  theme_blood() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7)
  )

# ----------------------------
# Panel C: TPS12 localization violin
# ----------------------------

cell_summary <- obj@meta.data %>%
  group_by(figure_cell_state) %>%
  summarise(
    n_cells = n(),
    median_TPS12 = median(TPS12_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_TPS12))

top_states <- cell_summary %>%
  filter(n_cells >= 20) %>%
  slice_head(n = 12) %>%
  pull(figure_cell_state)

p1c <- obj@meta.data %>%
  filter(figure_cell_state %in% top_states) %>%
  mutate(
    figure_cell_state = fct_reorder(
      figure_cell_state,
      TPS12_score,
      .fun = median
    )
  ) %>%
  ggplot(aes(x = figure_cell_state, y = TPS12_score)) +
  geom_violin(
    fill = "grey85",
    color = "grey30",
    linewidth = 0.2,
    scale = "width"
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    linewidth = 0.25
  ) +
  coord_flip() +
  labs(
    title = "TPS12 localization across cell states",
    x = NULL,
    y = "TPS12 module score"
  ) +
  theme_blood()

# ----------------------------
# Panel D: TPS12 refinement schematic
# ----------------------------

refinement_df <- data.frame(
  x = c(1, 2.5, 4, 5.5),
  y = c(1, 1, 1, 1),
  label = c(
    "Transcriptome-wide\nsingle-cell discovery",
    "Proliferative / stress-\nadaptive B-cell state",
    "Candidate persistence\nprogram",
    "Locked TPS12\npanel"
  )
)

arrow_df <- data.frame(
  x = c(1.45, 2.95, 4.45),
  xend = c(2.05, 3.55, 5.05),
  y = c(1, 1, 1),
  yend = c(1, 1, 1)
)

TPS12_label <- paste(
  "TPS12 genes:",
  "MKI67, RRM2, CHEK1, FOXP1, TCL1A,",
  "MCM2–6, PCNA, TYMS",
  sep = "\n"
)

p1d <- ggplot() +
  geom_label(
    data = refinement_df,
    aes(x = x, y = y, label = label),
    size = 2.7,
    label.size = 0.25,
    fill = "grey95"
  ) +
  geom_segment(
    data = arrow_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.18, "cm")),
    linewidth = 0.35
  ) +
  annotate(
    "text",
    x = 5.5,
    y = 0.55,
    label = TPS12_label,
    size = 2.4,
    hjust = 0.5
  ) +
  xlim(0.4, 6.1) +
  ylim(0.35, 1.3) +
  labs(title = "Refinement of the TPS12 panel") +
  theme_void(base_size = 8) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0)
  )

# ----------------------------
# Assemble 4-panel figure
# ----------------------------

fig1 <- (p1a | p1d) / (p1b | p1c) +
  plot_annotation(tag_levels = "A")

# ----------------------------
# Save TIFF
# ----------------------------

ggsave(
  filename = "results/main_figures/Figure_1_single_cell_ecosystem_TPS12_4panel.tiff",
  plot = fig1,
  width = 12,
  height = 9,
  units = "in",
  dpi = 600,
  compression = "lzw"
)