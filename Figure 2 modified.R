# ============================================================
# FIGURE 2 — CopyKAT malignancy validation and TPS12 enrichment
# Compact 2x2 Blood-style version
#
# Panels:
#   A. CopyKAT-inferred copy-number state UMAP
#   B. TPS12 by inferred copy-number state
#   C. Aneuploid enrichment by state
#   D. Composition of the malignant compartment
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

set.seed(123)

dir.create("results/main_figures", recursive = TRUE, showWarnings = FALSE)

theme_blood <- function(base_size = 9) {
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

safe_as_character <- function(x) {
  if (is.list(x)) x <- unlist(x)
  as.character(x)
}

copykat_path <- "data/processed/GSE182434_Bcell_CopyKAT.rds"

if (!file.exists(copykat_path)) {
  stop("Missing file: ", copykat_path)
}

b_copykat <- readRDS(copykat_path)
b_copykat$cell_annotation <- safe_as_character(Idents(b_copykat))

if (!"copykat_prediction" %in% colnames(b_copykat@meta.data)) {
  copykat_cols <- grep(
    "copykat|prediction|aneuploid|diploid",
    colnames(b_copykat@meta.data),
    value = TRUE,
    ignore.case = TRUE
  )
  
  if (length(copykat_cols) == 0) {
    stop("No CopyKAT prediction column found in metadata.")
  }
  
  b_copykat$copykat_prediction <- b_copykat@meta.data[[copykat_cols[1]]]
}

b_copykat$copykat_prediction <- safe_as_character(b_copykat$copykat_prediction)

b_copykat$copykat_prediction <- dplyr::case_when(
  grepl("aneuploid", b_copykat$copykat_prediction, ignore.case = TRUE) ~ "aneuploid",
  grepl("diploid", b_copykat$copykat_prediction, ignore.case = TRUE) ~ "diploid",
  TRUE ~ "not.defined"
)

b_copykat$copykat_prediction <- factor(
  b_copykat$copykat_prediction,
  levels = c("aneuploid", "diploid", "not.defined")
)

TPS12_genes <- c(
  "MKI67", "RRM2", "CHEK1", "FOXP1", "TCL1A",
  "MCM2", "MCM3", "MCM4", "MCM5", "MCM6",
  "PCNA", "TYMS"
)

TPS12_present <- intersect(TPS12_genes, rownames(b_copykat))

if (length(TPS12_present) < 2) {
  stop("Too few TPS12 genes found in object.")
}

if (!"TPS121" %in% colnames(b_copykat@meta.data)) {
  b_copykat <- AddModuleScore(
    b_copykat,
    features = list(TPS12_present),
    name = "TPS12"
  )
}

b_copykat$TPS12_score <- as.numeric(b_copykat$TPS121)

meta_df <- as.data.frame(b_copykat@meta.data, stringsAsFactors = FALSE)

meta_df$cell_annotation <- safe_as_character(meta_df$cell_annotation)
meta_df$copykat_prediction <- safe_as_character(meta_df$copykat_prediction)
meta_df$TPS12_score <- as.numeric(meta_df$TPS12_score)

meta_df$copykat_prediction <- factor(
  meta_df$copykat_prediction,
  levels = c("aneuploid", "diploid", "not.defined")
)

meta_df$manuscript_state <- dplyr::case_when(
  meta_df$cell_annotation == "Highly_Cycling_B" ~
    "Stress-Adaptive Proliferative B",
  meta_df$cell_annotation %in% c("Cycling_B", "Cycling_B_2") ~
    "Proliferative B",
  meta_df$cell_annotation == "B_cells_IFN" ~
    "IFN-Responsive B",
  meta_df$cell_annotation == "B_cells_STAT1" ~
    "STAT1-high B",
  meta_df$cell_annotation == "Plasma_CTSS" ~
    "Antigen-Presenting Plasma-like B",
  meta_df$cell_annotation == "Plasma_IFN" ~
    "IFN-Responsive Plasma-like B",
  meta_df$cell_annotation == "Plasma_like" ~
    "Plasma-like B",
  TRUE ~ NA_character_
)

state_order <- c(
  "Stress-Adaptive Proliferative B",
  "Proliferative B",
  "IFN-Responsive B",
  "Antigen-Presenting Plasma-like B",
  "Plasma-like B",
  "IFN-Responsive Plasma-like B",
  "STAT1-high B"
)

state_colors <- c(
  "Stress-Adaptive Proliferative B" = "#7CAE00",
  "Proliferative B" = "#F8766D",
  "IFN-Responsive B" = "#00BFC4",
  "STAT1-high B" = "#619CFF",
  "Antigen-Presenting Plasma-like B" = "#A3A500",
  "IFN-Responsive Plasma-like B" = "#C77CFF",
  "Plasma-like B" = "#E76BF3"
)

state_labels_short <- c(
  "Stress-Adaptive Proliferative B" = "Stress-Adaptive\nProlif. B",
  "Proliferative B" = "Proliferative B",
  "IFN-Responsive B" = "IFN-Resp. B",
  "Antigen-Presenting Plasma-like B" = "AP Plasma-like B",
  "Plasma-like B" = "Plasma-like B",
  "IFN-Responsive Plasma-like B" = "IFN-Resp.\nPlasma-like B",
  "STAT1-high B" = "STAT1-high B"
)

b_copykat$manuscript_state <- meta_df$manuscript_state

# ----------------------------
# Panel A
# ----------------------------

pal_copykat <- c(
  aneuploid = "#B2182B",
  diploid = "grey65",
  not.defined = "grey94"
)

pA <- DimPlot(
  b_copykat,
  reduction = "umap",
  group.by = "copykat_prediction",
  cols = pal_copykat,
  raster = TRUE,
  pt.size = 0.45
) +
  coord_cartesian(
    xlim = c(-15, 4),
    ylim = c(-13, 15)
  ) +
  labs(
    title = "CopyKAT-inferred copy-number state",
    x = "umap_1",
    y = "umap_2",
    color = "CopyKAT"
  ) +
  theme_blood(base_size = 9) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.title = element_text(face = "bold")
  )

# ----------------------------
# Panel B
# ----------------------------

copykat_clean <- meta_df %>%
  dplyr::filter(copykat_prediction %in% c("diploid", "aneuploid")) %>%
  dplyr::filter(!is.na(TPS12_score))

p_wilcox <- wilcox.test(
  TPS12_score ~ copykat_prediction,
  data = copykat_clean
)$p.value

p_label <- ifelse(
  p_wilcox < 1e-300,
  "Wilcoxon P < 1e-300",
  paste0("Wilcoxon P = ", format.pval(p_wilcox, digits = 3, eps = 1e-300))
)

pB <- ggplot(
  copykat_clean,
  aes(x = copykat_prediction, y = TPS12_score, fill = copykat_prediction)
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    color = "grey30",
    linewidth = 0.3
  ) +
  geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    linewidth = 0.3,
    fill = "white"
  ) +
  scale_fill_manual(
    values = c(diploid = "grey65", aneuploid = "#B2182B")
  ) +
  labs(
    title = "TPS12 by copy-number state",
    subtitle = p_label,
    x = NULL,
    y = "TPS12 module score"
  ) +
  theme_blood(base_size = 9) +
  theme(legend.position = "none")

# ----------------------------
# Panel C
# ----------------------------

aneuploid_tab <- meta_df %>%
  dplyr::filter(copykat_prediction %in% c("aneuploid", "diploid")) %>%
  dplyr::filter(!is.na(manuscript_state)) %>%
  dplyr::group_by(manuscript_state, copykat_prediction) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(manuscript_state) %>%
  dplyr::mutate(
    total = sum(n),
    fraction = n / total
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(copykat_prediction == "aneuploid") %>%
  dplyr::mutate(
    percent_aneuploid = 100 * fraction,
    manuscript_state = factor(
      manuscript_state,
      levels = rev(state_order)
    )
  )

pC <- ggplot(
  aneuploid_tab,
  aes(x = percent_aneuploid, y = manuscript_state)
) +
  geom_col(fill = "#B2182B", width = 0.7) +
  geom_text(
    aes(label = paste0(round(percent_aneuploid), "%")),
    hjust = -0.15,
    size = 2.7
  ) +
  scale_y_discrete(labels = state_labels_short) +
  scale_x_continuous(
    limits = c(0, 115),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  labs(
    title = "Aneuploid enrichment by state",
    x = "Aneuploid cells (%)",
    y = NULL
  ) +
  theme_blood(base_size = 9)

# ----------------------------
# Panel D
# ----------------------------

aneuploid_composition <- meta_df %>%
  dplyr::filter(copykat_prediction == "aneuploid") %>%
  dplyr::filter(!is.na(manuscript_state)) %>%
  dplyr::group_by(manuscript_state) %>%
  dplyr::summarise(
    n_aneuploid = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    percent_of_aneuploid = 100 * n_aneuploid / sum(n_aneuploid),
    manuscript_state = factor(
      manuscript_state,
      levels = rev(state_order)
    )
  )

max_d <- max(aneuploid_composition$percent_of_aneuploid, na.rm = TRUE)

pD <- ggplot(
  aneuploid_composition,
  aes(
    x = percent_of_aneuploid,
    y = manuscript_state,
    fill = manuscript_state
  )
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = paste0(round(percent_of_aneuploid), "%")),
    hjust = -0.15,
    size = 2.7
  ) +
  scale_y_discrete(labels = state_labels_short) +
  scale_fill_manual(values = state_colors) +
  scale_x_continuous(
    limits = c(0, max_d * 1.15),
    breaks = pretty_breaks(n = 5)
  ) +
  labs(
    title = "Composition of the malignant compartment",
    x = "Fraction of all aneuploid cells (%)",
    y = NULL
  ) +
  theme_blood(base_size = 9) +
  theme(legend.position = "none")

# ----------------------------
# Assemble compact 2x2 Figure 2
# ----------------------------

fig2 <- (pA | pB) /
  (pC | pD) +
  plot_layout(
    widths = c(1.35, 1),
    heights = c(1.25, 1)
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 11, face = "bold")
  )

# ----------------------------
# Save
# ----------------------------

ggsave(
  filename = "results/main_figures/Figure_2_CopyKAT_TPS12_malignant_enrichment_RECAP.pdf",
  plot = fig2,
  width = 10.5,
  height = 7.5,
  units = "in",
  limitsize = FALSE
)

ggsave(
  filename = "results/main_figures/Figure_2_CopyKAT_TPS12_malignant_enrichment_RECAP.tiff",
  plot = fig2,
  width = 10.5,
  height = 7.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  limitsize = FALSE
)

write.csv(
  aneuploid_tab,
  "results/main_figures/Table_Figure2_aneuploid_fraction_by_state.csv",
  row.names = FALSE
)

write.csv(
  aneuploid_composition,
  "results/main_figures/Table_Figure2_aneuploid_composition_by_state.csv",
  row.names = FALSE
)

write.csv(
  copykat_clean,
  "results/main_figures/Table_Figure2_TPS12_by_CopyKAT_state.csv",
  row.names = FALSE
)

message("Figure 2 recreated successfully as a compact 2x2 Blood-style version.")
message("Saved PDF and TIFF in results/main_figures/")