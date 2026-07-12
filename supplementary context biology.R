# ============================================================
# SUPPLEMENTARY FIGURE — Context-dependent TPS biology
#
# Combines:
#   A. GSE23501 TPS12 by molecular subtype
#   B. GSE23501 TPS12 by treatment context
#   C. GSE31312 TPS by GEP subgroup
#   D. GSE31312 TPS by IHC subgroup
#
# Outputs:
#   results/supplementary_figures/Supplementary_Figure_TPS_context_dependence.pdf
#   results/supplementary_figures/Supplementary_Figure_TPS_context_dependence.tiff
# ============================================================

library(dplyr)
library(ggplot2)
library(patchwork)

dir.create("results/supplementary_figures", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Theme
# ============================================================

theme_blood <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2),
      axis.title = element_text(size = base_size + 1),
      axis.text = element_text(size = base_size),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# ============================================================
# Load GSE23501 harmonized metadata
# ============================================================

gse23501_path <- "results/GSE23501_context_dependence/GSE23501_TPS12_context_metadata.csv"

if (!file.exists(gse23501_path)) {
  stop("Missing file: ", gse23501_path)
}

gse23501 <- read.csv(
  gse23501_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Clean factors if needed
gse23501$COO <- factor(gse23501$COO)
gse23501$Treatment <- factor(gse23501$Treatment)
gse23501$TPS12 <- as.numeric(gse23501$TPS12)

# ============================================================
# Load GSE31312 harmonized metadata
# ============================================================

gse31312_path <- "results/GSE31312_context_dependence/GSE31312_context_harmonized_metadata.csv"

if (!file.exists(gse31312_path)) {
  stop("Missing file: ", gse31312_path)
}

gse31312 <- read.csv(
  gse31312_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gse31312$COO_context <- factor(gse31312$COO_context)
gse31312$IHC_context <- factor(gse31312$IHC_context)
gse31312$TPS_context <- as.numeric(gse31312$TPS_context)

# ============================================================
# Panel A: GSE23501 TPS12 by molecular subtype
# ============================================================

pA <- ggplot(
  gse23501 %>%
    filter(!is.na(COO), !is.na(TPS12)),
  aes(x = COO, y = TPS12, fill = COO)
) +
  geom_violin(trim = TRUE, color = "grey30", linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  labs(
    title = "GSE23501: TPS12 by molecular subtype",
    x = "Molecular subtype",
    y = "TPS12"
  ) +
  theme_blood(base_size = 9)

# ============================================================
# Panel B: GSE23501 TPS12 by treatment context
# ============================================================

# Drop tiny/blank treatment groups only if needed
gse23501_tx <- gse23501 %>%
  filter(!is.na(Treatment), !is.na(TPS12)) %>%
  mutate(
    Treatment = as.character(Treatment),
    Treatment = ifelse(Treatment == "" | Treatment == "NA", NA, Treatment)
  ) %>%
  filter(!is.na(Treatment)) %>%
  group_by(Treatment) %>%
  mutate(n_treatment = n()) %>%
  ungroup()

# Keep all treatments but order by median TPS12
gse23501_tx$Treatment <- reorder(
  gse23501_tx$Treatment,
  gse23501_tx$TPS12,
  FUN = median,
  na.rm = TRUE
)

pB <- ggplot(
  gse23501_tx,
  aes(x = Treatment, y = TPS12, fill = Treatment)
) +
  geom_boxplot(outlier.size = 0.6, linewidth = 0.25) +
  labs(
    title = "GSE23501: TPS12 by treatment context",
    x = "Treatment",
    y = "TPS12"
  ) +
  theme_blood(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1, size = 6)
  )

# ============================================================
# Panel C: GSE31312 TPS by GEP subgroup
# ============================================================

pC <- ggplot(
  gse31312 %>%
    filter(!is.na(COO_context), !is.na(TPS_context)),
  aes(x = COO_context, y = TPS_context, fill = COO_context)
) +
  geom_violin(trim = TRUE, color = "grey30", linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  labs(
    title = "GSE31312: TPS by GEP subgroup",
    x = "GEP subgroup",
    y = "TPS"
  ) +
  theme_blood(base_size = 9)

# ============================================================
# Panel D: GSE31312 TPS by IHC subgroup
# ============================================================

pD <- ggplot(
  gse31312 %>%
    filter(!is.na(IHC_context), !is.na(TPS_context)),
  aes(x = IHC_context, y = TPS_context, fill = IHC_context)
) +
  geom_violin(trim = TRUE, color = "grey30", linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  labs(
    title = "GSE31312: TPS by IHC subgroup",
    x = "IHC subgroup",
    y = "TPS"
  ) +
  theme_blood(base_size = 9)

# ============================================================
# Assemble supplementary figure
# ============================================================

supp_context_fig <- (pA | pC) /
  (pB | pD) +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(
    tag_levels = "A",
    title = "Context-dependent TPS biology across external DLBCL cohorts"
  ) &
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0)
  )

# ============================================================
# Save PDF and TIFF
# ============================================================

ggsave(
  filename = "results/supplementary_figures/Supplementary_Figure_TPS_context_dependence.pdf",
  plot = supp_context_fig,
  width = 11,
  height = 8,
  units = "in"
)

ggsave(
  filename = "results/supplementary_figures/Supplementary_Figure_TPS_context_dependence.tiff",
  plot = supp_context_fig,
  width = 11,
  height = 8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

message("Supplementary context-dependence figure saved successfully.")