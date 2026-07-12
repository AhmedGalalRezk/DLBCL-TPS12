# ============================================================
# FIGURE 5 — Biological and clinical correlates of TPS in GSE31312
# ============================================================

library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(patchwork)
library(stringr)

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

cox_extract <- function(model, label, variable_name = "TPS_z") {
  s <- summary(model)
  row_name <- grep(variable_name, rownames(s$coefficients), value = TRUE)[1]
  
  data.frame(
    association = label,
    variable = row_name,
    estimate = s$coefficients[row_name, "exp(coef)"],
    lower95 = s$conf.int[row_name, "lower .95"],
    upper95 = s$conf.int[row_name, "upper .95"],
    p_value = s$coefficients[row_name, "Pr(>|z|)"],
    metric = "HR per SD TPS",
    row.names = NULL
  )
}

glm_extract <- function(model, label, variable_name = "TPS_z") {
  s <- summary(model)
  row_name <- grep(variable_name, rownames(s$coefficients), value = TRUE)[1]
  
  beta <- s$coefficients[row_name, "Estimate"]
  se <- s$coefficients[row_name, "Std. Error"]
  p <- s$coefficients[row_name, "Pr(>|z|)"]
  
  data.frame(
    association = label,
    variable = row_name,
    estimate = exp(beta),
    lower95 = exp(beta - 1.96 * se),
    upper95 = exp(beta + 1.96 * se),
    p_value = p,
    metric = "OR per SD TPS",
    row.names = NULL
  )
}

# ============================================================
# Load GSE31312 metadata
# ============================================================

metadata_path <- "results/GSE31312_context_dependence/GSE31312_context_harmonized_metadata.csv"

if (!file.exists(metadata_path)) {
  stop("Missing file: ", metadata_path)
}

df <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
colnames(df) <- make.names(colnames(df), unique = TRUE)

# ============================================================
# Harmonize TPS / PFS / IPI
# ============================================================

if ("TPS_context" %in% colnames(df)) {
  df$TPS <- as.numeric(df$TPS_context)
} else if ("TPS12_Score" %in% colnames(df)) {
  df$TPS <- as.numeric(df$TPS12_Score)
} else if ("TherapyPersistenceScore" %in% colnames(df)) {
  df$TPS <- as.numeric(df$TherapyPersistenceScore)
} else {
  stop("No TPS column found.")
}

df$TPS_z <- as.numeric(scale(df$TPS))

df$TPS_group <- factor(
  ifelse(df$TPS >= median(df$TPS, na.rm = TRUE), "High", "Low"),
  levels = c("Low", "High")
)

df$PFS_time_years <- if ("PFS_time_years" %in% colnames(df)) {
  as.numeric(df$PFS_time_years)
} else if ("Time_context" %in% colnames(df)) {
  as.numeric(df$Time_context)
} else {
  stop("No PFS time column found.")
}

df$PFS_event <- if ("PFS_event" %in% colnames(df)) {
  as.numeric(df$PFS_event)
} else if ("Event_context" %in% colnames(df)) {
  as.numeric(df$Event_context)
} else {
  stop("No PFS event column found.")
}

ipi_col <- grep("IPI|ipi", colnames(df), value = TRUE)[1]

if (is.na(ipi_col)) {
  stop("No IPI column found.")
}

df$IPI_num <- suppressWarnings(as.numeric(df[[ipi_col]]))

if (all(is.na(df$IPI_num))) {
  df$IPI_num <- suppressWarnings(
    as.numeric(str_extract(as.character(df[[ipi_col]]), "\\d+"))
  )
}

df$High_IPI <- factor(
  ifelse(df$IPI_num >= 3, "High IPI", "Low/intermediate IPI"),
  levels = c("Low/intermediate IPI", "High IPI")
)

# Optional response variable
response_col <- grep(
  "response|clinical_response|best_response|Response|CR|PR|SD|PD",
  colnames(df),
  value = TRUE,
  ignore.case = TRUE
)[1]

if (!is.na(response_col)) {
  df$Response <- toupper(trimws(as.character(df[[response_col]])))
  df$CR_binary <- ifelse(
    df$Response == "CR",
    1,
    ifelse(is.na(df$Response) | df$Response == "", NA, 0)
  )
} else {
  message("No response column found. CR vs non-CR will be skipped.")
  df$Response <- NA_character_
  df$CR_binary <- NA_real_
}

df$Early_progression <- ifelse(
  df$PFS_event == 1 & df$PFS_time_years <= 1,
  1,
  ifelse(!is.na(df$PFS_time_years) & !is.na(df$PFS_event), 0, NA)
)

# ============================================================
# Panel A: Hallmark correlates
# ============================================================

hallmark_path <- "results/main_figures/Table_main_GSE31312_pathway_correlations.csv"

if (!file.exists(hallmark_path)) {
  stop("Missing Hallmark/pathway correlation file: ", hallmark_path)
}

hallmark <- read.csv(hallmark_path, stringsAsFactors = FALSE, check.names = FALSE)
colnames(hallmark) <- make.names(colnames(hallmark), unique = TRUE)

pathway_col <- grep(
  "pathway|hallmark|gene_set|geneset|term|name",
  colnames(hallmark),
  value = TRUE,
  ignore.case = TRUE
)[1]

rho_col <- grep(
  "rho|spearman|correlation|cor",
  colnames(hallmark),
  value = TRUE,
  ignore.case = TRUE
)[1]

fdr_col <- grep(
  "FDR|adj|q.value|padj|p.adjust",
  colnames(hallmark),
  value = TRUE,
  ignore.case = TRUE
)[1]

if (is.na(pathway_col)) {
  stop("Could not identify pathway column. Columns are: ", paste(colnames(hallmark), collapse = ", "))
}

if (is.na(rho_col)) {
  stop("Could not identify rho/correlation column. Columns are: ", paste(colnames(hallmark), collapse = ", "))
}

message("Pathway column used: ", pathway_col)
message("Rho column used: ", rho_col)
message("FDR column used: ", fdr_col)

hallmark_plot <- hallmark %>%
  mutate(
    pathway_raw = as.character(.data[[pathway_col]]),
    pathway_clean = gsub("HALLMARK_", "", pathway_raw),
    pathway_clean = gsub("_", " ", pathway_clean),
    pathway_clean = toupper(pathway_clean),
    rho = as.numeric(.data[[rho_col]]),
    FDR_value = if (!is.na(fdr_col)) as.numeric(.data[[fdr_col]]) else NA_real_,
    FDR_sig = ifelse(!is.na(FDR_value) & FDR_value < 0.05, TRUE, FALSE)
  ) %>%
  filter(!is.na(rho)) %>%
  filter(
    grepl(
      "E2F|G2M|MYC|DNA REPAIR|P53|INTERFERON GAMMA|APOPTOSIS|TNFA|HYPOXIA|INFLAMMATORY",
      pathway_clean
    )
  ) %>%
  distinct(pathway_clean, .keep_all = TRUE)

if (nrow(hallmark_plot) == 0) {
  stop("No selected Hallmark pathways matched. Check pathway names in Table_main_GSE31312_pathway_correlation.csv.")
}

hallmark_plot <- hallmark_plot %>%
  arrange(rho) %>%
  mutate(
    pathway_clean = factor(pathway_clean, levels = unique(pathway_clean))
  )

pA <- ggplot(
  hallmark_plot,
  aes(x = rho, y = pathway_clean, fill = FDR_sig)
) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(
    values = c("FALSE" = "grey75", "TRUE" = "#B2182B"),
    name = "FDR < 0.05"
  ) +
  labs(
    title = "Hallmark correlates of TPS",
    x = "Spearman rho",
    y = NULL
  ) +
  theme_blood(base_size = 9) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 8)
  )

# ============================================================
# Panel B: High-IPI Kaplan-Meier
# ============================================================

km_df <- df %>%
  filter(
    High_IPI == "High IPI",
    !is.na(PFS_time_years),
    !is.na(PFS_event),
    !is.na(TPS_group),
    PFS_time_years > 0
  )

if (nrow(km_df) < 10 || length(unique(km_df$PFS_event)) < 2) {
  stop("Insufficient high-IPI survival data for Kaplan-Meier plot.")
}

fit_high_ipi <- survfit(
  Surv(PFS_time_years, PFS_event) ~ TPS_group,
  data = km_df
)

pB_surv <- ggsurvplot(
  fit_high_ipi,
  data = km_df,
  pval = TRUE,
  risk.table = FALSE,
  conf.int = FALSE,
  palette = c("steelblue", "firebrick"),
  legend.title = "TPS group",
  legend.labs = c("Low", "High"),
  xlab = "PFS time (years)",
  ylab = "PFS probability",
  title = "High-IPI PFS by TPS group",
  ggtheme = theme_blood(base_size = 9)
)

pB <- pB_surv$plot +
  theme(
    legend.position = c(0.82, 0.82),
    legend.background = element_blank()
  )


# ============================================================
# Assemble Figure 5
# ============================================================

fig5 <- pA / pB +
  plot_layout(heights = c(1.05, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 12, face = "bold")
  )
# ============================================================
# Save
# ============================================================

ggsave(
  filename = "results/main_figures/Figure_5_GSE31312_biologic_clinical_context_REVISED.pdf",
  plot = fig5,
  width = 9,
  height = 8,
  units = "in"
)

ggsave(
  filename = "results/main_figures/Figure_5_GSE31312_biologic_clinical_context_REVISED.tiff",
  plot = fig5,
  width = 9,
  height = 8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

message("Revised Figure 5 saved successfully.")