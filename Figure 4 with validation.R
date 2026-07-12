
# ============================================================
# FIGURE 4 — GSE10846 survival validation + multi-cohort forest
#
# Panels:
#   A. Kaplan–Meier OS by TPS group in GSE10846
#   B. Multivariable Cox model in GSE10846
#   C. TPS distribution by TPS group in GSE10846
#   D. External cohort Cox validation forest plot
#
# Outputs:
#   results/main_figures/Figure_4_GSE10846_survival_validation_RECAP.pdf
#   results/main_figures/Figure_4_GSE10846_survival_validation_RECAP.tiff
# ============================================================

library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(patchwork)
library(readr)

dir.create("results/main_figures", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Theme
# ============================================================

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

safe_read_csv <- function(path) {
  if (file.exists(path)) {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    NULL
  }
}

# ============================================================
# Load GSE10846 data
# ============================================================

merged_path <- "results/GSE10846/GSE10846_survival_merged_with_TPS.csv"
cox_path <- "results/GSE10846/GSE10846_multivariable_cox_TPS_age_stage.csv"

if (!file.exists(merged_path)) {
  stop("Missing file: ", merged_path)
}

if (!file.exists(cox_path)) {
  stop("Missing file: ", cox_path)
}

merged <- read.csv(
  merged_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cox_results <- read.csv(
  cox_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ============================================================
# Harmonize survival variables
# ============================================================

if (!"OS_time" %in% colnames(merged)) {
  if ("futime(year)" %in% colnames(merged)) {
    merged$OS_time <- as.numeric(merged$`futime(year)`)
  } else {
    stop("Could not find OS_time or futime(year).")
  }
}

if (!"OS_event" %in% colnames(merged)) {
  if ("fustaut" %in% colnames(merged)) {
    merged$OS_event <- as.numeric(merged$fustaut)
  } else {
    stop("Could not find OS_event or fustaut.")
  }
}

if (!"TherapyPersistenceScore" %in% colnames(merged)) {
  stop("TherapyPersistenceScore column not found.")
}

if (!"TherapyPersistenceGroup" %in% colnames(merged)) {
  merged$TherapyPersistenceGroup <- ifelse(
    merged$TherapyPersistenceScore >= median(
      merged$TherapyPersistenceScore,
      na.rm = TRUE
    ),
    "High",
    "Low"
  )
}

merged$OS_time <- as.numeric(merged$OS_time)
merged$OS_event <- as.numeric(merged$OS_event)

merged$TherapyPersistenceGroup <- factor(
  merged$TherapyPersistenceGroup,
  levels = c("Low", "High")
)

# ============================================================
# Panel A: Kaplan–Meier OS
# ============================================================

km_df <- merged %>%
  filter(
    !is.na(OS_time),
    !is.na(OS_event),
    !is.na(TherapyPersistenceGroup)
  )

fit <- survfit(
  Surv(OS_time, OS_event) ~ TherapyPersistenceGroup,
  data = km_df
)

pA_surv <- ggsurvplot(
  fit,
  data = km_df,
  pval = TRUE,
  risk.table = FALSE,
  conf.int = FALSE,
  palette = c("steelblue", "firebrick"),
  legend.title = "TPS group",
  legend.labs = c("Low", "High"),
  xlab = "Overall survival time (years)",
  ylab = "Overall survival probability",
  title = "GSE10846 overall survival by TPS",
  ggtheme = theme_blood(base_size = 9)
)

pA <- pA_surv$plot +
  theme(
    legend.position = c(0.78, 0.82),
    legend.background = element_blank()
  )

# ============================================================
# Panel B: Multivariable Cox forest plot
# ============================================================

cox_plot_df <- cox_results %>%
  mutate(
    label = case_when(
      variable == "TherapyPersistenceScore" ~ "TPS",
      variable == "Age_num" ~ "Age",
      variable == "Stage_cleanII" ~ "Stage II",
      variable == "Stage_cleanIII" ~ "Stage III",
      variable == "Stage_cleanIV" ~ "Stage IV",
      TRUE ~ variable
    ),
    label = factor(
      label,
      levels = rev(c("TPS", "Age", "Stage II", "Stage III", "Stage IV"))
    )
  )

pB <- ggplot(
  cox_plot_df,
  aes(x = HR, y = label)
) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35) +
  geom_errorbarh(
    aes(xmin = lower_95, xmax = upper_95),
    height = 0.18,
    linewidth = 0.45
  ) +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(
    title = "Multivariable Cox model",
    subtitle = "Adjusted for age and stage",
    x = "Hazard ratio, log scale",
    y = NULL
  ) +
  theme_blood(base_size = 9)

# ============================================================
# Panel C: TPS distribution
# ============================================================

pC <- merged %>%
  filter(
    !is.na(TherapyPersistenceScore),
    !is.na(TherapyPersistenceGroup)
  ) %>%
  ggplot(
    aes(
      x = TherapyPersistenceGroup,
      y = TherapyPersistenceScore,
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
    width = 0.12,
    outlier.shape = NA,
    linewidth = 0.25,
    fill = "white"
  ) +
  scale_fill_manual(
    values = c("Low" = "steelblue", "High" = "firebrick")
  ) +
  labs(
    title = "TPS distribution",
    x = "TPS group",
    y = "TranscriptionalPersistenceScore"
  ) +
  theme_blood(base_size = 9) +
  theme(legend.position = "none")

# ============================================================
# Panel D: Multi-cohort validation forest plot
# ============================================================

# Try to read cohort-specific Cox outputs if available.
# If a file is missing, the corresponding row remains NA and will be dropped.

gse10846_row <- data.frame(
  Cohort = "GSE10846",
  HR = cox_results$HR[cox_results$variable == "TherapyPersistenceScore"][1],
  lower95 = cox_results$lower_95[cox_results$variable == "TherapyPersistenceScore"][1],
  upper95 = cox_results$upper_95[cox_results$variable == "TherapyPersistenceScore"][1]
)

# GSE23501
gse23501_path <- "results/GSE23501_full_TPS/GSE23501_TPS12_vs_fullTPS_cox_results.csv"
gse23501_tab <- safe_read_csv(gse23501_path)

gse23501_row <- data.frame(
  Cohort = "GSE23501",
  HR = NA_real_,
  lower95 = NA_real_,
  upper95 = NA_real_
)

if (!is.null(gse23501_tab)) {
  gse23501_tps12_os <- gse23501_tab %>%
    filter(
      Endpoint == "OS",
      Panel == "TPS12"
    )
  
  if (nrow(gse23501_tps12_os) > 0) {
    gse23501_row$HR <- gse23501_tps12_os$HR[1]
    
    # If CIs are unavailable in this file, they remain NA.
    if ("lower95" %in% colnames(gse23501_tps12_os)) {
      gse23501_row$lower95 <- gse23501_tps12_os$lower95[1]
    }
    if ("upper95" %in% colnames(gse23501_tps12_os)) {
      gse23501_row$upper95 <- gse23501_tps12_os$upper95[1]
    }
  }
}

# GSE31312
gse31312_path <- "results/GSE31312/GSE31312_TPS12_survival_cox_results.csv"
gse31312_tab <- safe_read_csv(gse31312_path)

gse31312_row <- data.frame(
  Cohort = "GSE31312",
  HR = NA_real_,
  lower95 = NA_real_,
  upper95 = NA_real_
)

if (!is.null(gse31312_tab)) {
  candidate_hr <- grep("HR|hazard", colnames(gse31312_tab), value = TRUE, ignore.case = TRUE)
  candidate_lo <- grep("lower|low|L95|lower95", colnames(gse31312_tab), value = TRUE, ignore.case = TRUE)
  candidate_hi <- grep("upper|high|U95|upper95", colnames(gse31312_tab), value = TRUE, ignore.case = TRUE)
  
  if (length(candidate_hr) > 0) {
    gse31312_row$HR <- as.numeric(gse31312_tab[[candidate_hr[1]]][1])
  }
  if (length(candidate_lo) > 0) {
    gse31312_row$lower95 <- as.numeric(gse31312_tab[[candidate_lo[1]]][1])
  }
  if (length(candidate_hi) > 0) {
    gse31312_row$upper95 <- as.numeric(gse31312_tab[[candidate_hi[1]]][1])
  }
}

# GSE32918
gse32918_path <- "results/GSE32918/GSE32918_TPS12_survival_cox_results.csv"
gse32918_tab <- safe_read_csv(gse32918_path)

gse32918_row <- data.frame(
  Cohort = "GSE32918",
  HR = NA_real_,
  lower95 = NA_real_,
  upper95 = NA_real_
)

if (!is.null(gse32918_tab)) {
  candidate_hr <- grep("HR|hazard", colnames(gse32918_tab), value = TRUE, ignore.case = TRUE)
  candidate_lo <- grep("lower|low|L95|lower95", colnames(gse32918_tab), value = TRUE, ignore.case = TRUE)
  candidate_hi <- grep("upper|high|U95|upper95", colnames(gse32918_tab), value = TRUE, ignore.case = TRUE)
  
  if (length(candidate_hr) > 0) {
    gse32918_row$HR <- as.numeric(gse32918_tab[[candidate_hr[1]]][1])
  }
  if (length(candidate_lo) > 0) {
    gse32918_row$lower95 <- as.numeric(gse32918_tab[[candidate_lo[1]]][1])
  }
  if (length(candidate_hi) > 0) {
    gse32918_row$upper95 <- as.numeric(gse32918_tab[[candidate_hi[1]]][1])
  }
}

validation_forest <- bind_rows(
  gse10846_row,
  gse31312_row,
  gse23501_row,
  gse32918_row
)

# If CI is unavailable but HR exists, create NA-safe columns.
validation_forest <- validation_forest %>%
  mutate(
    HR = as.numeric(HR),
    lower95 = as.numeric(lower95),
    upper95 = as.numeric(upper95)
  )

write.csv(
  validation_forest,
  "results/main_figures/Table_Figure4_multicohort_validation_forest.csv",
  row.names = FALSE
)

validation_plot_df <- validation_forest %>%
  filter(
    !is.na(HR),
    !is.na(lower95),
    !is.na(upper95)
  ) %>%
  mutate(
    Cohort = factor(Cohort, levels = rev(Cohort))
  )

if (nrow(validation_plot_df) == 0) {
  
  pD <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "Add external cohort HR and 95% CI values\nto Table_Figure4_multicohort_validation_forest.csv",
      size = 3.2
    ) +
    xlim(-1, 1) +
    ylim(-1, 1) +
    theme_void() +
    labs(title = "External cohort validation")
  
} else {
  
  pD <- ggplot(
    validation_plot_df,
    aes(x = HR, y = Cohort)
  ) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35) +
    geom_errorbarh(
      aes(xmin = lower95, xmax = upper95),
      height = 0.18,
      linewidth = 0.45
    ) +
    geom_point(size = 2) +
    scale_x_log10() +
    labs(
      title = "External cohort validation",
      x = "Hazard ratio, log scale",
      y = NULL
    ) +
    theme_blood(base_size = 9)
}

# ============================================================
# Assemble figure
# ============================================================

fig4 <- pA /
  ((pB | pC) / pD) +
  plot_layout(
    heights = c(1.15, 1.25)
  ) +
  plot_annotation(tag_levels = "A")

# ============================================================
# Save PDF and TIFF
# ============================================================

ggsave(
  filename = "results/main_figures/Figure_4_GSE10846_survival_validation_RECAP.pdf",
  plot = fig4,
  width = 9,
  height = 9,
  units = "in"
)

ggsave(
  filename = "results/main_figures/Figure_4_GSE10846_survival_validation_RECAP.tiff",
  plot = fig4,
  width = 9,
  height = 9,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

message("Figure 4 recreated and saved as PDF and TIFF.")
message("Multi-cohort forest table saved to results/main_figures/Table_Figure4_multicohort_validation_forest.csv")