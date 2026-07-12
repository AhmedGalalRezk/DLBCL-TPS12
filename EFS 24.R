# ============================================================
# SUPPLEMENTARY ANALYSIS — TPS12 association with EFS24 failure
#
# Input:
#   results/GSE31312_context_dependence/GSE31312_context_harmonized_metadata.csv
#
# Required columns:
#   TPS_context
#   PFS_time_years
#   PFS_event
#
# Definition:
#   EFS24_failure = 1 if PFS_event == 1 and PFS_time_years < 2
#   EFS24_failure = 0 if PFS_time_years >= 2
#   censored patients with PFS_time_years < 2 are excluded
#
# Outputs:
#   results/supplementary_figures/Supplementary_TPS12_EFS24_failure.pdf
#   results/supplementary_figures/Supplementary_TPS12_EFS24_failure.tiff
#   results/supplementary_figures/Table_TPS12_EFS24_logistic_regression.csv
#   results/supplementary_figures/Table_TPS12_EFS24_group_test.csv
# ============================================================

library(dplyr)
library(ggplot2)
library(survival)
library(survminer)

set.seed(123)

dir.create(
  "results/supplementary_figures",
  recursive = TRUE,
  showWarnings = FALSE
)

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

# ============================================================
# Load data
# ============================================================

input_file <- "results/GSE31312_context_dependence/GSE31312_context_harmonized_metadata.csv"

if (!file.exists(input_file)) {
  stop("Missing file: ", input_file)
}

df <- read.csv(input_file, stringsAsFactors = FALSE)

required_cols <- c(
  "TPS_context",
  "PFS_time_years",
  "PFS_event"
)

missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ============================================================
# Clean and derive strict EFS24 endpoint
# ============================================================

dat0 <- df %>%
  dplyr::mutate(
    TPS12 = as.numeric(TPS_context),
    PFS_time_years = as.numeric(PFS_time_years),
    PFS_event = as.numeric(PFS_event)
  ) %>%
  dplyr::filter(
    !is.na(TPS12),
    !is.na(PFS_time_years),
    !is.na(PFS_event)
  )

# Strict EFS24 handling:
# - Event before 24 months = EFS24 failure
# - Event-free at or beyond 24 months = EFS24 achieved
# - Censored before 24 months = indeterminate; exclude

dat <- dat0 %>%
  dplyr::mutate(
    EFS24_status = dplyr::case_when(
      PFS_event == 1 & PFS_time_years < 2 ~ "EFS24 failure",
      PFS_time_years >= 2 ~ "EFS24 achieved",
      PFS_event == 0 & PFS_time_years < 2 ~ "Indeterminate",
      TRUE ~ NA_character_
    ),
    EFS24_failure = dplyr::case_when(
      EFS24_status == "EFS24 failure" ~ 1,
      EFS24_status == "EFS24 achieved" ~ 0,
      TRUE ~ NA_real_
    )
  )

efs_summary <- dat %>%
  dplyr::count(EFS24_status, name = "n")

print(efs_summary)

write.csv(
  efs_summary,
  "results/supplementary_figures/Table_TPS12_EFS24_endpoint_summary.csv",
  row.names = FALSE
)

dat_efs <- dat %>%
  dplyr::filter(!is.na(EFS24_failure)) %>%
  dplyr::mutate(
    EFS24_status = factor(
      EFS24_status,
      levels = c("EFS24 achieved", "EFS24 failure")
    ),
    TPS_group = ifelse(
      TPS12 >= median(TPS12, na.rm = TRUE),
      "TPS-high",
      "TPS-low"
    ),
    TPS_group = factor(
      TPS_group,
      levels = c("TPS-low", "TPS-high")
    )
  )

cat("\nStrict EFS24 analysis set:\n")
print(table(dat_efs$EFS24_status))
print(table(dat_efs$TPS_group, dat_efs$EFS24_status))

# ============================================================
# Statistical tests
# ============================================================

# Continuous TPS12 by EFS24 status
wilcox_res <- wilcox.test(
  TPS12 ~ EFS24_status,
  data = dat_efs
)

# Logistic regression, continuous TPS12
fit_cont <- glm(
  EFS24_failure ~ TPS12,
  family = binomial(),
  data = dat_efs
)

fit_group <- glm(
  EFS24_failure ~ TPS_group,
  family = binomial(),
  data = dat_efs
)

or_table_cont <- data.frame(
  model = "Continuous TPS12",
  term = names(coef(fit_cont)),
  OR = exp(coef(fit_cont)),
  lower95 = exp(confint.default(fit_cont)[, 1]),
  upper95 = exp(confint.default(fit_cont)[, 2]),
  p_value = summary(fit_cont)$coefficients[, 4],
  row.names = NULL
)

or_table_group <- data.frame(
  model = "TPS-high vs TPS-low",
  term = names(coef(fit_group)),
  OR = exp(coef(fit_group)),
  lower95 = exp(confint.default(fit_group)[, 1]),
  upper95 = exp(confint.default(fit_group)[, 2]),
  p_value = summary(fit_group)$coefficients[, 4],
  row.names = NULL
)

or_table <- dplyr::bind_rows(
  or_table_cont,
  or_table_group
)

print(or_table)

write.csv(
  or_table,
  "results/supplementary_figures/Table_TPS12_EFS24_logistic_regression.csv",
  row.names = FALSE
)

# Fisher test for TPS-high/low vs EFS24 failure
tab_group <- table(
  dat_efs$TPS_group,
  dat_efs$EFS24_status
)

fisher_res <- fisher.test(tab_group)

group_table <- as.data.frame(tab_group)
group_table$fisher_p <- fisher_res$p.value

write.csv(
  group_table,
  "results/supplementary_figures/Table_TPS12_EFS24_group_test.csv",
  row.names = FALSE
)

# ============================================================
# Plot 1: TPS12 by EFS24 status
# ============================================================

p1 <- ggplot(
  dat_efs,
  aes(
    x = EFS24_status,
    y = TPS12,
    fill = EFS24_status
  )
) +
  geom_violin(
    trim = TRUE,
    color = "grey30",
    linewidth = 0.3,
    scale = "width"
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = c(
      "EFS24 achieved" = "grey70",
      "EFS24 failure" = "#B2182B"
    )
  ) +
  labs(
    title = "TPS12 by EFS24 status",
    subtitle = paste0(
      "Wilcoxon P = ",
      format.pval(wilcox_res$p.value, digits = 3, eps = 1e-4)
    ),
    x = NULL,
    y = "TPS12 score"
  ) +
  theme_blood(base_size = 9) +
  theme(
    legend.position = "none"
  )

# ============================================================
# Plot 2: EFS24 failure proportion by TPS group
# ============================================================

prop_df <- dat_efs %>%
  dplyr::group_by(TPS_group, EFS24_status) %>%
  dplyr::summarise(
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::group_by(TPS_group) %>%
  dplyr::mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(EFS24_status == "EFS24 failure")

p2 <- ggplot(
  prop_df,
  aes(
    x = TPS_group,
    y = 100 * proportion,
    fill = TPS_group
  )
) +
  geom_col(
    width = 0.65,
    color = "grey30",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = paste0(round(100 * proportion), "%")),
    vjust = -0.4,
    size = 3
  ) +
  scale_fill_manual(
    values = c(
      "TPS-low" = "grey70",
      "TPS-high" = "#B2182B"
    )
  ) +
  scale_y_continuous(
    limits = c(0, max(100 * prop_df$proportion, na.rm = TRUE) * 1.25),
    breaks = pretty_breaks(n = 5)
  ) +
  labs(
    title = "EFS24 failure by TPS12 group",
    subtitle = paste0(
      "Fisher P = ",
      format.pval(fisher_res$p.value, digits = 3, eps = 1e-4)
    ),
    x = NULL,
    y = "EFS24 failure (%)"
  ) +
  theme_blood(base_size = 9) +
  theme(
    legend.position = "none"
  )

# ============================================================
# Optional PFS Kaplan-Meier analysis
# ============================================================

surv_dat <- dat0 %>%
  dplyr::mutate(
    TPS_group = ifelse(
      TPS12 >= median(TPS12, na.rm = TRUE),
      "TPS-high",
      "TPS-low"
    ),
    TPS_group = factor(
      TPS_group,
      levels = c("TPS-low", "TPS-high")
    )
  )

surv_obj <- survival::Surv(
  time = surv_dat$PFS_time_years,
  event = surv_dat$PFS_event
)

km_fit <- survival::survfit(
  surv_obj ~ TPS_group,
  data = surv_dat
)

cox_fit <- survival::coxph(
  surv_obj ~ TPS12,
  data = surv_dat
)

cox_table <- data.frame(
  term = names(coef(cox_fit)),
  HR = exp(coef(cox_fit)),
  lower95 = exp(confint(cox_fit)[, 1]),
  upper95 = exp(confint(cox_fit)[, 2]),
  p_value = summary(cox_fit)$coefficients[, 5],
  row.names = NULL
)

write.csv(
  cox_table,
  "results/supplementary_figures/Table_TPS12_PFS_Cox.csv",
  row.names = FALSE
)

km_plot <- survminer::ggsurvplot(
  km_fit,
  data = surv_dat,
  risk.table = FALSE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c("grey60", "#B2182B"),
  xlab = "Progression-free survival time (years)",
  ylab = "Progression-free survival probability",
  legend.title = NULL,
  legend.labs = c("TPS-low", "TPS-high"),
  ggtheme = theme_blood(base_size = 9)
)

p3 <- km_plot$plot +
  labs(
    title = "Progression-free survival by TPS12 group"
  ) +
  theme(
    plot.title = element_text(face = "bold")
  )

# ============================================================
# Assemble figure
# ============================================================

efs_fig <- (p1 | p2) / p3 +
  patchwork::plot_layout(
    heights = c(1, 1.1)
  ) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 10, face = "bold")
  )

# ============================================================
# Save
# ============================================================

ggsave(
  filename = "results/supplementary_figures/Supplementary_TPS12_EFS24_failure.pdf",
  plot = efs_fig,
  width = 9,
  height = 8,
  units = "in",
  limitsize = FALSE
)

ggsave(
  filename = "results/supplementary_figures/Supplementary_TPS12_EFS24_failure.tiff",
  plot = efs_fig,
  width = 9,
  height = 8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  limitsize = FALSE
)

message("TPS12-EFS24/PFS analysis complete.")
message("Results saved in results/supplementary_figures/")