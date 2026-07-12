# =========================================================
# GSE31312 treatment-context / response-dependence analysis
# =========================================================

library(dplyr)
library(ggplot2)
library(survival)
library(survminer)

dir.create(
  "results/GSE31312_treatment_context",
  recursive = TRUE,
  showWarnings = FALSE
)

bio_df <- readRDS(
  "results/GSE31312_biologic_validation/GSE31312_biologic_validation_table.rds"
)

# =========================================
# 1. Response association
# =========================================

table(bio_df$Respon, useNA = "ifany")

kruskal.test(
  TPS_full_Score ~ Respon,
  data = bio_df
)

pairwise.wilcox.test(
  bio_df$TPS_full_Score,
  bio_df$Respon,
  p.adjust.method = "BH"
)

p_response <- ggplot(
  bio_df,
  aes(x = Respon, y = TPS_full_Score, fill = Respon)
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  theme_classic(base_size = 14) +
  labs(
    title = "TPS_full score by treatment response",
    x = "Response",
    y = "TPS_full score"
  )

print(p_response)

ggsave(
  "results/GSE31312_treatment_context/TPS_by_response.pdf",
  p_response,
  width = 7,
  height = 5
)

# =========================================
# 2. Early progression association
# =========================================

bio_df$Early_progression <- ifelse(
  bio_df$PFS_event == 1 & bio_df$PFS_time_years < 1,
  "Early progression",
  "No early progression"
)

bio_df$Early_progression <- factor(
  bio_df$Early_progression,
  levels = c("No early progression", "Early progression")
)

table(bio_df$Early_progression)

wilcox.test(
  TPS_full_Score ~ Early_progression,
  data = bio_df
)

p_early <- ggplot(
  bio_df,
  aes(x = Early_progression, y = TPS_full_Score, fill = Early_progression)
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  theme_classic(base_size = 14) +
  labs(
    title = "TPS_full score and early progression",
    x = "",
    y = "TPS_full score"
  )

print(p_early)

ggsave(
  "results/GSE31312_treatment_context/TPS_early_progression.pdf",
  p_early,
  width = 6,
  height = 5
)

# =========================================
# 3. Clinical-risk adjusted progression model
# =========================================

bio_df$Age_num <- as.numeric(bio_df$Age)
bio_df$IPI_num <- as.numeric(bio_df$`IPI sc`)

early_model <- glm(
  Early_progression ~ TPS_full_Score + Age_num + IPI_num,
  data = bio_df,
  family = binomial
)

summary(early_model)

exp(cbind(
  OR = coef(early_model),
  confint(early_model)
))

# =========================================
# 4. PFS Cox model with interaction by IPI risk group
# =========================================

bio_df$IPI_group <- factor(bio_df$`IPI code`)

table(bio_df$IPI_group, useNA = "ifany")

cox_ipi_interaction <- coxph(
  Surv(PFS_time_years, PFS_event) ~
    TPS_full_Score * IPI_group,
  data = bio_df
)

summary(cox_ipi_interaction)

cox_ipi_no_interaction <- coxph(
  Surv(PFS_time_years, PFS_event) ~
    TPS_full_Score + IPI_group,
  data = bio_df
)

anova(
  cox_ipi_no_interaction,
  cox_ipi_interaction,
  test = "LRT"
)

# =========================================
# 5. TPS effect separately within IPI groups
# =========================================

ipi_specific_results <- lapply(
  levels(bio_df$IPI_group),
  function(g) {
    
    df_g <- bio_df %>%
      filter(IPI_group == g)
    
    if (nrow(df_g) < 20 || length(unique(df_g$PFS_event)) < 2) {
      return(data.frame(
        IPI_group = g,
        n = nrow(df_g),
        events = sum(df_g$PFS_event, na.rm = TRUE),
        HR = NA,
        p = NA
      ))
    }
    
    fit <- coxph(
      Surv(PFS_time_years, PFS_event) ~ TPS_full_Score,
      data = df_g
    )
    
    s <- summary(fit)
    
    data.frame(
      IPI_group = g,
      n = nrow(df_g),
      events = sum(df_g$PFS_event, na.rm = TRUE),
      HR = exp(coef(fit)),
      p = s$coefficients[, "Pr(>|z|)"]
    )
  }
) %>%
  bind_rows()

ipi_specific_results

# =========================================
# 6. Biologic conservation within response groups
# =========================================

selected_sets <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_MYC_TARGETS_V1"
)

response_pathway_cor <- data.frame()

for (resp in unique(bio_df$Respon)) {
  
  df_r <- bio_df %>% filter(Respon == resp)
  
  if (nrow(df_r) < 10) next
  
  for (pathway in selected_sets) {
    
    ct <- cor.test(
      df_r$TPS_full_Score,
      df_r[[pathway]],
      method = "spearman"
    )
    
    response_pathway_cor <- rbind(
      response_pathway_cor,
      data.frame(
        Response = resp,
        Pathway = pathway,
        n = nrow(df_r),
        rho = unname(ct$estimate),
        p = ct$p.value
      )
    )
  }
}

response_pathway_cor$FDR <- p.adjust(
  response_pathway_cor$p,
  method = "fdr"
)

response_pathway_cor

# =========================================================
# IPI-specific Kaplan-Meier analysis
# =========================================================

library(survival)
library(survminer)
library(dplyr)

# ---------------------------------------------
# Create simplified IPI groups
# ---------------------------------------------

bio_df$IPI_group_simple <- ifelse(
  bio_df$`IPI code` %in% c("high, int-high"),
  "High IPI",
  "Low/Intermediate IPI"
)

table(bio_df$IPI_group_simple, useNA = "ifany")

# =========================================================
# KM WITHIN HIGH IPI
# =========================================================

df_high_ipi <- bio_df %>%
  filter(IPI_group_simple == "High IPI")

fit_high_ipi <- survfit(
  Surv(PFS_time_years, PFS_event) ~ TPS_full_Group,
  data = df_high_ipi
)

p_high_ipi <- ggsurvplot(
  fit_high_ipi,
  data = df_high_ipi,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("steelblue", "firebrick"),
  title = "High IPI subgroup",
  xlab = "Progression-free survival time (years)",
  ylab = "Progression-free survival probability"
)

print(p_high_ipi)

# Cox model within high IPI

cox_high_ipi <- coxph(
  Surv(PFS_time_years, PFS_event) ~ TPS_full_Score,
  data = df_high_ipi
)

summary(cox_high_ipi)

# =========================================================
# KM WITHIN LOW/INTERMEDIATE IPI
# =========================================================

df_low_ipi <- bio_df %>%
  filter(IPI_group_simple == "Low/Intermediate IPI")

fit_low_ipi <- survfit(
  Surv(PFS_time_years, PFS_event) ~ TPS_full_Group,
  data = df_low_ipi
)

p_low_ipi <- ggsurvplot(
  fit_low_ipi,
  data = df_low_ipi,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("steelblue", "firebrick"),
  title = "Low/Intermediate IPI subgroup",
  xlab = "Progression-free survival time (years)",
  ylab = "Progression-free survival probability"
)

print(p_low_ipi)

# Cox model within low/intermediate IPI

cox_low_ipi <- coxph(
  Surv(PFS_time_years, PFS_event) ~ TPS_full_Score,
  data = df_low_ipi
)

summary(cox_low_ipi)

# =========================================================
# Save figures
# =========================================================

dir.create(
  "results/GSE31312_treatment_context",
  recursive = TRUE,
  showWarnings = FALSE
)

pdf(
  "results/GSE31312_treatment_context/KM_high_IPI.pdf",
  width = 7,
  height = 6
)
print(p_high_ipi)
dev.off()

pdf(
  "results/GSE31312_treatment_context/KM_low_intermediate_IPI.pdf",
  width = 7,
  height = 6
)
print(p_low_ipi)
dev.off()

# =========================================================
# Extract HR tables
# =========================================================

high_sum <- summary(cox_high_ipi)
low_sum  <- summary(cox_low_ipi)

ipi_km_results <- data.frame(
  Subgroup = c("High IPI", "Low/Intermediate IPI"),
  N = c(nrow(df_high_ipi), nrow(df_low_ipi)),
  Events = c(
    sum(df_high_ipi$PFS_event, na.rm = TRUE),
    sum(df_low_ipi$PFS_event, na.rm = TRUE)
  ),
  HR = c(
    exp(coef(cox_high_ipi)),
    exp(coef(cox_low_ipi))
  ),
  Lower95 = c(
    exp(confint(cox_high_ipi))[1],
    exp(confint(cox_low_ipi))[1]
  ),
  Upper95 = c(
    exp(confint(cox_high_ipi))[2],
    exp(confint(cox_low_ipi))[2]
  ),
  P_value = c(
    high_sum$coefficients[,"Pr(>|z|)"],
    low_sum$coefficients[,"Pr(>|z|)"]
  )
)

ipi_km_results

write.csv(
  ipi_km_results,
  "results/GSE31312_treatment_context/IPI_specific_KM_results.csv",
  row.names = FALSE
)

# =========================================
# Save results
# =========================================

write.csv(
  ipi_specific_results,
  "results/GSE31312_treatment_context/IPI_specific_TPS_PFS_results.csv",
  row.names = FALSE
)

write.csv(
  response_pathway_cor,
  "results/GSE31312_treatment_context/response_specific_TPS_pathway_correlations.csv",
  row.names = FALSE
)

saveRDS(
  bio_df,
  "results/GSE31312_treatment_context/GSE31312_treatment_context_TPS.rds"
)