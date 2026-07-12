# =========================================
# GSE10846 survival validation
# Therapy-proliferation signature
# =========================================

library(GEOquery)
library(limma)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

# =========================================
# Download GEO cohort
# =========================================

gse <- getGEO("GSE10846", GSEMatrix = TRUE)
gse <- gse[[1]]

expr <- exprs(gse)
pheno <- pData(gse)
features <- fData(gse)

# =========================================
# Inspect annotation columns
# =========================================

colnames(features)

# =========================================
# Find gene symbol column
# =========================================

symbol_col <- grep(
  "symbol|gene.symbol|gene symbol",
  colnames(features),
  ignore.case = TRUE,
  value = TRUE
)[1]

symbol_col

features$gene_symbol <- features[[symbol_col]]

valid <- !is.na(features$gene_symbol) &
  features$gene_symbol != ""

expr_valid <- expr[valid, ]
features_valid <- features[valid, ]

features_valid$gene_symbol <- sapply(
  strsplit(as.character(features_valid$gene_symbol), " /// "),
  function(x) trimws(x[1])
)

expr_gene <- avereps(
  expr_valid,
  ID = features_valid$gene_symbol
)

# =========================================
# Therapy-persistence genes
# =========================================

therapy_genes <- c(
  "MKI67","RRM2","CHEK1","FOXP1",
  "TCL1A","TCL1B","HMGB2","FBL",
  "PGAM1","SPIB","MCM2","MCM3",
  "MCM4","MCM5","MCM6","PCNA","TYMS"
)

sig_present <- therapy_genes[
  therapy_genes %in% rownames(expr_gene)
]

sig_missing <- setdiff(
  therapy_genes,
  sig_present
)

sig_present
sig_missing

# =========================================
# Build proliferation score
# =========================================

expr_sig <- expr_gene[
  sig_present,
  ,
  drop = FALSE
]

expr_sig_z <- t(scale(t(expr_sig)))

pheno$TherapyPersistenceScore <- colMeans(
  expr_sig_z,
  na.rm = TRUE
)

pheno$TherapyPersistenceGroup <- ifelse(
  pheno$TherapyPersistenceScore >= median(
    pheno$TherapyPersistenceScore,
    na.rm = TRUE
  ),
  "High",
  "Low"
)

table(pheno$TherapyPersistenceGroup)

# =========================================
# Inspect phenotype columns
# =========================================

colnames(pheno)

grep(
  "surv|death|event|status|follow|time|os|pfs|efs|dfs|progress",
  colnames(pheno),
  ignore.case = TRUE,
  value = TRUE
)

# =========================================
# Load clinical survival annotations
# =========================================

library(readxl)

clinical <- read_excel(
  "data/external_validation/mmc6.xlsx"
)

colnames(clinical)
head(clinical)

# =========================================
# Merge GEO expression + clinical survival
# =========================================

merged <- merge(
  pheno,
  clinical,
  by.x = "geo_accession",
  by.y = "Sample"
)

dim(merged)

table(merged$TherapyPersistenceGroup)

merged$OS_time <- as.numeric(merged$`futime(year)`)
merged$OS_event <- as.numeric(merged$fustaut)

fit <- survfit(
  Surv(OS_time, OS_event) ~ TherapyPersistenceGroup,
  data = merged
)

ggsurvplot(
  fit,
  data = merged,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("steelblue", "firebrick"),
  title = "Overall survival by TherapyPersistenceScore"
)

cox_model <- coxph(
  Surv(OS_time, OS_event) ~ TherapyPersistenceScore,
  data = merged
)

summary(cox_model)


# =========================================
# Prepare covariates
# =========================================

merged$Age_num <- as.numeric(merged$Age)

merged$Stage_clean <- gsub("Stage ", "", merged$Stage)
merged$Stage_clean <- factor(
  merged$Stage_clean,
  levels = c("I", "II", "III", "IV")
)

# COO subtype from GEO annotation
merged$molecular_COO <- gsub(
  "gene expression profiling subgroup: ",
  "",
  merged$characteristics_ch1
)

merged$molecular_COO <- factor(
  merged$molecular_COO,
  levels = c("GCB", "ABC", "UC")
)

# =========================================
# Check covariates
# =========================================

summary(merged$Age_num)

table(merged$Stage_clean, useNA = "ifany")

table(merged$molecular_COO, useNA = "ifany")

colSums(is.na(
  merged[, c(
    "OS_time",
    "OS_event",
    "TherapyPersistenceScore",
    "Age_num",
    "Stage_clean",
    "molecular_COO"
  )]
))

sum(complete.cases(
  merged[, c(
    "OS_time",
    "OS_event",
    "TherapyPersistenceScore",
    "Age_num",
    "Stage_clean",
    "molecular_COO"
  )]
))
# =========================================
# Multivariable Cox model: age + stage
# =========================================

cox_multivariable_age_stage <- coxph(
  Surv(OS_time, OS_event) ~
    TherapyPersistenceScore +
    Age_num +
    Stage_clean,
  data = merged
)

summary(cox_multivariable_age_stage)


# =========================================
# Save locked panel genes
# =========================================

locked_panel <- data.frame(
  gene = sig_present,
  panel = "DLBCL_TherapyPersistenceScore_v1"
)

write.csv(
  locked_panel,
  "results/GSE10846/DLBCL_TherapyPersistenceScore_v1_genes.csv",
  row.names = FALSE
)

write.csv(
  merged,
  "results/GSE10846/GSE10846_survival_merged_with_TPS.csv",
  row.names = FALSE
)


# =========================================
# Forest plot for multivariable Cox model
# =========================================

ggforest(
  cox_multivariable_age_stage,
  data = merged,
  main = "Multivariable Cox model: Overall survival",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 1.0
)

# =========================================
# Save Cox model results
# =========================================

cox_summary <- summary(cox_multivariable_age_stage)

cox_results <- data.frame(
  variable = rownames(cox_summary$coefficients),
  coef = cox_summary$coefficients[, "coef"],
  HR = cox_summary$coefficients[, "exp(coef)"],
  p_value = cox_summary$coefficients[, "Pr(>|z|)"],
  lower_95 = cox_summary$conf.int[, "lower .95"],
  upper_95 = cox_summary$conf.int[, "upper .95"]
)

cox_results

write.csv(
  cox_results,
  "results/GSE10846/GSE10846_multivariable_cox_TPS_age_stage.csv",
  row.names = FALSE
)


# =========================================
# Candidate operational panel
# =========================================

panel_genes <- c(
  "MKI67",
  "RRM2",
  "CHEK1",
  "FOXP1",
  "TCL1A",
  "MCM2",
  "MCM3",
  "MCM4",
  "MCM5",
  "MCM6",
  "PCNA",
  "TYMS"
)

panel_present <- intersect(
  panel_genes,
  rownames(expr_gene)
)

panel_present
length(panel_present)

# =========================================
# Build reduced panel score
# =========================================

expr_panel <- expr_gene[panel_present, , drop = FALSE]

expr_panel_z <- t(scale(t(expr_panel)))

reducedPanelScore <- colMeans(
  expr_panel_z,
  na.rm = TRUE
)

merged$ReducedPanelScore <- reducedPanelScore[
  match(merged$geo_accession, names(reducedPanelScore))
]

summary(merged$ReducedPanelScore)

# =========================================
# Reduced panel survival model
# =========================================

cox_reduced <- coxph(
  Surv(OS_time, OS_event) ~ ReducedPanelScore,
  data = merged
)

summary(cox_reduced)

# =========================================
# Reduced panel KM plot
# =========================================

merged$ReducedPanelGroup <- ifelse(
  merged$ReducedPanelScore >
    median(merged$ReducedPanelScore, na.rm = TRUE),
  "High",
  "Low"
)

fit_reduced <- survfit(
  Surv(OS_time, OS_event) ~ ReducedPanelGroup,
  data = merged
)

ggsurvplot(
  fit_reduced,
  data = merged,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("steelblue", "firebrick"),
  title = "Overall survival by reduced resistance panel"
)

# =========================================
# RT-relevant mechanism scores
# =========================================

rt_checkpoint_genes <- c("CHEK1", "RRM2", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "PCNA")
rt_dna_synthesis_genes <- c("TYMS", "PCNA", "RRM2")
rt_proliferation_genes <- c("MKI67", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6")

score_bulk_module <- function(expr_gene, genes) {
  present <- intersect(genes, rownames(expr_gene))
  message("Genes used: ", paste(present, collapse = ", "))
  z <- t(scale(t(expr_gene[present, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

rt_checkpoint_score <- score_bulk_module(expr_gene, rt_checkpoint_genes)
rt_dna_synthesis_score <- score_bulk_module(expr_gene, rt_dna_synthesis_genes)
rt_proliferation_score <- score_bulk_module(expr_gene, rt_proliferation_genes)

merged$RT_CheckpointScore <- rt_checkpoint_score[match(merged$geo_accession, names(rt_checkpoint_score))]
merged$RT_DNASynthesisScore <- rt_dna_synthesis_score[match(merged$geo_accession, names(rt_dna_synthesis_score))]
merged$RT_ProliferationScore <- rt_proliferation_score[match(merged$geo_accession, names(rt_proliferation_score))]

# =========================================
# Compare RT-relevant biology by TPS group
# =========================================

wilcox.test(RT_CheckpointScore ~ TherapyPersistenceGroup, data = merged)
wilcox.test(RT_DNASynthesisScore ~ TherapyPersistenceGroup, data = merged)
wilcox.test(RT_ProliferationScore ~ TherapyPersistenceGroup, data = merged)

boxplot(
  RT_CheckpointScore ~ TherapyPersistenceGroup,
  data = merged,
  main = "RT checkpoint biology by TherapyPersistenceGroup",
  ylab = "RT checkpoint score"
)

boxplot(
  RT_DNASynthesisScore ~ TherapyPersistenceGroup,
  data = merged,
  main = "DNA synthesis / repair metabolism by TherapyPersistenceGroup",
  ylab = "DNA synthesis score"
)

boxplot(
  RT_ProliferationScore ~ TherapyPersistenceGroup,
  data = merged,
  main = "Proliferation-linked RT biology by TherapyPersistenceGroup",
  ylab = "RT proliferation score"
)

# =========================================
# Survival models for RT-relevant biology
# =========================================

cox_rt_checkpoint <- coxph(
  Surv(OS_time, OS_event) ~ RT_CheckpointScore + Age_num + Stage_clean,
  data = merged
)

cox_rt_dna <- coxph(
  Surv(OS_time, OS_event) ~ RT_DNASynthesisScore + Age_num + Stage_clean,
  data = merged
)

cox_rt_prolif <- coxph(
  Surv(OS_time, OS_event) ~ RT_ProliferationScore + Age_num + Stage_clean,
  data = merged
)

summary(cox_rt_checkpoint)
summary(cox_rt_dna)
summary(cox_rt_prolif)


# =========================================
# Treatment-modality biology modules
# =========================================

rt_module <- c("CHEK1", "RRM2", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "PCNA", "TYMS", "MKI67")

chemo_resistance_module <- c("BCL2", "BCL2L1", "MCL1", "BIRC5", "MYC", "MKI67", "TOP2A", "TYMS")

car_t_bispecific_module <- c("CD19", "MS4A1", "CD22", "CD79A", "CD79B", "B2M", "HLA-A", "HLA-B", "HLA-C", "CD58")

adc_polivy_module <- c("CD79A", "CD79B", "MS4A1", "CD19", "CD22")

targeted_abc_module <- c("FOXP1", "IRF4", "MYD88", "CD79B", "NFKB1", "NFKB2", "REL", "BCL2", "BTK", "CARD11")

score_bulk_module <- function(expr_gene, genes) {
  present <- intersect(genes, rownames(expr_gene))
  message("Genes used: ", paste(present, collapse = ", "))
  z <- t(scale(t(expr_gene[present, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

modality_scores <- data.frame(
  geo_accession = colnames(expr_gene),
  RT_DNA_Damage = score_bulk_module(expr_gene, rt_module),
  Chemo_Resistance = score_bulk_module(expr_gene, chemo_resistance_module),
  CAR_T_Bispecific_Visibility = score_bulk_module(expr_gene, car_t_bispecific_module),
  Polivy_ADC_Target = score_bulk_module(expr_gene, adc_polivy_module),
  Targeted_ABC_Biology = score_bulk_module(expr_gene, targeted_abc_module)
)

merged <- merge(
  merged,
  modality_scores,
  by = "geo_accession",
  all.x = TRUE
)

# Compare each modality module by TherapyPersistenceGroup
wilcox.test(RT_DNA_Damage ~ TherapyPersistenceGroup, data = merged)
wilcox.test(Chemo_Resistance ~ TherapyPersistenceGroup, data = merged)
wilcox.test(CAR_T_Bispecific_Visibility ~ TherapyPersistenceGroup, data = merged)
wilcox.test(Polivy_ADC_Target ~ TherapyPersistenceGroup, data = merged)
wilcox.test(Targeted_ABC_Biology ~ TherapyPersistenceGroup, data = merged)


# =========================================
# TPS-high subtype discovery
# =========================================

high_tps <- subset(
  merged,
  TherapyPersistenceGroup == "High"
)

modality_matrix <- high_tps[, c(
  "RT_DNA_Damage",
  "Chemo_Resistance",
  "CAR_T_Bispecific_Visibility",
  "Polivy_ADC_Target",
  "Targeted_ABC_Biology"
)]

rownames(modality_matrix) <- high_tps$geo_accession

heatmap(
  as.matrix(modality_matrix),
  scale = "column",
  margins = c(8, 8),
  main = "TPS-high therapeutic vulnerability states"
)

# =========================================
# Cluster TPS-high tumors
# =========================================

set.seed(123)

k <- 3

clusters <- kmeans(
  scale(modality_matrix),
  centers = k
)

high_tps$TPS_subtype <- factor(clusters$cluster)

table(high_tps$TPS_subtype)

aggregate(
  modality_matrix,
  by = list(Subtype = high_tps$TPS_subtype),
  mean
)


# =========================================
# Survival by TPS-high subtype
# =========================================

fit_subtypes <- survfit(
  Surv(OS_time, OS_event) ~ TPS_subtype,
  data = high_tps
)

ggsurvplot(
  fit_subtypes,
  data = high_tps,
  pval = TRUE,
  risk.table = TRUE,
  title = "Overall survival by TPS-high biologic subtype"
)

aggregate(
  modality_matrix,
  by = list(Subtype = high_tps$TPS_subtype),
  median
)


# =========================================
# Publication heatmap: TPS-high vulnerability states
# =========================================

library(pheatmap)

subtype_annotation <- data.frame(
  TPS_subtype = high_tps$TPS_subtype
)

rownames(subtype_annotation) <- high_tps$geo_accession

pheatmap(
  t(scale(t(as.matrix(modality_matrix)))),
  annotation_row = subtype_annotation,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "TPS-high therapeutic vulnerability states",
  filename = "results/GSE10846/TPS_high_vulnerability_heatmap.pdf",
  width = 8,
  height = 10
)

# =========================================
# Subtype biology summary table
# =========================================

subtype_summary_mean <- aggregate(
  modality_matrix,
  by = list(Subtype = high_tps$TPS_subtype),
  mean
)

subtype_summary_median <- aggregate(
  modality_matrix,
  by = list(Subtype = high_tps$TPS_subtype),
  median
)

write.csv(
  subtype_summary_mean,
  "results/GSE10846/TPS_high_subtype_mean_scores.csv",
  row.names = FALSE
)

write.csv(
  subtype_summary_median,
  "results/GSE10846/TPS_high_subtype_median_scores.csv",
  row.names = FALSE
)

summary(cox_multivariable_age_stage)$concordance
summary(cox_reduced)$concordance

cox_stage_only <- coxph(Surv(OS_time, OS_event) ~ Age_num + Stage_clean, data = merged)

anova(cox_stage_only, cox_multivariable_age_stage, test = "LRT")
AIC(cox_stage_only, cox_multivariable_age_stage)

install.packages("timeROC")
library(timeROC)

roc3 <- timeROC(
  T = merged$OS_time,
  delta = merged$OS_event,
  marker = merged$TherapyPersistenceScore,
  cause = 1,
  times = c(2, 5),
  iid = TRUE
)

roc3$AUC
plot(roc3, time = 2)
plot(roc3, time = 5, add = TRUE)

# =========================================
# Multimodality vulnerability map
# =========================================

vulnerability_map <- data.frame(
  TPS_subtype = c("Subtype 1", "Subtype 2", "Subtype 3"),
  Dominant_biology = c(
    "ADC-retained / immune-visible persistence",
    "Replication-stressed / immune-low persistence",
    "Aggressive hybrid: RT-stress + chemo-resistance + ABC biology"
  ),
  Potential_vulnerability = c(
    "Polivy-like ADC / bispecific consideration",
    "RT + radiosensitization / CHK1-ATR axis",
    "Combination strategy: RT + ADC/immune/targeted therapy"
  )
)

write.csv(
  vulnerability_map,
  "results/GSE10846/TPS_high_vulnerability_map.csv",
  row.names = FALSE
)

vulnerability_map

library(gridExtra)
library(grid)

pdf(
  "results/GSE10846/TPS_high_vulnerability_map_table.pdf",
  width = 11,
  height = 4
)

grid.table(vulnerability_map)

dev.off()

# =========================================
# Continuous vulnerability ranking per patient
# =========================================

vulnerability_scores <- merged[, c(
  "geo_accession",
  "TherapyPersistenceScore",
  "TherapyPersistenceGroup",
  "RT_DNA_Damage",
  "Chemo_Resistance",
  "CAR_T_Bispecific_Visibility",
  "Polivy_ADC_Target",
  "Targeted_ABC_Biology"
)]

score_cols <- c(
  "RT_DNA_Damage",
  "Chemo_Resistance",
  "CAR_T_Bispecific_Visibility",
  "Polivy_ADC_Target",
  "Targeted_ABC_Biology"
)

vulnerability_scores$Dominant_axis <- score_cols[
  max.col(vulnerability_scores[, score_cols], ties.method = "first")
]

table(vulnerability_scores$Dominant_axis)

write.csv(
  vulnerability_scores,
  "results/GSE10846/patient_level_continuous_vulnerability_axes.csv",
  row.names = FALSE
)

library(ggplot2)

ggplot(
  vulnerability_scores,
  aes(x = Dominant_axis, fill = TherapyPersistenceGroup)
) +
  geom_bar(position = "dodge") +
  theme_classic(base_size = 14) +
  labs(
    title = "Dominant therapeutic vulnerability axis per patient",
    x = "Dominant vulnerability axis",
    y = "Number of patients"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# =========================================
# Radar / spider plot
# =========================================

install.packages("fmsb")
library(fmsb)

radar_df <- subtype_summary_median
rownames(radar_df) <- paste0("Subtype_", radar_df$Subtype)
radar_df$Subtype <- NULL

radar_plot_df <- rbind(
  apply(radar_df, 2, max),
  apply(radar_df, 2, min),
  radar_df
)

pdf("results/GSE10846/TPS_high_subtype_radar_plot.pdf", width = 8, height = 8)

radarchart(
  radar_plot_df,
  axistype = 1,
  title = "TPS-high therapeutic vulnerability subtypes"
)

legend(
  "topright",
  legend = rownames(radar_df),
  bty = "n",
  lty = 1
)

dev.off()

# =========================================
# Survival plot save
# =========================================

pdf("results/GSE10846/TPS_high_subtype_survival.pdf", width = 8, height = 7)

ggsurvplot(
  fit_subtypes,
  data = high_tps,
  pval = TRUE,
  risk.table = TRUE,
  title = "Overall survival by TPS-high biologic subtype"
)

dev.off()