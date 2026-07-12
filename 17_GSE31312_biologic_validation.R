# =========================================================
# GSE31312 biologic external validation
# Conserved TPS biology despite heterogeneous survival
# =========================================================

library(msigdbr)
library(GSVA)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape2)

dir.create(
  "results/GSE31312_biologic_validation",
  recursive = TRUE,
  showWarnings = FALSE
)

# =========================================================
# REQUIRED OBJECTS
# =========================================================
# Assumes already loaded from prior script:
#
# expr_gene
# analysis31312
# TPS12_Score
# TPS_full_Score
#
# =========================================================

dim(expr_gene)
dim(analysis31312)

# =========================================================
# 1. Hallmark gene sets
# =========================================================

hallmark_df <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

hallmark_list <- split(
  hallmark_df$gene_symbol,
  hallmark_df$gs_name
)

names(hallmark_list)

# =========================================================
# 2. Pathways of interest
# =========================================================

selected_sets <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_HYPOXIA",
  "HALLMARK_P53_PATHWAY"
)

selected_hallmarks <- hallmark_list[selected_sets]

selected_sets

# =========================================================
# 3. GSVA / ssGSEA scoring - new GSVA syntax
# =========================================================

library(GSVA)

ssgsea_param <- ssgseaParam(
  exprData = as.matrix(expr_gene),
  geneSets = selected_hallmarks,
  normalize = TRUE
)

hallmark_scores <- gsva(ssgsea_param)

dim(hallmark_scores)

# =========================================================
# 4. Match samples to analysis cohort
# =========================================================

hallmark_scores_df <- as.data.frame(t(hallmark_scores))

hallmark_scores_df$geo_accession <- rownames(hallmark_scores_df)

analysis31312$geo_accession <- as.character(
  analysis31312$geo_accession
)

bio_df <- merge(
  analysis31312,
  hallmark_scores_df,
  by = "geo_accession"
)

dim(bio_df)

# =========================================================
# 5. Correlation analysis with TPS
# =========================================================

pathway_correlations <- data.frame()

for (pathway in selected_sets) {
  
  cor_test <- cor.test(
    bio_df$TPS_full_Score,
    bio_df[[pathway]],
    method = "spearman"
  )
  
  tmp <- data.frame(
    Pathway = pathway,
    Spearman_rho = cor_test$estimate,
    P_value = cor_test$p.value
  )
  
  pathway_correlations <- rbind(
    pathway_correlations,
    tmp
  )
}

pathway_correlations$FDR <- p.adjust(
  pathway_correlations$P_value,
  method = "fdr"
)

pathway_correlations <- pathway_correlations %>%
  arrange(desc(Spearman_rho))

pathway_correlations

# =========================================================
# 6. Plot pathway correlations
# =========================================================

p_cor <- ggplot(
  pathway_correlations,
  aes(
    x = reorder(Pathway, Spearman_rho),
    y = Spearman_rho
  )
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "GSE31312 biologic validation of TPS",
    x = "",
    y = "Spearman correlation with TPS_full"
  ) +
  theme_bw(base_size = 12)

print(p_cor)

pdf(
  "results/GSE31312_biologic_validation/TPS_pathway_correlations.pdf",
  width = 8,
  height = 5
)

print(p_cor)

dev.off()

# =========================================================
# 7. TPS-high vs TPS-low pathway comparison
# =========================================================

pathway_group_stats <- data.frame()

for (pathway in selected_sets) {
  
  test <- wilcox.test(
    bio_df[[pathway]] ~ bio_df$TPS_full_Group
  )
  
  med_high <- median(
    bio_df[[pathway]][bio_df$TPS_full_Group == "High"],
    na.rm = TRUE
  )
  
  med_low <- median(
    bio_df[[pathway]][bio_df$TPS_full_Group == "Low"],
    na.rm = TRUE
  )
  
  tmp <- data.frame(
    Pathway = pathway,
    Median_High = med_high,
    Median_Low = med_low,
    Delta = med_high - med_low,
    P_value = test$p.value
  )
  
  pathway_group_stats <- rbind(
    pathway_group_stats,
    tmp
  )
}

pathway_group_stats$FDR <- p.adjust(
  pathway_group_stats$P_value,
  method = "fdr"
)

pathway_group_stats <- pathway_group_stats %>%
  arrange(desc(Delta))

pathway_group_stats

# =========================================================
# 8. Heatmap of pathway activity
# =========================================================

heatmap_mat <- bio_df[, selected_sets]

rownames(heatmap_mat) <- bio_df$geo_accession

heatmap_mat <- scale(heatmap_mat)

annotation_df <- data.frame(
  TPS = bio_df$TPS_full_Group
)

rownames(annotation_df) <- bio_df$geo_accession

pdf(
  "results/GSE31312_biologic_validation/TPS_hallmark_heatmap.pdf",
  width = 8,
  height = 10
)

pheatmap(
  t(heatmap_mat),
  annotation_col = annotation_df,
  show_colnames = FALSE,
  fontsize_row = 10,
  clustering_method = "ward.D2"
)

dev.off()

# =========================================================
# 9. COO subtype analysis
# =========================================================

grep(
  "molecular|subtype|coo|gcb|abc",
  colnames(bio_df),
  ignore.case = TRUE,
  value = TRUE
)

# modify this line if needed depending on column name
if ("molecular_subtype_by_gene_expression:ch1" %in% colnames(bio_df)) {
  
  bio_df$COO <- bio_df$`molecular_subtype_by_gene_expression:ch1`
  
} else {
  
  bio_df$COO <- NA
}

table(bio_df$COO, useNA = "ifany")

# =========================================================
# 10. TPS distribution by COO
# =========================================================

p_coo <- ggplot(
  bio_df,
  aes(
    x = COO,
    y = TPS_full_Score,
    fill = COO
  )
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA
  ) +
  labs(
    title = "TPS by COO subtype",
    x = "COO subtype",
    y = "TPS_full score"
  ) +
  theme_bw(base_size = 12)

print(p_coo)

pdf(
  "results/GSE31312_biologic_validation/TPS_by_COO.pdf",
  width = 6,
  height = 5
)

print(p_coo)

dev.off()

# =========================================================
# 11. Response association
# =========================================================

table(bio_df$Respon, useNA = "ifany")

response_stats <- bio_df %>%
  group_by(Respon) %>%
  summarise(
    n = n(),
    median_TPS = median(TPS_full_Score, na.rm = TRUE),
    mean_TPS = mean(TPS_full_Score, na.rm = TRUE)
  )

response_stats

p_response <- ggplot(
  bio_df,
  aes(
    x = Respon,
    y = TPS_full_Score,
    fill = Respon
  )
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA
  ) +
  labs(
    title = "TPS by response category",
    x = "Clinical response",
    y = "TPS_full score"
  ) +
  theme_bw(base_size = 12)

print(p_response)

pdf(
  "results/GSE31312_biologic_validation/TPS_by_response.pdf",
  width = 7,
  height = 5
)

print(p_response)

dev.off()

# =========================================================
# 12. Early progression analysis
# =========================================================

bio_df$Early_progression <- ifelse(
  bio_df$PFS_event == 1 &
    bio_df$PFS_time_years < 1,
  "Early progression",
  "No early progression"
)

table(bio_df$Early_progression)

wilcox.test(
  TPS_full_Score ~ Early_progression,
  data = bio_df
)

p_early <- ggplot(
  bio_df,
  aes(
    x = Early_progression,
    y = TPS_full_Score,
    fill = Early_progression
  )
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA
  ) +
  labs(
    title = "TPS and early progression",
    x = "",
    y = "TPS_full score"
  ) +
  theme_bw(base_size = 12)

print(p_early)

pdf(
  "results/GSE31312_biologic_validation/TPS_early_progression.pdf",
  width = 6,
  height = 5
)

print(p_early)

dev.off()

# =========================================================
# 13. Save outputs
# =========================================================

write.csv(
  pathway_correlations,
  "results/GSE31312_biologic_validation/pathway_correlations.csv",
  row.names = FALSE
)

write.csv(
  pathway_group_stats,
  "results/GSE31312_biologic_validation/pathway_group_comparisons.csv",
  row.names = FALSE
)

write.csv(
  response_stats,
  "results/GSE31312_biologic_validation/response_associations.csv",
  row.names = FALSE
)

write.csv(
  bio_df,
  "results/GSE31312_biologic_validation/GSE31312_biologic_validation_table.csv",
  row.names = FALSE
)

saveRDS(
  bio_df,
  "results/GSE31312_biologic_validation/GSE31312_biologic_validation_table.rds"
)

# =========================================================
# 14. Final console outputs
# =========================================================

pathway_correlations

pathway_group_stats

response_stats

summary(bio_df$TPS_full_Score)

table(bio_df$TPS_full_Group)

table(bio_df$Early_progression)