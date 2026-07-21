# ============================================================
# FIGURE 8 — TPS12 in established DLBCL molecular contexts
#
# Panel A: TPS12 across GEP-defined COO — GSE31312
# Panel B: Independent COO validation — GSE23501
# Panel C: TPS12 vs DZ/DHIT-like transcriptional variation
#          within GCB-DLBCL — GSE31312
# Panel D: TPS12 across published genetic subtypes
#          MCD, BN2, N1, EZB, Other — NCI DLBCL cohort
#
# Final output:
# results/main_figures/Figure_8_TPS12_molecular_context_FINAL.tiff
# ============================================================


# ============================================================
# 0. Libraries
# ============================================================

library(tidyverse)
library(patchwork)
library(GEOquery)
library(Biobase)
library(readxl)

set.seed(123)


# ============================================================
# 1. Output directories
# ============================================================

outdir_context <-
  "results/TPS12_molecular_classification_context"

outdir_genetic <-
  "results/TPS12_genetic_subgroups_NCI"

outdir_final <-
  "results/main_figures"

dir.create(
  outdir_context,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  outdir_genetic,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  outdir_final,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# PART I
# COO CONTEXT — GSE31312 AND GSE23501
# ============================================================


# ============================================================
# 2. Input metadata
# ============================================================

file31312 <-
  "./results/GSE31312_context_dependence/GSE31312_context_harmonized_metadata.csv"

file23501 <-
  "./results/GSE23501_context_dependence/GSE23501_TPS12_context_metadata.csv"

stopifnot(
  file.exists(file31312),
  file.exists(file23501)
)

dat31312 <- read_csv(
  file31312,
  show_col_types = FALSE
)

dat23501 <- read_csv(
  file23501,
  show_col_types = FALSE
)

cat("\n===== GSE31312 DIMENSIONS =====\n")
print(dim(dat31312))

cat("\n===== GSE23501 DIMENSIONS =====\n")
print(dim(dat23501))


# ============================================================
# 3. Clean GSE31312
# ============================================================

d31312 <- dat31312 %>%
  transmute(
    sample_id = as.character(geo_accession),
    TPS12 = as.numeric(TPS12_Score),
    COO = as.character(COO_context)
  ) %>%
  mutate(
    COO = case_when(
      str_detect(
        str_to_upper(COO),
        "ABC"
      ) ~ "ABC",
      
      str_detect(
        str_to_upper(COO),
        "GCB"
      ) ~ "GCB",
      
      TRUE ~ "Unclassified"
    ),
    
    COO = factor(
      COO,
      levels = c(
        "GCB",
        "ABC",
        "Unclassified"
      )
    )
  ) %>%
  filter(
    is.finite(TPS12)
  )


# ============================================================
# 4. Clean GSE23501
# ============================================================

d23501 <- dat23501 %>%
  transmute(
    TPS12 = as.numeric(TPS12),
    COO = as.character(COO)
  ) %>%
  mutate(
    COO = case_when(
      str_detect(
        str_to_upper(COO),
        "ABC"
      ) ~ "ABC",
      
      str_detect(
        str_to_upper(COO),
        "GCB"
      ) ~ "GCB",
      
      TRUE ~ "Unclassified"
    ),
    
    COO = factor(
      COO,
      levels = c(
        "GCB",
        "ABC",
        "Unclassified"
      )
    )
  ) %>%
  filter(
    is.finite(TPS12)
  )


# ============================================================
# 5. COO counts
# ============================================================

cat("\n===== GSE31312 COO COUNTS =====\n")
print(table(d31312$COO))

cat("\n===== GSE23501 COO COUNTS =====\n")
print(table(d23501$COO))


# ============================================================
# 6. COO statistical testing
# ============================================================

kw31312 <- kruskal.test(
  TPS12 ~ COO,
  data = d31312
)

pair31312 <- pairwise.wilcox.test(
  d31312$TPS12,
  d31312$COO,
  p.adjust.method = "BH",
  exact = FALSE
)

kw23501 <- kruskal.test(
  TPS12 ~ COO,
  data = d23501
)

pair23501 <- pairwise.wilcox.test(
  d23501$TPS12,
  d23501$COO,
  p.adjust.method = "BH",
  exact = FALSE
)

cat("\n===== GSE31312 KRUSKAL-WALLIS =====\n")
print(kw31312)

cat("\n===== GSE31312 PAIRWISE WILCOXON =====\n")
print(pair31312)

cat("\n===== GSE23501 KRUSKAL-WALLIS =====\n")
print(kw23501)

cat("\n===== GSE23501 PAIRWISE WILCOXON =====\n")
print(pair23501)


# ============================================================
# PART II
# DZ / DHIT-LIKE TRANSCRIPTIONAL CONTEXT
# ============================================================


# ============================================================
# 7. Load GSE31312 expression
# ============================================================

gse_list <- getGEO(
  "GSE31312",
  GSEMatrix = TRUE
)

if (is.list(gse_list)) {
  eset <- gse_list[[1]]
} else {
  eset <- gse_list
}

expr_gse31312 <- exprs(eset)

cat("\n===== GSE31312 EXPRESSION DIMENSIONS =====\n")
print(dim(expr_gse31312))

cat("\n===== GSE31312 PLATFORM =====\n")
print(annotation(eset))


# ============================================================
# 8. DZ / DHIT-like transcriptional genes
# ============================================================

dz_up <- c(
  "OR13A1",
  "MYC",
  "SLC25A27",
  "ALOX5",
  "TNFSF8",
  "PEG10",
  "GAMT",
  "SNHG19",
  "QRSL1",
  "RGCC",
  "JCHAIN",
  "CD24",
  "AFMID",
  "SMIM14",
  "SYBU"
)

dz_down <- c(
  "GPR137B",
  "CDK5R1",
  "LY75",
  "VASP",
  "RFFL",
  "MIR155HG",
  "VOPP1",
  "BATF",
  "STAT3",
  "IRF4",
  "SGPP2",
  "CD80",
  "SEMA7A",
  "EBI3",
  "IL21R"
)


# ============================================================
# 9. Map GPL570 probes to unique gene symbols
# ============================================================

feature_annot <- fData(eset) %>%
  as.data.frame() %>%
  rownames_to_column(
    "probe_id"
  ) %>%
  transmute(
    probe_id,
    gene_symbol_raw =
      as.character(
        `Gene Symbol`
      )
  ) %>%
  filter(
    !is.na(gene_symbol_raw),
    gene_symbol_raw != ""
  )

feature_annot_unique <- feature_annot %>%
  mutate(
    gene_symbol =
      str_trim(
        gene_symbol_raw
      )
  ) %>%
  filter(
    !str_detect(
      gene_symbol,
      "///"
    )
  ) %>%
  select(
    probe_id,
    gene_symbol
  )

cat(
  "\nUnique probe-to-gene mappings retained:",
  nrow(feature_annot_unique),
  "\n"
)


# ============================================================
# 10. Collapse probes to gene-level expression
# ============================================================

expr_probe_df <- as.data.frame(
  expr_gse31312
) %>%
  rownames_to_column(
    "probe_id"
  ) %>%
  inner_join(
    feature_annot_unique,
    by = "probe_id"
  )

expr_gene_long <- expr_probe_df %>%
  select(
    -probe_id
  ) %>%
  pivot_longer(
    cols = -gene_symbol,
    names_to = "sample_id",
    values_to = "expression"
  ) %>%
  group_by(
    gene_symbol,
    sample_id
  ) %>%
  summarise(
    expression =
      median(
        expression,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

expr_gene_mat <- expr_gene_long %>%
  pivot_wider(
    names_from = sample_id,
    values_from = expression
  ) %>%
  column_to_rownames(
    "gene_symbol"
  ) %>%
  as.matrix()

cat("\n===== GENE-LEVEL MATRIX DIMENSIONS =====\n")
print(dim(expr_gene_mat))


# ============================================================
# 11. DZ gene availability
# ============================================================

dz_up_present <- intersect(
  dz_up,
  rownames(expr_gene_mat)
)

dz_down_present <- intersect(
  dz_down,
  rownames(expr_gene_mat)
)

dz_up_missing <- setdiff(
  dz_up,
  rownames(expr_gene_mat)
)

dz_down_missing <- setdiff(
  dz_down,
  rownames(expr_gene_mat)
)

cat(
  "\nDZ upregulated genes present:",
  length(dz_up_present),
  "of",
  length(dz_up),
  "\n"
)

print(dz_up_present)

cat("\nDZ upregulated genes missing:\n")
print(dz_up_missing)

cat(
  "\nDZ downregulated genes present:",
  length(dz_down_present),
  "of",
  length(dz_down),
  "\n"
)

print(dz_down_present)

cat("\nDZ downregulated genes missing:\n")
print(dz_down_missing)

if (
  length(dz_up_present) < 8 ||
  length(dz_down_present) < 8
) {
  stop(
    "Insufficient DZ/DHIT-like gene representation."
  )
}


# ============================================================
# 12. Calculate exploratory DZ/DHIT-like score
# ============================================================

dz_genes_present <- unique(
  c(
    dz_up_present,
    dz_down_present
  )
)

dz_expr <- expr_gene_mat[
  dz_genes_present,
  ,
  drop = FALSE
]

# Gene-wise z-scoring across samples
dz_expr_z <- t(
  scale(
    t(dz_expr)
  )
)

# Remove genes with non-finite values
dz_expr_z <- dz_expr_z[
  apply(
    dz_expr_z,
    1,
    function(x) {
      all(
        is.finite(x)
      )
    }
  ),
  ,
  drop = FALSE
]

dz_up_final <- intersect(
  dz_up_present,
  rownames(dz_expr_z)
)

dz_down_final <- intersect(
  dz_down_present,
  rownames(dz_expr_z)
)

dz_score <-
  colMeans(
    dz_expr_z[
      dz_up_final,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  ) -
  colMeans(
    dz_expr_z[
      dz_down_final,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )

dz_scores_df <- tibble(
  sample_id =
    names(dz_score),
  
  DZ_DHIT_like_score =
    as.numeric(
      dz_score
    )
)


# ============================================================
# 13. Merge DZ score with TPS12 and COO
# ============================================================

molecular_context <- d31312 %>%
  select(
    sample_id,
    TPS12,
    COO
  ) %>%
  inner_join(
    dz_scores_df,
    by = "sample_id"
  ) %>%
  filter(
    is.finite(TPS12),
    is.finite(
      DZ_DHIT_like_score
    )
  )

cat("\n===== MERGED TPS12 + DZ DIMENSIONS =====\n")
print(dim(molecular_context))

cat("\n===== MERGED COO COUNTS =====\n")
print(table(molecular_context$COO))


# ============================================================
# 14. Restrict DZ analysis to GCB
# ============================================================

gcb_dz <- molecular_context %>%
  filter(
    COO == "GCB"
  )

cat(
  "\nGCB tumors available for DZ/DHIT-like analysis:",
  nrow(gcb_dz),
  "\n"
)

cor_dz_tps <- cor.test(
  gcb_dz$TPS12,
  gcb_dz$DZ_DHIT_like_score,
  method = "spearman",
  exact = FALSE
)

rho_dz <- unname(
  cor_dz_tps$estimate
)

p_dz <- cor_dz_tps$p.value

cat(
  "\n===== TPS12 VS DZ/DHIT-LIKE SCORE IN GCB =====\n",
  "Spearman rho = ",
  sprintf(
    "%.3f",
    rho_dz
  ),
  "\nP = ",
  sprintf(
    "%.3f",
    p_dz
  ),
  "\n"
)


# ============================================================
# PART III
# NCI GENETIC SUBTYPE ANALYSIS
# ============================================================


# ============================================================
# 15. NCI input files
# ============================================================

nci_expr_file <-
  "data/NCI_DLBCL/RNAseq_gene_expression_562.txt"

nci_subtype_file <-
  "data/NCI_DLBCL/Supplementary_Appendix_2.xlsx"

stopifnot(
  file.exists(nci_expr_file),
  file.exists(nci_subtype_file)
)


# ============================================================
# 16. TPS12 genes
# ============================================================

tps12_genes <- c(
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


# ============================================================
# 17. Read NCI RNA-seq matrix
# ============================================================

nci_expr_raw <- read.delim(
  nci_expr_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n===== NCI EXPRESSION DIMENSIONS =====\n")
print(dim(nci_expr_raw))

cat("\n===== NCI FIRST COLUMNS =====\n")
print(
  head(
    names(nci_expr_raw),
    10
  )
)


# ============================================================
# 18. Build NCI long-format gene expression
# ============================================================

if (
  !"Gene" %in%
  names(nci_expr_raw)
) {
  stop(
    "Expected gene-symbol column 'Gene' not found in NCI expression file."
  )
}

nci_expr_clean <- nci_expr_raw %>%
  rename(
    gene_symbol = Gene
  ) %>%
  filter(
    !is.na(gene_symbol),
    gene_symbol != ""
  )

nci_expr_long <- nci_expr_clean %>%
  pivot_longer(
    cols = -c(
      gene_symbol,
      Accession,
      Gene_ID
    ),
    names_to = "sample_id",
    values_to = "expr_value"
  ) %>%
  mutate(
    expr_value =
      as.numeric(
        expr_value
      )
  ) %>%
  filter(
    is.finite(
      expr_value
    )
  )


# ============================================================
# 19. Check TPS12 gene availability in NCI cohort
# ============================================================

nci_available_genes <- unique(
  nci_expr_long$gene_symbol
)

tps_present <- intersect(
  tps12_genes,
  nci_available_genes
)

tps_missing <- setdiff(
  tps12_genes,
  nci_available_genes
)

cat(
  "\nTPS12 genes present:",
  length(tps_present),
  "of 12\n"
)

print(tps_present)

cat("\nTPS12 genes missing:\n")
print(tps_missing)

if (
  length(tps_present) < 10
) {
  stop(
    "Too few TPS12 genes available in NCI cohort."
  )
}


# ============================================================
# 20. Calculate TPS12 in NCI cohort
# ============================================================

nci_tps_expr <- nci_expr_long %>%
  filter(
    gene_symbol %in%
      tps_present
  ) %>%
  group_by(
    gene_symbol
  ) %>%
  mutate(
    z =
      as.numeric(
        scale(
          expr_value
        )
      )
  ) %>%
  ungroup()

nci_tps_scores <- nci_tps_expr %>%
  group_by(
    sample_id
  ) %>%
  summarise(
    TPS12 =
      mean(
        z,
        na.rm = TRUE
      ),
    .groups = "drop"
  )


# ============================================================
# 21. Read official NCI subtype metadata
# ============================================================

clinical_subtypes <- read_excel(
  nci_subtype_file,
  sheet =
    "Tab S9 Characteristics DLBCL"
)

cat("\n===== NCI TAB S9 DIMENSIONS =====\n")
print(dim(clinical_subtypes))


# ============================================================
# 22. Merge TPS12 with official genetic subtype assignments
# ============================================================

genetic_meta <- clinical_subtypes %>%
  transmute(
    sample_id =
      `dbGaP submitted subject ID`,
    
    COO =
      `Gene Expression Subgroup`,
    
    Genetic_Subtype =
      `Genetic Subtype`
  ) %>%
  filter(
    !is.na(sample_id),
    !is.na(Genetic_Subtype)
  ) %>%
  distinct(
    sample_id,
    .keep_all = TRUE
  )

genetic_tps <- nci_tps_scores %>%
  inner_join(
    genetic_meta,
    by = "sample_id"
  ) %>%
  filter(
    Genetic_Subtype %in%
      c(
        "MCD",
        "BN2",
        "N1",
        "EZB",
        "Other"
      )
  ) %>%
  mutate(
    Genetic_Subtype =
      factor(
        Genetic_Subtype,
        levels =
          c(
            "MCD",
            "BN2",
            "N1",
            "EZB",
            "Other"
          )
      )
  )

cat(
  "\n===== NCI GENETIC SUBTYPE COUNTS AFTER MERGE =====\n"
)

print(
  table(
    genetic_tps$Genetic_Subtype
  )
)


# ============================================================
# 23. NCI genetic subtype statistical tests
# ============================================================

kw_genetic <- kruskal.test(
  TPS12 ~ Genetic_Subtype,
  data = genetic_tps
)

pairwise_genetic <- pairwise.wilcox.test(
  genetic_tps$TPS12,
  genetic_tps$Genetic_Subtype,
  p.adjust.method = "BH",
  exact = FALSE
)

cat(
  "\n===== KRUSKAL-WALLIS: TPS12 BY GENETIC SUBTYPE =====\n"
)

print(
  kw_genetic
)

cat(
  "\n===== PAIRWISE WILCOXON, BH-ADJUSTED =====\n"
)

print(
  pairwise_genetic
)


# ============================================================
# 24. Genetic subtype summary
# ============================================================

genetic_summary <- genetic_tps %>%
  group_by(
    Genetic_Subtype
  ) %>%
  summarise(
    n = n(),
    
    median_TPS12 =
      median(
        TPS12,
        na.rm = TRUE
      ),
    
    IQR_TPS12 =
      IQR(
        TPS12,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

print(
  genetic_summary
)


# ============================================================
# PART IV
# FIGURE 8 PANELS
# ============================================================


# ============================================================
# 25. Shared theme
# ============================================================

theme_context <- theme_classic(
  base_size = 11
) +
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 11
      ),
    
    plot.subtitle =
      element_text(
        size = 9.5
      ),
    
    axis.title =
      element_text(
        size = 10
      ),
    
    axis.text =
      element_text(
        size = 9
      ),
    
    plot.tag =
      element_text(
        face = "bold",
        size = 14
      ),
    
    legend.position =
      "none"
  )


# ============================================================
# 26. Panel A — GSE31312 COO
# ============================================================

pA <- ggplot(
  d31312,
  aes(
    x = COO,
    y = TPS12,
    fill = COO
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.35
  ) +
  geom_boxplot(
    width = 0.20,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.10,
    size = 0.8,
    alpha = 0.25
  ) +
  annotate(
    "text",
    x = 2,
    y =
      max(
        d31312$TPS12,
        na.rm = TRUE
      ),
    label =
      "Kruskal-Wallis P < 0.001",
    vjust = -0.5,
    size = 3.3
  ) +
  theme_context +
  labs(
    title =
      "TPS12 across GEP-defined COO",
    subtitle =
      "GSE31312",
    x = NULL,
    y =
      "TPS12 score"
  )


# ============================================================
# 27. Panel B — Independent COO validation
# ============================================================

pB <- ggplot(
  d23501,
  aes(
    x = COO,
    y = TPS12,
    fill = COO
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.35
  ) +
  geom_boxplot(
    width = 0.20,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.10,
    size = 0.8,
    alpha = 0.25
  ) +
  annotate(
    "text",
    x = 2,
    y =
      max(
        d23501$TPS12,
        na.rm = TRUE
      ),
    label = paste0(
      "Kruskal-Wallis P = ",
      sprintf(
        "%.4f",
        kw23501$p.value
      )
    ),
    vjust = -0.5,
    size = 3.3
  ) +
  theme_context +
  labs(
    title =
      "Independent COO validation",
    subtitle =
      "GSE23501",
    x = NULL,
    y =
      "TPS12 score"
  )


# ============================================================
# 28. Panel C — DZ/DHIT-like biology
# ============================================================

pC <- ggplot(
  gcb_dz,
  aes(
    x =
      DZ_DHIT_like_score,
    y =
      TPS12
  )
) +
  geom_point(
    size = 1.8,
    alpha = 0.70
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x =
      min(
        gcb_dz$DZ_DHIT_like_score,
        na.rm = TRUE
      ),
    y =
      max(
        gcb_dz$TPS12,
        na.rm = TRUE
      ),
    hjust = 0,
    vjust = 1,
    label = paste0(
      "Spearman \u03c1 = ",
      sprintf(
        "%.3f",
        rho_dz
      ),
      "\nP = ",
      sprintf(
        "%.3f",
        p_dz
      )
    ),
    size = 3.3
  ) +
  theme_context +
  labs(
    title =
      "TPS12 relationship with DZ/DHIT-like biology",
    subtitle =
      "GSE31312 GCB-DLBCL",
    x =
      "DZ/DHIT-like expression score",
    y =
      "TPS12 score"
  )


# ============================================================
# 29. Panel D — NCI genetic subtypes
# ============================================================

pD <- ggplot(
  genetic_tps,
  aes(
    x =
      Genetic_Subtype,
    y =
      TPS12,
    fill =
      Genetic_Subtype
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.45
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.80
  ) +
  geom_jitter(
    width = 0.12,
    size = 0.65,
    alpha = 0.30
  ) +
  theme_context +
  labs(
    title =
      "TPS12 across genetic subtypes of DLBCL",
    subtitle =
      "NCI DLBCL cohort; Kruskal-Wallis P < 0.0001",
    x =
      "Genetic subtype",
    y =
      "TPS12 score"
  )


# ============================================================
# 30. Assemble final Figure 8
# ============================================================

Figure8_final <- (
  pA | pB
) / (
  pC | pD
) +
  plot_annotation(
    tag_levels = "A"
  )


# ============================================================
# 31. Preview
# ============================================================

print(
  Figure8_final
)


# ============================================================
# 32. Save final Figure 8
# ============================================================

figure8_tiff <- file.path(
  outdir_final,
  "Figure_8_TPS12_molecular_context_FINAL.tiff"
)

figure8_png <- file.path(
  outdir_final,
  "Figure_8_TPS12_molecular_context_FINAL_preview.png"
)

ggsave(
  filename =
    figure8_tiff,
  plot =
    Figure8_final,
  device =
    "tiff",
  width =
    12,
  height =
    9,
  units =
    "in",
  dpi =
    600,
  compression =
    "lzw",
  bg =
    "white"
)

ggsave(
  filename =
    figure8_png,
  plot =
    Figure8_final,
  width =
    12,
  height =
    9,
  units =
    "in",
  dpi =
    300,
  bg =
    "white"
)


# ============================================================
# 33. Export key results
# ============================================================

write_csv(
  molecular_context,
  file.path(
    outdir_context,
    "GSE31312_TPS12_DZ_DHIT_like_scores.csv"
  )
)

write_csv(
  tibble(
    analysis =
      "TPS12 vs DZ/DHIT-like score in GCB-DLBCL",
    
    n =
      nrow(
        gcb_dz
      ),
    
    spearman_rho =
      rho_dz,
    
    p_value =
      p_dz,
    
    n_up_genes =
      length(
        dz_up_final
      ),
    
    n_down_genes =
      length(
        dz_down_final
      )
  ),
  file.path(
    outdir_context,
    "GSE31312_TPS12_DZ_DHIT_like_correlation.csv"
  )
)

write_csv(
  genetic_summary,
  file.path(
    outdir_genetic,
    "TPS12_genetic_subtype_summary.csv"
  )
)

write_csv(
  genetic_tps,
  file.path(
    outdir_genetic,
    "TPS12_NCI_genetic_subtype_merged_data.csv"
  )
)


# ============================================================
# 34. Final console summary
# ============================================================

cat(
  "\n============================================\n",
  "FIGURE 8 COMPLETE\n",
  "============================================\n"
)

cat(
  "\nPanel A:\n",
  "GSE31312 COO Kruskal-Wallis P = ",
  format(
    kw31312$p.value,
    scientific = TRUE,
    digits = 3
  ),
  "\n"
)

cat(
  "\nPanel B:\n",
  "GSE23501 COO Kruskal-Wallis P = ",
  format(
    kw23501$p.value,
    scientific = TRUE,
    digits = 3
  ),
  "\n"
)

cat(
  "\nPanel C:\n",
  "TPS12 vs DZ/DHIT-like score in GCB-DLBCL\n",
  "Spearman rho = ",
  sprintf(
    "%.3f",
    rho_dz
  ),
  "\nP = ",
  sprintf(
    "%.3f",
    p_dz
  ),
  "\n"
)

cat(
  "\nPanel D:\n",
  "NCI genetic subtype Kruskal-Wallis P = ",
  format(
    kw_genetic$p.value,
    scientific = TRUE,
    digits = 3
  ),
  "\n"
)

cat(
  "\nFinal TIFF saved to:\n",
  figure8_tiff,
  "\n"
)

cat(
  "\nFinal PNG preview saved to:\n",
  figure8_png,
  "\n"
)