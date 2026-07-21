# ============================================================
# NCI DLBCL — TPS12 across published genetic subtypes
#
# Genetic subtypes:
# MCD, BN2, N1, EZB, Other
#
# Source:
# Schmitz et al., NEJM 2018
# ============================================================

library(tidyverse)
library(readxl)
library(patchwork)

set.seed(123)

# ============================================================
# 1. Paths
# ============================================================

expr_file <- "data/NCI_DLBCL/RNAseq_gene_expression_562.txt"

subtype_file <- "data/NCI_DLBCL/Supplementary_Appendix_2.xlsx"

outdir <- "results/TPS12_genetic_subgroups_NCI"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(file.exists(expr_file))
stopifnot(file.exists(subtype_file))

# ============================================================
# 2. TPS12 genes
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
# 3. Read expression matrix
# ============================================================

expr_raw <- read.delim(
  expr_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\nExpression dimensions:\n")
print(dim(expr_raw))

cat("\nFirst columns:\n")
print(head(names(expr_raw), 10))

cat("\nFirst rows:\n")
print(head(expr_raw[, 1:min(6, ncol(expr_raw))]))

# ============================================================
# 4. Detect gene-symbol column
# ============================================================

gene_col_candidates <- names(expr_raw)[
  grepl(
    "gene|symbol",
    names(expr_raw),
    ignore.case = TRUE
  )
]

cat("\nPossible gene columns:\n")
print(gene_col_candidates)

if (length(gene_col_candidates) == 0) {
  stop(
    "Could not identify gene-symbol column. ",
    "Inspect names(expr_raw)."
  )
}

gene_col <- gene_col_candidates[1]

# ============================================================
# 5. Build gene-level expression matrix
# ============================================================

expr_clean <- expr_raw %>%
  rename(
    gene_symbol = all_of(gene_col)
  ) %>%
  filter(
    !is.na(gene_symbol),
    gene_symbol != ""
  )

expr_long <- expr_clean %>%
  pivot_longer(
    cols = -c(gene_symbol, Accession, Gene_ID),
    names_to = "sample_id",
    values_to = "expr_value"
  ) %>%
  mutate(
    expr_value = as.numeric(expr_value)
  ) %>%
  filter(
    is.finite(expr_value)
  )

# ============================================================
# 6. Check TPS12 gene availability
# ============================================================

available_genes <- unique(
  expr_gene$gene_symbol
)

tps_present <- intersect(
  tps12_genes,
  available_genes
)

tps_missing <- setdiff(
  tps12_genes,
  available_genes
)

cat(
  "\nTPS12 genes present:",
  length(tps_present),
  "of 12\n"
)

print(tps_present)

cat("\nTPS12 genes missing:\n")
print(tps_missing)

if (length(tps_present) < 10) {
  stop(
    "Too few TPS12 genes available."
  )
}

tps_expr <- expr_long %>%
  filter(
    gene_symbol %in% tps_present
  ) %>%
  group_by(
    gene_symbol
  ) %>%
  mutate(
    z = as.numeric(
      scale(expr_value)
    )
  ) %>%
  ungroup()

tps_scores <- tps_expr %>%
  group_by(
    sample_id
  ) %>%
  summarise(
    TPS12 = mean(
      z,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# ============================================================
# 7. Extract official NCI genetic subtype assignments
# ============================================================

library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)

clinical_subtypes <- read_excel(
  subtype_file,
  sheet = "Tab S9 Characteristics DLBCL"
)

cat("\n===== TAB S9 DIMENSIONS =====\n")
print(dim(clinical_subtypes))

cat("\n===== TAB S9 COLUMNS =====\n")
print(names(clinical_subtypes))

cat("\n===== FIRST 10 ROWS =====\n")
print(
  clinical_subtypes %>%
    select(
      `dbGaP submitted subject ID`,
      `dbGaP accession`,
      `Gene Expression Subgroup`,
      `Genetic Subtype`
    ) %>%
    head(10)
)

cat("\n===== GENETIC SUBTYPE COUNTS =====\n")
print(
  clinical_subtypes %>%
    count(`Genetic Subtype`, sort = TRUE)
)

cat("\n===== GENE EXPRESSION SUBGROUP COUNTS =====\n")
print(
  clinical_subtypes %>%
    count(`Gene Expression Subgroup`, sort = TRUE)
)

# ============================================================
# 8. Check which identifier matches RNA-seq sample IDs
# ============================================================

expr_ids <- unique(tps_scores$sample_id)

id_check <- tibble(
  identifier = c(
    "dbGaP submitted subject ID",
    "dbGaP accession"
  ),
  n_matches = c(
    sum(
      clinical_subtypes$`dbGaP submitted subject ID`
      %in% expr_ids,
      na.rm = TRUE
    ),
    sum(
      clinical_subtypes$`dbGaP accession`
      %in% expr_ids,
      na.rm = TRUE
    )
  )
)

cat("\n===== SAMPLE ID MATCHING =====\n")
print(id_check)

cat("\n===== FIRST TPS12 SAMPLE IDs =====\n")
print(head(expr_ids, 20))

cat("\n===== FIRST TAB S9 SUBJECT IDs =====\n")
print(
  head(
    clinical_subtypes$`dbGaP submitted subject ID`,
    20
  )
)

cat("\n===== FIRST TAB S9 ACCESSION IDs =====\n")
print(
  head(
    clinical_subtypes$`dbGaP accession`,
    20
  )
)


# ============================================================
# 9. Merge TPS12 with official genetic subtype assignments
# ============================================================

genetic_meta <- clinical_subtypes %>%
  transmute(
    sample_id = `dbGaP submitted subject ID`,
    COO = `Gene Expression Subgroup`,
    Genetic_Subtype = `Genetic Subtype`
  ) %>%
  filter(
    !is.na(sample_id),
    !is.na(Genetic_Subtype)
  ) %>%
  distinct(sample_id, .keep_all = TRUE)

genetic_tps <- tps_scores %>%
  inner_join(
    genetic_meta,
    by = "sample_id"
  )

cat("\n===== MERGED TPS12 + GENETIC SUBTYPE =====\n")
print(dim(genetic_tps))

cat("\n===== GENETIC SUBTYPE COUNTS AFTER MERGE =====\n")
print(
  genetic_tps %>%
    count(Genetic_Subtype, sort = TRUE)
)

cat("\n===== COO × GENETIC SUBTYPE =====\n")
print(
  table(
    genetic_tps$COO,
    genetic_tps$Genetic_Subtype,
    useNA = "ifany"
  )
)

# ============================================================
# 10. Set genetic subtype order
# ============================================================

genetic_tps <- genetic_tps %>%
  filter(
    Genetic_Subtype %in% c(
      "MCD", "BN2", "N1", "EZB", "Other"
    )
  ) %>%
  mutate(
    Genetic_Subtype = factor(
      Genetic_Subtype,
      levels = c(
        "MCD", "BN2", "N1", "EZB", "Other"
      )
    )
  )

# ============================================================
# 11. Overall comparison
# ============================================================

kw_genetic <- kruskal.test(
  TPS12 ~ Genetic_Subtype,
  data = genetic_tps
)

cat("\n===== KRUSKAL-WALLIS: TPS12 BY GENETIC SUBTYPE =====\n")
print(kw_genetic)

pairwise_genetic <- pairwise.wilcox.test(
  genetic_tps$TPS12,
  genetic_tps$Genetic_Subtype,
  p.adjust.method = "BH",
  exact = FALSE
)

cat("\n===== PAIRWISE WILCOXON, BH-ADJUSTED =====\n")
print(pairwise_genetic)

# ============================================================
# 12. Summary statistics
# ============================================================

genetic_summary <- genetic_tps %>%
  group_by(Genetic_Subtype) %>%
  summarise(
    n = n(),
    median_TPS12 = median(TPS12, na.rm = TRUE),
    IQR_TPS12 = IQR(TPS12, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n===== GENETIC SUBTYPE SUMMARY =====\n")
print(genetic_summary)

write.csv(
  genetic_summary,
  file.path(
    outdir,
    "TPS12_genetic_subtype_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 13. Figure
# ============================================================

p_genetic <- ggplot(
  genetic_tps,
  aes(
    x = Genetic_Subtype,
    y = TPS12,
    fill = Genetic_Subtype
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.45
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    width = 0.12,
    size = 0.8,
    alpha = 0.35
  ) +
  labs(
    title = "TPS12 across genetic subtypes of DLBCL",
    subtitle = "NCI DLBCL cohort; Kruskal-Wallis P < 0.0001",
    x = "Genetic subtype",
    y = "TPS12 score"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

print(p_genetic)

# ============================================================
# 14. Save figure
# ============================================================

ggsave(
  file.path(
    outdir,
    "Figure_TPS12_NCI_genetic_subtypes.tiff"
  ),
  p_genetic,
  width = 7,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(
    outdir,
    "Figure_TPS12_NCI_genetic_subtypes_preview.png"
  ),
  p_genetic,
  width = 7,
  height = 5,
  units = "in",
  dpi = 200
)

# Save full merged data
write.csv(
  genetic_tps,
  file.path(
    outdir,
    "TPS12_NCI_genetic_subtype_merged_data.csv"
  ),
  row.names = FALSE
)

cat(
  "\nAnalysis complete.\n",
  "Output directory:\n",
  outdir,
  "\n"
)