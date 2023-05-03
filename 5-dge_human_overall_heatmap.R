# Extract L1000 from DGEs for Overall Comparisons for each brain region and generate heat map

library(tidyverse)
library(pheatmap)
library(ggplot2)


files <- list.files("results/2. DGEs (Galaxy) Post SVA", "*.csv", full.names = TRUE, recursive = TRUE) |>
  keep(~ str_detect(.x, "Overall"))


# Anterior Insula ---------------------------------------------------------

dge_ai <- read_csv(files[1], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC, PValue) |>
  drugfindR::prepare_signature(gene_column = "GeneID") |>
  dplyr::select(Name_GeneSymbol, AI = Value_LogDiffExp, AI_Pval = Significance_pvalue )


# Cingulate Gyrus ---------------------------------------------------------

dge_cg <- read_csv(files[2], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC, PValue) |>
  drugfindR::prepare_signature(gene_column = "GeneID") |>
  dplyr::select(Name_GeneSymbol, CG = Value_LogDiffExp, CG_Pval = Significance_pvalue )

# DLPFC ------------------------------------------------------------------

dge_dlpfc <- read_csv(files[3], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC, PValue) |>
  drugfindR::prepare_signature(gene_column = "GeneID") |>
  dplyr::select(Name_GeneSymbol, DLPFC = Value_LogDiffExp, DLPFC_Pval = Significance_pvalue )

# NA ---------------------------------------------------------------------

dge_na <- read_csv(files[4], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC, PValue) |>
  drugfindR::prepare_signature(gene_column = "GeneID") |>
  dplyr::select(Name_GeneSymbol, na = Value_LogDiffExp, NA_Pval = Significance_pvalue )
# OFC ---------------------------------------------------------------------

dge_ofc <- read_csv(files[5], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC, PValue) |>
  drugfindR::prepare_signature(gene_column = "GeneID") |>
  dplyr::select(Name_GeneSymbol, OFC = Value_LogDiffExp, OFC_Pval = Significance_pvalue )

# Sub ---------------------------------------------------------------------

dge_sub <- read_csv(files[6], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC, PValue) |>
  drugfindR::prepare_signature(gene_column = "GeneID") |>
  dplyr::select(Name_GeneSymbol, SUB = Value_LogDiffExp, SUB_Pval = Significance_pvalue )


# Combine Dataframes ------------------------------------------------------

diffexp <- dge_ai |>
  inner_join(dge_cg) |>
  inner_join(dge_dlpfc) |>
  inner_join(dge_na) |>
  inner_join(dge_ofc) |>
  inner_join(dge_sub) |>

  write_csv("results/3. L1000 Heat Maps Post SVA/Human_DGE-L1000_Matrix.csv")

# Heatmap Generation ------------------------------------------------------

diff_expr_modified <- read_csv("results/3. L1000 Heat Maps Post SVA/Human_DGE-L1000_Matrix.csv") |>
  select(-ends_with("Pval")) |>
  column_to_rownames("Name_GeneSymbol") |>
  as.matrix()
heatmap(diff_expr_modified)
#need to add scale

# Run correlation analysis -------------------------------------------------
cor(diff_expr_modified, method = "spearman")
