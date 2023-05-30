# Extract top bottom 10% of genes from L1000 for each brain region and generate heat map (Overall Comparison)

library(tidyverse)
library(pheatmap)
library(ggplot2)


files <- list.files("results/2. DGEs (Galaxy) Post SVA", "*.csv", full.names = TRUE, recursive = TRUE) |>
  keep(~ str_detect(.x, "Overall"))


# Anterior Insula ---------------------------------------------------------

#Read anterior insula CSV file and obtain the genes, LFC values, and p values
dge_ai <- read_csv(files[1], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC, PValue) |>
  drugfindR::prepare_signature(gene_column = "GeneID") |>
  dplyr::select(Name_GeneSymbol, AI = Value_LogDiffExp, AI_Pval = Significance_pvalue )

#Obtain the top 10% DEGs from our dataframe
dge_ai_top <- dge_ai |>
  slice_max(AI, prop = 0.1)

#Obtain the bottom 10% DEGs from our dataframe
dge_ai_bottom <- dge_ai |>
  slice_min(AI, prop = 0.1)

#Combine the top and bottom DEGs into a single dataframe
dge_ai_selected <- bind_rows(dge_ai_top, dge_ai_bottom)

# Cingulate Gyrus ---------------------------------------------------------

dge_cg <- read_csv(files[2], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC, PValue) |>
  drugfindR::prepare_signature(gene_column = "GeneID") |>
  dplyr::select(Name_GeneSymbol, CG = Value_LogDiffExp, CG_Pval = Significance_pvalue )

dge_cg_top <- dge_cg |>
  slice_max(CG, prop = 0.1)

dge_cg_bottom <- dge_cg |>
  slice_min(CG, prop = 0.1)

dge_cg_selected <- bind_rows(dge_cg_top, dge_cg_bottom)

# DLPFC ------------------------------------------------------------------

dge_dlpfc <- read_csv(files[3], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC, PValue) |>
  drugfindR::prepare_signature(gene_column = "GeneID") |>
  dplyr::select(Name_GeneSymbol, DLPFC = Value_LogDiffExp, DLPFC_Pval = Significance_pvalue )

dge_dlpfc_top <- dge_dlpfc |>
  slice_max(DLPFC, prop = 0.1)

dge_dlpfc_bottom <- dge_dlpfc |>
  slice_min(DLPFC, prop = 0.1)

dge_dlpfc_selected <- bind_rows(dge_dlpfc_top, dge_dlpfc_bottom)

# NA ---------------------------------------------------------------------

dge_na <- read_csv(files[4], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC, PValue) |>
  drugfindR::prepare_signature(gene_column = "GeneID") |>
  dplyr::select(Name_GeneSymbol, NAC = Value_LogDiffExp, NA_Pval = Significance_pvalue )

dge_na_top <- dge_na |>
  slice_max(NAC, prop = 0.1)

dge_na_bottom <- dge_na |>
  slice_min(NAC, prop = 0.1)

dge_na_selected <- bind_rows(dge_na_top, dge_na_bottom)

# OFC ---------------------------------------------------------------------

dge_ofc <- read_csv(files[5], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC, PValue) |>
  drugfindR::prepare_signature(gene_column = "GeneID") |>
  dplyr::select(Name_GeneSymbol, OFC = Value_LogDiffExp, OFC_Pval = Significance_pvalue )

dge_ofc_top <- dge_ofc |>
  slice_max(OFC, prop = 0.1)

dge_ofc_bottom <- dge_ofc |>
  slice_min(OFC, prop = 0.1)

dge_ofc_selected <- bind_rows(dge_ofc_top, dge_ofc_bottom)

# Sub ---------------------------------------------------------------------

dge_sub <- read_csv(files[6], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC, PValue) |>
  drugfindR::prepare_signature(gene_column = "GeneID") |>
  dplyr::select(Name_GeneSymbol, SUB = Value_LogDiffExp, SUB_Pval = Significance_pvalue)

dge_sub_top <- dge_sub |>
  slice_max(SUB, prop = 0.1)

dge_sub_bottom <- dge_sub |>
  slice_min(SUB, prop = 0.1)

dge_sub_selected <- bind_rows(dge_sub_top, dge_sub_bottom)

# Combine Dataframes ------------------------------------------------------

diffexp <- dge_ai_selected |>
  full_join(dge_cg_selected) |>
  full_join(dge_dlpfc_selected) |>
  full_join(dge_na_selected) |>
  full_join(dge_ofc_selected) |>
  full_join(dge_sub_selected)

#change all NA's in dataframe to '0'
diffexp <- replace(diffexp, is.na(diffexp), 0) |>
  write_csv("results/5. Top 10% Heat Maps Post SVA/Overall/Human_DGE_L1000_10%_Matrix.csv")

# Heatmap Generation ------------------------------------------------------

diff_expr_modified <- read_csv("results/5. Top 10% Heat Maps Post SVA/Overall/Human_DGE_L1000_10%_Matrix.csv") |>
  select(-ends_with("Pval")) |>
  column_to_rownames("Name_GeneSymbol") |>
  as.matrix()
heatmap(diff_expr_modified)
#need to add scale

# Run correlation analysis -------------------------------------------------
cor(diff_expr_modified, method = "spearman")
