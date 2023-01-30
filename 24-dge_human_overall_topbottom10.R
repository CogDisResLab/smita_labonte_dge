# DGE Heatmap Plot


library(tidyverse)
library(org.Hs.eg.db)
library(pheatmap)
library(ggplot2)


files <- list.files("Galaxy Output/", "*.csv", full.names = TRUE, recursive = TRUE) |>
  keep(~ str_detect(.x, "Overall"))

# Anterior Insula ---------------------------------------------------------

dge_ai <- read_csv(files[1], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC) |>
  dplyr::rename(AI = logFC)

mapped_ids <- select(org.Hs.eg.db, as.character(dge_ai$GeneID), columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID")

dge_ai_symboled <- dge_ai |>
  inner_join(mapped_ids, by = c(GeneID = "ENTREZID")) |>
  dplyr::select(-GeneID, SYMBOL, AI) |>
  mutate(PValue = rep(0, nrow(mapped_ids))) |>
  drugfindR::prepare_signature(gene_column = "SYMBOL", logfc_column = "AI") |>
  dplyr::select(Name_GeneSymbol, AI = Value_LogDiffExp)

dge_ai_symboled_top <- dge_ai_symboled |>
  slice_max(AI, prop = 0.1)

dge_ai_symboled_bottom <- dge_ai_symboled |>
  slice_min(AI, prop = 0.1)

dge_ai_selected <- bind_rows(dge_ai_symboled_top, dge_ai_symboled_bottom)

# Cingulate Gyrus ---------------------------------------------------------

dge_cg <- read_csv(files[2], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC) |>
  dplyr::rename(CG = logFC)

mapped_ids <- select(org.Hs.eg.db, as.character(dge_cg$GeneID), columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID")

dge_cg_symboled <- dge_cg |>
  inner_join(mapped_ids, by = c(GeneID = "ENTREZID")) |>
  dplyr::select(-GeneID, SYMBOL, CG) |>
  mutate(PValue = rep(0, nrow(mapped_ids))) |>
  drugfindR::prepare_signature(gene_column = "SYMBOL", logfc_column = "CG") |>
  dplyr::select(Name_GeneSymbol, CG = Value_LogDiffExp)

dge_cg_symboled_top <- dge_cg_symboled |>
  slice_max(CG, prop = 0.1)

dge_cg_symboled_bottom <- dge_cg_symboled |>
  slice_min(CG, prop = 0.1)

dge_cg_selected <- bind_rows(dge_cg_symboled_top, dge_cg_symboled_bottom)

# DLPFC ------------------------------------------------------------------

dge_dlpfc <- read_csv(files[3], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC) |>
  dplyr::rename(DLPFC = logFC)

mapped_ids <- select(org.Hs.eg.db, as.character(dge_dlpfc$GeneID), columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID")

dge_dlpfc_symboled <- dge_dlpfc |>
  inner_join(mapped_ids, by = c(GeneID = "ENTREZID")) |>
  dplyr::select(-GeneID, SYMBOL, DLPFC) |>
  mutate(PValue = rep(0, nrow(mapped_ids))) |>
  drugfindR::prepare_signature(gene_column = "SYMBOL", logfc_column = "DLPFC") |>
  dplyr::select(Name_GeneSymbol, DLPFC = Value_LogDiffExp)

dge_dlpfc_symboled_top <- dge_dlpfc_symboled |>
  slice_max(DLPFC, prop = 0.1)

dge_dlpfc_symboled_bottom <- dge_dlpfc_symboled |>
  slice_min(DLPFC, prop = 0.1)

dge_dlpfc_selected <- bind_rows(dge_dlpfc_symboled_top, dge_dlpfc_symboled_bottom)

# NA ---------------------------------------------------------------------

dge_na <- read_csv(files[6], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC) |>
  dplyr::rename(NAA = logFC)

mapped_ids <- select(org.Hs.eg.db, as.character(dge_na$GeneID), columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID")

dge_na_symboled <- dge_na |>
  inner_join(mapped_ids, by = c(GeneID = "ENTREZID")) |>
  dplyr::select(-GeneID, SYMBOL, NAA) |>
  mutate(PValue = rep(0, nrow(mapped_ids))) |>
  drugfindR::prepare_signature(gene_column = "SYMBOL", logfc_column = "NAA") |>
  dplyr::select(Name_GeneSymbol, NAA = Value_LogDiffExp)

dge_na_symboled_top <- dge_na_symboled |>
  slice_max(NAA, prop = 0.1)

dge_na_symboled_bottom <- dge_na_symboled |>
  slice_min(NAA, prop = 0.1)

dge_na_selected <- bind_rows(dge_na_symboled_top, dge_na_symboled_bottom)

# OFC ---------------------------------------------------------------------

dge_ofc <- read_csv(files[7], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC) |>
  dplyr::rename(OFC = logFC)

mapped_ids <- select(org.Hs.eg.db, as.character(dge_ofc$GeneID), columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID")

dge_ofc_symboled <- dge_ofc |>
  inner_join(mapped_ids, by = c(GeneID = "ENTREZID")) |>
  dplyr::select(-GeneID, SYMBOL, OFC) |>
  mutate(PValue = rep(0, nrow(mapped_ids))) |>
  drugfindR::prepare_signature(gene_column = "SYMBOL", logfc_column = "OFC") |>
  dplyr::select(Name_GeneSymbol, OFC = Value_LogDiffExp)

dge_ofc_symboled_top <- dge_ofc_symboled |>
  slice_max(OFC, prop = 0.1)

dge_ofc_symboled_bottom <- dge_ofc_symboled |>
  slice_min(OFC, prop = 0.1)

dge_ofc_selected <- bind_rows(dge_ofc_symboled_top, dge_ofc_symboled_bottom)


# Sub ---------------------------------------------------------------------

dge_sub <- read_csv(files[8], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC) |>
  dplyr::rename(SUB = logFC)

mapped_ids <- select(org.Hs.eg.db, as.character(dge_sub$GeneID), columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID")

dge_sub_symboled <- dge_sub |>
  inner_join(mapped_ids, by = c(GeneID = "ENTREZID")) |>
  dplyr::select(-GeneID, SYMBOL, SUB) |>
  mutate(PValue = rep(0, nrow(mapped_ids))) |>
  drugfindR::prepare_signature(gene_column = "SYMBOL", logfc_column = "SUB") |>
  dplyr::select(Name_GeneSymbol, SUB = Value_LogDiffExp)

dge_sub_symboled_top <- dge_sub_symboled |>
  slice_max(SUB, prop = 0.1)

dge_sub_symboled_bottom <- dge_sub_symboled |>
  slice_min(SUB, prop = 0.1)

dge_sub_selected <- bind_rows(dge_sub_symboled_top, dge_sub_symboled_bottom)

# Combine Dataframes ------------------------------------------------------

diffexp <- dge_ai_selected |>
  full_join(dge_cg_selected) |>
  full_join(dge_dlpfc_selected) |>
  full_join(dge_na_selected) |>
  full_join(dge_ofc_selected) |>
  full_join(dge_sub_selected)

#change all NA's in dataframe to '0'
diffexp <- replace(diffexp, is.na(diffexp), 0) |>

  write_csv("results/Human_DGE-L1000_Matrix_Overall_topbottom10.csv")


# Heatmap Generation ------------------------------------------------------

diff_expr_modified <- read_csv("results/Human_DGE-L1000_Matrix_Overall_topbottom10.csv") |>
  column_to_rownames("Name_GeneSymbol") |>
  as.matrix()

heatmap(diff_expr_modified)
#need to add scale
