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


# Combine Dataframes ------------------------------------------------------

 diffexp <- dge_ai_symboled |>
  inner_join(dge_cg_symboled) |>
  inner_join(dge_dlpfc_symboled) |>
  inner_join(dge_na_symboled) |>
  inner_join(dge_ofc_symboled) |>
  inner_join(dge_sub_symboled) |>


write_csv("results/Human_DGE-L1000_Matrix.csv")


# Heatmap Generation ------------------------------------------------------

diff_expr_modified <- read_csv("results/Human_DGE-L1000_Matrix.csv") |>
column_to_rownames("Name_GeneSymbol") |>
as.matrix()
heatmap(diff_expr_modified)
#need to add scale

