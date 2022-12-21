# DGE Heatmap Plot


library(tidyverse)
library(org.Hs.eg.db)

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

mapped_ids <- select(org.Hs.eg.db, as.character(dge_ai$GeneID), columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID")

dge_cg_symboled <- dge_cg |>
  inner_join(mapped_ids, by = c(GeneID = "ENTREZID")) |>
  dplyr::select(-GeneID, SYMBOL, CG) |>
  mutate(PValue = rep(0, nrow(mapped_ids))) |>
  drugfindR::prepare_signature(gene_column = "SYMBOL", logfc_column = "CG") |>
  dplyr::select(Name_GeneSymbol, CG = Value_LogDiffExp)



# Combine Dataframes ------------------------------------------------------

diffexp <- dge_ai_symboled |>
  inner_join(dge_cg_symboled) |>
  inner_join(blah)
