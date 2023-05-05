# Extract top bottom 10% of genes from L1000 for each brain region and generate heat map (Mouse Overall Comparison)

library(tidyverse)
library(Orthology.eg.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(pheatmap)
library(ggplot2)

files <- list.files("results/2. Mouse DGEs (Galaxy) No SVA", "*.csv", full.names = TRUE, recursive = TRUE) |>
  keep(~ str_detect(.x, "Overall"))

# NA ---------------------------------------------------------------------

dge_na <- read_csv(files[1], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC, PValue) |>
  dplyr::rename(NAC = logFC) |>
  mutate(GeneID = as.character(GeneID))

mouse_symbol <- select(org.Mm.eg.db, as.character(dge_na$GeneID), c("ENTREZID","SYMBOL"),"ENTREZID")

mapped <- select(Orthology.eg.db, mouse_symbol$ENTREZID, "Homo_sapiens","Mus_musculus") %>%
  filter(!is.na(Homo_sapiens)) %>%
  mutate(Homo_sapiens = as.character(Homo_sapiens))

human_genes <- mapIds(org.Hs.eg.db, mapped$Homo_sapiens, "SYMBOL", "ENTREZID") %>%
  enframe(name = "ENTREZID", value = "SYMBOL")

mapped_mouse <- mapped |>
  inner_join(human_genes, by = c("Homo_sapiens" = "ENTREZID"))

dge_na <- dge_na %>% inner_join(mapped_mouse, by = c("GeneID"="Mus_musculus"))

mapped_ids <- select(org.Mm.eg.db, as.character(dge_na$GeneID), columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID")

dge_na_symboled <- dge_na |>
  dplyr::select(-GeneID, SYMBOL, NAC, PValue) |>
  drugfindR::prepare_signature(gene_column = "SYMBOL", logfc_column = "NAC") |>
  dplyr::select(Name_GeneSymbol, NAC = Value_LogDiffExp, NAC_Pval = Significance_pvalue)

dge_na_symboled_top <- dge_na_symboled |>
  slice_max(NAC, prop = 0.1)

dge_na_symboled_bottom <- dge_na_symboled |>
  slice_min(NAC, prop = 0.1)

dge_na_selected <- bind_rows(dge_na_symboled_top, dge_na_symboled_bottom)

# PFC ---------------------------------------------------------------------

dge_pfc <- read_csv(files[2], col_types = cols(.default = col_character())) |>
  dplyr::select(GeneID, logFC,PValue) |>
  dplyr::rename(PFC = logFC) |>
  mutate(GeneID = as.character(GeneID))

mouse_symbol <- select(org.Mm.eg.db, as.character(dge_na$GeneID), c("ENTREZID","SYMBOL"),"ENTREZID")

mapped <- select(Orthology.eg.db, mouse_symbol$ENTREZID, "Homo_sapiens","Mus_musculus") %>%
  filter(!is.na(Homo_sapiens)) %>%
  mutate(Homo_sapiens = as.character(Homo_sapiens))

human_genes <- mapIds(org.Hs.eg.db, mapped$Homo_sapiens, "SYMBOL", "ENTREZID") %>%
  enframe(name = "ENTREZID", value = "SYMBOL")

mapped_mouse <- mapped |>
  inner_join(human_genes, by = c("Homo_sapiens" = "ENTREZID"))

dge_pfc <- dge_pfc %>% inner_join(mapped_mouse, by = c("GeneID"="Mus_musculus"))

mapped_ids <- select(org.Mm.eg.db, as.character(dge_na$GeneID), columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID")

dge_pfc_symboled <- dge_pfc |>
  dplyr::select(-GeneID, SYMBOL, PFC, PValue) |>
  drugfindR::prepare_signature(gene_column = "SYMBOL", logfc_column = "PFC") |>
  dplyr::select(Name_GeneSymbol, PFC = Value_LogDiffExp, PFC_Pval = Significance_pvalue)

dge_pfc_symboled_top <- dge_pfc_symboled |>
  slice_max(PFC, prop = 0.1)

dge_pfc_symboled_bottom <- dge_pfc_symboled |>
  slice_min(PFC, prop = 0.1)

dge_pfc_selected <- bind_rows(dge_pfc_symboled_top, dge_pfc_symboled_bottom)

# Combine Dataframes ------------------------------------------------------

diffexp <- dge_na_selected |>
  full_join(dge_pfc_selected)

#change all NA's in dataframe to '0'
diffexp <- replace(diffexp, is.na(diffexp), 0) |>

  write_csv("results/5. Mouse Top 10% Heat Maps No SVA/Overall/Mouse_DGE_L1000_10%_Matrix.csv")

# Heatmap Generation ------------------------------------------------------

diff_expr_modified <- read_csv("results/5. Mouse Top 10% Heat Maps No SVA/Overall/Mouse_DGE_L1000_10%_Matrix.csv") |>
  dplyr::select(-ends_with("Pval")) |>
  column_to_rownames("Name_GeneSymbol") |>
  as.matrix()

heatmap(diff_expr_modified)
#need to add scale

# Run correlation analysis -------------------------------------------------
cor(diff_expr_modified, method = "spearman")
