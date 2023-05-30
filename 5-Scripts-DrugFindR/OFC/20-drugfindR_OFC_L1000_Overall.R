# drugfindR: Manual procurement of all the drug similarity scores
# Comparison: OFC L1000 Overall

library(tidyverse)
library(drugfindR)

#Input DGE matrix
dge <- read_csv("results/3. L1000 Heat Maps Post SVA/Overall/Human_DGE_L1000_Overall_Matrix.csv")

#Extract gene signature (gene name and LFC values)
ofc_dge_sig <- dge |>
  prepare_signature(gene_column = "Name_GeneSymbol", logfc_column = "OFC",
                    pval_column = "OFC_Pval")

#Extract concordant chemical perturbagens (CPs) for the L1000 genes
ofc_up_concordants <- get_concordants(ofc_dge_sig)

#Extract discordant CPs for the L1000 genes
ofc_dn_concordants <- get_concordants(ofc_dge_sig)

#Combine all CPs and filter
ofc_concordants_all <- ofc_up_concordants |>
  mutate(similarity_type = if_else(similarity < 0, "Discordant", "Concordant")) |>
  group_by(similarity_type) |>
  nest()

#Filters CPs by top concordant, top discordant, bottom concordant, bottom discordant
x <- ofc_concordants_all %>%
  mutate(sheet_name = str_c(similarity_type)) %>%
  ungroup() %>%
  select(-similarity_type) %>%
  select(sheet_name, data) %>%
  deframe() %>%
  writexl::write_xlsx("results/6. DrugFindR Output/OFC/L1000/Overall/OFC_L1000_Overall_drugfindr.xlsx")
