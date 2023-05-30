# drugfindR: Manual procurement of all the drug similarity scores
# Comparison: Mouse PFC L1000 Females

library(tidyverse)
library(drugfindR)

#Input DGE matrix
dge <- read_csv("results/3. Mouse L1000 Heat Maps No SVA/Females/Mouse_DGE_L1000_Females_Matrix.csv")

#Extract gene signature (gene name and LFC values)
pfc_dge_sig <- dge |>
  prepare_signature(gene_column = "Name_GeneSymbol", logfc_column = "PFC",
                    pval_column = "PFC_Pval")

#Extract concordant chemical perturbagens (CPs) for the L1000 genes
pfc_up_concordants <- get_concordants(pfc_dge_sig)

#Extract discordant CPs for the L1000 genes
pfc_dn_concordants <- get_concordants(pfc_dge_sig)

#Combine all CPs and filter
pfc_concordants_all <- pfc_up_concordants |>
  mutate(similarity_type = if_else(similarity < 0, "Discordant", "Concordant")) |>
  group_by(similarity_type) |>
  nest()

#Filters CPs by top concordant, top discordant, bottom concordant, bottom discordant
x <- pfc_concordants_all %>%
  mutate(sheet_name = str_c(similarity_type)) %>%
  ungroup() %>%
  select(-similarity_type) %>%
  select(sheet_name, data) %>%
  deframe() %>%
  writexl::write_xlsx("results/6. Mouse DrugFindR Output/Mouse PFC/L1000/Females/Mouse_PFC_L1000_Females_drugfindr.xlsx")
