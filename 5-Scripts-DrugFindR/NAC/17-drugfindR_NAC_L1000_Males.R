# drugfindR: Manual procurement of all the drug similarity scores
# Comparison: NAC L1000 Males

library(tidyverse)
library(drugfindR)

#Input DGE matrix
dge <- read_csv("results/3. L1000 Heat Maps Post SVA/Males/Human_DGE_L1000_Males_Matrix.csv")

#Extract gene signature (gene name and LFC values)
nac_dge_sig <- dge |>
  prepare_signature(gene_column = "Name_GeneSymbol", logfc_column = "NAC",
                    pval_column = "NA_Pval")

#Extract concordant chemical perturbagens (CPs) for the L1000 genes
nac_up_concordants <- get_concordants(nac_dge_sig)

#Extract discordant CPs for the L1000 genes
nac_dn_concordants <- get_concordants(nac_dge_sig)

#Combine all CPs and filter
nac_concordants_all <- nac_up_concordants |>
  mutate(similarity_type = if_else(similarity < 0, "Discordant", "Concordant")) |>
  group_by(similarity_type) |>
  nest()

#Filters CPs by top concordant, top discordant, bottom concordant, bottom discordant
x <- nac_concordants_all %>%
  mutate(sheet_name = str_c(similarity_type)) %>%
  ungroup() %>%
  select(-similarity_type) %>%
  select(sheet_name, data) %>%
  deframe() %>%
  writexl::write_xlsx("results/6. DrugFindR Output/NAC/L1000/Males/NAC_L1000_Males_drugfindr.xlsx")
