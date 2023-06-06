# drugfindR: Manual procurement of all the drug similarity scores
# Comparison: DLPFC top/bottom 5% from L1000 Males Meds

library(tidyverse)
library(drugfindR)

#Input DGE matrix with genes for our comparison of interest
dge <- read_csv("results/4. Top 5% Heat Maps Post SVA/Males Meds/Human_DGE_L10005%_Males_Meds_Matrix.csv") |>
  select(Name_GeneSymbol, starts_with("DLPFC")) %>%
  filter(DLPFC != 0)

#Extract gene signature (gene name and LFC values)
dlpfc_dge_sig <- dge |>
  prepare_signature(gene_column = "Name_GeneSymbol", logfc_column = "DLPFC",
                    pval_column = "DLPFC_Pval")

#Create a file listing top X% gene signature
dlpfc_up_sig <- dlpfc_dge_sig |>
  filter_signature("up", threshold = 0) |>
  write_csv("results/6. DrugFindR Output/DLPFC/5%fromL1000/Males Meds/DLPFC_L1000_5%_males_meds_top_genes_sig.csv")

#Create a file listing bottom X% gene signature
dlpfc_dn_sig <- dlpfc_dge_sig |>
  filter_signature("down", threshold = 0) |>
  write_csv ("results/6. DrugFindR Output/DLPFC/5%fromL1000/Males Meds/DLPFC_L1000_5%_males_meds_bottom_genes_sig.csv")

#Extract concordant chemical perturbagens (CPs) for the top/bottom X% of genes
dlpfc_up_concordants <- get_concordants(dlpfc_up_sig, sig_direction = "Up")
#write_csv ("results/8. DrugFindR Output/AI/AI_concordant_drugs_for_top_and_bot.csv")

#Extract discordant CPs for the top/bottom X% of genes
dlpfc_dn_concordants <- get_concordants(dlpfc_dn_sig, sig_direction = "Down")
#write_csv ("results/8. DrugFindR Output/AI/AI_discordant_drugs_for_top_and_bot.csv")

#Combine all CPs and filter
dlpfc_concordants_all <- bind_rows(dlpfc_up_concordants, dlpfc_dn_concordants) |>
  mutate(similarity_type = if_else(similarity < 0, "Discordant", "Concordant")) |>
  group_by(similarity_type, sig_direction) |>
  nest()

#Filters CPs by top concordant, top discordant, bottom concordant, bottom discordant
x <- dlpfc_concordants_all %>%
  mutate(sheet_name = str_c(sig_direction, similarity_type, collapse = "-")) %>%
  ungroup() %>%
  select(-sig_direction, -similarity_type) %>%
  select(sheet_name, data) %>%
  deframe() %>%
  writexl::write_xlsx("results/6. DrugFindR Output/DLPFC/5%fromL1000/Males Meds/DLPFC_L1000_5%_Males_Meds_drugfindr.xlsx")
