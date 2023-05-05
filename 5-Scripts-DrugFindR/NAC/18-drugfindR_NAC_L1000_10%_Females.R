# drugfindR: Manual procurement of all the drug similarity scores
# Comparison: NAC top and bottom 10% from L1000 Females

library(tidyverse)
library(drugfindR)

#Input DGE matrix
dge <- read_csv("results/5. Top 10% Heat Maps Post SVA/Females/Human_DGE_L1000_10%_Females_Matrix.csv")

#Extract gene signature (gene name and LFC values)
nac_dge_sig <- dge |>
  prepare_signature(gene_column = "Name_GeneSymbol", logfc_column = "NAC",
                    pval_column = "NA_Pval")

#Extract top X% gene signature
nac_up_sig <- nac_dge_sig |>
  filter_signature("up", prop = 0.05) |>
  write_csv("results/6. DrugFindR Output/NAC/10%fromL1000/Females/NAC_L1000_10%_females_top_genes_sig.csv")

#Extract bottom X% gene signature
nac_dn_sig <- nac_dge_sig |>
  filter_signature("down", prop = 0.05) |>
  write_csv ("results/6. DrugFindR Output/NAC/10%fromL1000/Females/NAC_L1000_10%_females_bottom_genes_sig.csv")

#Extract concordant chemical perturbagens (CPs) for the top/bottom X% of genes
nac_up_concordants <- get_concordants(nac_up_sig, sig_direction = "Up")
#write_csv ("results/8. DrugFindR Output/AI/AI_concordant_drugs_for_top_and_bot.csv")

#Extract discordant CPs for the top/bottom X% of genes
nac_dn_concordants <- get_concordants(nac_dn_sig, sig_direction = "Down")
#write_csv ("results/8. DrugFindR Output/AI/AI_discordant_drugs_for_top_and_bot.csv")

#Combine all CPs and filter
nac_concordants_all <- bind_rows(nac_up_concordants, nac_dn_concordants) |>
  #write_csv("results/8. DrugFindR Output/AI/AI_all_drugs.csv") |>
  mutate(similarity_type = if_else(similarity < 0, "Discordant", "Concordant")) |>
  group_by(similarity_type, sig_direction) |>
  nest()

#Filters CPs by top concordant, top discordant, bottom concordant, bottom discordant
x <- nac_concordants_all %>%
  mutate(sheet_name = str_c(sig_direction, similarity_type, collapse = "-")) %>%
  ungroup() %>%
  select(-sig_direction, -similarity_type) %>%
  select(sheet_name, data) %>%
  deframe() %>%
  writexl::write_xlsx("results/6. DrugFindR Output/NAC/10%fromL1000/Females/NAC_L1000_10%_Females_drugfindr.xlsx")
