# drugfindR: Manual procurement of all the drug similarity scores
# Comparison: CG top and bottom 5% from L1000 Males Meds

library(tidyverse)
library(drugfindR)

#Input DGE matrix
dge <- read_csv("results/4. Top 5% Heat Maps Post SVA/Males Meds/Human_DGE_L10005%_Males_Meds_Matrix.csv")

#Extract gene signature (gene name and LFC values)
cg_dge_sig <- dge |>
  prepare_signature(gene_column = "Name_GeneSymbol", logfc_column = "CG",
                    pval_column = "CG_Pval")

#Extract top X% gene signature
cg_up_sig <- cg_dge_sig |>
  filter_signature("up", prop = 0.05) |>
  write_csv("results/6. DrugFindR Output/CG/5%fromL1000/Males Meds/CG_L1000_5%_males_meds_top_genes_sig.csv")

#Extract bottom X% gene signature
cg_dn_sig <- cg_dge_sig |>
  filter_signature("down", prop = 0.05) |>
  write_csv ("results/6. DrugFindR Output/CG/5%fromL1000/Males Meds/CG_L1000_5%_males_meds_bottom_genes_sig.csv")

#Extract concordant chemical perturbagens (CPs) for the top/bottom X% of genes
cg_up_concordants <- get_concordants(cg_up_sig, sig_direction = "Up")
#write_csv ("results/8. DrugFindR Output/AI/AI_concordant_drugs_for_top_and_bot.csv")

#Extract discordant CPs for the top/bottom X% of genes
cg_dn_concordants <- get_concordants(cg_dn_sig, sig_direction = "Down")
#write_csv ("results/8. DrugFindR Output/AI/AI_discordant_drugs_for_top_and_bot.csv")

#Combine all CPs and filter
cg_concordants_all <- bind_rows(cg_up_concordants, cg_dn_concordants) |>
  #write_csv("results/8. DrugFindR Output/AI/AI_all_drugs.csv") |>
  mutate(similarity_type = if_else(similarity < 0, "Discordant", "Concordant")) |>
  group_by(similarity_type, sig_direction) |>
  nest()

#Filters CPs by top concordant, top discordant, bottom concordant, bottom discordant
x <- cg_concordants_all %>%
  mutate(sheet_name = str_c(sig_direction, similarity_type, collapse = "-")) %>%
  ungroup() %>%
  select(-sig_direction, -similarity_type) %>%
  select(sheet_name, data) %>%
  deframe() %>%
  writexl::write_xlsx("results/6. DrugFindR Output/CG/5%fromL1000/Males Meds/CG_L1000_5%_Males_Meds_drugfindr.xlsx")
