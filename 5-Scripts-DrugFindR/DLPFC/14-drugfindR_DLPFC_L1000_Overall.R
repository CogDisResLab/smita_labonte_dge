# drugfindR: Manual procurement of all the drug similarity scores
# Comparison: DLPFC L1000 Overall

library(tidyverse)
library(drugfindR)

#Input DGE matrix
dge <- read_csv("results/3. L1000 Heat Maps Post SVA/Overall/Human_DGE_L1000_Overall_Matrix.csv")

#Extract gene signature (gene name and LFC values)
dlpfc_dge_sig <- dge |>
  prepare_signature(gene_column = "Name_GeneSymbol", logfc_column = "CG",
                    pval_column = "CG_Pval")

#Extract top X% gene signature
dlpfc_up_sig <- dlpfc_dge_sig |>
  filter_signature("up", prop = 0.05) |>
  write_csv("results/6. DrugFindR Output/DLPFC/L1000/Overall/DLPFC_overall_top_genes_sig.csv")

#Extract bottom X% gene signature
dlpfc_dn_sig <- dlpfc_dge_sig |>
  filter_signature("down", prop = 0.05) |>
  write_csv ("results/6. DrugFindR Output/DLPFC/L1000/Overall/DLPFC_overall_bottom_genes_sig.csv")

#Extract concordant chemical perturbagens (CPs) for the top/bottom X% of genes
dlpfc_up_concordants <- get_concordants(dlpfc_up_sig, sig_direction = "Up")
#write_csv ("results/8. DrugFindR Output/AI/AI_concordant_drugs_for_top_and_bot.csv")

#Extract discordant CPs for the top/bottom X% of genes
dlpfc_dn_concordants <- get_concordants(dlpfc_dn_sig, sig_direction = "Down")
#write_csv ("results/8. DrugFindR Output/AI/AI_discordant_drugs_for_top_and_bot.csv")

#Combine all CPs and filter
dlpfc_concordants_all <- bind_rows(dlpfc_up_concordants, dlpfc_dn_concordants) |>
  #write_csv("results/8. DrugFindR Output/AI/AI_all_drugs.csv") |>
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
  writexl::write_xlsx("results/6. DrugFindR Output/DLPFC/L1000/Overall/DLPFC_L1000_Overall_drugfindr.xlsx")
