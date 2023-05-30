# drugfindR: Manual procurement of all the drug similarity scores
# Comparison: AI top/bottom 10% from L1000 Males

library(tidyverse)
library(drugfindR)

#Input DGE matrix with genes for our comparison of interest
dge <- read_csv("results/5. Top 10% Heat Maps Post SVA/Males/Human_DGE_L1000_10%_Males_Matrix.csv") |>
  select(Name_GeneSymbol, starts_with("AI")) %>%
  filter(AI != 0)

#Extract gene signature (gene name and LFC values)
ai_dge_sig <- dge |>
  prepare_signature(gene_column = "Name_GeneSymbol", logfc_column = "AI",
                    pval_column = "AI_Pval")

#Create a file listing top X% gene signature
ai_up_sig <- ai_dge_sig |>
  filter_signature("up", threshold = 0) |>
  write_csv("results/6. DrugFindR Output/AI/10%fromL1000/Males/AI_L1000_10%_males_top_genes_sig.csv")

#Create a file listing bottom X% gene signature
ai_dn_sig <- ai_dge_sig |>
  filter_signature("down", threshold = 0) |>
  write_csv ("results/6. DrugFindR Output/AI/10%fromL1000/Males/AI_L1000_10%_males_bottom_genes_sig.csv")

#Extract concordant chemical perturbagens (CPs) for the top/bottom X% of genes
ai_up_concordants <- get_concordants(ai_up_sig, sig_direction = "Up")

#Extract discordant CPs for the top/bottom X% of genes
ai_dn_concordants <- get_concordants(ai_dn_sig, sig_direction = "Down")

#Combine all CPs and filter
ai_concordants_all <- bind_rows(ai_up_concordants, ai_dn_concordants) |>
  #write_csv("results/8. DrugFindR Output/AI/AI_all_drugs.csv") |>
  mutate(similarity_type = if_else(similarity < 0, "Discordant", "Concordant")) |>
  group_by(similarity_type, sig_direction) |>
  nest()

#Filters CPs by top concordant, top discordant, bottom concordant, bottom discordant
x <- ai_concordants_all %>%
  mutate(sheet_name = str_c(sig_direction, similarity_type, collapse = "-")) %>%
  ungroup() %>%
  select(-sig_direction, -similarity_type) %>%
  select(sheet_name, data) %>%
  deframe() %>%
  writexl::write_xlsx("results/6. DrugFindR Output/AI/10%fromL1000/Males/AI_L1000_10%_Males_drugfindr.xlsx")
