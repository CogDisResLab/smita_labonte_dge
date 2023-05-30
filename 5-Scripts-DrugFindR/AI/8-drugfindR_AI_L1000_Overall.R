# drugfindR: Manual procurement of all the drug similarity scores
# Comparison: AI L1000 Overall

library(tidyverse)
library(drugfindR)

#Input DGE matrix
dge <- read_csv("results/3. L1000 Heat Maps Post SVA/Overall/Human_DGE_L1000_Overall_Matrix.csv")

#Extract gene signature (gene name and LFC values)
ai_dge_sig <- dge |>
  prepare_signature(gene_column = "Name_GeneSymbol", logfc_column = "AI",
                    pval_column = "AI_Pval")

#Extract concordant chemical perturbagens (CPs) for the L1000 genes
ai_up_concordants <- get_concordants(ai_dge_sig)
#write_csv ("results/6. DrugFindR Output/AI/L1000/Overall/AI_concordant_drugs_for_L1000.csv")

#Extract discordant CPs for the L1000 genes
ai_dn_concordants <- get_concordants(ai_dge_sig)
#write_csv ("results/6. DrugFindR Output/AI/L1000/Overall/AI_discordant_drugs_for_L1000.csv")

#Combine all CPs and filter
ai_concordants_all <- ai_up_concordants |>
  #write_csv("results/6. DrugFindR Output/AI/L1000/Overall/AI_all_drugs.csv") |>
  mutate(similarity_type = if_else(similarity < 0, "Discordant", "Concordant")) |>
  group_by(similarity_type) |>
  nest()

#Filters CPs by top concordant, top discordant, bottom concordant, bottom discordant
x <- ai_concordants_all %>%
  mutate(sheet_name = str_c(similarity_type)) %>%
  ungroup() %>%
  select(-similarity_type) %>%
  select(sheet_name, data) %>%
  deframe() %>%
  writexl::write_xlsx("results/6. DrugFindR Output/AI/L1000/Overall/AI_L1000_Overall_drugfindr.xlsx")
