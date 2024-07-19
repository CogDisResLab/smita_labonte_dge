# drugfindR: Manual procurement of all the drug similarity scores
# Comparison: CG top/bottom 50% from L1000 Females

library(tidyverse)
library(drugfindR)

#Input DGE matrix with genes for our comparison of interest
dge <- read_csv("results/5. Top 50% Heat Maps Post SVA/Females/Human_DGE_L1000_50%_Females_Matrix.csv") |>
  select(Name_GeneSymbol, starts_with("CG")) %>%
  filter(CG != 0)

#Extract gene signature (gene name and LFC values)
cg_dge_sig <- dge |>
  prepareSignature(geneColumn = "Name_GeneSymbol", logfcColumn = "CG",
                    pvalColumn = "CG_Pval")

#Create a file listing top X% gene signature
cg_up_sig <- cg_dge_sig |>
  filterSignature("up", threshold = 0) |>
  write_csv("results/6. DrugFindR Output/CG/50%fromL1000/Females/CG_L1000_50%_females_top_genes_sig.csv")

#Create a file listing bottom X% gene signature
cg_dn_sig <- cg_dge_sig |>
  filterSignature("down", threshold = 0) |>
  write_csv ("results/6. DrugFindR Output/CG/50%fromL1000/Females/CG_L1000_50%_females_bottom_genes_sig.csv")

#Extract concordant chemical perturbagens (CPs) for the top/bottom X% of genes
cg_up_concordants <- getConcordants(cg_up_sig)
#write_csv ("results/8. DrugFindR Output/AI/AI_concordant_drugs_for_top_and_bot.csv")

#Extract discordant CPs for the top/bottom X% of genes
cg_dn_concordants <- getConcordants(cg_dn_sig)
#write_csv ("results/8. DrugFindR Output/AI/AI_discordant_drugs_for_top_and_bot.csv")

#Combine all CPs and filter
cg_concordants_all <- bind_rows(cg_up_concordants, cg_dn_concordants) |>
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
  writexl::write_xlsx("results/6. DrugFindR Output/CG/50%fromL1000/Females/CG_L1000_50%_Females_drugfindr.xlsx")
