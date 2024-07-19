# drugfindR: Manual procurement of all the drug similarity scores
# Comparison: CG top/bottom 25% from L1000 Overall

library(tidyverse)
library(drugfindR)

#Input DGE matrix with genes for our comparison of interest
dge <- read_csv("results/4. Top 25% Heat Maps Post SVA/Overall/Human_DGE_L100025%_Matrix.csv") |>
  select(Name_GeneSymbol, starts_with("CG")) %>%
  filter(CG != 0)

#Extract gene signature (gene name and LFC values)
cg_dge_sig <- dge |>
  prepareSignature(geneColumn = "Name_GeneSymbol", logfcColumn = "CG",
                    pvalColumn = "CG_Pval")

#Create a file listing top X% gene signature
cg_up_sig <- cg_dge_sig |>
  filterSignature("up", threshold = 0) |>
  write_csv("results/6. DrugFindR Output/CG/25%fromL1000/Overall/CG_L1000_25%_overall_top_genes_sig.csv")

#Create a file listing bottom X% gene signature
cg_dn_sig <- cg_dge_sig |>
  filterSignature("down", threshold = 0) |>
  write_csv ("results/6. DrugFindR Output/CG/25%fromL1000/Overall/CG_L1000_25%_overall_bottom_genes_sig.csv")

#Extract concordant chemical perturbagens (CPs) for the top/bottom X% of genes
cg_up_concordants <- getConcordants(cg_up_sig)

#Extract discordant CPs for the top/bottom X% of genes
cg_dn_concordants <- getConcordants(cg_dn_sig)

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
  writexl::write_xlsx("results/6. DrugFindR Output/CG/25%fromL1000/Overall/CG_L1000_25%_Overall_drugfindr.xlsx")
