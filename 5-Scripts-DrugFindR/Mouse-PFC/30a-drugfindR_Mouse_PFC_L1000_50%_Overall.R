# drugfindR: Manual procurement of all the drug similarity scores
# Comparison: PFC top/bottom 50% from L1000 Overall

library(tidyverse)
library(drugfindR)

#Input DGE matrix with genes for our comparison of interest
dge <- read_csv("results/5. Mouse Top 50% Heat Maps No SVA/Overall/Mouse_DGE_L1000_50%_Matrix.csv") |>
  select(Name_GeneSymbol, starts_with("PFC")) |>
  filter(PFC != 0)

#Extract gene signature (gene name and LFC values)
pfc_dge_sig <- dge |>
  prepareSignature(geneColumn = "Name_GeneSymbol", logfcColumn = "PFC",
                    pvalColumn = "PFC_Pval")

#Create a file listing top X% gene signature
pfc_up_sig <- pfc_dge_sig |>
  filterSignature("up", threshold = 0) |>
  write_csv("results/6. Mouse DrugFindR Output/Mouse-PFC/50%fromL1000/Overall/Mouse_PFC_L1000_50%_overall_top_genes_sig.csv")

#Create a file listing bottom X% gene signature
pfc_dn_sig <- pfc_dge_sig |>
  filterSignature("down", threshold = 0) |>
  write_csv("results/6. Mouse DrugFindR Output/Mouse-PFC/50%fromL1000/Overall/Mouse_PFC_L1000_50%_overall_bottom_genes_sig.csv")

#Extract concordant chemical perturbagens (CPs) for the top/bottom X% of genes
pfc_up_concordants <- getConcordants(pfc_up_sig)

#Extract discordant CPs for the top/bottom X% of genes
pfc_dn_concordants <- getConcordants(pfc_dn_sig)

#Combine all CPs and filter
pfc_concordants_all <- bind_rows(pfc_up_concordants, pfc_dn_concordants) |>
  mutate(similarity_type = if_else(similarity < 0, "Discordant", "Concordant")) |>
  group_by(similarity_type, sig_direction) |>
  nest()

#Filters CPs by top concordant, top discordant, bottom concordant, bottom discordant
x <- pfc_concordants_all %>%
  mutate(sheet_name = str_c(sig_direction, similarity_type, collapse = "-")) %>%
  ungroup() %>%
  select(-sig_direction, -similarity_type) %>%
  select(sheet_name, data) %>%
  deframe() %>%
  writexl::write_xlsx("results/6. Mouse DrugFindR Output/Mouse-PFC/50%fromL1000/Overall/Mouse_PFC_L1000_50%_Overall_drugfindr.xlsx")
