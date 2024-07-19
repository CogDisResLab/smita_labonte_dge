# drugfindR: Manual procurement of all the drug similarity scores
# Comparison: OFC top/bottom 25% from L1000 Females

library(tidyverse)
library(drugfindR)

#Input DGE matrix with genes for our comparison of interest
dge <- read_csv("results/4. Top 25% Heat Maps Post SVA/Females/Human_DGE_L100025%_Females_Matrix.csv") |>
  select(Name_GeneSymbol, starts_with("OFC")) %>%
  filter(OFC != 0)

#Extract gene signature (gene name and LFC values)
ofc_dge_sig <- dge |>
  prepareSignature(geneColumn = "Name_GeneSymbol", logfcColumn = "OFC",
                    pvalColumn = "OFC_Pval")

#Create a file listing top X% gene signature
ofc_up_sig <- ofc_dge_sig |>
  filterSignature("up", threshold = 0) |>
  write_csv("results/6. DrugFindR Output/OFC/25%fromL1000/Females/OFC_L1000_25%_females_top_genes_sig.csv")

#Create a file listing bottom X% gene signature
ofc_dn_sig <- ofc_dge_sig |>
  filterSignature("down", threshold = 0) |>
  write_csv ("results/6. DrugFindR Output/OFC/25%fromL1000/Females/OFC_L1000_25%_females_bottom_genes_sig.csv")

#Extract concordant chemical perturbagens (CPs) for the top/bottom X% of genes
ofc_up_concordants <- getConcordants(ofc_up_sig)

#Extract discordant CPs for the top/bottom X% of genes
ofc_dn_concordants <- getConcordants(ofc_dn_sig)

#Combine all CPs and filter
ofc_concordants_all <- bind_rows(ofc_up_concordants, ofc_dn_concordants) |>
  mutate(similarity_type = if_else(similarity < 0, "Discordant", "Concordant")) |>
  group_by(similarity_type, sig_direction) |>
  nest()

#Filters CPs by top concordant, top discordant, bottom concordant, bottom discordant
x <- ofc_concordants_all %>%
  mutate(sheet_name = str_c(sig_direction, similarity_type, collapse = "-")) %>%
  ungroup() %>%
  select(-sig_direction, -similarity_type) %>%
  select(sheet_name, data) %>%
  deframe() %>%
  writexl::write_xlsx("results/6. DrugFindR Output/OFC/25%fromL1000/Females/OFC_L1000_25%_Females_drugfindr.xlsx")
