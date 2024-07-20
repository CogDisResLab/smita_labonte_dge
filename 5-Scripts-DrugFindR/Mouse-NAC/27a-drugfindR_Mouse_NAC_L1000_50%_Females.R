# drugfindR: Manual procurement of all the drug similarity scores
# Comparison: NAC top/bottom 50% from L1000 Females

library(tidyverse)
library(drugfindR)

#Input DGE matrix with genes for our comparison of interest
dge <- read_csv("results/5. Mouse Top 50% Heat Maps No SVA/Females/Mouse_DGE_L1000_50%_Females_Matrix.csv") |>
  select(Name_GeneSymbol, starts_with("NAC")) |>
  filter(NAC != 0)

#Extract gene signature (gene name and LFC values)
nac_dge_sig <- dge |>
  prepareSignature(geneColumn = "Name_GeneSymbol", logfcColumn = "NAC",
                    pvalColumn = "NAC_Pval")

#Create a file listing top X% gene signature
nac_up_sig <- nac_dge_sig |>
  filterSignature("up", threshold = 0) |>
  write_csv("results/6. Mouse DrugFindR Output/Mouse-NAC/50%fromL1000/Females/Mouse_NAC_L1000_50%_females_top_genes_sig.csv")

#Create a file listing bottom X% gene signature
nac_dn_sig <- nac_dge_sig |>
  filterSignature("down", threshold = 0) |>
  write_csv("results/6. Mouse DrugFindR Output/Mouse-NAC/50%fromL1000/Females/Mouse_NAC_L1000_50%_females_bottom_genes_sig.csv")

#Extract concordant chemical perturbagens (CPs) for the top/bottom X% of genes
nac_up_concordants <- getConcordants(nac_up_sig)

#Extract discordant CPs for the top/bottom X% of genes
nac_dn_concordants <- getConcordants(nac_dn_sig)

#Combine all CPs and filter
nac_concordants_all <- bind_rows(nac_up_concordants, nac_dn_concordants) |>
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
  writexl::write_xlsx("results/6. Mouse DrugFindR Output/Mouse-NAC/50%fromL1000/Females/Mouse_NAC_L1000_50%_Females_drugfindr.xlsx")
