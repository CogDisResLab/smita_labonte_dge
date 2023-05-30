# drugfindR: Manual procurement of all the drug similarity scores
# Comparison: SUB L1000 Females

library(tidyverse)
library(drugfindR)

#Input DGE matrix
dge <- read_csv("results/3. L1000 Heat Maps Post SVA/Females/Human_DGE_L1000_Females_Matrix.csv")

#Extract gene signature (gene name and LFC values)
sub_dge_sig <- dge |>
  prepare_signature(gene_column = "Name_GeneSymbol", logfc_column = "SUB",
                    pval_column = "SUB_Pval")

#Extract concordant chemical perturbagens (CPs) for the L1000 genes
sub_up_concordants <- get_concordants(sub_dge_sig)

#Extract discordant CPs for the L1000 genes
sub_dn_concordants <- get_concordants(sub_dge_sig)

#Combine all CPs and filter
sub_concordants_all <- sub_up_concordants |>
  mutate(similarity_type = if_else(similarity < 0, "Discordant", "Concordant")) |>
  group_by(similarity_type) |>
  nest()

#Filters CPs by top concordant, top discordant, bottom concordant, bottom discordant
x <- sub_concordants_all %>%
  mutate(sheet_name = str_c(similarity_type)) %>%
  ungroup() %>%
  select(-similarity_type) %>%
  select(sheet_name, data) %>%
  deframe() %>%
  writexl::write_xlsx("results/6. DrugFindR Output/SUB/L1000/Females/SUB_L1000_Females_drugfindr.xlsx")
