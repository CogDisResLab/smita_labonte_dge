# drugfindR: Manual procurement of all the drug similarity scores
# Comparison: SUB top/bottom 5% from L1000 Females

library(tidyverse)
library(drugfindR)

#Input DGE matrix with genes for our comparison of interest
dge <- read_csv("results/4. Top 5% Heat Maps Post SVA/Females/Human_DGE_L10005%_Females_Matrix.csv") |>
  select(Name_GeneSymbol, starts_with("SUB")) %>%
  filter(SUB != 0)

#Extract gene signature (gene name and LFC values)
sub_dge_sig <- dge |>
  prepare_signature(gene_column = "Name_GeneSymbol", logfc_column = "SUB",
                    pval_column = "SUB_Pval")

#Create a file listing top X% gene signature
sub_up_sig <- sub_dge_sig |>
  filter_signature("up", threshold = 0) |>
  write_csv("results/6. DrugFindR Output/SUB/5%fromL1000/Females/SUB_L1000_5%_females_top_genes_sig.csv")

#Create a file listing bottom X% gene signature
sub_dn_sig <- sub_dge_sig |>
  filter_signature("down", threshold = 0) |>
  write_csv ("results/6. DrugFindR Output/SUB/5%fromL1000/Females/SUB_L1000_5%_females_bottom_genes_sig.csv")

#Extract concordant chemical perturbagens (CPs) for the top/bottom X% of genes
sub_up_concordants <- get_concordants(sub_up_sig, sig_direction = "Up")

#Extract discordant CPs for the top/bottom X% of genes
sub_dn_concordants <- get_concordants(sub_dn_sig, sig_direction = "Down")

#Combine all CPs and filter
sub_concordants_all <- bind_rows(sub_up_concordants, sub_dn_concordants) |>
  mutate(similarity_type = if_else(similarity < 0, "Discordant", "Concordant")) |>
  group_by(similarity_type, sig_direction) |>
  nest()

#Filters CPs by top concordant, top discordant, bottom concordant, bottom discordant
x <- sub_concordants_all %>%
  mutate(sheet_name = str_c(sig_direction, similarity_type, collapse = "-")) %>%
  ungroup() %>%
  select(-sig_direction, -similarity_type) %>%
  select(sheet_name, data) %>%
  deframe() %>%
  writexl::write_xlsx("results/6. DrugFindR Output/SUB/5%fromL1000/Females/SUB_L1000_5%_Females_drugfindr.xlsx")
