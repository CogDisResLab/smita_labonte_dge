# Manual procurement of the drug similarity scores

library(tidyverse)
library(drugfindR)


dge <- read_csv("results/Human_DGE-L1000_Matrix_Overall.csv")

ai_dge_sig <- dge |>
  prepare_signature(gene_column = "Name_GeneSymbol", logfc_column = "AI",
                    pval_column = NA)

ai_up_sig <- ai_dge_sig |>
  filter_signature("up", prop = 0.05)

ai_dn_sig <- ai_dge_sig |>
  filter_signature("down", prop = 0.05)

ai_up_concordants <- get_concordants(ai_up_sig)

ai_dn_concordants <- get_concordants(ai_dn_sig)

ai_concordants_all <- bind_rows(ai_up_concordants, ai_dn_concordants) |>
  write_csv("/path/to/file.csv") |>
  mutate(similarity_type = if_else(similarity < 0, "Discordant", "Concordant")) |>
  group_by(similarity_type, sig_direction) |>
  nest()

