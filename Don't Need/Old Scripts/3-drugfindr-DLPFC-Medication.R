# Single drugfindR
#
#

library(tidyverse)
library(drugfindR)

dge <- read.csv("results/kallisto/DLPFC_medication_dge.csv")

drugfind_results <- investigate_signature(dge, "CP",
                                          similarity_threshold = 0.2,
                                          filter_threshold = 0.5,
                                          gene_column = "SYMBOL",
                                          logfc_column = "logFC",
                                          pval_column = "PValue",
                                          source_name = "DLPFC_medication") |>

  write_csv("figures/DrugFindR_Old_Output/DLPFC_Med1.csv")

drugfind_results_2 <- investigate_signature(dge, "CP",
                                            similarity_threshold = 0.2,
                                            filter_threshold = 0.5,
                                            discordant = TRUE,
                                            gene_column = "SYMBOL",
                                            logfc_column = "logFC",
                                            pval_column = "PValue",
                                            source_name = "DLPFC_medication") |>

  write_csv("figures/DrugFindR_Old_Output/DLPFC_Med2.csv")
