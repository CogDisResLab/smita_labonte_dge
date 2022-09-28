# Single drugfindR
#
#

library(tidyverse)
library(drugfindR)

dge <- read.csv("results/kallisto/Orbitofrontal_Cortex_overall_dge.csv")

drugfind_results <- investigate_signature(dge, "CP",
                                          similarity_threshold = 0.2,
                                          gene_column = "SYMBOL",
                                          logfc_column = "logFC",
                                          pval_column = "PValue",
                                          source_name = "OFC_overall")

drugfind_results_2 <- investigate_signature(dge, "CP",
                                          similarity_threshold = 0.2,
                                          discordant = TRUE,
                                          gene_column = "SYMBOL",
                                          logfc_column = "logFC",
                                          pval_column = "PValue",
                                          source_name = "OFC_overall")

