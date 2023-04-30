# Single drugfindR
#
#

library(tidyverse)
library(drugfindR)
library(org.Hs.eg.db)

dge <- read.csv("Galaxy Output/OFC/OFC_Female_DGE.csv") %>% mutate(GeneID = as.character(GeneID))
gns<-select(org.Hs.eg.db, as.character(dge$GeneID), c("ENTREZID","SYMBOL"),"ENTREZID")
dge <- dge %>% inner_join(gns, by = c("GeneID"="ENTREZID"))

drugfind_results <- investigate_signature(dge, "CP",
                                          similarity_threshold = 0.2,
                                          filter_prop = 0.95,
                                          gene_column = "SYMBOL",
                                          logfc_column = "logFC",
                                          pval_column = "PValue",
                                          source_name = "OFC_Female") |>

  write_csv("figures/DrugFindR_Output/OFC/OFC_Female1.csv")

drugfind_results_2 <- investigate_signature(dge, "CP",
                                            similarity_threshold = 0.2,
                                            filter_prop = 0.95,
                                            discordant = TRUE,
                                            gene_column = "SYMBOL",
                                            logfc_column = "logFC",
                                            pval_column = "PValue",
                                            source_name = "OFC_Female") |>
  write_csv("figures/DrugFindR_Output/OFC/OFC_Female2.csv")
