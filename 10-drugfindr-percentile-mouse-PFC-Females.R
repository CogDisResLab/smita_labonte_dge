# Single drugfindR
#
#

library(tidyverse)
library(drugfindR)
library(Orthology.eg.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

dge <- read.csv("Galaxy Output/Mouse PFC/Mouse_PFC_Females_DGE.csv") %>% mutate(GeneID = as.character(GeneID))

gns <- select(org.Mm.eg.db, as.character(dge$GeneID), c("ENTREZID","SYMBOL"),"ENTREZID")

mapped <- select(Orthology.eg.db, gns$ENTREZID, "Homo_sapiens","Mus_musculus") %>%
  filter(!is.na(Homo_sapiens)) %>%
  mutate(Homo_sapiens = as.character(Homo_sapiens))

human_genes <- mapIds(org.Hs.eg.db, mapped$Homo_sapiens, "SYMBOL", "ENTREZID") %>%
  enframe(name = "ENTREZID", value = "SYMBOL")

mapped_mouse <- mapped |>
  inner_join(human_genes, by = c("Homo_sapiens" = "ENTREZID"))

dge <- dge %>% inner_join(mapped_mouse, by = c("GeneID"="Mus_musculus"))

drugfind_results <- investigate_signature(dge, "CP",
                                          similarity_threshold = 0.2,
                                          filter_prop = 0.95,
                                          gene_column = "SYMBOL",
                                          logfc_column = "logFC",
                                          pval_column = "PValue",
                                          source_name = "Mouse_PFC_Female") |>

  write_csv("figures/DrugFindR_Output/Mouse PFC/Mouse_PFC_Female1.csv")

drugfind_results_2 <- investigate_signature(dge, "CP",
                                            similarity_threshold = 0.2,
                                            filter_prop = 0.95,
                                            discordant = TRUE,
                                            gene_column = "SYMBOL",
                                            logfc_column = "logFC",
                                            pval_column = "PValue",
                                            source_name = "Mouse_PFC_Female") |>
  write_csv("figures/DrugFindR_Output/Mouse PFC/Mouse_PFC_Female2.csv")
