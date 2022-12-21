#Heatmap Generation

library(tidyverse)
library(ComplexHeatmap)
library(pheatmap)
library(org.Hs.eg.db)
library(Orthology.eg.db)
library(org.Mm.eg.db)


#Change GeneIDs to Gene Symbols in DGE
dge <- read.csv("Galaxy Output/Mouse PFC/Mouse_PFC_Overall_DGE.csv") %>%
  mutate(GeneID = as.character(GeneID))
gns <- select(org.Mm.eg.db, as.character(dge$GeneID), c("ENTREZID","SYMBOL"),"ENTREZID")
mapped <- select(Orthology.eg.db, gns$ENTREZID, "Homo_sapiens","Mus_musculus") %>%
  filter(!is.na(Homo_sapiens)) %>%
  mutate(Homo_sapiens = as.character(Homo_sapiens))
human_genes <- mapIds(org.Hs.eg.db, mapped$Homo_sapiens, "SYMBOL", "ENTREZID") %>%
  enframe(name = "ENTREZID", value = "SYMBOL")
mapped_mouse <- mapped |>
  inner_join(human_genes, by = c("Homo_sapiens" = "ENTREZID"))
dge <- dge %>% inner_join(mapped_mouse, by = c("GeneID"="Mus_musculus"))

#Change GeneIDs to Gene Symbols in Count Matrix
count <- read.csv("Galaxy Count Matrix Output/Mouse PFC/Mouse_PFC_Overall_Counts.csv") %>%
  mutate(Geneid = as.character(Geneid))
gns1 <- select(org.Mm.eg.db, as.character(count$Geneid), c("ENTREZID","SYMBOL"),"ENTREZID")
count_mapped <- select(Orthology.eg.db, gns1$ENTREZID, "Homo_sapiens","Mus_musculus") %>%
  filter(!is.na(Homo_sapiens)) %>%
  mutate(Homo_sapiens = as.character(Homo_sapiens))
human_genes1 <- mapIds(org.Hs.eg.db, mapped$Homo_sapiens, "SYMBOL", "ENTREZID") %>%
  enframe(name = "ENTREZID", value = "SYMBOL")
mapped_mouse1 <- mapped |>
  inner_join(human_genes1, by = c("Homo_sapiens" = "ENTREZID"))
count <- count %>% inner_join(mapped_mouse1, by = c("Geneid"="Mus_musculus"))


#Get top 10 genes from DGE
ofc_top_dge <- dge %>%
  arrange(FDR) %>%
  slice_head(n=10)

#Get bottom 10 genes from DGE
ofc_bot_dge <- dge %>%
  arrange(FDR) %>%
  slice_tail(n=10)


filtered_count <- count %>%
  filter (SYMBOL %in% c(ofc_top_dge$SYMBOL, ofc_bot_dge$SYMBOL))

count_matrix <- filtered_count %>%
  column_to_rownames ("SYMBOL") %>%
  dplyr::select(-Geneid) %>%
  mutate(across(everything(), ~ as.numeric(.x))) %>%
  as.matrix()

ctmx <- count_matrix + 1

lctmtx <- log(ctmx)

colnames(lctmtx) <- paste0("Sample", 1:39)

# heatmap(lctmtx)

pheatmap(lctmtx)
