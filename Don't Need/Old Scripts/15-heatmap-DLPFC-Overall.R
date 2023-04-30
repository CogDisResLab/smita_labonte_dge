#Heatmap Generation

library(tidyverse)
library(ComplexHeatmap)
library(pheatmap)
library(org.Hs.eg.db)


#Change GeneIDs to Gene Symbols in DGE
dge <- read.csv("Galaxy Output/DLPFC/DLPFC_Overall_DGE.csv") %>% mutate(GeneID = as.character(GeneID))
gns <-select(org.Hs.eg.db, as.character(dge$GeneID), c("ENTREZID","SYMBOL"),"ENTREZID")
dge <- dge %>% inner_join(gns, by = c("GeneID"="ENTREZID")) %>%
  filter(!is.na(SYMBOL))

#Change GeneIDs to Gene Symbols in Count Matrix
count <- read.csv("Galaxy Count Matrix Output/DLPFC/DLPFC_Overall_Counts.csv") %>% mutate(Geneid = as.character(Geneid))
gns1 <-select(org.Hs.eg.db, as.character(count$Geneid), c("ENTREZID","SYMBOL"),"ENTREZID")
count <- count %>% inner_join(gns, by = c("Geneid"="ENTREZID"))%>%
  filter(!is.na(SYMBOL))

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
  as.matrix()

ctmx <- count_matrix + 1

lctmtx <- log(ctmx)

colnames(lctmtx) <- paste0("Sample", 1:48)

# heatmap(lctmtx)

pheatmap(lctmtx)
