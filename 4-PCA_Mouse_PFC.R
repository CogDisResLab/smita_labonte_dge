library(tidyverse)

#Load Mouse NA metadata:--------------------------------------------------------
PFC_metadata <- read_csv("data/GSE102556-SRA-Mouse-PFC-Metadata.csv")
PFC_count <-
  read_csv("Galaxy Count Matrix Mapped/Mouse PFC/PFC_Mouse_Overall_Mapped_Final.csv") %>%
  column_to_rownames("Geneid")

mask <- rowSums(PFC_count) > 48

PFC_count_filtered <- PFC_count[mask,]

#Select row and gender columns for further analysis
PFC_metadata_truncated <- PFC_metadata %>%
  dplyr::select(Run, gender, Phenotype) %>%
  column_to_rownames("Run") %>%
  as.matrix()

#Extract principal components
pca = prcomp(t(PFC_count_filtered), scale = TRUE)

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

## now make a fancy looking plot that shows the PCs and the variation:
library(ggplot2)

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Z=pca$x[,3],
                       Sex=PFC_metadata_truncated[rownames(pca$x), 1],
                       Diagnosis=PFC_metadata_truncated[rownames(pca$x), 2])
pca.data

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample, shape=Sex, color=Diagnosis, size=Z)) +
  geom_point() +
  geom_label() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  scale_x_continuous(limits = c(-400, 400)) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  scale_y_continuous(limits = c(-400, 400)) +
  scale_size_continuous(name = paste("PC3 - ", pca.var.per[3], "%", sep="")) +
  theme_bw() +
  ggtitle("Mouse PFC PCA Graph") +
  scale_color_manual(values = c("darkgreen", "red"), breaks = c("Stress", "CTRL"))

## get the name of the top 10 measurements (genes) that contribute
## most to pc1.
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes ## show the names of the top 10 genes

pca$rotation[top_10_genes,1] ## show the scores (and +/- sign)

