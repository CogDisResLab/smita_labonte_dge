library(tidyverse)

#Load DLPFC metadata: Comparison MDD vs. CTRL ------------------------------------
DLPFC_metadata <- read_csv("data/GSE102556-DLPFC-metadata.csv")
DLPFC_count <-
  read_csv("Galaxy Count Matrix Mapped/DLPFC/DLPFC_Overall_Mapped_Final.csv") %>%
  column_to_rownames("Geneid")

mask <- rowSums(DLPFC_count) > 48

DLPFC_count_filtered <- DLPFC_count[mask,]

#Select row and gender columns for further analysis
DLPFC_metadata_truncated <- DLPFC_metadata %>%
  dplyr::select(Run, gender, phenotype, medication, Cause_of_death) %>%
  column_to_rownames("Run") %>%
  as.matrix()

#Extract principal components
pca = prcomp(t(DLPFC_count_filtered), scale = TRUE)

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

## now make a fancy looking plot that shows the PCs and the variation:
library(ggplot2)

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Sex=DLPFC_metadata$gender,
                       Diagnosis=DLPFC_metadata$phenotype,
                       Death=DLPFC_metadata$Cause_of_death)
pca.data

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample, shape = Death, color = Diagnosis)) +
  geom_point() +
  #geom_label() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("DLPFC MDD vs CTRL (Overall) PCA Graph") +
  scale_color_manual(values = c("darkgreen", "red"), breaks = c("MDD", "CTRL"))

## get the name of the top 10 measurements (genes) that contribute
## most to pc1.
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes ## show the names of the top 10 genes

pca$rotation[top_10_genes,1] ## show the scores (and +/- sign)

