library(tidyverse)

#Load CG metadata: Comparison MDD vs. CTRL -------------------------------------
CG_metadata <- read_csv("data/GSE102556-CG-metadata.csv")
CG_count <-
  read_csv("data/GSE102556-CG-normalized-varianced-counts.csv") %>%
  column_to_rownames("Gene")

mask <- rowSums(CG_count) > 48

CG_count_filtered <- CG_count[mask,]

#Select row and gender columns for further analysis
CG_metadata_truncated <- CG_metadata %>%
  dplyr::select(Run, gender, phenotype, medication, Cause_of_death) %>%
  column_to_rownames("Run") %>%
  as.matrix()

#Extract principal components
pca = prcomp(t(CG_count_filtered), scale = TRUE)

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

## now make a fancy looking plot that shows the PCs and the variation:
library(ggplot2)

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Z=pca$x[,3],
                       Sex=CG_metadata_truncated[rownames(pca$x), 1],
                       Diagnosis=CG_metadata_truncated[rownames(pca$x), 2],
                       Death=CG_metadata_truncated[rownames(pca$x), 4])
pca.data

#remove color in line 43 and scale on line 49 when filtering by gender
ggplot(data=pca.data, aes(x=X, y=Y, label=Sample, shape=Death, color=Diagnosis, size=Z)) +
  geom_point() +
  #geom_label() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  scale_x_continuous(limits = c(-400, 400)) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  scale_y_continuous(limits = c(-400, 400)) +
  scale_size_continuous(name = paste("PC3 - ", pca.var.per[3], "%", sep="")) +
  theme_bw() +
  ggtitle("CG MDD vs CTRL (Overall) PCA Graph") +
  scale_color_manual(values = c("darkgreen", "red"), breaks = c("MDD", "CTRL"))

## get the name of the top 10 measurements (genes) that contribute
## most to pc1.
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes ## show the names of the top 10 genes

pca$rotation[top_10_genes,1] ## show the scores (and +/- sign)
