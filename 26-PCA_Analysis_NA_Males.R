library(tidyverse)

#Load NA metadata: Comparison MDD vs. CTRL for Males --------------------------
NA_metadata <- read_csv("data/GSE102556-NA-metadata.csv")
NA_count <-
  read_csv("Galaxy Count Matrix Mapped/NA/NA_Overall_Mapped_Final.csv") %>%
  column_to_rownames("Geneid")

mask <- rowSums(NA_count) > 48

NA_count_filtered <- NA_count[mask,]

#Select row and gender columns for further analysis
NA_metadata_truncated <- NA_metadata %>%
  dplyr::select(Run, gender, phenotype, medication, Cause_of_death) %>%
  column_to_rownames("Run") %>%
  as.matrix()

#Add this line if we want to filter by gender
na_count_male <- NA_count_filtered[, NA_metadata_truncated[,1] == "male"]

#Extract principal components
pca = prcomp(t(na_count_male), scale = TRUE)

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

## now make a fancy looking plot that shows the PCs and the variation:
library(ggplot2)

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Sex=NA_metadata$gender[NA_metadata$gender == "male"],
                       Diagnosis=NA_metadata$phenotype[NA_metadata$gender == "male"])
pca.data

#remove color in line 43 and scale on line 49 when filtering by gender
ggplot(data=pca.data, aes(x=X, y=Y, label=Sample, shape = Sex, color = Diagnosis)) +
  geom_point() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("NA MDD vs CTRL (Males) PCA Graph") +
  scale_color_manual(values = c("darkgreen", "red"), breaks = c("MDD", "CTRL"))

## get the name of the top 10 measurements (genes) that contribute
## most to pc1.
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes ## show the names of the top 10 genes

pca$rotation[top_10_genes,1] ## show the scores (and +/- sign)

