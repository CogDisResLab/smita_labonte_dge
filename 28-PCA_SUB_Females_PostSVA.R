library(tidyverse)

#Load SUB metadata: Comparison MDD vs. CTRL for Females ------------------------
SUB_metadata <- read_csv("data/GSE102556-SUB-metadata.csv")
SUB_count <-
  read_csv("data/GSE102556-Sub-normalized-varianced-counts.csv") %>%
  column_to_rownames("Gene")

mask <- rowSums(SUB_count) > 48

SUB_count_filtered <- SUB_count[mask,]

#Select row and gender columns for further analysis
SUB_metadata_truncated <- SUB_metadata %>%
  dplyr::select(Run, gender, phenotype, medication, Cause_of_death) %>%
  column_to_rownames("Run") %>%
  as.matrix()

#Add this line if we want to filter by gender
sub_count_female <- SUB_count_filtered[, SUB_metadata_truncated[,1] == "female"]

#for male samples, filter out rows that have low/0 variance
sub_count_female <- sub_count_female[
  rowSums(sub_count_female) > sum(SUB_metadata_truncated[,1] == "female"),
]

#Extract principal components
pca = prcomp(t(sub_count_female), scale = TRUE)

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

## now make a fancy looking plot that shows the PCs and the variation:
library(ggplot2)

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Sex=SUB_metadata$gender[SUB_metadata$gender == "female"],
                       Diagnosis=SUB_metadata$phenotype[SUB_metadata$gender == "female"],
                       Death=SUB_metadata$Cause_of_death[SUB_metadata$gender == "female"])
pca.data

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample, shape = Death, color = Diagnosis)) +
  geom_point() +
  geom_label() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  scale_x_continuous(limits = c(-400, 400)) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  scale_y_continuous(limits = c(-400, 400)) +
  theme_bw() +
  ggtitle("SUB MDD vs CTRL (Females) PCA Graph") +
  scale_color_manual(values = c("darkgreen", "red"), breaks = c("MDD", "CTRL"))

## get the name of the top 10 measurements (genes) that contribute
## most to pc1.
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes ## show the names of the top 10 genes

pca$rotation[top_10_genes,1] ## show the scores (and +/- sign)
