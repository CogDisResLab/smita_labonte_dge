#BiocManager::install("DEFormats")
library(DEFormats)
library(dplyr)
library(tidyverse)
setwd("C:/Users/xzhang11/OneDrive - University of Toledo/Desktop/Harry Project/DEG")
#setwd("C:/Users/daisy/Desktop/Harry Project/DEG")
df <- read.csv("Control_vs_Bipolar.genes.counts.csv", row.names =1, check.names = FALSE)
head(df)
colnames(df)[1] <- "SYMBOL"
colnames(df)
list.files()
# remove"1005_s" form sample names
(colnames(df) <- gsub(x = colnames(df), pattern = "1005_S", replacement = ""))


#library(readxl)
#phenotype <- read_excel("Stanley Foundation cohort data.xlsx")
#head(phenotype)
#colnames(phenotype)
#phenotype <- mutate(phenotype, sample = paste(`S code`, DDx, sep = '_'))
#phenotype$sample <- gsub(x = phenotype$sample, pattern = "SZ", replacement = "Schizophrenia")
#phenotype$sample <- gsub(x = phenotype$sample, pattern = "C", replacement = "Control")
#phenotype$sample <- gsub(x = phenotype$sample, pattern = "BD", replacement = "BipolarDisorder")
#phenotype$sample <- gsub(x = phenotype$sample, pattern = "MDD", replacement = "MajorDepression")
#write.csv(phenotype, "phenotype.csv", row.names = F)
phenotype <- read.csv("phenotype.csv", check.names = FALSE)

library(edgeR)
library(DESeq2)
library('variancePartition')
library(BiocParallel)
library (scales)
cData <- df
cData$ENSEMBL <- rownames(df)
gData <- cData[, 32, drop = FALSE]
cData <- cData[,-c(1,32)]
head(cData)
pData <- phenotype
pData <- pData %>% 
  rename("lifetime antipsychotics (fluphenazine eq)" = "antipsychotics",
         "smoking history" = "smoking_history",
         "substance abuse severity" = "substance_abuse_severity",
         "ethanol severity" = "ethanol_severity",
         "brain weight" = "brain_weight")
head(pData)
# replace missing values with appropriate values
pData$ethanol_severity[is.na(pData$ethanol_severity)] <- -1
pData$onset[is.na(pData$onset)] <- 0
pData$duration[is.na(pData$duration)] <- 0
pData$smoking_history[is.na(pData$smoking_history)] <- 0
pData$substance_abuse_severity <- as.factor(pData$substance_abuse_severity)
pData$ethanol_severity <- as.factor(pData$ethanol_severity)
pData <- pData[c(52,51,45,19,54,47,49,41,13,26,21,60,59,40,42,27,28,43,15,25,48,6,18,11,35,5,22,20,12,29),]
head(pData)
rownames(pData) <- pData$sample
colnames(cData) == rownames(pData)
group<-factor(pData$DDx)
t <- pData


# creat DGEList object
x <- DGEList(counts = cData,lib.size = colSums(cData), 
             norm.factors = rep(1,ncol(cData)),samples = pData, 
             genes = gData, group = group, remove.zeros = F)

dds = as.DESeqDataSet(x)
dds
dds <- estimateSizeFactors( dds )
ddsNormalized <- counts(dds,normalized=TRUE)
write.csv(ddsNormalized, file = "Control_vs_Bipolar.genes.nCounts.csv")
# filter counts
keep.exprs <- filterByExpr(x, group = group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
dds = as.DESeqDataSet(x)
dds
keep <- rowSums(fpm(dds,robust = T)>1) >= 0.5 * ncol(dds)
sum(keep)
dds <- dds[keep,]
dds
dds <- estimateSizeFactors( dds )
ddsNormalized <- counts(dds,normalized=TRUE)
write.csv(ddsNormalized, file = "Control_vs_Bipolar.genes.filtered.nCounts.csv")

#select which variable to use
# First find canonical correlation between the variables highly correlated will be removed
colnames(colData(dds))
form <- ~  ZT  + age + onset + duration + pH + brain_weight + PMI + antipsychotics + 
  suicide + psychosis + Sex +  group + smoking_history + substance_abuse_severity + 
  ethanol_severity + hemisphere + DDx
coldata <- as.data.frame(colData(dds))
C = canCorPairs(form, coldata)
plotCorrMatrix( C )

#Variance explained

pData$ZT <- rescale(pData$ZT, to = c(0, 100))
pData$age <- rescale(pData$age, to = c(0, 100))
pData$onset <- rescale(pData$onset, to = c(0, 100))
pData$duration <- rescale(pData$duration, to = c(0, 100))
pData$pH <- rescale(pData$pH, to = c(0, 100))
pData$brain_weight <- rescale(pData$brain_weight, to = c(0, 100))
pData$PMI <- rescale(pData$PMI, to = c(0, 100))
pData$antipsychotics <- rescale(pData$antipsychotics, to = c(0, 100))





cData_edgeR <- DGEList(cData, group = pData$group)
cData_edgeR <- cData_edgeR[filterByExpr(cData_edgeR, group = 'group'), ,keep.lib.sizes=FALSE]
cData_edgeR <- calcNormFactors(cData_edgeR)
cData_norm <- cpm(cData_edgeR, normalized.lib.sizes=TRUE, log = FALSE)




head(pData)
varFormula <- ~ ZT  + age + onset + duration + pH + brain_weight + PMI + antipsychotics + 
  (1|DDx) + (1|suicide) + (1|psychosis) + (1|Sex) +  (1|group) + (1|smoking_history) + (1|substance_abuse_severity) + 
  (1|ethanol_severity) + (1|hemisphere) 
varPart <- fitExtractVarPartModel(ddsNormalized, varFormula, pData)
#core_num <- 20
#param <- SnowParam(workers = core_num, type = "SOCK")
# varPart <- fitExtractVarPartModel(cData_norm, varFormula, pData) # 1167 s

vp <- sortCols( varPart )
head(vp)
plotPercentBars( vp[1:10,] )
plotVarPart( vp )
abline(h = 5)
head(varPart[order(varPart$group, decreasing=TRUE),])
                         # NA to 0    # ethanol_severity NA to -1, onset&duration NA to 0, smoking_history NA to 0
mean(varPart$group) *100 #0.6035091   #0.8498122 
#ZT: 1.931757   1.702718
#age: 3.522274    3.598414
#onset: 6.156074    6.096972
#duration: 9.109074   9.905056
#pH: 6.12779    7.2824
#brain_weight: 1.727815   1.535046
#PMI: 3.329859    3.070166
#antipsychotics: 2.651601   2.402741
#DDx: 0.6752026   0.8408914
#suicide: 6.284175    7.114361
#psychosis: 1.35117   1.135828
#Sex: 1.908446    2.024489
#smoking_history: 1.397521    1.095201
#substance_abuse_severity: 7.937876   7.819332
#ethanol_severity: 3.377884   3.833099
#hemisphere: 1.25608    1.039426
#Residuals:     38.65405
save(varPart, file = "Control_vs_Bipolar_varPart_Unbinned_RemoveNA.RData")
write.csv(vp, file = "Control_vs_Bipolar_varPart_Unbinned_RemoveNA.csv")
