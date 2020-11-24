suppressWarnings(suppressMessages(library("tidyverse", quietly = T)))
library(ggrepel)
library(factoextra)
library(viridis)
library(DESeq2)
library(EnhancedVolcano)
library(pca3d)
library(rgl)
library(NbClust)
library(pheatmap)
library(RColorBrewer)

# Read in Data
#### obtain the data ####
# load the dataset from input directory
df <- read.delim(paste("C:/Users/Desktop/ST537/count_matrix.txt"))

# Make a matrix of the counts. Any column with the word "SRR" in its
# title contains count data.
datacols<-grep("SRR",names(df))
m<-as.matrix(ceiling(df[datacols]))

# Assign row names to the matrix. Combine the values in the EnsemblId and
# GeneName columns into a row label.
rownames(m)=paste(df$ensgene)

# Read the data frame of sample information, including sample name, 
# subject, tissue type, and treatment.
sampleinfo = read.delim(paste("C:/Users/Desktop/ST537/sample_info.txt"))
sinfo<- sampleinfo
sinfo2<- sampleinfo

#Start DE analysis with DESeq2 (default)
#Rename columns of sample info 
colnames(sinfo) <- c("Sample","Factor2","Factor1")
colnames(sinfo2) <- c("Sample","Factor2","Factor1")

Factor2 = factor(sinfo[,ncol(sinfo)])
Factor1 = factor(sinfo[,ncol(sinfo)-1])

# The DE analysis is focus on Factor1 by default.
dds=DESeqDataSetFromMatrix(countData=m,
                           colData=sinfo,
                           design= ~ Factor2 + Factor1)

# Heatmaps and PCA
## Prepare count matrix for heatmap and PCA plot
#Function to convert input file to log2(n+1) for exploratory plots
log2<- function(x,y){
  #calculates the row average for log2
  avg<- rowMeans(x,na.rm = TRUE)
  
  #calculates the log2(n+1)
  logfunc<- log((y+1)/(avg+1),2)
}

newmat<- df
for(i in c(2:ncol(df)) ){
  newmat[,i]<- log2(df[,2:ncol(df)],df[,i])
}

# Run DESeq
dds=DESeq(dds)

# Get the results.
results=results(dds, alpha=0.05)

resultsNames(dds)
dds_results=lfcShrink(dds,
                      coef=5,
                      res=results,type='normal')
# Generate 
summary(results)

results(dds, alpha = 0.05, lfcThreshold = 0.58)

# Save the results to a file.
write.csv(results,file.path("C:/Users/Desktop/ST537/deseq2_results.csv"))

#Volcano Plots
png(file.path("C:/Users/Desktop/ST537/VolcanoPlot.png"))
cols <- densCols(results$log2FoldChange, -log10(results$pvalue))
plot(results$log2FoldChange, -log10(results$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(0.05), col="brown")

gn.selected <- abs(results$log2FoldChange) > 2.5 & results$padj < 0.05 
text(results$log2FoldChange[gn.selected],
     -log10(results$padj)[gn.selected],
     lab=rownames(results)[gn.selected ], cex=0.4)
dev.off()

# calls on the heatmap and PCA script to generate those plots
source("Z:/data/eahome/Multivariate_DE_analysis/Heatmaps.R")
source("Z:/data/eahome/Multivariate_DE_analysis/PCA.R")
