# Make a matrix of the counts. Any column with the word "SRR in its
# title contains count data.
datacols<-grep("SRR",names(newmat))
m<-as.data.frame((newmat[datacols]))

#Heatmaps
# List with colors for each annotation.
colors<- c("cadetblue1","firebrick")
names(colors) = c(strsplit(unique(sinfo$Factor1), ','))

ann_colors =list(
  dex=colors)

# Assign row names to the matrix. Combine the values in the EnsemblId and
# GeneName columns into a row label.
rownames(sinfo)=colnames(m)
colnames(sinfo) <- colnames(sampleinfo) 

sinfo<- sinfo[,-1]

png(file.path("C:/Users/Desktop/ST537/HeatmapAllGenes.png"))
pheatmap(
  mat               = m,
  color             = inferno(10),
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  annotation_col    = sinfo,
  annotation_colors = ann_colors,
  drop_levels       = TRUE,
  cellwidth         = 35, 
  main              = "All Genes in Comparison"
)
dev.off()

###################### Create heatmap for significant genes
sig<- as.data.frame(results) 

sig<- sig %>% mutate(ensgene=row.names(sig))

#filter out non significant genes
sigd<- subset(sig,padj < 0.05) 

#join tables
sgene<- inner_join(df,sigd,by="ensgene")
sgene<- sgene %>% select(-c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))

newmat1<- sgene
for(i in c(2:ncol(sgene)) ){
  newmat1[,i]<- log2(sgene[,2:ncol(sgene)],sgene[,i])
}

datacols<-grep("SRR",names(newmat1))
msig<-as.data.frame((newmat1[datacols]))

# List with colors for each annotation.
colors<- c("cadetblue1","firebrick")
names(colors) = c(strsplit(unique(sinfo2$Factor1), ','))

ann_colors =list(
  dex=colors)

# Assign row names to the matrix. Combine the values in the EnsemblId and
# GeneName columns into a row label.
rownames(sinfo2)=colnames(msig)
colnames(sinfo2) <- colnames(sampleinfo) 

sinfo2<- sinfo2[,-1]

png(file.path("C:/Users/Desktop/ST537/HeatmapDEGenes.png"))
print(pheatmap(
  mat               = msig,
  color             = inferno(10),
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  annotation_col    = sinfo2,
  annotation_colors = ann_colors,
  drop_levels       = TRUE,
  cellwidth         = 35, 
  main              = "Differentially Expressed Genes"
))
dev.off()
