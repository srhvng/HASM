# Make a matrix of the counts. Any column with the word "Sample" in its
# title contains count data.
datacols<-grep("Sample",names(newmat))
m<-as.data.frame((newmat[datacols]))

#Heatmaps
# List with colors for each annotation.
colors<- c("cadetblue1","firebrick")
names(colors) = c(strsplit(unique(sinfo$Factor1), ','))

ann_colors =list(
  Treatment=colors)

if(is.null(params$Breaks)){
  params$Breaks = seq(-9,0,1)
  brks=params$Breaks
} else {
  brks=params$Breaks
}

if(is.null(params$Breaks)){
  params$Lbreaks = c(-9,-6,-3,0)
  lgnd_brks=params$Lbreaks
} else {
  lgnd_brks=params$Lbreaks
}

# Assign row names to the matrix. Combine the values in the EnsemblId and
# GeneName columns into a row label.
rownames(sinfo)=colnames(m)
colnames(sinfo) <- colnames(sampleinfo) 

sinfo<- sinfo[,-1]

png(file.path(params$outdir, "HeatmapAllGenes.png"))
pheatmap(
  mat               = m,
  breaks            = brks,
  legend_breaks     = lgnd_brks,
  color             = inferno(10),
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  annotation_col    = sinfo,
  annotation_colors = ann_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "All Genes in Comparison"
)
dev.off()

###################### Create heatmap for significant genes
sig<- as.data.frame(results) 

sig<- sig %>% mutate(GeneName=row.names(sig))

if(is.null(params$tool)){
  sigd<- subset(sig,padj < 0.05) 
} else {
  sigd<- subset(sig,adj.P.Val < 0.05) 
}

sgene<- inner_join(df,sigd,by="GeneName")

newmat1<- sgene
for(i in c(2:ncol(sgene)) ){
  newmat1[,i]<- log2(sgene[,2:ncol(sgene)],sgene[,i])
}

datacols<-grep("Sample",names(newmat1))
msig<-as.data.frame((newmat1[datacols]))

# List with colors for each annotation.
colors<- c("cadetblue1","firebrick")
names(colors) = c(strsplit(unique(sinfo2$Factor1), ','))

ann_colors =list(
  Treatment=colors)

# Assign row names to the matrix. Combine the values in the EnsemblId and
# GeneName columns into a row label.
rownames(sinfo2)=colnames(msig)
colnames(sinfo2) <- colnames(sampleinfo) 

sinfo2<- sinfo2[,-1]

png(file.path(params$outdir, "HeatmapDEGenes.png"))
pheatmap(
  mat               = msig,
  breaks            = brks,
  legend_breaks     = lgnd_brks,
  color             = inferno(10),
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  annotation_col    = sinfo2,
  annotation_colors = ann_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "Differentially Expressed Genes"
)
dev.off()
