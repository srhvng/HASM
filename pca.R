###PCA plots
# Make a matrix of the counts. Any column with the word "Sample" in its
# title contains count data.
datacols<-grep("Sample",names(newmat))
m<-as.data.frame((newmat[datacols]))

#transpose data
model<- prcomp(t(m))

#rename sample info columns
colnames(sinfo) <- c("Factor2","Factor1")
colnames(sinfo2) <- c("Factor2","Factor1")

#3d PCA plot
#gr<- factor(sinfo$Factor1)
#summary(gr)
#pca3d(model,group=gr,legend="topleft")
#snapshotPCA3d(file=file.path(params$outdir, "PCA-allGenes.png"))
#print(snapshotPCA3d(file=file.path(params$outdir, "PCA-allGenes.png")))

#2D pca plot
df_out <- as.data.frame(model$x)
df_out$group <- sinfo[,ncol(sinfo)]

percentage <- round(model$sdev / sum(model$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

p<- ggplot(df_out,aes(x=PC1,y=PC2,color=group ))
p<- p+geom_point()+theme+ xlab(percentage[1]) + ylab(percentage[2])
png(file.path(params$outdir, "PCA_ALLGenes.png"))
print(p)
dev.off()


######## Create PCA plot for significant genes
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

#transpose data
models<- prcomp(t(msig), retx = TRUE, center = TRUE, scale. = FALSE)

#3d PCA plot
#gr<- factor(sinfo2$Factor1)
#summary(gr)
#pca3d(models,group=gr,legend="topleft")
#snapshotPCA3d(file=file.path(params$outdir, "PCA-DEGenes.png"))
#print(snapshotPCA3d(file=file.path(params$outdir, "PCA-DEGenes.png")))

#2D pca plot
df_out <- as.data.frame(models$x)
df_out$group <- sinfo2[,ncol(sinfo2)]

percentage <- round(models$sdev / sum(models$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

p<- ggplot(df_out,aes(x=PC1,y=PC2,color=group ))
p<- p+geom_point()+theme+ xlab(percentage[1]) + ylab(percentage[2])
png(file.path(params$outdir,"PCA_DEGenes.png"))
print(p)
dev.off()
