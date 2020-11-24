###PCA plots
# Make a matrix of the counts. Any column with the word "SRR" in its
# title contains count data.
datacols<-grep("SRR",names(newmat))
m<-as.data.frame((newmat[datacols]))

#transpose data
model<- prcomp(t(m))
#summary(model)

#rename sample info columns
colnames(sinfo) <- c("Factor2","Factor1")
colnames(sinfo2) <- c("Factor2","Factor1")

#3d PCA plot
gr<- factor(sinfo$Factor1)
summary(gr)
pca3d(model,group=gr,legend="topleft")
snapshotPCA3d(file=file.path("C:/Users/Desktop/ST537/PCA-allGenes.png"))
print(snapshotPCA3d(file=file.path("C:/Users/Desktop/ST537/PCA-allGenes.png")))

######## Create PCA plot for significant genes
sig<- as.data.frame(results) 

sig<- sig %>% mutate(ensgene=row.names(sig))

sigd<- subset(sig,padj < 0.05) 

sgene<- inner_join(df,sigd,by="ensgene")
sgene<- sgene %>% select(-c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))

newmat1<- sgene
for(i in c(2:ncol(sgene)) ){
  newmat1[,i]<- log2(sgene[,2:ncol(sgene)],sgene[,i])
}

datacols<-grep("SRR",names(newmat1))
msig<-as.data.frame((newmat1[datacols]))

#transpose data
models<- prcomp(t(msig), retx = TRUE, center = TRUE, scale. = FALSE)

#3d PCA plot
gr<- factor(sinfo2$Factor1)
summary(gr)
pca3d(models,group=gr,legend="topleft")
snapshotPCA3d(file=file.path("C:/Users/Desktop/ST537/PCA-DEGenes.png"))
print(snapshotPCA3d(file=file.path("C:/Users/Desktop/ST537/PCA-DEGenes.png")))
