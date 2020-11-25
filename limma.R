# Create DGElist object
dobj<- DGEList(m)

# Create separate factor objects for each experiment factor.
if(is.null(params$NFactors)){
  params$NFactors = 2
  
  colnames(sinfo) <- c("Sample","Factor2","Factor1")
  colnames(sinfo2) <- c("Sample","Factor2","Factor1")
  Factor1 = factor(sinfo[,ncol(sinfo)])
  Factor2 = factor(sinfo[,ncol(sinfo)-1])  
  
  # Filter low expressed genes in at least 3 sample, by default
  keep<- filterByExpr(dobj,group=Factor1)
  d<- dobj[keep,,keep.lib.sizes=FALSE]
  
  #calculate normalization factors for downstream
  #Trimmed mean of M-values - scaling factor for library sizes
  dobj<- calcNormFactors(d, method="TMM")
  dobj$samples$norm.factors
  
  #Specify model to fit. 
  ##This is for voom as it uses variances of the model residuals
  mm<- model.matrix(~ 0 + Factor1 + Factor2)
  
} else if(params$NFactors==3){
  
    colnames(sinfo) <- c("Sample","Factor3","Factor2","Factor1")
    colnames(sinfo2) <- c("Sample","Factor3","Factor2","Factor1")
    
    Factor1 = factor(sinfo[,ncol(sinfo)])
    Factor2 = factor(sinfo[,ncol(sinfo)-1])
    Factor3 = factor(sinfo[,ncol(sinfo)-2])
    
    # Filter low expressed genes in at least 3 sample, by default
    keep<- filterByExpr(dobj,group=Factor1)
    d<- dobj[keep,,keep.lib.sizes=FALSE]
    
    #calculate normalization factors for downstream
    #Trimmed mean of M-values - scaling factor for library sizes
    dobj<- calcNormFactors(d, method="TMM")
    dobj$samples$norm.factors
    
    #Specify model to fit. 
    ##This is for voom as it uses variances of the model residuals
    mm<- model.matrix(~ 0 + Factor1 + Factor2 +Factor3)
  
} else if(params$NFactors==4){
  
    colnames(sinfo) <- c("Sample","Factor4","Factor3","Factor2","Factor1")
    colnames(sinfo2) <- c("Sample","Factor4","Factor3","Factor2","Factor1")
    
    Factor1 = factor(sinfo[,ncol(sinfo)])
    Factor2 = factor(sinfo[,ncol(sinfo)-1])
    Factor3 = factor(sinfo[,ncol(sinfo)-2])
    Factor4 = factor(sinfo[,ncol(sinfo)-3])
    
    # Filter low expressed genes in at least 3 sample, by default
    keep<- filterByExpr(dobj,group=Factor1)
    d<- dobj[keep,,keep.lib.sizes=FALSE]
    
    #calculate normalization factors for downstream
    #Trimmed mean of M-values - scaling factor for library sizes
    dobj<- calcNormFactors(d, method="TMM")
    dobj$samples$norm.factors
    
    #Specify model to fit. 
    ##This is for voom as it uses variances of the model residuals
    mm<- model.matrix(~ 0 + Factor1 + Factor2 + Factor3 + Factor4)
  
} else if(params$NFactors==5){
  
    colnames(sinfo) <- c("Sample","Factor5","Factor4","Factor3","Factor2","Factor1")
    colnames(sinfo2) <- c("Sample","Factor5","Factor4","Factor3","Factor2","Factor1")
    
    Factor1 = factor(sinfo[,ncol(sinfo)])
    Factor2 = factor(sinfo[,ncol(sinfo)-1])
    Factor3 = factor(sinfo[,ncol(sinfo)-2])
    Factor4 = factor(sinfo[,ncol(sinfo)-3])
    Factor5 = factor(sinfo[,ncol(sinfo)-4])
    
    # Filter low expressed genes in at least 3 sample, by default
    keep<- filterByExpr(dobj,group=Factor1)
    d<- dobj[keep,,keep.lib.sizes=FALSE]
    
    #calculate normalization factors for downstream
    #Trimmed mean of M-values - scaling factor for library sizes
    dobj<- calcNormFactors(d, method="TMM")
    dobj$samples$norm.factors
    
    #Specify model to fit. 
    ##This is for voom as it uses variances of the model residuals
    mm<- model.matrix(~ 0 + Factor1 + Factor2 +Factor3 + Factor4 + Factor5)
  
} else if(params$NFactors==6){
  
    colnames(sinfo) <- c("Sample","Factor6","Factor5","Factor4","Factor3","Factor2","Factor1")
    colnames(sinfo2) <- c("Sample","Factor6","Factor5","Factor4","Factor3","Factor2","Factor1")
    
    Factor1 = factor(sinfo[,ncol(sinfo)])
    Factor2 = factor(sinfo[,ncol(sinfo)-1])
    Factor3 = factor(sinfo[,ncol(sinfo)-2])
    Factor4 = factor(sinfo[,ncol(sinfo)-3])
    Factor5 = factor(sinfo[,ncol(sinfo)-4])
    Factor6 = factor(sinfo[,ncol(sinfo)-5])
    
    # Filter low expressed genes in at least 3 sample, by default
    keep<- filterByExpr(dobj,group=Factor1)
    d<- dobj[keep,,keep.lib.sizes=FALSE]
    
    #calculate normalization factors for downstream
    #Trimmed mean of M-values - scaling factor for library sizes
    dobj<- calcNormFactors(d, method="TMM")
    dobj$samples$norm.factors
    
    #Specify model to fit. 
    ##This is for voom as it uses variances of the model residuals
    mm<- model.matrix(~ 0 + Factor1 + Factor2 +Factor3 + Factor4 + Factor5 + Factor6)
  
} else if(params$NFactors==7){
  
    colnames(sinfo) <- c("Sample","Factor7","Factor6","Factor5","Factor4","Factor3","Factor2","Factor1")
    colnames(sinfo2) <- c("Sample","Factor7","Factor6","Factor5","Factor4","Factor3","Factor2","Factor1")
    
    Factor1 = factor(sinfo[,ncol(sinfo)])
    Factor2 = factor(sinfo[,ncol(sinfo)-1])
    Factor3 = factor(sinfo[,ncol(sinfo)-2])
    Factor4 = factor(sinfo[,ncol(sinfo)-3])
    Factor5 = factor(sinfo[,ncol(sinfo)-4])
    Factor6 = factor(sinfo[,ncol(sinfo)-5])
    Factor7 = factor(sinfo[,ncol(sinfo)-6])
     
    # Filter low expressed genes in at least 3 sample, by default
    keep<- filterByExpr(dobj,group=Factor1)
    d<- dobj[keep,,keep.lib.sizes=FALSE]
    
    #calculate normalization factors for downstream
    #Trimmed mean of M-values - scaling factor for library sizes
    dobj<- calcNormFactors(d, method="TMM")
    dobj$samples$norm.factors
    
    #Specify model to fit. 
    ##This is for voom as it uses variances of the model residuals
    mm<- model.matrix(~ 0 + Factor1 + Factor2 +Factor3 + Factor4 + Factor5 + Factor6 + Factor7)
  
} else {
  print("NFactor - min:2 and max:7 Please provide a number within the range.")
}

#Specify contrast
contr<- makeContrasts(
  paste(colnames(mm)[1],"-",colnames(mm)[2], sep=""),
  levels=colnames(mm)
)
#Check if our contrast is right
contr

#Apply voom function
y<-voom(d,mm,plot=F)

#Fit model
fit <- lmFit(y, mm)
vfit<- contrasts.fit(fit,contrasts = contr)
head(coef(fit))

#Apply DE analysis
efit<- eBayes(vfit)

#What genes are DE?
results<- topTable(efit, sort.by = "P", n=Inf)

#To compare with voom: mean-var plot upstream
#plotSA(efit, main="Final Model: Mean-variance trend")

de<-decideTests(efit)

#export results
write.csv(results,file.path(params$outdir,"limma_results.csv"))

#Volcano Plots
png(file.path(params$outdir,"VolcanoPlot.png"))
print(EnhancedVolcano(results,lab=rownames(results),
                x='logFC',
                y='P.Value',
                pCutoff = params$Pcutoff,
                FCcutoff = params$FCcutoff))
dev.off()
