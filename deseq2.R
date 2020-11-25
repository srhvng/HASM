# Create separate factor objects for each experiment factor.
if(is.null(params$NFactors)){
  #Rename columns of sample info 
  colnames(sinfo) <- c("Sample","Factor2","Factor1")
  colnames(sinfo2) <- c("Sample","Factor2","Factor1")
  
  Factor2 = factor(sinfo[,ncol(sinfo)])
  Factor1 = factor(sinfo[,ncol(sinfo)-1])
  
  # Create a DESeqDataSet object. Make sure that Site is the last
  # factor in the experiment design expression.
  dds=DESeqDataSetFromMatrix(countData=m,
                             colData=sinfo,
                             design= ~ Factor2 + Factor1)
} else if(params$NFactors==3) {
  
    #Rename columns of sample info 
    colnames(sinfo) <- c("Sample","Factor3","Factor2","Factor1")
    colnames(sinfo2) <- c("Sample","Factor3","Factor2","Factor1")
    
    Factor3 = factor(sinfo[,ncol(sinfo)])
    Factor2 = factor(sinfo[,ncol(sinfo)-1])
    Factor1 = factor(sinfo[,ncol(sinfo)-2])
    
    dds=DESeqDataSetFromMatrix(countData=m,
                               colData=sinfo,
                               design= ~ Factor3 + Factor2 + Factor1)
  
} else if(params$NFactors==4){
  
    #Rename columns of sample info 
    colnames(sinfo) <- c("Sample","Factor4","Factor3","Factor2","Factor1")
    colnames(sinfo2) <- c("Sample","Factor4","Factor3","Factor2","Factor1")
    
    Factor4 = factor(sinfo[,ncol(sinfo)-3])
    Factor3 = factor(sinfo[,ncol(sinfo)-2])
    Factor2 = factor(sinfo[,ncol(sinfo)-1])
    Factor1 = factor(sinfo[,ncol(sinfo)])
    
    dds=DESeqDataSetFromMatrix(countData=m,
                               colData=sinfo,
                               design= ~ Factor1 + Factor2 + Factor3 + Factor4)
  
} else if(params$NFactors==5){
  
    #Rename columns of sample info 
    colnames(sinfo) <- c("Sample","Factor5","Factor4","Factor3","Factor2","Factor1")
    colnames(sinfo2) <- c("Sample","Factor5","Factor4","Factor3","Factor2","Factor1")
    
    Factor5 = factor(sinfo[,ncol(sinfo)-4])
    Factor4 = factor(sinfo[,ncol(sinfo)-3])
    Factor3 = factor(sinfo[,ncol(sinfo)-2])
    Factor2 = factor(sinfo[,ncol(sinfo)-1])
    Factor1 = factor(sinfo[,ncol(sinfo)])
    
    dds=DESeqDataSetFromMatrix(countData=m,
                               colData=sinfo,
                               design= ~ Factor1 + Factor2 + Factor3 + Factor4 + Factor5)
  
} else if(params$NFactors==6){
  
    #Rename columns of sample info 
    colnames(sinfo) <- c("Sample","Factor6","Factor5","Factor4","Factor3","Factor2","Factor1")
    colnames(sinfo2) <- c("Sample","Factor6","Factor5","Factor4","Factor3","Factor2","Factor1")
    
    Factor6 = factor(sinfo[,ncol(sinfo)-5])
    Factor5 = factor(sinfo[,ncol(sinfo)-4])
    Factor4 = factor(sinfo[,ncol(sinfo)-3])
    Factor3 = factor(sinfo[,ncol(sinfo)-2])
    Factor2 = factor(sinfo[,ncol(sinfo)-1])
    Factor1 = factor(sinfo[,ncol(sinfo)])
    
    dds=DESeqDataSetFromMatrix(countData=m,
                               colData=sinfo,
                               design= ~ Factor1 + Factor2 + Factor3 + Factor4 + Factor5 + Factor6)
  
} else if(params$NFactors==7){
  
    #Rename columns of sample info 
    colnames(sinfo) <- c("Sample","Factor7","Factor6","Factor5","Factor4","Factor3","Factor2","Factor1")
    colnames(sinfo2) <- c("Sample","Factor7","Factor6","Factor5","Factor4","Factor3","Factor2","Factor1")
    
    Factor7 = factor(sinfo[,ncol(sinfo)-6])
    Factor6 = factor(sinfo[,ncol(sinfo)-5])
    Factor5 = factor(sinfo[,ncol(sinfo)-4])
    Factor4 = factor(sinfo[,ncol(sinfo)-3])
    Factor3 = factor(sinfo[,ncol(sinfo)-2])
    Factor2 = factor(sinfo[,ncol(sinfo)-1])
    Factor1 = factor(sinfo[,ncol(sinfo)])
    
    dds=DESeqDataSetFromMatrix(countData=m,
                               colData=sinfo,
                               design= ~ Factor1 + Factor2 + Factor3 + Factor4 + Factor5 + Factor6 + Factor7)
  
} else {
  print("NFactor - min:2 and max:7 Please provide a number within the range.")
}

# Run DESeq
dds=DESeq(dds)

# Get the results.
results=results(dds)

#resultsNames(dds)
#dds_results=lfcShrink(dds,
#                      coef="Treatment_Vehicle_vs_Drug",
#                      #contrast=c(sapply(strsplit(unique(Factor1), ','), function(x) toString(dQuote(x)))),
#                     res=dds_results,type='normal')
summary(results)

# Save the results to a file.
write.csv(results,file.path(params$outdir,"deseq2_results.csv"))

#Volcano Plots
#png("VolcanoPlot.png")
png(file.path(params$outdir,"VolcanoPlot.png"))
print(EnhancedVolcano(results(dds),lab=rownames(results(dds)),
                      x='log2FoldChange',
                      y='pvalue',
                      pCutoff = params$Pcutoff,
                      FCcutoff = params$FCcutoff))
dev.off()
