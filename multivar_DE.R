suppressWarnings(suppressMessages(library("tidyverse", quietly = T)))
library(ggrepel)
library(getopt)
library(factoextra)
library(viridis)

#### set up the arguments to get from docker call in multivariate_de_analysis.R ####
spec <- matrix(c(
  'Countdir'             , 'c', 1, "character", "Directory and file that contains the rsem_count.txt files. Required",
  'SampleInfo'           , 's', 1, "character", "Directory that contains the sample information file (sampleinfo.txt). Required.",
  'outdir'               , 'o', 2, "character", "alternate output directory, default is current working directory",
  'tool'                 , 't', 2, "character", "Specify DESeq2 or limma tool. Default set to DESeq2",
  'FCcutoff'             , 'f', 2, "integer"  , "Fold change cutoff for volcano plot. Default: 0.585",
  'Pcutoff'              , 'p', 2, "integer"  , "P-value cutoff for volcano plot. Default: 0.05",
  'NFactors'             , 'n', 2, "interger" , "Number for factors going into design model. Default: 2. Maximum: 7.",
  'Breaks'               , 'b', 2, "character", "Set breaks for heatmap. Default: seq(-9,0,1). Please provide in same format.",
  'Lbreaks'              , 'l', 2, "character", "Set legend breaks for heatmap. Default: c(-9,-6,-3,0). Please provide in same format.",
  'help'                 , 'h', 0, "logical"  , "Help request"
),ncol=5,byrow=T)

args <- commandArgs(trailingOnly = TRUE)
hh <- paste(unlist(args),collapse=' ')
listoptions <- unlist(strsplit(hh,'--'))[-1]
options.args <- sapply(listoptions,function(x){
  unlist(strsplit(x, ' '))[-1]
}, simplify=FALSE)
options.names <- sapply(listoptions,function(x){
  option <-  unlist(strsplit(x, ' '))[1]
})
names(options.args) <- unlist(options.names)

params = options.args

######Testing
#params$Countdir = "Z:/User/count_matrix.txt"
#params$SampleInfo = "Z:/User/sample_info.txt"

# if help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(params$help)){
  cat(getopt(spec, usage =T))
  q(status=1)
}
if(is.null(params$Countdir)){
  cat('\n--Countdir and --SampleInfo are both required, use the --help flag for more information.\n')
  q(status=1)
}
if(is.null(params$SampleInfo)){
  cat('\n--Countdir and --SampleInfo are both required, use the --help flag for more information.\n')
  q(status=1)
}
if(is.null(params$outdir)){
  params$outdir = getwd()
} else {
  if(substr(x = params$outdir,start = nchar(params$outdir),stop = nchar(params$outdir))=="/"){
    params$outdir = substr(x = params$outdir,start = 1,stop = nchar(params$outdir)-1)
  }
}
if(is.null(params$FCcutoff)){
  params$FCcutoff = 0.585
} 
if(is.null(params$Pcutoff)){
  params$Pcutoff = 0.05
} 

library(DESeq2)
#library(SummarizedExperiment)
library(edgeR)
library(limma)
library(pheatmap)
library(RColorBrewer)
library(pca3d)
library(NbClust)
library(EnhancedVolcano)

# Read in Data
#### obtain the data ####
# load the dataset from input directory
df <- read.delim(paste(params$Countdir))
#df <- read.delim(sapply(params$Countdir, function(x) toString(dQuote(x))))

# Make a matrix of the counts. Any column with the word "Sample" in its
# title contains count data.
datacols<-grep("Sample",names(df))
m<-as.matrix(ceiling(df[datacols]))

# Assign row names to the matrix. Combine the values in the EnsemblId and
# GeneName columns into a row label.
rownames(m)=paste(df$GeneName)

# Read the data frame of sample information, including sample name, 
# subject, tissue type, and treatment.
sampleinfo = read.delim(paste(params$SampleInfo))
#sampleinfo = read.delim(sapply(params$SampleInfo, function(x) toString(dQuote(x)))) 
sinfo<- sampleinfo
sinfo2<- sampleinfo

#Start DE analysis with DESeq2 (default)
if(is.null(params$tool)){
  source("deseq2.R")
} else {
  source("limma.R")
}

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

source("heatmaps.R")
source("pca.R")
