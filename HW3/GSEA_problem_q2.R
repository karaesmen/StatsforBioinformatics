#########################
######  GSEA  ######
#########################  
rm(list=ls())
graphics.off()
require(GEOquery)
require(affy)
require("biomaRt")
install.packages("piano")
library(piano)

#### GSA Example ####
## 1. Load the gse19439 data and perform the Kruskal Wallis scan mentioned in HW\#1:

# The Kruskal-Wallis H test (sometimes also called the "one-way ANOVA on ranks") 
# is a rank-based nonparametric test that can be used to determine if there are 
# statistically significant differences between two or more groups of 
# an independent variable on a continuous or ordinal dependent variable.


LoadScanFlag=F  # flip this flag to T to load the data and perform the KW scan
#takes long time
if(LoadScanFlag){
  
  load("~/Google Drive/Stats for Bioinformatics/HW3/gse19439.RData")
#   gse19439=getGEO("GSE19439",GSEMatrix=T)
#   gse19439=gse19439[[1]]
  
  tmp=as.character(pData(phenoData(gse19439))[,1])
  J=length(tmp) # J=number of samples
  TBgroup=rep("",J)
  for(j in 1:J) TBgroup[j]=substring(tmp[j],1,3)
  # make a factor for TBgroup
  FTB=factor(TBgroup,levels=c("CON","LTB","PTB"))
  
  # get our expression set
  X=exprs(gse19439)
  
  #
  # Let's do a simple Kruskal-Wallis Scan across all features 
  #
  
  # do a quick kruskal-wallis scan ?
  myKrusk=function(i){
    cat(i,"...",fill=F)
    kruskal.test(x=X[i,],g=FTB)$p.value
  }
  
  myPvals=mapply(myKrusk,1:(dim(X)[1])) ;save(file="myPvals.RData",myPvals)  
  save(gse19439,myPvals,FTB,TBgroup,X,file="HW3loadscan.RData")
}

if(!LoadScanFlag) load("HW3loadscan.RData")

# 2. Now, we have to look at the "features" of our ExpressionSet. 
#    What EntrezIDs do they map to?

myfeat=featureData(gse19439)
# what annotation information do we have?
varLabels(myfeat)

# o.k., but what do those varLabels means?
myvarmdat=varMetadata(myfeat)

# Now we know what info we have, how can we extract it and have a look?
mypdat=pData(myfeat)

# now get the EntrezIDs
Entrezs <- mypdat$Entrez_Gene_ID

# 3. Get minimum p-value (the one we got from the 1st part with Kriskal-Wallis) 
#    across the set of features that map to the same EntrezID.

UnqFlag=F # change this flag to T to run this code the first time..

if(UnqFlag){
  # This is the EntrezIDs:
  
  unqEntrez=unique(Entrezs)
  mapEntrz=rep(NA,length(unqEntrez))
  
  # for loop takes a long time
  for(i in 1:length(unqEntrez)){
    edx=which(Entrezs==unqEntrez[i])
    mapEntrz[i]=edx[order(myPvals[edx])[1]]
  }
  
  mapEntrz=mapEntrz[unqEntrez!=""]
  names(mapEntrz)=unqEntrez[unqEntrez!=""]
  myEntrezs=names(mapEntrz)
  myX=X[mapEntrz,]
  myP=myPvals[mapEntrz]
  
  save(mapEntrz,myEntrezs,myX,myP,file="UnqEntrezs.RData")
}

if(!UnqFlag) load("UnqEntrezs.RData")

# 4. We want to run the runGSA command in piano. 
# In order to do that, we will need to map each EntrezID to it's GO ID(s). 
# Note that EntrezIDs may map to multiple GO IDs. We will use biomaRt to accomplish this..

require("biomaRt")
require("piano")

listEnsembl()

# Trying to find the attribute I am interested in
# Exported and opened with a text editor and found the attribute i'm looking for 
atrr <- listAttributes(ensembl)
write.table(atrr, file="./attr.txt", sep="\t")
listAttributes(ensembl)[36,]

filter <- listFilters(ensembl)
write.table(filter, file="./filt.txt", sep="\t")



# first, we need an ensembl.. Since our data is on human subjects..
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

GOIDSflag=F # trip this flag to  run getBM - this one takes some time, however.

if(GOIDSflag){
  goids = getBM(attributes=c("entrezgene","go_id"), filters="entrezgene", 
                values=myEntrezs, mart=ensembl)  
  myGsc <- loadGSC(goids)
  save(goids,myGsc,file="goids.RData")
}

if(!GOIDSflag) load("goids.RData")

head(goids)