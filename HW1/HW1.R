rm(list=ls())
library(knitr)
library(GEOquery)

load(file="gse19439.RData")

# make a factor object for Control, latent TB and Positive TB
tmp=as.character(pData(phenoData(gse19439))[,1]) 

J=length(tmp) # J=number of samples 
TBgroup=rep("",J)

for(j in 1:J) TBgroup[j]=substring(tmp[j],1,3)

# make a factor for TBgroup
FTB=factor(TBgroup,levels=c("CON","LTB","PTB")) # get our expression set
X=exprs(gse19439)

# do a quick kruskal-wallis scan
myKrusk=function(i){ 
  cat(i,"...",fill=F) 
  kruskal.test(x=X[i,],g=FTB)$p.value
}


# originally, I ran this in the HW1 directory:
# myPvals=mapply(myKrusk,1:(dim(X)[1])) ;save(file="myPvals.RData",myPvals)
# so that now I can save time and just load:
load("myPvals.RData")

# populate vector with last names of the groups in the class.
# note: the code is written this way, becasue it used to be student nameas and not Group names
GroupLabels=c("Group I","Group II","Group III","Group IV")

# pick the best 4 p-values and assign them to the students.
best4=order(myPvals)[1:4]
# print out list of best 4
print(best4)



###################
#### HW1 ####
myF <- featureData(gse19439)
myPh <- pData(myF)
Our_ID <- rownames(X)[26058]
Our_gene_info <- myPh[which(rownames(myPh) == Our_ID),]
our_gene <- Our_gene_info[,"Symbol"] 
G3 <- X[26058,])

group_data=X[26058,]
probe_name=rownames(X)[26058]

probe_data_feature=pData(featureData(gse19439))
index=match(probe_name,probe_data_feature[,"ID"])
gene_name=probe_data_feature[index,]$Symbol
gene_name$Symbol


probe_data_pheno=pData(phenoData(gse19439))[,1]
probe_data_type=as.factor(substring(probe_data_pheno,1,3))
boxplot(group_data~probe_data_type)

install.packages('vioplot')
library(vioplot)
vioplot(group_data[as.numeric(probe_data_type)==1],
        group_data[as.numeric(probe_data_type)==2],
        group_data[as.numeric(probe_data_type)==3],
        names=c("CON", "LTB", "PTB"),
        col="tomato1")

best20=order(myPvals)[1:20]

library(beanplot)
beanplot(group_data[as.numeric(probe_data_type)==1],
         group_data[as.numeric(probe_data_type)==2],
         group_data[as.numeric(probe_data_type)==3],
         names=c("CON", "LTB", "PTB"), col=c( "turquoise", "purple","royalblue", "tomato"))

library(matrixStats)
library(gplots)
featureVars=rowVars(X)
best20=order(myPvals)[1:20]
best20.X=X[best20,]
colnames(best20.X) <- TBgroup 

par(mar=c(8,10,10,4))
heatmap.2(best20.X, trace="none", margins=c(5,8), scale="row")
