# 
library("GEOquery")
# 
# # Berry UK Test (array)
# gse19444=getGEO("GSE19444",GSEMatrix=T)
# show(gse19444)
pdat19444=pData(phenoData(gse19444[[1]]))[,c(1,11:16)]
# dim(pdat19444)

# # Cleaning the phenotype data while keeping the original unchanged
pdat <- pdat19444[,-1]
#
# 
foo <- function(x){
  temp <- strsplit(as.character(x), ": ")
  sapply(temp,"[", 2)}
pdat <- apply(pdat, 2, foo)
colnames(pdat) <- c("gender", "ethnicity", "illness", "geographical_region", "bcg_vaccinated", "region_of_birth")
rownames(pdat) <- rownames(pdat19444)

save(file="pdat.RData", pdat)
save(file="gse19444.RData", gse19444)

rm(list=ls())

load("gse19444.RData")
load("pdat.RData")

## gene expression
exp.data <- exprs(gse19444[[1]])
dim(exp.data)

## c ###

pc <- prcomp(exp.data, center = TRUE, scale = TRUE)
summary(pc)
plot(pc, type="lines")

par(mfrow=c(2,3), mar=c(2,2,2,2))
palette <- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02')
var<-NULL
for(i in 1:ncol(pdat)){
  var <-as.factor(pdat[,i])
  clr <- palette[1:length(levels(var))]
  plot(pc$rotation, col=clr[var], pch=19, main=colnames(pdat)[i], cex.main=1.8)
  legend("topright", pch=19, levels(var), col=clr)
}

### d

library(diptest)

exp.dip <- apply(exp.data, 1, dip.test)
dip.pval<-NULL
for(i in 1:nrow(exp.data)){dip.pval[i] <- dip.test(x=exp.data[i,])$p.value}

feat.select <- function(cutoff){

feats <- exp.data[which(dip.pval<sort(dip.pval)[cutoff+1]),]
pc.f <- prcomp(feats, center = TRUE, scale = F)

par(mfrow=c(2,3), mar=c(2,2,2,2))
palette <- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02')
var<-NULL
for(i in 1:ncol(pdat)){
  var <-as.factor(pdat[,i])
  clr <- palette[1:length(levels(var))]
  plot(pc.f$rotation, col=clr[var], pch=19, main=paste(colnames(pdat)[i], "-", cutoff, "features", sep=" "), cex.main=1.8)
  legend("topright", pch=19, levels(var), col=clr)
}
}

for(i in 1:5){feat.select(cutoff = c(50, 200, 500, 1000, 1500)[i])}


#### Q2 ####

### Clustering - Dendrogram ###
rm(list=ls())
library("GEOquery")
library("dendextend")

load("gse19444.RData")
load("pdat.RData")#cleaned phenotype data saved as in the first problem
exp.data <- exprs(gse19444[[1]])
pdat <- as.data.frame(pdat)
summary(pdat$illness)
 
groupCodes <- c(rep("Control", 12), rep("Latent",21), rep("PTB", 21))
colnames(exp.data) <- groupCodes
colorCodes <- c(Control="#008837", Latent='#c2a5cf', PTB="#7b3294")


d.euc <- dist(t(exp.data), method = "euclidean")
d.man <- dist(t(exp.data), method = "manhattan")
linkage <- c( "single", "complete", "centroid")

myclust <-function(dis){
par(mfrow=c(2,2), mar=c(3,2,3,2))
for(i in 1:length(linkage)){
  hc <- hclust(dis, method=linkage[i])  # apply hirarchical clustering with different linkages
  hcd<-as.dendrogram(hc) #, dLeaf=0.1, h=0.3)
  
  # Assigning the labels of dendrogram object with new colors:
  labels_colors(hcd) <- colorCodes[groupCodes][order.dendrogram(hcd)]
  # Plotting the new dendrogram
  plot(hcd, main=paste(linkage[i], "linkage", sep=" "))
  }
}
myclust(d.euc)
myclust(d.man)

# c #

cut100 <-exp.data[which(dip.pval<sort(dip.pval)[100+1]),]
cut1k <- exp.data[which(dip.pval<sort(dip.pval)[1000+1]),]

d.euc100 <- dist(t(cut100), method = "euclidean")
d.man <- dist(t(cut1k), method = "manhattan")

myclust(d.euc100)
