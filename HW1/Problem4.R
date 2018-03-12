rm(list=ls())

source("http://bioconductor.org/biocLite.R")
biocLite("multtest")				
library(multtest)	
library(affy)

setwd("~/Google Drive/Stats for Bioinformatics/HW1/")
# this provides an AffyBatch object
load("PSpikeAffyBatch.RData")

sampleNames(affydata)

# affy_deg=AffyRNAdeg(affydata)
# plotAffyRNAdeg(affy_deg)
# summaryAffyRNAdeg(affy_deg)

#read in the csv file containing information on 
spikeDF=read.table(file="AffyProbeSpikeValues.csv",sep="\t")
# levels(spikeDF[,2])
# summary(spikeDF[,2])

# let's clean up the object and replce MC and MF with NAs and get
# a real valued vector of FoldChanges..

SpikeFC=as.numeric(levels(spikeDF[,2])[spikeDF[,2]])
names(SpikeFC)=spikeDF[,1]

#removing zeros and na
indexNA = which(( SpikeFC != 0) & (!is.na( SpikeFC))) 

#expression values of above selected probes, class is factor
v =  SpikeFC[indexNA]

# probes that have the same concentration in both C and S
indexEC = which(v == 1)

# probes that have different concentrations between C and S
indexDC = which(v != 1)


routes <- read.table("route table.txt", header = FALSE, sep="/", 
                     strip.white = TRUE,na.strings = "EMPTY", stringsAsFactors = F)
colnames(routes) <- c("Route", "Background Correction", "Probe Normalization", "PM correction", "Summarization")
#write.table(routes, "./Routes.txt", quote = F, sep="\t", col.names = colnames(routes))

#for loop generates different expression values for different 
#routes specified in the "route table" or data frame "routes"
exp_vals <- list()
for(i in 1:nrow(routes)){
  temp = NULL
  temp <- expresso(affydata,bgcorrect.method = routes$`Background Correction`[i], 
                            normalize.method = routes$`Probe Normalization`[i],				
                            pmcorrect.method = routes$`PM correction`[i], 
                            summary.method = routes$Summarization[i])
  exp_vals[[i]] <- exprs(temp)[indexNA,]
}

save(exp_vals, file="./expression_values.RData")

load("./expression_values.RData")

# setting the probes with same conc. to 0
# vs different to 1				
ln = dim(exp_vals[[1]])[1]
myresponse = rep(NA, ln)
myresponse[indexEC] = 0	
myresponse[indexDC] = 1	
group = factor(c(rep("A", 9), rep("B", 9)))	

myfunc <- function(x){
  test <- mt.maxT(x, group, B = 10000)		
  #vector of test statistics in the original data order
  adjpv <- test$adjp[order(test$index)]
  roc(response = myresponse, predictor = adjpv)
}

rocc <- lapply(exp_vals, myfunc)
save(rocc, file="./ROC_results.RData")

## plotting the curves
quartz()
legend_list <- NULL
mycol <- rainbow(8)
plot(rocc[[1]], main= "ROC for Different Routes", col=mycol[1])
for(i in 1:8){
  plot(rocc[[i]], add=TRUE, col=mycol[i])
  legend_list[i] <- paste("Route", i, "- AUC:", 
                          round(as.numeric(rocc[[i]]$auc),2),sep=" ")
}

legend("bottomright", legend_list, lty=rep(1, 8), 				
         lwd=rep(4, 8) , col=mycol) 				
  

###### check how good the methods are #####

par(mfrow=c(2,2), mar=c(2,2,2,2))
for(i in 1:8){boxplot(exp_vals[[i]], main=paste("Route", i, sep=" "))}
# check for minimums
mins <- lapply(exp_vals, min)
save(mins, file="./mins")
# check for NAs
foo <- function(x){sum(1*is.na(x))}
nas<-lapply(exp_vals, foo) #Route 6 has more than 5000 NAs
save(nas, file="./nas")

### Clustering - Dendrogram ###

# vector of colors
labelColors = c("#CDB380", "#036564", "#EB6841", "#EDC951", 
                "cornflowerblue", "firebrick")
par(mfrow=c(4,2), mar=c(3,2,3,2))

for(i in 1:8){
d <- dist(t(exp_vals[[i]]))   # find distance matrix uses euclidean for dist calc
hc <- hclust(d)                # apply hirarchical clustering 
hcd<-as.dendrogram(hc)

# cut dendrogram in 6 clusters
clusMember = cutree(hc, 6)

# function to get color labels
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}

# using dendrapply
clusDendro = dendrapply(hcd, colLab)

# make plot
plot(clusDendro, main = paste("Route", i, sep=" "), type = "triangle",
     cex = 1.1, label.offset = 1)

}

### testing the clustering with same concentrations for all groups


for(i in 1:8){
  d <- dist(t(exp_vals[[i]][indexEC,]))   # find distance matrix uses euclidean for dist calc
  hc <- hclust(d)                # apply hirarchical clustering 
  hcd<-as.dendrogram(hc)
  
  # cut dendrogram in 6 clusters
  clusMember = cutree(hc, 6)
  
  # function to get color labels
  colLab <- function(n) {
    if (is.leaf(n)) {
      a <- attributes(n)
      labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
      attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    n
  }
  
  # using dendrapply
  clusDendro = dendrapply(hcd, colLab)
  
  # make plot
  plot(clusDendro, main = paste("Route", i, sep=" "), type = "triangle",
       cex = 1.1, label.offset = 1)
  
}

### PCA ###

source("http://www.stat.usu.edu/jrstevens/stat5570/pc2.R")
mfun <- function(x){
  a <- t(na.omit(x))
  pc<-pc2(a, scale=T)
  pc
}
pc <- lapply(exp_vals, mfun)
par(mfrow=c(2,4), mar=c(2,2,2,2))
for(i in 1:8){plot(pc[[i]], main="Screeplot")}
labelColors = c("#CDB380", "#036564", "#EB6841", "#EDC951", 
                "cornflowerblue", "firebrick")
mycols = c(rep(labelColors[1], 3),
          rep(labelColors[2], 3),
          rep(labelColors[3], 3),
          rep(labelColors[4], 3),
          rep(labelColors[5], 3),
          rep(labelColors[6], 3))

par(mfrow=c(2,5), mar=c(2,2,2,2))
for(i in 1:8){
plot(pc[[i]]$scores[,1],pc[[i]]$scores[,2], pch=19, col=mycols, 
     xlab="Component 1", ylab="Component 2", cex=1, main=paste("Route", i, sep=" "))
}

# legend("center", legend=rownames(pc$scores), col=mycols, pch=19, cex=0.6)

#homemade legend
x.p <- rep( c(seq(0.5, 7, by=2.8)) , 6)
y <- c(rep(6,3), rep(5,3), rep(4,3), rep(3,3), rep(2,3), rep(1,3) )
x.t <- rep(c(seq(1.8, 9, by=2.8)), 6)
plot(x.p, y, axes=FALSE, frame.plot=TRUE, col=mycols, 
     xlab="", ylab="", main = "Legend", xlim=c(0,8), ylim=c(0.8,6.2), pch=19, cex=1.8)
text(x.t, y, labels=rownames(pc[[1]]$scores))

# check for same conc ####
EC_exp <- lapply(exp_vals,"[",indexEC,,drop=FALSE)
pc <- lapply(EC_exp, mfun)
par(mfrow=c(2,4), mar=c(2,2,2,2))
for(i in 1:8){plot(pc[[i]], main="Screeplot")}
labelColors = c("#CDB380", "#036564", "#EB6841", "#EDC951", 
                "cornflowerblue", "firebrick")
mycols = c(rep(labelColors[1], 3),
           rep(labelColors[2], 3),
           rep(labelColors[3], 3),
           rep(labelColors[4], 3),
           rep(labelColors[5], 3),
           rep(labelColors[6], 3))

par(mfrow=c(2,5), mar=c(2,2,2,2))
for(i in 1:8){
  plot(pc[[i]]$scores[,1],pc[[i]]$scores[,2], pch=19, col=mycols, 
       xlab="Component 1", ylab="Component 2", cex=1, main=paste("Route", i, sep=" "))
}

# legend("center", legend=rownames(pc$scores), col=mycols, pch=19, cex=0.6)

#homemade legend
x.p <- rep( c(seq(0.5, 7, by=2.8)) , 6)
y <- c(rep(6,3), rep(5,3), rep(4,3), rep(3,3), rep(2,3), rep(1,3) )
x.t <- rep(c(seq(1.8, 9, by=2.8)), 6)
plot(x.p, y, axes=FALSE, frame.plot=TRUE, col=mycols, 
     xlab="", ylab="", main = "Legend", xlim=c(0,8), ylim=c(0.8,6.2), pch=19, cex=1.8)
text(x.t, y, labels=rownames(pc[[1]]$scores))