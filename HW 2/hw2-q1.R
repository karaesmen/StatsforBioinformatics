
rm(list=ls())
setwd("/Users/aarizvi/Dropbox/Group2/HW2/ARizvi/")
library(GEOquery)
## Berry UK Training (array)
gse19439 <- getGEO("GSE19439",GSEMatrix=T)

pdat19439 <- pData(phenoData(gse19439[[1]]))[,c(1,11:16)]
pdat19 <- pdat19439
## altering 'factor' columns so they are more interpretable
colnames(pdat19439) <- c("group", "gender", "ethnicity", "illness", "geographical_region", "bcg_vaccinated", "region_of_birth")
rename <- function(x){
        temp <- strsplit(as.character(x), ": ")
        mat <- matrix(unlist(temp), ncol=2, byrow=TRUE)        
}
pdat19439$group <- as.factor(substr(pdat19439$group, 1, 3))
for(i in 2:ncol(pdat19439)){pdat19439[,i] <- as.factor(rename(pdat19439[,i])[,2])}
#save(file="pdat19439.RData", pdat19439)
#save(file="gse19439.RData", gse19439)
load("gse19439.RData")
load("pdat19439.RData")

## grab assay data
assay.data <- exprs(gse19439[[1]])
dim(assay.data)

## part c
pc <- prcomp(assay.data, center = TRUE, scale = TRUE)
summary(pc)
plot(pc, type="lines")

titles <- colnames(pdat19439)

factorplot <- function(prcomp.object, group.factor, titles){
        colpal <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33')
        colpal <- colpal[1:length(levels(group.factor))]
        plot(prcomp.object$rotation, col=colpal[group.factor], pch=19)
        legend("topleft", pch=19, levels(group.factor), col=colpal)
        title(titles)
}

par(mfrow=c(2,4))
for(i in 1:ncol(pdat19439)){
        factorplot(pc, pdat19439[,i], titles[i])
}

## part d
feature.subset <- function(data, cutoff){
        MyDip <- function(i) {dip.test(x=data[i,])$p.value}
        I <- dim(data)[1] # number of features
        mfilt <- as.vector(mapply(MyDip,1:I))
        data[which(mfilt<sort(mfilt)[cutoff+1]),]
}

our.feature <- feature.subset(assay.data, 50)
# So, we have a histogram of p-values of the dip test. 
# which ones should we keep?
# which is to say, we will keep all the probesets with a corresponding
# p-value below some cut-off that we specify.
# But what cut-off?
# Choices:
# 1. We could decide to keep the best M  (i.e., the lowest M p-values)
# 2. We could control some kind of error rate, like FDR.
# we will keep best M p-values
feature.selection.pca <- function(feature.subset, phenolabel){
        fit <- prcomp(t(feature.subset),center=T, scale=T)
        pcx1 <- t(feature.subset)%*%fit$rotation[,1]
        pcx2 <- t(feature.subset)%*%fit$rotation[,2]
        plot(pcx1,pcx2,pch=" ",main=paste("Filtered Dataset (p = ", cutoff, ")"))
        text(pcx1,pcx2,label=phenolabel,col=as.integer(phenolabel))
}


#example of how you would plot a feature group of interest
#first argument is any column ("Feature") from pdat19439 (e.g. pdat19439$gender)
#second input argument is any level of that column (e.g. levels(pdat19439$gender)[1])
feature.selection <- function(feature, subgroup){assay.data[,colnames(assay.data) %in% rownames(pdat19439[feature == subgroup,])]}

feature.selection.plot <- function(feature, subgroup){
        select <- assay.data[,colnames(assay.data) %in% rownames(pdat19439[feature == subgroup,])]
        pca <- prcomp(select, center = TRUE, scale = TRUE)
        plot(pca$rotation, pch=19)
        title(paste("PCA of", deparse(substitute(feature)),  'of', subgroup, "patients only", sep=" "))
}

female.assay <- feature.selection(pdat19439$gender, "Female")
dim(female.assay) #23 females
female.assay[1:5,1:5] #assay data
feature.selection.plot(pdat19439$gender, "Female")
