rm(list=ls())
setwd("~/Google Drive/Stats for Bioinformatics/HW3/HW3_Q3/")
load("~/Google Drive/Stats for Bioinformatics/HW3/Run1/Array.RData") #ArrayData
load("~/Google Drive/Stats for Bioinformatics/HW3/GeneMap.RData")    #GeneMap
# load("~/Google Drive/Stats for Bioinformatics/HW3/Run1/Valid.RData") #ValidData

# Our ArrayData consists of 20 arrays, 10 for Group A and 10 for Group B 
# Rownames of the ArrayData are the feature names

ArrayData[1:3,c(1,20,21,40)]
n.gr <- length(unique(GeneMap)) # We have 50 gene sets
hist(GeneMap) # Each set have different number of features
              # But each feature only belongs to one set

#### 3b #####
data= ArrayData; dist.m="euclidean"
clus.m="complete"; K=3

ftest <-function(data= ArrayData, dist.m="euclidean", clus.m="complete", K=3){

d <- dist(data, method=dist.m)   # generate distance matrix
hc <- hclust(d, method = clus.m) # apply hirarchical clustering

# cut dendrogram in K clusters
clusMember = cutree(hc, K)

# put all group information in a data frame
feat.groups <- data.frame(Features=rownames(ArrayData), 
                          Sets=GeneMap,
                          Clusters = clusMember, row.names=1)
group.table <- t(table(Clusters=feat.groups$Clusters,  Sets=feat.groups$Sets))

ftest <- fisher.test(group.table)
}

ftest()

# my.dist <- c("euclidean", "manhattan")
# my.clust <- c("single", "complete", "average", "centroid")
# dis = NULL
# clus = NULL
# pval= NULL
# pvtab <- data.frame(distance=dis, clustering=clus, pval=NULL)
# for(i in 1:length(my.dist)){
#   dis <- my.dist[i]
#   for(j in 1:length(my.clust)){
#     clus <- my.clust[j]
#     pval <- ftest(dis, clus)
#     pvtab[i,] <- data.frame(distance=dis, clustering=clus, pvalue=pval)
#   }
# }



mo.plot <- function(data=ArrayData, dist.m="euclidean", clus.m="complete", K=3){
  # dist.m: "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
  # clus.m: "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), 
  #         "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
  # K: cut dendrogram in K clusters, K must be between 1-6, and integer
  
  d <- dist(data, method=dist.m)   # generate distance matrix
  hc <- hclust(d, method = clus.m) # apply hirarchical clustering 
  hcd<-as.dendrogram(hc)
  
  # cut dendrogram in K clusters
  clusMember = cutree(hc, K)
  
  # put all group information in a data frame
  feat.groups <- data.frame(Features=rownames(ArrayData), 
                            Sets=GeneMap,
                            Clusters = clusMember, row.names=1)
  
  
  # I created 60 rows of dummy features to have the first set 0 to show 
  # number of clusters properly
  a = 60
  labl <- data.frame(Features=paste("dummy", seq(1:a), sep=""),
                     Sets = rep(0, a),
                     Clusters = c(rep.int(1:K, a/K)), row.names=1)
             
  feat.groups <- rbind(feat.groups,labl)
  
  group.table <- table( Sets=feat.groups$Sets,Clusters=feat.groups$Clusters)
  mosaicplot(group.table, cex.axis=0.6,
             main=paste(dist.m, "distance ,", clus.m, "linkage",  "for K =", K), 
             color=c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f'))
}

j=1
pdf(file=paste("./mosaic", j ,".pdf",  sep="" ))
par(mfrow=c(4, 1), mar=c(2,2,2,2), cex.main=0.8)
for(i in 3:6) mo.plot(K=i)
dev.off()

j=2
pdf(file=paste("./mosaic", j ,".pdf",  sep="" ))
par(mfrow=c(4, 1), mar=c(2,2,2,2), cex.main=0.8)
for(i in 3:6) mo.plot(K=i, dist.m="manhattan", clus.m="centroid")
dev.off()

j=3
pdf(file=paste("./mosaic", j ,".pdf",  sep="" ))
par(mfrow=c(4, 1), mar=c(2,2,2,2), cex.main=0.8)
for(i in 3:6) mo.plot(dist.m = "manhattan", clus.m = "complete" , K =i)
dev.off()

for(i in 3:6) mo.plot(K=i, dist.m="manhattan", clus.m="centroid")
for(i in 2:5) mo.plot(dist.m = "manhattan", clus.m = "complete" , K =i)

#### 3c #####

# Extract only non-NULL elements from ValidData list
valD <- ValidData[vapply(ValidData, Negate(is.null), NA)]

# Extract all elements of the list as a matrix
# now we have a matrix with 130 A & 130 B samples, and 10 features
valD <- do.call("rbind", valD)
valSet <- GeneMap[which(rownames(ArrayData) %in% rownames(valD))]

# Create a pie chart with actual sets
pie(c(5,rep(1,5)), labels=paste("Set", unique(valSet[order(valSet)])),
    col=c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f'),
    radius = 1, cex=2)

