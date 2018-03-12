rm(list=ls())

## Change this to your directory!
load("~/Google Drive/Stats for Bioinformatics/HW3/Run1/Array.RData") #ArrayData
load("~/Google Drive/Stats for Bioinformatics/HW3/GeneMap.RData")    #GeneMap
load("~/Google Drive/Stats for Bioinformatics/HW3/Run1/Valid.RData") #ValidData

# GeneMap is a numeric vector, specifying the set each feature belongs to
n.gr <- length(unique(GeneMap)) # We have 50 gene sets
hist(GeneMap) # Each set have different number of features
              # But each feature only belongs to one set

#### 3b #####

## mo.plot is a function that creates mosaic plots for cluster groups and actual gene sets.
## mo.plot accepts different distance or clustering methods as defined in the function as comments.
## it also accepts number of clusters to be created, K must be between 1-6.
# Function also also creates 15 rows of dummy features to have the first set 0 to show 
# number of clusters properly on the mosaic plot


mo.plot <- function(data=ArrayData, dist.m="euclidean", clus.m="complete", K=3){
  # data input is the ArrayData (a data.frame) your group received during the biomarker challenge.
  # choose distance metric - dist.m: "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
  # choose clustering method - clus.m: "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), 
  #                                   "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
  # number of clusters - K: cut dendrogram in K clusters, K must be between 1-6, and integer
  
  d <- dist(data, method=dist.m)   # generate distance matrix
  hc <- hclust(d, method = clus.m) # apply hirarchical clustering 
  hcd<-as.dendrogram(hc)
  
  # cut dendrogram in K clusters
  clusMember = cutree(hc, K)
  
  # put all group information in a data frame
  feat.groups <- data.frame(Features=rownames(ArrayData), 
                            Sets=GeneMap,
                            Clusters = clusMember, row.names=1)
  
  
  # I created 15 rows of dummy features to have the first set 0 to show 
  # number of clusters properly
  a = 60
  labl <- data.frame(Features=paste("dummy", seq(1:a), sep=""),
                     Sets = rep(0, a),
                     Clusters = c(rep.int(1:K, a/K)), row.names=1)
             
  feat.groups <- rbind(feat.groups,labl)
  
  group.table <- table( Sets=feat.groups$Sets,Clusters=feat.groups$Clusters)
  mosaicplot(group.table, cex.axis=0.6,
             main=paste(dist.m, "distance ,", clus.m, "linkage",  "for 50 Sets and K =", K), 
             color=c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f'))
}

par(mfrow=c(4, 1), mar=c(2,2,2,2))
for(i in 3:6){mo.plot(K=i)}

par(mfrow=c(4, 1), mar=c(2,2,2,2))
for(i in 3:6){mo.plot(dist.m = "euclidean", clus.m = "ward" , K =i)}

#### 3c #####
# ValidData is the list you got from the Biomarker Challange
# It is a list with the features you decided to validate during the Biomarker challenge
# All other features that weren't selected are NULL

# Extract only non-NULL elements from ValidData list
valD <- ValidData[vapply(ValidData, Negate(is.null), NA)]

# Extract all elements of the list as a matrix
# now we have a matrix with the number of samples your group decided to validate 
# and features you selected
valD <- do.call("rbind", valD)
valSet <- GeneMap[which(rownames(ArrayData) %in% rownames(valD))]

# Create a pie chart with actual sets from the MapGene
pie(c(5,rep(1,5)), labels=paste("Set", unique(valSet[order(valSet)])),
    col=c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f'),
    radius = 1, cex=2)

