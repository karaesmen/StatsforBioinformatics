#Insert code from problem one here. Need the gse19439 data
#Problem 2
Problem 2 Using student 1 as an example:
  Info = pData(phenoData(gse19439))[,1:2]

stu1 <- X[6874,];stu2 <- X[10685,];stu3 <- X[26058,]
Info[,1] <- FTB;Info[,3] <- stu1;Info[,4] <- stu2;Info[,5] <- stu3

stu1_CON <- Info$V3[which(Info$title=="CON")]

stu1_LTB <- Info$V3[which(Info$title=="LTB")]

stu1_PTB <- Info$V3[which(Info$title=="PTB")]

par(mfrow = c(1,2))

boxplot(stu1_PTB,stu1_LTB,stu1_CON)

vioplot(stu1_PTB,stu1_LTB,stu1_CON)

#Problem 3 
library(gplots)
best20 =order(myPvals)[1:20]
heatmap = heatmap.2(X[best20,])
library(NMF)
Info = pData(phenoData(gse19439))[,1:2]
Info[,1] <- FTB
group = Info[1]
heatmap2 = aheatmap(X[best20,], annCol=group, border=True)

#Below here code from 1-3 is not necesary

rm(list=ls())

#This provides an AffyBatch object
spikeDF = read.table(file = "AffyProbeSpikeValues.csv", sep = "\t")				
levels(spikeDF[, 2])				
summary(spikeDF[, 2])		

#Converting factors to numeric 
SpikeFC = as.numeric(levels(spikeDF[, 2])[spikeDF[, 2]])				
names(SpikeFC) = spikeDF[, 1]	

#removing zeros and na
data2 = which(( spikeDF[, 2]	 != 0) & (!is.na( spikeDF[, 2]	)))
v =  spikeDF[, 2]	[data2]				
v2 = SpikeFC[data2]				
data3 = which(v == 1)				
data4 = which(v != 1)	

#loading the data
load("PSpikeAffyBatch.RData")
require(affydata)		

#ROUTE1:
#Analysis route with background correction method as rma,
#probe normalization method as constant, PM correction method as pmonly and 
#summarization method as mas
route1 <- expresso(affydata, bgcorrect.method = "rma", normalize.method = "constant",				
                   pmcorrect.method = "pmonly", summary.method = "mas")	

#matrix of expression values
rt1 = exprs(route1)[data2, ]				
ln = dim(rt1)[1]
myresponse = rep(NA, ln)
myresponse[data3] = 0	
myresponse[data4] = 1	
source("http://bioconductor.org/biocLite.R")
biocLite("multtest")				
require(multtest)				
group = factor(c(rep("A", 9), rep("B", 9)))	

#compute permutation adjusted p-values for step-down multiple testing procedures
resT1 <- mt.maxT(rt1, group, B = 10000)		
#vector of test statistics in the original data order
teststat1 <- resT1$teststat[order(resT1$index)]				
require(pROC)				
roc1 = roc(response = myresponse, predictor = abs(teststat1))				
x11()				
plot(roc1, main = "ROC for route1")


#ROUTE2
#Analysis route with background correction method as mas,
#probe normalization method as quantiles, PM correction method as pmonly and 
#summarization method as avgdiff
route2 <- expresso(affydata, bgcorrect.method = "mas", normalize.method = "quantiles",				
                   pmcorrect.method = "pmonly", summary.method = "avgdiff")				
rt2 = exprs(route2)[data2, ]				
ln2 = dim(rt2)[1]	
group = factor(c(rep("A", 9), rep("B", 9)))				
resT2 <- mt.maxT(rt2, group, B = 10000)	
teststat2 <- resT2$teststat[order(resT2$index)]				
myresponse = rep(NA, ln)	
myresponse[data3] = 0		
myresponse[data4] = 1
roc2 = roc(response = myresponse, predictor = abs(teststat2))				
x11()				
#plotting roc curves				
plot(roc2, add=TRUE, colorize=TRUE, main = "ROC for route2")	


#ROUTE3
#Analysis route with background correction method as rma,
#probe normalization method as invariantset, PM correction method as pmonly and 
#summarization method as mas
route3 <- expresso(affydata, bgcorrect.method = "rma", normalize.method = "invariantset",				
                   pmcorrect.method = "pmonly", summary.method = "mas")				
rt3 = exprs(route3)[data2, ]				
ln3 = dim(rt3)[1]	
group = factor(c(rep("A", 9), rep("B", 9)))				
resT3 <- mt.maxT(rt3, group, B = 10000)	
teststat3 <- resT3$teststat[order(resT3$index)]				
myresponse = rep(NA, ln)	
myresponse[data3] = 0		
myresponse[data4] = 1
roc3 = roc(response = myresponse, predictor = abs(teststat3))				
x11()				
#plotting roc curves				
plot(roc3, add=TRUE, colorize=TRUE, main = "ROC for route3")	



#ROUTE4
#Analysis route with background correction method as rma,
#probe normalization method as loess, PM correction method as mas and 
#summarization method as avgdiff
route4 <- expresso(affydata, bgcorrect.method = "rma", normalize.method = "loess",				
                   pmcorrect.method = "mas", summary.method = "avgdiff")				
rt4 = exprs(route4)[data2, ]				
ln4 = dim(rt4)[1]	
group = factor(c(rep("A", 9), rep("B", 9)))				
resT4 <- mt.maxT(rt4, group, B = 10000)	
teststat4 <- resT4$teststat[order(resT4$index)]				
myresponse = rep(NA, ln)	
myresponse[data3] = 0		
myresponse[data4] = 1
roc4 = roc(response = myresponse, predictor = abs(teststat4))				
x11()				
#plotting roc curves				
plot(roc4, add=TRUE, colorize=TRUE, main = "ROC for route4")	



#ROUTE5
#Analysis route with background correction method as mas,
#probe normalization method as constant, PM correction method as subtractmm and 
#summarization method as mas
route5 <- expresso(affydata, bgcorrect.method = "mas", normalize.method = "constant",				
                   pmcorrect.method = "subtractmm", summary.method = "mas")				
rt5 = exprs(route5)[data2, ]				
ln5 = dim(rt5)[1]	
group = factor(c(rep("A", 9), rep("B", 9)))				
resT5 <- mt.maxT(rt5, group, B = 10000)	
teststat5 <- resT5$teststat[order(resT5$index)]				
myresponse = rep(NA, ln)	
myresponse[data3] = 0		
myresponse[data4] = 1
roc5 = roc(response = myresponse, predictor = abs(teststat5))				
x11()				
#plotting roc curves				
plot(roc5, add=TRUE, colorize=TRUE, main = "ROC for route5")	




#ROUTE6
#Analysis route with background correction method as rma,
#probe normalization method as quantiles, PM correction method as mas and 
#summarization method as medianpolish
route6 <- expresso(affydata, bgcorrect.method = "rma", normalize.method = "quantiles",				
                   pmcorrect.method = "mas", summary.method = "medianpolish")				
rt6 = exprs(route6)[data2, ]				
ln6 = dim(rt6)[1]	
group = factor(c(rep("A", 9), rep("B", 9)))				
resT6 <- mt.maxT(rt6, group, B = 10000)	
teststat6 <- resT6$teststat[order(resT6$index)]				
myresponse = rep(NA, ln)	
myresponse[data3] = 0		
myresponse[data4] = 1
roc6 = roc(response = myresponse, predictor = abs(teststat6))				
x11()				
#plotting roc curves				
plot(roc6, add=TRUE, colorize=TRUE, main = "ROC for route6")	



#ROUTE7
#Analysis route with background correction method as mas,
#probe normalization method as invariantset, PM correction method as subtractmm and 
#summarization method as medianpolish
route7 <- expresso(affydata, bgcorrect.method = "mas", normalize.method = "invariantset",				
                   pmcorrect.method = "subtractmm", summary.method = "medianpolish")				
rt7 = exprs(route7)[data2, ]				
ln7 = dim(rt7)[1]	
group = factor(c(rep("A", 9), rep("B", 9)))				
resT7 <- mt.maxT(rt7, group, B = 10000)	
teststat7 <- resT7$teststat[order(resT7$index)]				
myresponse = rep(NA, ln)	
myresponse[data3] = 0		
myresponse[data4] = 1
roc7 = roc(response = myresponse, predictor = abs(teststat7))				
x11()				
#plotting roc curves				
plot(roc7, add=TRUE, colorize=TRUE, main = "ROC for route7")



#ROUTE8
#Analysis route with background correction method as rma,
#probe normalization method as qspline, PM correction method as subtractmm and 
#summarization method as avgdiff
route8 <- expresso(affydata, bgcorrect.method = "rma", normalize.method = "qspline",				
                   pmcorrect.method = "subtractmm", summary.method = "avgdiff")				
rt8 = exprs(route8)[data2, ]				
ln8 = dim(rt8)[1]	
group = factor(c(rep("A", 9), rep("B", 9)))				
resT8 <- mt.maxT(rt8, group, B = 10000)	
teststat8 <- resT8$teststat[order(resT8$index)]				
myresponse = rep(NA, ln)	
myresponse[data3] = 0		
myresponse[data4] = 1
roc8 = roc(response = myresponse, predictor = abs(teststat8))				
x11()				
#plotting roc curves				
plot(roc8, add=TRUE, colorize=TRUE, main = "ROC for route8")


#plotting roc curves
X11()
plot(roc1, main = "ROC for routes", col='BLUE')				
plot(roc2, add=TRUE, col='GREEN' )				
plot(roc3, add=TRUE, col='RED')				
plot(roc4, add=TRUE, col='YELLOW')				
plot(roc5, add=TRUE, col='darkorange')				
plot(roc6, add=TRUE, col='darkorchid')				
plot(roc7, add=TRUE, col='cyan')				
plot(roc8, add=TRUE, col='magenta')				
legend("bottomright", c("Route 1","Route 2","Route 3","Route 4","Route 5",
                        "Route 6","Route 7","Route 8"), lty=c(1,1,1,1,1,1,1,1), 				
       lwd=c(1,1,1,1,1,1,1,1),col=c("BLUE","GREEN","RED","YELLOW",
                                    "DARKORANGE","DARKORCHID","CYAN","MAGENTA")) 				

