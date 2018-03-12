dvals=seq(0,3.5,length=101)
pwr=rep(0,length(dvals))
for(i in 1:(length(dvals))) pwr[i]=pwr.t.test(n=10,d=dvals[i],sig.level=0.05)$power
plot(dvals,pwr,type="l",xlab=expression(delta),ylab="power",lwd=2,ylim=c(0,1))
abline(h=0.05,lty=2)
pwrBon=rep(0,length(dvals))
for(i in 1:(length(dvals))) pwrBon[i]=pwr.t.test(n=10,d=dvals[i],sig.level=0.05/200)$power






SimDatHW2=function(p.alt=10, # the number of differentially expressed features ... 
                   p.null=90, # the number of non-differentially expressed features 
                   n=20, # the number of samples in each of the treatment and control groups
                   rho.alt=0.2, # correlation of the alt hyp variables
                   rho.null=0.1, # correlation of the null hyp variables
                   delta=2,# the mean of the p.alt features in the "treatment" group
                   sdC=0 # the variance of the sample specific centering error
){
  p=p.alt+p.null
  Sigma=array(rep(0,p^2),dim=c(p,p))
  Sigma[1:p.alt,1:p.alt]=rho.alt
  Sigma[(p.alt+(1:p.null)),(p.alt+(1:p.null))]=rho.null
  diag(Sigma)=1
  Xc=mvrnorm(n,mu=rep(0,p),Sigma=Sigma)
  Xt=mvrnorm(n,mu=c(rep(delta,p.alt),rep(0,p.null)),Sigma=Sigma)
  x=t(rbind(Xc,Xt))
  colnames(x)=c(paste("cntrl",1:n,sep=""),paste("trt",1:n,sep=""))
  rownames(x)=c(paste("DEgene",1:(p.alt)),paste("nonDEgene",1:p.null,sep=""))
  if(sdC>0){
    CentError=rnorm(2*n,sd=sdC)
    for(j in 1:(2*n)) x[,j]=x[,j]+CentError[j]
  }
  return(x)
}

library(MASS)
library(multtest)
nreps <- 1000
f <- factor(c(rep("cntrl", 10), rep("trt", 10)))
PhiVec <- rep(NA,nreps)
for(k in 1:nreps){
  X=SimDatHW2(p.alt=1,p.null=199,n=10,rho.alt=0.0,rho.null=0.75,delta=2)
  # calculate minP adjusted p-values
  resT <- mt.maxT(X, f, B = 2500)
  PhiVec[k]=as.integer(resT$adjp[resT$index==1]<=0.05)
}

plot(dvals,pwr,type="l",xlab=expression(delta),ylab="power",ylim=c(0,1))
abline(h=0.05,lty=2)
lines(dvals,pwrBon,col=4,lwd=2)
abline(h=0.05/200,lty=3,col=4)
# add our simulation point
points(2,mean(PhiVec),pch="+",col="darkgreen", cex=2)
