rm(list=ls())
source("Gaussian processes/GPfunc.R")

#draw<-sampleGP(noDraw=20,xVal=seq(-6,6,length.out=100),covFunc=SqrExpCov,lengthScale=2,sigma=1,func=sin)
#draw<-sampleGP(noDraw=20,xVal=seq(-6,6,length.out=100),covFunc=SqrExpCov,lengthScale=2,sigma=1,func=sin,onlyPara=T)
#round(draw$covMat,4)
#image(draw$covMat)
#plotGP(drawMatrix=draw$draws)

# data:
wages<-read.delim(file="Gaussian processes/CanadianWages2.csv",sep=" ")
# standardize data:
wages$age<-scale(wages$age)
plot(wages$age,wages$logWage,pch=20,ylim=c(10,17),col="blue")
# calculate posterior

# prior:
# let l=0.2 and sigma=1 (change those later)
# assume sigmaError to be sd of data
priorSigma<-1
priorLengthScale<-0.2

sigmaError<-sd(wages$logWage)
xGrid<-seq(-1.5,2.2,length=300)
draws<-sampleGP(noDraw=30,xVal=xGrid,covFunc=SqrExpCov,func=meanWages,sigma=priorSigma,lengthScale=priorLengthScale)
# prior mean function:
plotGP(draws$draws)
curve(expr=meanWages,from=-1.7,to=2.3,ylim=c(10,15),add=TRUE,lwd=4,col="red")


#plot(wages$age,wages$logWage,pch=20,ylim=c(10,17),,xlim=c(-1.7,2.3),col="blue",xlab="age",ylab="log(wages)")
#curve(expr=meanWages,from=-1.7,to=2.3,ylim=c(10,15),add=TRUE,lwd=3)


test<-posteriorDist(xGrid=xGrid,priorMean=meanWages,obsData=wages,covFunc=SqrExpCov,sigmaError=sigmaError,sigma=priorSigma,lengthScale=priorLengthScale)

plotTheoreticalProbBand(postDist=test,predictiveMean=FALSE,plotData=TRUE)
plotTheoreticalProbBand(postDist=test,predictiveMean=TRUE,plotData=TRUE)




library(MASS)
drawsPosterior<-mvrnorm(n=2000,mu=fBar,Sigma=posteriorCov)
colnames(drawsPosterior)<-round(xGrid,3)
# prior mean function:
plotGP(drawsPosterior)
lines(x=xGrid,y=fBar,col="red",lwd=3)
plotProbBandsGP(drawMatrix=drawsPosterior,postMean=fBar,xGrid=xGrid)

#plot(x=xGrid,y=fBar,type="l",main=as.character(now()))
conf<-apply(X=drawsPosterior,MARGIN=2,FUN=quantile,prob=c(0.05,0.95))








