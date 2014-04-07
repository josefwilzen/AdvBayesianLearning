rm(list=ls())
graphics.off()
source("Gaussian processes/GPfunc.R")

# data:
wages<-read.delim(file="Gaussian processes/CanadianWages2.csv",sep=" ")
# standardize data:
wages$age<-scale(wages$age)
plot(wages$age,wages$logWage,pch=20,ylim=c(10,15.5),col="blue")
# calculate posterior

# prior:
# let l=0.2 and sigma=1 (change those later)
# assume sigmaError to be sd of data
priorSigma<-c(0.5,1.5)
priorLengthScale<-c(0.01,0.5,5)

sigmaError<-sd(wages$logWage)
xGrid<-seq(-2,3,length=300)
draws<-sampleGP(noDraw=30,xVal=xGrid,covFunc=SqrExpCov,func=meanWages,sigma=priorSigma[1],lengthScale=priorLengthScale[1])
# prior mean function:
plotGP(draws$draws)
curve(expr=meanWages,from=min(xGrid),to=max(xGrid),ylim=c(10,15),add=TRUE,lwd=4,col="red")


#plot(wages$age,wages$logWage,pch=20,ylim=c(10,17),,xlim=c(-1.7,2.3),col="blue",xlab="age",ylab="log(wages)")
#curve(expr=meanWages,from=-1.7,to=2.3,ylim=c(10,15),add=TRUE,lwd=3)

modelNo<-1
par(mfrow=c(3,2))
for(i in 1:length(priorLengthScale)){  #over lengthscale
  for(j in 1:length(priorSigma)){  # over sigma
    assign(x=paste("model",modelNo,sep=""),value=posteriorDist(xGrid=xGrid,
        priorMean=meanWages,obsData=wages,covFunc=SqrExpCov,sigmaError=sigmaError,
        sigma=priorSigma[j],lengthScale=priorLengthScale[i]))
    plotTheoreticalProbBand(postDist=get(x=paste("model",modelNo,sep="")),plotLim=c(8,17),predictiveMean=TRUE,plotData=TRUE,
            main=paste("model ",modelNo,":  ","sigma=",priorSigma[j]," l=",priorLengthScale[i] ,sep=""),xlab="Std. age",ylab="log(wages)")
    modelNo<-modelNo+1
  }    
} 
par(mfrow=c(1,1))



marginalLikelihood(obsData=wages,meanFunc=meanWages,covFunc=SqrExpCov,sigmaError=sigmaError,sigma=2,lengthScale=0.01)

evaluateMarginalLikelihood<-function(sigmaVect){
  
  
}

sigmaVect<-seq(0.01,20,by=0.01)
marginalLikelihood(obsData=wages,meanFunc=meanWages,covFunc=SqrExpCov,sigmaError=sigmaError,sigma=0.5,lengthScale=0.5)

temp<-list(obsData=wages,meanFunc=meanWages,covFunc=SqrExpCov,sigmaError=sigmaError,sigma=0.5,lengthScale=0.5)
do.call(what=marginalLikelihood,args=temp)


names(unlist(formals(SqrExpCov)[-(1:2)]))
a<-"sigma"
b<-10
list((as.name(a)))
list(as.character(a)=b)
optim()

optimMarginalLikelihood<-function(para,covFunc,...){
  paraName<-names(unlist(formals(covFunc)[-(1:2)]))
  otherPara<-as.list(match.call())[-(1:2)]
  paraList<-vector("list",length(paraName))
  for(i in 1:length(paraName)){
    paraList[[i]]<-para[i]
    names(paraList)[i]<-paraName[i]
  }
  allPara<-c(otherPara,paraList)
  print(allPara)
  res<-do.call(marginalLikelihood,allPara)
  return(res)
}

optimMarginalLikelihood(para=c(1,0.01),obsData=wages,meanFunc=meanWages,covFunc=SqrExpCov,sigmaError=sigmaError,sigma=0.5,lengthScale=0.5)





blurFunc<-function(...){
  return(match.call())
}

blurFunc<-function(...){
  return(as.list(environment()))
}
tempf <- function(a, b = 2, ...) {
  argg <- c(as.list(environment()), list(...))
  print(argg)
}


tempf <- function(a, b = 2, ...) {
  argg <- as.list(args(tmpfun))
  print(argg)
}
tempf(a=1, b=2, d=4,c=3)

tmpfun <- function(a,b,...) {
  hej<-as.list(match.call())[-(1:2)]
  print(hej)
  print("--------------")
  h<-hej[-1]
  print(h)
  return(hej)
  #print(as.list(match.call(expand.dots=FALSE)))
}
a<-tmpfun(a=1, b=2, d=4,c=3)
class(a)
str(a)
a[-1]
tempf(1, c = 3)
c(ex1,ex2)
ex1<-list(a=1,b=3,c=34)
ex2<-list(g="lkd",h=1:4)

list(unlist(ex1),unlist(ex2))
d<-blurFunc(a=1,b=2,c=3)
append
match.call()
as.list(environment())
as.list(args(tmpfun))
dput
Reduce


