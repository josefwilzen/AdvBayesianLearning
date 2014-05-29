rm(list=ls())
graphics.off()
setwd("/home/joswi05/Dropbox/Josef/Advanced Bayesian Learning/AdvBayesianLearning")
source("Gaussian processes/GPfunc.R")
# data:
wages<-read.delim(file="Gaussian processes/CanadianWages2.csv",sep=" ")
# standardize data:
wages$age<-scale(wages$age)
plot(wages$age,wages$logWage,pch=20,ylim=c(10,15.5),col="blue",xlab="age",ylab="log(wages)")
# calculate posterior

# prior:
# let l=0.2 and sigma=1 (change those later)
# assume sigmaError to be sd of data
priorSigma<-c(2,1.5)
priorLengthScale<-c(0.1,0.5,5)

sigmaError<-sd(wages$logWage)
xGrid<-seq(-2,3,length=300)
draws<-sampleGP(noDraw=30,xVal=xGrid,covFunc=SqrExpCov,func=meanWages,sigma=priorSigma[1],lengthScale=priorLengthScale[1])
# prior mean function:
plotGP(draws$draws,xlab="age",ylab="log(wages)")
curve(expr=meanWages,from=min(xGrid),to=max(xGrid),ylim=c(10,15),add=TRUE,lwd=4,col="red")


#plot(wages$age,wages$logWage,pch=20,ylim=c(10,17),,xlim=c(-1.7,2.3),col="blue",xlab="age",ylab="log(wages)")
#curve(expr=meanWages,from=-1.7,to=2.3,ylim=c(10,15),add=TRUE,lwd=3)

# modelNo<-1
# par(mfrow=c(3,2))
# for(i in 1:length(priorLengthScale)){  #over lengthscale
#   for(j in 1:length(priorSigma)){  # over sigma
#     assign(x=paste("model",modelNo,sep=""),value=posteriorDist(xGrid=xGrid,
#         priorMean=meanWages,obsData=wages,covFunc=SqrExpCov,sigmaError=sigmaError,
#         sigma=priorSigma[j],lengthScale=priorLengthScale[i]))
#     plotTheoreticalProbBand(postDist=get(x=paste("model",modelNo,sep="")),plotLim=c(8,17),predictiveMean=TRUE,plotData=TRUE,
#             main=paste("model ",modelNo,":  ","sigma=",priorSigma[j]," l=",priorLengthScale[i] ,sep=""),xlab="Std. age",ylab="log(wages)")
#     modelNo<-modelNo+1
#   }    
# } 
# par(mfrow=c(1,1))

#-------------------------------------------------------------
# Kernel 1:  Squared Exponential kernel
# optimering
#-------------------------------------------------------------
lowerVal<-0.001
startValues<-c(sigmaError,priorSigma[1],priorLengthScale[1])
kernelOptim1<-optim(par=startValues,fn=optimMarginalLikelihood,method="L-BFGS-B",
            obsData=wages,meanFunc=meanWages,covFunc=SqrExpCov,
              lower=c(lowerVal,lowerVal,lowerVal),control=list(fnscale=-1),hessian=TRUE)

SqrExpSigmaError<-kernelOptim1$par[1]
SqrExpSigma<-kernelOptim1$par[2]
SqrExpLengthScale<-kernelOptim1$par[3]
  

postSqrExp<-posteriorDist(xGrid=xGrid,
              priorMean=meanWages,obsData=wages,covFunc=SqrExpCov,sigmaError=SqrExpSigmaError,
              sigma=SqrExpSigma,lengthScale=SqrExpLengthScale)
plotTheoreticalProbBand(postDist=postSqrExp,plotLim=c(8,17),predictiveMean=TRUE,plotData=TRUE,
                        main=paste("sigma model=",round(SqrExpSigmaError,3),
                          "  sigma kernel=",round(SqrExpSigma,3),"  l=",round(SqrExpLengthScale,3) ,sep=" "),
                        xlab="Std. age",ylab="log(wages)")

#-------------------------------------------------------------
# Other kernels
# ej optimering
#-------------------------------------------------------------

# exp gamma:

gammaSigmaError<-kernelOptim1$par[1]
gammaSigma<-kernelOptim1$par[2]
gammaGamma<-c(0.5,1.5)
gammaLengthScale<-kernelOptim1$par[3]


for(i in 1:2){
  postGamma<-posteriorDist(xGrid=xGrid,
                           priorMean=meanWages,obsData=wages,covFunc=gammaExponential,sigmaError=gammaSigmaError,
                           sigma=gammaSigma,lengthScale=gammaLengthScale,gamma=gammaGamma[i])
  plotTheoreticalProbBand(postDist=postGamma,plotLim=c(8,17),predictiveMean=TRUE,plotData=TRUE,
                          main=paste("gamma=",round(gammaGamma[i],3) ,sep=" "),
                          xlab="Std. age",ylab="log(wages)")
}

# matern kernel:


SqrExpSigmaError<-kernelOptim1$par[1]
SqrExpSigma<-kernelOptim1$par[2]
SqrExpLengthScale<-kernelOptim1$par[3]



#-------------------------------------------------------------
# Kernel 2: Gamma exponential
# optimering
#-------------------------------------------------------------

formals(gammaExponential)
lowerVal<-0.01
gamma<-0.5
startValues2<-c(sigmaError,priorSigma[1],gamma,priorLengthScale[1])
startValues2<-c(kernelOptim1$par[1:2],gamma,kernelOptim1$par[3])
kernelOptim2<-optim(par=startValues2,fn=optimMarginalLikelihood,method="L-BFGS-B",
                 obsData=wages,meanFunc=meanWages,covFunc=gammaExponential,
                 lower=c(lowerVal,lowerVal,lowerVal,lowerVal),upper=c(Inf,Inf,1.99,Inf),control=list(fnscale=-1),hessian=TRUE)


kernelOptim2<-optim(par=c(1,1,1,1),fn=optimMarginalLikelihood,method="BFGS",
                    obsData=wages,meanFunc=meanWages,covFunc=maternKernel,control=list(fnscale=-1),hessian=TRUE)

gammaSigmaError<-kernelOptim2$par[1]
gammaSigma<-kernelOptim2$par[2]
gammaGamma<-kernelOptim2$par[3]
gammaLengthScale<-kernelOptim2$par[4]

postGamma<-posteriorDist(xGrid=xGrid,
                          priorMean=meanWages,obsData=wages,covFunc=gammaExponential,sigmaError=gammaSigmaError,
                          sigma=gammaSigma,lengthScale=gammaLengthScale,gamma=gammaGamma)
plotTheoreticalProbBand(postDist=postGamma,plotLim=c(8,17),predictiveMean=TRUE,plotData=TRUE,
                        main=paste("sigma model=",round(gammaSigmaError,3),
                                   "sigma kernel=",round(gammaSigma,3),"gamma=",round(gammaGamma,3),"l=",round(gammaLengthScale,3) ,sep=" "),
                        xlab="Std. age",ylab="log(wages)")

#-------------------------------------------------------------
# Kernel 3: maternKernel
# optimering
#-------------------------------------------------------------

formals(maternKernel)
lowerVal<-0.01
startValues3<-c(sigmaError,priorSigma[1],priorLengthScale[1])
startValues3<-c(1,1,2)  # funkar tillsammans med v=1.5 och v=2.5
startValues3<-c(0.5,0.5,1)
kernelOptim2<-optim(par=startValues3,fn=optimMarginalLikelihood,method="L-BFGS-B",
                    obsData=wages,meanFunc=meanWages,covFunc=maternKernel,lower=c(lowerVal,lowerVal,lowerVal,lowerVal),control=list(fnscale=-1),hessian=TRUE)

gammaSigmaError<-kernelOptim2$par[1]
gammaSigma<-kernelOptim2$par[2]
gammaLengthScale<-kernelOptim2$par[3]

postGamma<-posteriorDist(xGrid=xGrid,priorMean=meanWages,obsData=wages,covFunc=maternKernel,sigmaError=gammaSigmaError,
                         sigma=gammaSigma,lengthScale=gammaLengthScale)
plotTheoreticalProbBand(postDist=postGamma,plotLim=c(8,17),predictiveMean=TRUE,plotData=TRUE,
                        main=paste("sigma model=",round(gammaSigmaError,3),"sigma kernel=",round(gammaSigma,3)),"l=",round(gammaLengthScale,3) ,sep=" ",
                        xlab="Std. age",ylab="log(wages)",quan=c(0.025,0.975))

plotTheoreticalProbBand(postDist=postGamma,plotLim=c(8,17),predictiveMean=TRUE,plotData=TRUE)

str(postGamma)
#-------------------------------------------------------------



