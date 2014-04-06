
meanFunc<-function(xValues,func=NULL){
  if(is.null(func)) return(rep(0,length(xValues)))
  yValues<-func(xValues)
  return(yValues)
}

SqrExpCov<-function(x1,x2,sigma=1,lengthScale=1){
  if(lengthScale<=0) stop("Values of lengthScale not allowed")
  if(sigma<=0) stop("Values of sigma not allowed")
  covVal<-sigma^2*exp(((x1-x2)^2)/(-2*lengthScale))
  return(covVal)
}


createCovMatrix<-function(xVal,yVal=NULL,covFunc,...){
  index<-1:length(xVal)
#   indexMat<-as.matrix(expand.grid(index,index))
#   covValues<-covFunc(x1=indexMat[,1],x2=indexMat[,2],...)
#   covMat<-matrix(0,length(xVal),length(xVal))
#   covMat[indexMat]<-covValues

  if(is.null(yVal)){
    covMat<-matrix(0,length(xVal),length(xVal))
    for(i in 1:length(xVal)){
      for(j in 1:length(xVal)){
        covMat[i,j]<-covFunc(x1=xVal[i],x2=xVal[j],...)
      } 
    }
  }else {
    covMat<-matrix(0,nrow=length(yVal),ncol=length(xVal))
    for(i in 1:length(yVal)){  # rows 
      for(j in 1:length(xVal)){ # cols
        covMat[i,j]<-covFunc(x1=yVal[i],x2=xVal[j],...)
      } 
    }
    colnames(covMat)<-paste("xVal",round(xVal,3))
    rownames(covMat)<-paste("yVal",round(yVal,3))
  }
  return(covMat)
}

sampleGP<-function(noDraw,xVal,covFunc,func=NULL,onlyPara=FALSE,...){
  meanVal<-meanFunc(xVal,func=func)
  covMat<-createCovMatrix(xVal=xVal,covFunc=covFunc,...)
  if(onlyPara){
    return(list(covMat=covMat,means=meanVal))
  } 
  library(MASS)
  draws<-mvrnorm(n=noDraw,mu=meanVal,Sigma=covMat)
  colnames(draws)<-xVal
  rownames(draws)<-1:noDraw
  res<-list(draws=draws,covMat=covMat,means=meanVal)
  return(res)
}

plotGP<-function(drawMatrix){
  maxVal<-max(drawMatrix)
  minVal<-min(drawMatrix)
  noDraws<-dim(drawMatrix)[1]
  if(noDraws==1) {
    plot(x=as.numeric(colnames(drawMatrix)),y=drawMatrix[1,], type="l", xlab="t",ylab="y")
  }else{
    plot(x=as.numeric(colnames(drawMatrix)),y=drawMatrix[1,], type="l", 
         xlab="t",ylab="y",ylim=c(minVal,maxVal))
    for(i in 2:noDraws)
    lines(x=as.numeric(colnames(drawMatrix)),y=drawMatrix[i,])
  }
}

meanWages<-function(x){
  y<-14*tanh(x+2.8)
  return(y)
}



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
# let l=1 and sigma=1 (change those later)
# assume sigmaError to be sd of data
sigmaError<-sd(wages$logWage)

# prior mean function:
plot(wages$age,wages$logWage,pch=20,ylim=c(10,17),,xlim=c(-1.7,2.3),col="blue",xlab="age",ylab="log(wages)")
curve(expr=meanWages,from=-1.7,to=2.3,ylim=c(10,15),add=TRUE,lwd=3)

# posterior mean
# grid of xValues
xGrid<-seq(-1.5,2.2,length=30)

meanXgrid<-as.matrix(meanFunc(xValues=xGrid,func=meanWages))
xKxGrid<-createCovMatrix(xVal=wages$age,yVal=xGrid,covFunc=SqrExpCov,sigma=1)
sqrtKy<-sqrt(createCovMatrix(xVal=wages$age,covFunc=SqrExpCov,sigma=1)+sigmaError*diag(1,nrow=dim(wages)[1]))
yDiff<-wages$logWage-meanFunc(xValues=wages$age,func=meanWages)
#fBar<-xKxGrid%*%sqrtKy%*%wages$logWage
fBar<-meanXgrid+xKxGrid%*%sqrtKy%*%wages$logWage
#plot(x=wages$age,y=meanFunc(xValues=wages$age,func=meanWages),type="l")
#plot(x=xGrid,y=meanXgrid,type="l")
graphics.off()
plot(x=xGrid,y=fBar,type="o")




temp3<-createCovMatrix(xVal=xGrid,yVal=wages$age,covFunc=SqrExpCov)
temp2<-createCovMatrix(xVal=wages$age,yVal=xGrid,covFunc=SqrExpCov)
temp1<-createCovMatrix(xVal=wages$age,covFunc=SqrExpCov)
temp4<-createCovMatrix(xVal=xGrid,covFunc=SqrExpCov)
dim(temp2)
dim(temp3)
dim(temp1)
dim(temp4)







