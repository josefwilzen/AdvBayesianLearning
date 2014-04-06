
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


createCovMatrix<-function(xVal,covFunc,...){
  index<-1:length(xVal)
#   indexMat<-as.matrix(expand.grid(index,index))
#   covValues<-covFunc(x1=indexMat[,1],x2=indexMat[,2],...)
#   covMat<-matrix(0,length(xVal),length(xVal))
#   covMat[indexMat]<-covValues
  covMat<-matrix(0,length(xVal),length(xVal))
  for(i in 1:length(xVal)){
    for(j in 1:length(xVal)){
      covMat[i,j]<-covFunc(x1=xVal[i],x2=xVal[j],...)
    } 
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

draw<-sampleGP(noDraw=20,xVal=seq(-6,6,length.out=100),covFunc=SqrExpCov,lengthScale=2,sigma=1,func=sin)
draw<-sampleGP(noDraw=20,xVal=seq(-6,6,length.out=100),covFunc=SqrExpCov,lengthScale=2,sigma=1,func=sin,onlyPara=T)
#round(draw$covMat,4)
image(draw$covMat)
plotGP(drawMatrix=draw$draws)

draw<-sampleGP(noDraw=20,xVal=seq(-6,6,length.out=100),covFunc=SqrExpCov)
