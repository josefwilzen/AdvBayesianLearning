
meanFunc<-function(xValues,func=NULL){
  if(is.null(func)) return(rep(0,length(xValues)))
  yValues<-func(xValues)
  return(yValues)
}

SqrExpCov<-function(x1,x2,sigma=1,lengthScale=1){
  if(lengthScale<=0) stop("Values of lengthScale not allowed")
  if(sigma<=0) stop("Values of sigma not allowed")
  covVal<-(sigma^2)*exp(((x1-x2)^2)/(-2*lengthScale))
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
  maxVal<-max(drawMatrix)+0.5
  minVal<-min(drawMatrix)-0.5
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


plotEmpiricalProbBand<-function(drawMatrix,postMean,xGrid,dataPoints=NULL){
  maxVal<-max(drawMatrix)+0.5
  minVal<-min(drawMatrix)-0.5
  noDraws<-dim(drawMatrix)[1]
  confBands<-apply(X=drawMatrix,MARGIN=2,FUN=quantile,prob=c(0.05,0.95))
  plot(x=xGrid,y=postMean, type="l", 
       xlab="t",ylab="y",ylim=c(minVal,maxVal),lwd=3)
  lines(x=xGrid,y=confBands[1,],col="blue")  # lower
  lines(x=xGrid,y=confBands[2,],col="blue")  # upper
}



meanWages<-function(x){
  y<-14*tanh(x+2.8)
  return(y)
}


posteriorDist<-function(xGrid,priorMean,obsData,covFunc,sigmaError,...){
  y<-obsData[,1]
  x<-obsData[,2] 
  # mean function for new x
  meanXgrid<-as.matrix(meanFunc(xValues=xGrid,func=priorMean))
  
  # the covariance matrix K(X_star,X):
  xKxGrid<-createCovMatrix(xVal=x,yVal=xGrid,covFunc=covFunc,...)
  
  # the inverse K_y = K(X,X)+simga^2*I:
  sqrtKy<-qr.solve(createCovMatrix(xVal=x,covFunc=covFunc,...)
                   +((sigmaError)^2)*diag(1,nrow=dim(obsData)[1]))
  
  # differnece between obs y and mean of X:
  yDiff<-y-meanFunc(xValues=x,func=priorMean)
  
  # posterior mean:
  fBar<-meanXgrid+xKxGrid%*%sqrtKy%*%yDiff
  
  # posterior covariance:
  temp1<-qr.solve(createCovMatrix(xVal=wages$age,covFunc=SqrExpCov)+((sigmaError)^2)*diag(1,nrow=dim(wages)[1]))
  temp2<-createCovMatrix(xVal=wages$age,yVal=xGrid,covFunc=SqrExpCov)
  temp3<-createCovMatrix(xVal=xGrid,yVal=wages$age,covFunc=SqrExpCov)
  temp4<-createCovMatrix(xVal=xGrid,covFunc=SqrExpCov)
  posteriorCov<-temp4-temp2%*%temp1%*%temp3
  res<-list(meanVal=fBar,covMat=posteriorCov,xVal=xGrid,sigmaError=sigmaError,obsData=obsData)
  return(res)
}


plotTheoreticalProbBand<-function(postDist,predictiveMean=FALSE,quan=c(0.025,0.975),plotLim=c(10,15),plotData=FALSE,...){
  varVect<-diag(postDist$covMat)
  meanVect<-postDist$meanVal
  obsData<-postDist$obsData
  x<-postDist$xVal
  sigmaError<-postDist$sigmaError
  borders<-qnorm(quan)
  bands<-matrix(0,2,length(meanVect))
  if(predictiveMean){
    varVect<-varVect+sigmaError^2
  }
  bands[1,]<-postDist$meanVal+borders[1]*sqrt(varVect)  # lower
  bands[2,]<-postDist$meanVal+borders[2]*sqrt(varVect)  # upper
  
  if(plotData){
    plot(x=obsData[,2],y=obsData[,1],ylim=plotLim,col="blue",...)
    lines(x=x,y=meanVect,lwd=3)
    lines(x=x,y=bands[1,],col="red",lwd=2)
    lines(x=x,y=bands[2,],col="red",lwd=2)
  }else{
    plot(x=x,y=meanVect,ylim=plotLim,type="l",lwd=3,...)
    lines(x=x,y=bands[1,],col="red")
    lines(x=x,y=bands[2,],col="red")
  }
}



