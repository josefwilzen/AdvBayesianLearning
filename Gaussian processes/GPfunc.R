
meanFunc<-function(xValues,func=NULL){
  if(is.null(func)) return(rep(0,length(xValues)))
  yValues<-func(xValues)
  return(yValues)
}

SqrExpCov<-function(x1,x2,sigma=1,lengthScale=1){
  if(lengthScale<=0) stop("Values of lengthScale not allowed")
  if(sigma<=0) stop("Values of sigma not allowed")
  covVal<-(sigma^2)*exp(((x1-x2)^2)/(-2*lengthScale^2))
  return(covVal)
}

maternKernel<-function(x1,x2,sigma=1,lengthScale=1){
  if(lengthScale<=0) stop("Values of lengthScale not allowed")
  if(sigma<=0) stop("Values of sigma not allowed")
  #if(nu<=0) stop("Values of nu not allowed")
  
  nu<-1.5 # 0.5, 1.5 2.5
  
  r<-abs(x1-x2)
  x<-(sqrt(2*nu)*r)/lengthScale
  y<-(2^(1-nu)/gamma(nu))*(((sqrt(2*nu)*r)/lengthScale)^nu)*besselI(x=x,nu=nu)
  return(y)
}

rationalQuaratic<-function(x1,x2,sigma=1,alpha=1,lengthScale=1){
  if(lengthScale<=0) stop("Values of lengthScale not allowed")
  if(sigma<=0) stop("Values of sigma not allowed")
  if(alpha<=0) stop("Values of alpha not allowed")
  
  r<-abs(x1-x2)
  y<-sigma*(1+r^2/(2*alpha*lengthScale^2))^(-alpha)
  return(y)
}

gammaExponential<-function(x1,x2,sigma=1,gamma=1,lengthScale=1){
  if(lengthScale<=0) stop("Values of lengthScale not allowed")
  if(sigma<=0) stop("Values of sigma not allowed")
  if(gamma<=0|gamma>2) stop("Values of gamma not allowed")
  r<-abs(x1-x2)
  y<- sigma*exp(-(r/lengthScale)^gamma)
  return(y)
}




createCovMatrix<-function(xVal,yVal=NULL,covFunc,...){
  # index<-1:length(xVal)
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


# covFunc- covariance function
# para - list with all parameters for covFunc

createCovMatrix2<-function(covFunc,para){
  # index<-1:length(xVal)
  #   indexMat<-as.matrix(expand.grid(index,index))
  #   covValues<-covFunc(x1=indexMat[,1],x2=indexMat[,2],...)
  #   covMat<-matrix(0,length(xVal),length(xVal))
  #   covMat[indexMat]<-covValues
  
  if(is.null(yVal)){
    covMat<-matrix(NA,length(xVal),length(xVal))
    for(i in 1:length(xVal)){
      for(j in 1:length(xVal)){
        if(j<=i){
          covMat[i,j]<-covFunc(x1=xVal[i],x2=xVal[j],...)
        }
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

plotGP<-function(drawMatrix,...){
  maxVal<-max(drawMatrix)+0.5
  minVal<-min(drawMatrix)-0.5
  noDraws<-dim(drawMatrix)[1]
  if(noDraws==1) {
    plot(x=as.numeric(colnames(drawMatrix)),y=drawMatrix[1,], type="l",...)
  }else{
    plot(x=as.numeric(colnames(drawMatrix)),y=drawMatrix[1,], type="l",ylim=c(minVal,maxVal),...)
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
  temp1<-qr.solve(createCovMatrix(xVal=x,covFunc=covFunc,...)+((sigmaError)^2)*diag(1,nrow=dim(obsData)[1]))
  temp2<-createCovMatrix(xVal=x,yVal=xGrid,covFunc=covFunc,...)
  temp3<-createCovMatrix(xVal=xGrid,yVal=x,covFunc=covFunc,...)
  temp4<-createCovMatrix(xVal=xGrid,covFunc=covFunc,...)
  posteriorCov<-temp4-temp2%*%temp1%*%temp3
  res<-list(meanVal=fBar,covMat=posteriorCov,xVal=xGrid,sigmaError=sigmaError,obsData=obsData)
  return(res)
}


plotTheoreticalProbBand<-function(postDist,predictiveMean=TRUE,quan=c(0.025,0.975),plotLim=c(10,15.5),plotData=TRUE,overlay=TRUE,...){
  varVect<-diag(postDist$covMat)
  meanVect<-postDist$meanVal
  obsData<-postDist$obsData
  x<-postDist$xVal
  sigmaError<-postDist$sigmaError
  borders<-qnorm(quan)
  bands<-matrix(0,2,length(meanVect))
  predBands<-matrix(0,2,length(meanVect))
  if(predictiveMean){
    varVectPred<-varVect+sigmaError^2
    predBands[1,]<-postDist$meanVal+borders[1]*sqrt(varVectPred)  # lower
    predBands[2,]<-postDist$meanVal+borders[2]*sqrt(varVectPred)  # upper
  }
  bands[1,]<-postDist$meanVal+borders[1]*sqrt(varVect)  # lower
  bands[2,]<-postDist$meanVal+borders[2]*sqrt(varVect)  # upper
  
  if(plotData&overlay){
    plot(x=obsData[,2],y=obsData[,1],ylim=plotLim,xlim=c(min(x),max(x)),col="grey",...)
    lines(x=x,y=meanVect,lwd=3)
    lines(x=x,y=bands[1,],col="red",lwd=2)
    lines(x=x,y=bands[2,],col="red",lwd=2)
    lines(x=x,y=predBands[1,],col="blue",lwd=2)
    lines(x=x,y=predBands[2,],col="blue",lwd=2)
    legend("bottomright",
           c("mean","prob. intervals","predictive intervals"), 
           col=c("black","red","blue"),horiz=FALSE,lty = c(1, 1, 1),lwd=c(3,2,2))
#     legend(-1, 1.9, c("sin", "cos", "tan"), col = c(3, 4, 6),
#            , lty = c(1, 1, 1),
#             bg = "gray90")
    
  }else if(plotData){
    plot(x=obsData[,2],y=obsData[,1],ylim=plotLim,xlim=c(min(x),max(x)),col="blue",...)
    lines(x=x,y=meanVect,lwd=3)
    lines(x=x,y=bands[1,],col="red",lwd=2)
    lines(x=x,y=bands[2,],col="red",lwd=2)
  }else{
    plot(x=x,y=meanVect,ylim=plotLim,xlim=c(min(x),max(x)),type="l",lwd=3,...)
    lines(x=x,y=bands[1,],col="red")
    lines(x=x,y=bands[2,],col="red")
  }
}


marginalLikelihood<-function(obsData,meanFunc,covFunc,sigmaError,...){
  if(sigmaError<=0) stop("Value for sigmaError is not allowed!")
  #browser()
  y<-obsData[,1]  
  x<-obsData[,2]
  meanDiff<-as.matrix(y-meanFunc(x))
  covMat<-createCovMatrix(xVal=x,covFunc=covFunc,...)+(sigmaError^2)*diag(length(x))
  #if(det(covMat)<0) return(-1e15)
  margLike<- -0.5*t(meanDiff)%*%solve(covMat)%*%meanDiff-0.5*log(det(covMat))-0.5*length(x)*log(2*pi)
  
  #if(!is.finite(-0.5*log(det(covMat)))) browser()
  
  #cat("\n part1: ",-0.5*t(meanDiff)%*%qr.solve(covMat)%*%meanDiff, " part2: ",-0.5*log(det(covMat))," part3: ",-0.5*length(x)*log(2*pi),"margLike: ",margLike,"\n")
  #if(is.nan(x=log(det(covMat)))) browser()
  if(det(covMat)<1e-15) return(-1e10) #print("det(covMat=0")#
  return(margLike)
}


optimMarginalLikelihood<-function(para,covFunc,...){
  paraName<-c("sigmaError",names(unlist(formals(covFunc)[-(1:2)])))
  otherPara<-as.list(match.call())[-(1:3)]
  paraList<-vector("list",length(paraName))
  for(i in 1:length(paraName)){
    paraList[[i]]<-para[i]
    names(paraList)[i]<-paraName[i]
  }
  allPara<-c(covFunc=covFunc,otherPara,paraList)
  #print(allPara)
  print(allPara[4:length(allPara)])
  res<-do.call(marginalLikelihood,allPara)
  # we want to maximize
  
  # controll only for -inf
#   if(res==-Inf){
#     res<- -1e10
#   }
  # controll both for inf and -inf
#   if(res==Inf){
#     res<-100000000
#   }else if(res==-Inf){
#     res<- -100000000
#   }
  cat("Value: ",as.vector(res),"\n")
  return(as.vector(res))
  #return(allPara)
}



# marginalLikelihood<-function(...,obsData,meanFunc,covFunc,sigmaError){
#   y<-obsData[,1]  
#   x<-obsData[,2]
#   meanDiff<-as.matrix(meanFunc(x)-y)
#   covMat<-createCovMatrix(xVal=x,covFunc=covFunc,...)+(sigmaError^2)*diag(length(x))
#   margLike<- -0.5*t(meanDiff)%*%covMat%*%meanDiff-0.5*log(det(covMat))-0.5*length(x)*log(2*pi)
#   return(margLike)
# }




SqrExpCovVect<-function(theta,x1,x2){
  sigma<-theta[2]
  sigmaError<-theta[1]
  l<-theta[3:length(theta)]
  x1<-as.matrix(x1)
  x2<-as.matrix(x2)
  M<-diag(l^(-2))
  y<-sigma^2*exp(-0.5*t(x1-x2)%*%M%*%(x1-x2))
  if(all(x1==x2)){
   y<-y+sigmaError^2
  }
  return(y)
}


covMatrixVect<-function(theta,myData){
  noObs<-dim(myData)[1]
  covMat<-matrix(0,noObs,noObs)
  for(i in 1:noObs){  # over rows 
    for(j in 1:noObs){  # over cols
      x1<-t(myData[i,])
      x2<-t(myData[j,])
      covMat[i,j]<-SqrExpCovVect(theta=theta,x1=x1,x2=x2)
      
    }
  }
  return(covMat)
}

# try to speed up with Matrix package

# SqrExpCovVect3<-function(x,noXvar,sigma,M){
#   #browser()
#   x1<-x[1:noXvar]
#   x2<-x[-(1:noXvar)]
#   y<-sigma^2*exp(-0.5*t(x1-x2)%*%M%*%(x1-x2))
#   return(y)
# }
# 
# covMatrixVect2<-function(theta,myData){
#   #browser()
#   myData<-Matrix(as.matrix(myData))
#   noObs<-dim(myData)[1]
#   noCols<-dim(myData)[2]
#   covMat<-Matrix(0,noObs,noObs)
#   index<-1:noObs
#   matrixIndex<-Matrix(as.matrix(expand.grid(index,index)))
#   allCombin<-cbind(myData[matrixIndex[,1],],myData[matrixIndex[,2],])
#   sigmaError<-theta[1]
#   sigma<-theta[2]
#   l<-theta[3:length(theta)]
#   M<-Diagonal(x=l^(-2))
#   Y<-apply(X=allCombin[[1]],MARGIN=1,FUN=SqrExpCovVect3,noXvar=noCols,sigma=sigma,M=as.matrix(M))
#   covMat[matrixIndex]<-Y
#   covMat<-covMat+sigmaError^2*diag(noObs)
#   return(covMat)
# }

# back up:

SqrExpCovVect3<-function(x,noXvar,sigma,M){
  x1<-x[1:noXvar]
  x2<-x[-(1:noXvar)]
  y<-sigma^2*exp(-0.5*t(x1-x2)%*%M%*%(x1-x2))
  return(y)
}

covMatrixVect2<-function(theta,myData){
  noObs<-dim(myData)[1]
  noCols<-dim(myData)[2]
  covMat<-matrix(0,noObs,noObs)
  index<-1:noObs
  matrixIndex<-as.matrix(expand.grid(index,index))
  allCombin<-cbind(myData[matrixIndex[,1],],myData[matrixIndex[,2],])
  sigmaError<-theta[1]
  sigma<-theta[2]
  l<-theta[3:length(theta)]
  M<-diag(l^(-2))
  Y<-apply(X=allCombin,MARGIN=1,FUN=SqrExpCovVect3,noXvar=noCols,sigma=sigma,M=M)
  covMat[matrixIndex]<-Y
  covMat<-covMat+sigmaError^2*diag(noObs)
  return(covMat)
}

library(Matrix)

# likelihood:
likelihoodGPclassify<-function(f,y,LOG=TRUE,...){
  f<-as.vector(f)
  y<-as.vector(y)
  if(LOG){
    likelihood<- ifelse(y==1,-log(1+exp(-f)),-f-log(1+exp(-f)))
  }else{
    likelihood<-ifelse(y==1,inv.logit(x=f),(1-inv.logit(x=f)))
  }
  return(likelihood)
}



logisticHessian<-function(f,y){
  W<--Diagonal(x=-likelihoodGPclassify(f=f,y=y,LOG=TRUE)*likelihoodGPclassify(f=f,y=-y))
  return(W)
}

logisticGradient<-function(f,y){
  tVect<-(y+1)/2
  G<-tVect-likelihoodGPclassify(f=f,y=y,LOG=TRUE)*likelihoodGPclassify(f=f,y=-y,LOG=TRUE)
  return(G)
}



















