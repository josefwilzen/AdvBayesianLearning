#-------------------------------------------------------------------------
# 1) VB
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# b) Gibbs sampler


gibbsHierarchialPoisson<-function(lambdaStart,betaStart,obsData,mySeed=NULL,nObs=1000,alpha=1,psi=1,theta=5,burnIn=300,thin=5){
  noDraws<-nObs*thin+burnIn
  set.seed(mySeed)
  # matrix for draws
  draws<-matrix(0,nrow=(noDraws+1),ncol=2)
  draws[1,]<-c(lambdaStart,betaStart)
  colnames(draws)<-c("lambda","beta")
  sumStat<-sum(obsData)
  n<-length(obsData)
  aLambda<-sumStat+alpha
  aBeta<-psi
  for(i in 2:dim(draws)[1]){
    # cond. for lambda:
    bLabmda<-n+draws[(i-1),2]
    draws[i,1]<-rgamma(n=1,shape=aLambda,rate=bLabmda)
    # cond. for beta:
    bBeta<-draws[(i-1),1]+theta
    draws[i,2]<-rgamma(n=1,shape=aBeta,rate=bBeta)
  }
  draws<-draws[-1,]
  draws<-draws[-(1:burnIn),]
  draws<-draws[seq(from=1,to=dim(draws)[1],by=thin),]
  return(draws)  
}

#-------------------------------------------------------------------------
# c) VB



vbParameterEstimate<-function(obsData,myB0=2,nIter=50,psi=1,alpha=1,theta=5){
  sumX<-sum(obsData)
  n<-length(obsData)
  aLambda<-sumX+alpha
  aBeta<-psi
  myB<-myB0
  for(i in 1:nIter){
    bLambda<-n+myB
    myLambda<-aLambda/bLambda
    bBeta<-theta+myLambda
    myB<-aBeta/bBeta
  }
  result<-matrix(0,2,2)
  colnames(result)<-c("lambda","beta")
  rownames(result)<-c("a","b")
  result[,1]<-c(aLambda,bLambda)
  result[,2]<-c(aBeta,bBeta)
  return(result)
}


vbDraws<-function(paraMatrix,noDraws=1000,mySeed=NULL){
  result<-matrix(0,noDraws,2)
  colnames(result)<-c("lambda","beta")
  set.seed(mySeed)
  result[,1]<-rgamma(n=noDraws,shape=paraMatrix[1,1],rate=paraMatrix[2,1])
  result[,2]<-rgamma(n=noDraws,shape=paraMatrix[1,2],rate=paraMatrix[2,2])
  return(result)
}

vb<-function(obsData,noDraws=1000,mySeed=NULL,...){
  para<-vbParameterEstimate(obsData=obsData,...)
  return(vbDraws(paraMatrix=para,noDraws=noDraws,mySeed=mySeed))
}




#-------------------------------------------------------------------------
# 2) ABC
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# a) 



euclideanDist<-function(x,y){
  z<-sqrt(sum((x-y)^2))
  return(z)
}



likeFreeRejectionSample2<-function(obsData,statistic=sum,nObs=1000,tol=1,psi=1,beta=1,theta=5,alpha=1,mySeed=NULL,storeData=FALSE,...){
   set.seed(mySeed)
  result<-matrix(0,nrow=nObs,ncol=3)
  colnames(result)<-c("lambda","beta","count")
  yStat<-statistic(obsData,...)
  n<-length(obsData)
  resultData<-matrix(0,nrow=nObs,ncol=n)
  for(i in 1:nObs){
    # generate the first proposal:
    beta<-rgamma(n=1,shape=psi,rate=theta)
    lambda<-rgamma(n=1,shape=alpha,rate=beta)
    z<-rpois(n=10,lambda=lambda)
    zStat<-statistic(z,...)
    count<-1
    while(euclideanDist(x=zStat,y=yStat)>tol){
      beta<-rgamma(n=1,shape=psi,rate=theta)
      lambda<-rgamma(n=1,shape=alpha,rate=beta)
      z<-rpois(n=n,lambda=lambda)
      zStat<-statistic(z,...)
      count<-count+1
    }
    tempDist<-euclideanDist(x=zStat,y=yStat)
    result[i,]<-c(lambda,beta,count)
    resultData[i,]<-z
    print(paste("i: ",i," count: ",count,"tempDist: ",tempDist))
    count<-1 
  }
  if(storeData){
    res<-list(para=result,dataSample=resultData)
    return(res)
  }
  return(result)
}

library(MASS)
library(mvtnorm)

posterior<-function(paraVect,obsData,psi=1,theta=5,alpha=1){
  # lambda first elemnt of paraVect and beta second element of paraVect
  lambda<-paraVect[1]
  beta<-paraVect[2]
  n<-length(obsData)
  postDens<-lambda^(sum(obsData))*exp(-lambda*n)*lambda^(alpha-1)*exp(-lambda*beta)*beta^(psi-1)*exp(-beta*theta)
  return(postDens)
}



likeFreeMCMC<-function(obsData,tuning,statistic=sum,nObs=1000,tol=1,psi=1,theta=5,alpha=1,mySeed=NULL,...){
  startObj<-likeFreeRejectionSample2(obsData=obsData,nObs=1,storeData=TRUE,mySeed=mySeed)
  
  startPara<-startObj$para[1,1:2]
  varProp<-tuning*diag(2)
  set.seed(mySeed)
  resultPara<-matrix(0,nrow=(nObs+1),ncol=4)
  resultPara[1,1:2]<-startPara
  colnames(resultPara)<-c("lambda","beta","all","acceted")
  yStat<-statistic(obsData,...)
  n<-length(obsData)
  resultData<-matrix(0,nrow=(nObs+1),ncol=n)
  resultData[1,]<-startObj$dataSample[1,]
  
  accept<-0
  
  for(i in 2:(nObs+1)){
    #if(i==2) browser()
    
    parameterProp<-rmvnorm(n=1,mean=resultPara[(i-1),1:2],sigma=varProp)
    parameterProp<-ifelse(parameterProp<=0,0.000001,parameterProp)
    
    lambda<-rgamma(n=1,shape=alpha,rate=parameterProp[,2])
    z<-rpois(n=10,lambda=lambda)
    zStat<-statistic(z,...)
    u<-runif(1)
    posteriorOld<-posterior(paraVect=resultPara[(i-1),1:2],obsData=obsData,psi=psi,theta=theta,alpha=alpha)
    posteriorNew<-posterior(paraVect=parameterProp,obsData=obsData,psi=psi,theta=theta,alpha=alpha)
    
    propGivenOld<-dmvnorm(x=parameterProp,mean=resultPara[(i-1),1:2],sigma=varProp)
    oldGivenProp<-dmvnorm(x=resultPara[(i-1),1:2],mean=parameterProp,sigma=varProp)
    
    posteriorRatio<-posteriorNew/posteriorOld
    proposalRatio<-oldGivenProp/propGivenOld
    
    cat("posteriorRatio: ",posteriorRatio,"\n")
    cond1<-posteriorRatio*proposalRatio>=u
    cond2<-euclideanDist(x=zStat,y=yStat)<=tol
    
    cat("cond1: ",cond1,"cond2: ",cond2,"\n")
    if(cond1&cond2){
      resultPara[i,1:2]<-parameterProp
      resultData[i,]<-z
    }else{
      resultPara[i,1:2]<-resultPara[(i-1),1:2]
      resultData[i,]<-resultData[(i-1),]
    }
    
    # generate the first proposal:
#     beta<-rgamma(n=1,shape=psi,rate=theta)
#     lambda<-rgamma(n=1,shape=alpha,rate=beta)
#     z<-rpois(n=10,lambda=lambda)
#     zStat<-statistic(z,...)
#     count<-1
#     while(euclideanDist(x=zStat,y=yStat)>tol){
#       beta<-rgamma(n=1,shape=psi,rate=theta)
#       lambda<-rgamma(n=1,shape=alpha,rate=beta)
#       z<-rpois(n=n,lambda=lambda)
#       zStat<-statistic(z,...)
#       count<-count+1
#     }
    tempDist<-euclideanDist(x=zStat,y=yStat)
    print(paste("i: ",i,"tempDist: ",tempDist))
    count<-1 
  }

  res<-list(parameters=resultPara,dataSample=resultData)

  return(res)
}













