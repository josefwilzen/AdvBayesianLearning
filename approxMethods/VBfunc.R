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
  aLambda<-sumStat+n*alpha+n+1
  aBeta<-n*psi-n+1
  for(i in 2:dim(draws)[1]){
    # cond. for lambda:
    bLabmda<-n+draws[(i-1),2]*n
    draws[i,1]<-rgamma(n=1,shape=aLambda,rate=bLabmda)
    # cond. for beta:
    bBeta<-n*draws[(i-1),1]+n*theta
    draws[i,2]<-rgamma(n=1,shape=aBeta,rate=bBeta)
  }
  draws<-draws[-1,]
  draws<-draws[-(1:burnIn),]
  draws<-draws[seq(from=1,to=dim(draws)[1],by=thin),]
  return(draws)  
}

#-------------------------------------------------------------------------
# c) VB



vbParameterEstimate<-function(obsData,myB0=2,nIter=20,psi=1,alpha=1,theta=5){
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

likeFreeRejectionSample2<-function(obsData,statistic=sum,distFunc=abs,nObs=1000,tol=0.0001,psi=1,beta=1,theta=5,...){
  result<-matrix(0,nrow=nObs,ncol=2)
  colnames(result)<-c("lambda","beta")
  yStat<-statistic(obsData,...)
  for(i in 1:nObs){
    
    # generate the first proposal:
    beta<-rgamma(n=1,shape=psi,rate=theta)
    lambda<-rgamma(n=1,shape=alpha,rate=beta)
    z<-rpois(n=10,lambda=lambda)
    zStat<-statistic(z,...)
    count<-0
    while(distFunc(yStat-zStat)>tol){
      #print(paste("Dist: ",distFunc(yStat-zStat)))
      beta<-rgamma(n=1,shape=psi,rate=5)
      lambda<-rgamma(n=1,shape=alpha,rate=beta)
      z<-rpois(n=10,lambda=lambda)
      zStat<-statistic(z,...)
      count<-count+1
    }
    tempDist<-distFunc(yStat-zStat)
    result[i,]<-c(lambda,beta)
    print(paste("i: ",i," count: ",count,"tempDist: ",tempDist))
    count<-1
    
  }

  return(result)
}





















