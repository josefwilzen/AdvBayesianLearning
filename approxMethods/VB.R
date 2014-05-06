#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
# Variational bayes 
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------



#-------------------------------------------------------------------------
# b)
rm(list=ls())
theta<-5
betaStart<-1
lambdaStart<-1
dataVect<-c(1,0,3,4,2,6,2,3,1,1)

source("/home/joswi05/Dropbox/Josef/Advanced Bayesian Learning/AdvBayesianLearning/approxMethods/VBfunc.R")

test<-gibbsHierarchialPoisson(lambdaStart=1,betaStart=1,obsData=dataVect,mySeed=4567,nObs=5000)

plot(test)
plot(test[,1],type="l")
plot(test[,2],type="l")
plot(density(test[,1]),main="lambda")
plot(density(test[,2]),main="beta")
cov(x=test[,1],y=test[,2])

library(ggplot2)
pl<-ggplot(data=as.data.frame(test),aes(x=lambda,y=beta))+geom_point()+theme_bw()+
  geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
print(pl)

#pl+scale_fill_gradient(limits=c(1e-5,8e-4))
#pl+stat_density2d(geom="point", aes(size = ..density..), contour = FALSE)

par(mfrow=c(2,2))
plot(density(test[,1]),main="lambda")
plot(density(test[,2]),main="beta")
#plot(test)
par(mfrow=c(1,1))
example(stat_density2d)



#  ABC
function(x){return(sum(log(x)))

test2<-likeFreeRejectionSample2(obsData=dataVect,tol=1e-5,nObs=1000,statistic=mad,distFunc=tempFunc)
plot(density(test2[,1]),main="lambda")
plot(density(test2[,2]),main="beta")
pl<-ggplot(data=as.data.frame(test2),aes(x=lambda,y=beta))+geom_point()+theme_bw()+
  geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
print(pl)

tempFunc<-function(x){
  y<-sum(abs(x))
  return(y)
}


par(mfrow=c(2,2))
plot(density(test[,1]),main="lambda")
plot(density(test[,2]),main="beta")
plot(density(test2[,1]),main="lambda")
plot(density(test2[,2]),main="beta")
par(mfrow=c(1,1))
dist(x=cbind(2:3,20:21))







