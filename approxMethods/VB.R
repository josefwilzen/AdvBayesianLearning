#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
# Variational bayes 
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------

#install.packages("gridExtra")
library(gridExtra)
library(ggplot2)

#-------------------------------------------------------------------------
# b)
rm(list=ls())
theta<-5
betaStart<-1
lambdaStart<-1
dataVect<-c(1,0,3,4,2,6,2,3,1,1)

source("/home/joswi05/Dropbox/Josef/Advanced Bayesian Learning/AdvBayesianLearning/approxMethods/VBfunc.R")

test<-gibbsHierarchialPoisson(lambdaStart=1,betaStart=1,obsData=dataVect,mySeed=49067,nObs=1000)
pGibbs1<-ggplot(data=as.data.frame(test),aes(x=lambda,y=beta))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
pGibbs2<-ggplot(as.data.frame(test), aes(lambda)) +geom_density(alpha = 0.2)+theme_bw()+xlab(expression(lambda)) 
pGibbs3<-ggplot(as.data.frame(test), aes(beta)) +geom_density(alpha = 0.2)+theme_bw()+xlab(expression(beta))
grid.arrange(pGibbs2,pGibbs3,pGibbs1, ncol=2)


#-------------------------------------------------------------------------
# c)

test4<-vb(obsData=dataVect,mySeed=20939)
pl<-ggplot(data=as.data.frame(test4),aes(x=lambda,y=beta))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
plVB1<-ggplot(as.data.frame(test4), aes(lambda)) +geom_density(alpha = 0.2)+theme_bw()+xlab(expression(lambda)) 
plVB2<-ggplot(as.data.frame(test4), aes(beta)) +geom_density(alpha = 0.2)+theme_bw()+xlab(expression(beta))
grid.arrange(plVB1,plVB2,pl, ncol=2)


#--------------------------------------------
# computing times:

system.time(test1<-gibbsHierarchialPoisson(lambdaStart=1,betaStart=1,obsData=dataVect,mySeed=49067,nObs=1000))
system.time(test2<-vb(obsData=dataVect,mySeed=20939))

#--------------------------------------------
# comparing:

testVB<-data.frame(test4,Method="VB")
testGibbs<-data.frame(test,Method="Gibbs")
testTotal<-rbind(testGibbs,testVB)

totLambda<-ggplot(testTotal, aes(x=lambda, color=Method)) + geom_density()+theme_bw()+xlab(expression(lambda))
totBeta<-ggplot(testTotal, aes(x=beta, color=Method)) + geom_density()+theme_bw()+xlab(expression(lambda))
grid.arrange(totLambda,totBeta,ncol=2)

#ggplot(df, aes(x=rating, fill=cond)) + geom_density(alpha=.3)

pl<-ggplot(data=as.data.frame(testTotal),aes(x=lambda,y=beta,col=Method))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
print(pl)
testTotal

#-------------------------------------------------------------------------
#  ABC


# a)

test2<-likeFreeRejectionSample2(obsData=dataVect,tol=1,nObs=1000,statistic=sum,mySeed=2344)
pReject<-ggplot(data=as.data.frame(test2),aes(x=lambda,y=beta))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
pReject1<-ggplot(as.data.frame(test2), aes(lambda)) +geom_density(alpha = 0.2)+theme_bw()+xlab(expression(lambda)) 
pReject2<-ggplot(as.data.frame(test2), aes(beta)) +geom_density(alpha = 0.2)+theme_bw()+xlab(expression(beta))
grid.arrange(pReject1,pReject2,pReject, ncol=2)
summary(test2[,3])


test3<-likeFreeRejectionSample2(obsData=dataVect,tol=1,nObs=1000,statistic=mean,mySeed=NULL)
pReject<-ggplot(data=as.data.frame(test4),aes(x=lambda,y=beta))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
pReject1<-ggplot(as.data.frame(test4), aes(lambda)) +geom_density(alpha = 0.2)+theme_bw()+xlab(expression(lambda)) 
pReject2<-ggplot(as.data.frame(test4), aes(beta)) +geom_density(alpha = 0.2)+theme_bw()+xlab(expression(beta))
grid.arrange(pReject1,pReject2,pReject, ncol=2)
summary(test2[,3])


rejection1<-likeFreeRejectionSample2(obsData=dataVect,tol=0,nObs=1000,statistic=mean,mySeed=2344)
rejection2<-likeFreeRejectionSample2(obsData=dataVect,tol=0,nObs=1000,statistic=sum,mySeed=3444)
rejection3<-likeFreeRejectionSample2(obsData=dataVect,tol=0,nObs=1000,statistic=median,mySeed=3554)
rejection4<-likeFreeRejectionSample2(obsData=dataVect,tol=0,nObs=1000,statistic=range,mySeed=34006)

plRejection1<-ggplot(data=as.data.frame(rejection1),aes(x=lambda,y=beta))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
plRejection2<-ggplot(data=as.data.frame(rejection2),aes(x=lambda,y=beta))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
plRejection3<-ggplot(data=as.data.frame(rejection3),aes(x=lambda,y=beta))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
plRejection4<-ggplot(data=as.data.frame(rejection4),aes(x=lambda,y=beta))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
grid.arrange(plRejection1,plRejection2,plRejection3,plRejection4, ncol=2)



rejection5<-likeFreeRejectionSample2(obsData=dataVect,tol=0,nObs=1000,statistic=mean,mySeed=2144)
rejection6<-likeFreeRejectionSample2(obsData=dataVect,tol=3,nObs=1000,statistic=mean,mySeed=22444)
rejection7<-likeFreeRejectionSample2(obsData=dataVect,tol=6,nObs=1000,statistic=mean,mySeed=355411)
rejection8<-likeFreeRejectionSample2(obsData=dataVect,tol=9,nObs=1000,statistic=mean,mySeed=34756)

plRejection5<-ggplot(data=as.data.frame(rejection5),aes(x=lambda,y=beta))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
plRejection6<-ggplot(data=as.data.frame(rejection6),aes(x=lambda,y=beta))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
plRejection7<-ggplot(data=as.data.frame(rejection7),aes(x=lambda,y=beta))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
plRejection8<-ggplot(data=as.data.frame(rejection8),aes(x=lambda,y=beta))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
grid.arrange(plRejection5,plRejection6,plRejection7,plRejection8, ncol=2)


testSample<-likeFreeRejectionSample2(obsData=dataVect,tol=0,nObs=1000,statistic=mean,mySeed=22144,storeData=TRUE)

str(testSample)

#--------------------------------------------------------------------------------
# b)
a<-likeFreeMCMC(obsData=dataVect,tuning=0.1,statistic=mean,nObs=1000,tol=1,mySeed=2873)

plot(a$parameters[,1],type="l")
plot(a$parameters[,2],type="l")
aPlot<-ggplot(data=as.data.frame(a$parameters[,1:2]),aes(x=lambda,y=beta))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
print(aPlot)