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
totBeta<-ggplot(testTotal, aes(x=beta, color=Method)) + geom_density()+theme_bw()+xlab(expression(beta))
grid.arrange(totLambda,totBeta,ncol=2)

#ggplot(df, aes(x=rating, fill=cond)) + geom_density(alpha=.3)

pl<-ggplot(data=as.data.frame(testTotal),aes(x=lambda,y=beta,col=Method))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
print(pl)


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


# compare methods:
testVB<-data.frame(test4,Method="VB")
testGibbs<-data.frame(test,Method="Gibbs")
testRejectionSampler<-data.frame(test2[,1:2],Method="Rejection")
testTotal<-rbind(testVB,testGibbs,testRejectionSampler)

aggregate(x=testTotal[,1],by=list(testTotal[,3]),FUN=summary)
aggregate(x=testTotal[,2],by=list(testTotal[,3]),FUN=summary)

totLambda<-ggplot(testTotal, aes(x=lambda, color=Method)) + geom_density()+theme_bw()+xlab(expression(lambda))
totBeta<-ggplot(testTotal, aes(x=beta, color=Method)) + geom_density()+theme_bw()+xlab(expression(beta))
grid.arrange(totLambda,totBeta,ncol=2)

pl<-ggplot(data=as.data.frame(testTotal),aes(x=lambda,y=beta,col=Method))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
print(pl)



# compare statistic

rejection1<-likeFreeRejectionSample2(obsData=dataVect,tol=1,nObs=1000,statistic=mean,mySeed=2344)
rejection2<-likeFreeRejectionSample2(obsData=dataVect,tol=1,nObs=1000,statistic=var,mySeed=3444)
rejection3<-likeFreeRejectionSample2(obsData=dataVect,tol=1,nObs=1000,statistic=median,mySeed=3554)
rejection4<-likeFreeRejectionSample2(obsData=dataVect,tol=1,nObs=1000,statistic=range,mySeed=34006)

plRejection1<-ggplot(data=as.data.frame(rejection1),aes(x=lambda,y=beta))+geom_point()+theme_bw()+
geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))+ggtitle("Mean")

plRejection2<-ggplot(data=as.data.frame(rejection2),aes(x=lambda,y=beta))+geom_point()+theme_bw()+
geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))+ggtitle("Variance")

plRejection3<-ggplot(data=as.data.frame(rejection3),aes(x=lambda,y=beta))+geom_point()+theme_bw()+
geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))+ggtitle("Median")

plRejection4<-ggplot(data=as.data.frame(rejection4),aes(x=lambda,y=beta))+geom_point()+theme_bw()+
geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))+ggtitle("Range")

grid.arrange(plRejection1,plRejection3,plRejection2,plRejection4, ncol=2)


# compare error tolarence:

rejection5<-likeFreeRejectionSample2(obsData=dataVect,tol=0,nObs=1000,statistic=mean,mySeed=2144)
rejection6<-likeFreeRejectionSample2(obsData=dataVect,tol=3,nObs=1000,statistic=mean,mySeed=22444)
rejection7<-likeFreeRejectionSample2(obsData=dataVect,tol=6,nObs=1000,statistic=mean,mySeed=355411)
rejection8<-likeFreeRejectionSample2(obsData=dataVect,tol=9,nObs=1000,statistic=mean,mySeed=34756)

plRejection5<-ggplot(data=as.data.frame(rejection5),aes(x=lambda,y=beta))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))+ggtitle("Tol: 0")
plRejection6<-ggplot(data=as.data.frame(rejection6),aes(x=lambda,y=beta))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))+ggtitle("Tol: 3")
plRejection7<-ggplot(data=as.data.frame(rejection7),aes(x=lambda,y=beta))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))+ggtitle("Tol: 6")
plRejection8<-ggplot(data=as.data.frame(rejection8),aes(x=lambda,y=beta))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))+ggtitle("Tol: 9")
grid.arrange(plRejection5,plRejection6,plRejection7,plRejection8, ncol=2)


testSample<-likeFreeRejectionSample2(obsData=dataVect,tol=0,nObs=1000,statistic=mean,mySeed=22144,storeData=TRUE)

str(testSample)

#--------------------------------------------------------------------------------
# b)

# bäst hitintills: tuning=0.05 och tol=1
#a<-likeFreeMCMC(obsData=dataVect,tuning=0.05,statistic=mean,nObs=1000,tol=1,mySeed=28765)
b<-likeFreeMCMC(obsData=dataVect,tuning=0.05,statistic=mean,nObs=1000,tol=1,mySeed=68285)
#d<-likeFreeMCMC(obsData=dataVect,tuning=0.05,statistic=mean,nObs=1000,tol=1,mySeed=227655)
#e<-likeFreeMCMC(obsData=dataVect,tuning=0.05,statistic=mean,nObs=1000,tol=1,mySeed=65334)
#cat("\n\n\n\n\n")

plot(b$parameters[,1],type="l")
plot(b$parameters[,2],type="l")
bPlot<-ggplot(data=as.data.frame(b$parameters[,1:2]),aes(x=lambda,y=beta))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
print(bPlot)

plot(e$parameters[,1],type="l")
plot(e$parameters[,2],type="l")
ePlot<-ggplot(data=as.data.frame(e$parameters[,1:2]),aes(x=lambda,y=beta))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
print(ePlot)


str(b$parameters)

testGibbs<-data.frame(test,Method="Gibbs")
testABCMCMC<-data.frame(b$parameters[,1:2],Method="ABC MCMC")
testRejectionSampler<-data.frame(test2[,1:2],Method="Rejection")
testTotal<-rbind(testGibbs,testABCMCMC,testRejectionSampler)

totLambda<-ggplot(testTotal, aes(x=lambda, color=Method)) + geom_density()+theme_bw()+xlab(expression(lambda))
totBeta<-ggplot(testTotal, aes(x=beta, color=Method)) + geom_density()+theme_bw()+xlab(expression(beta))
grid.arrange(totLambda,totBeta,ncol=2)



testGibbs<-data.frame(test,Method="Gibbs")
testABCMCMC<-data.frame(b$parameters[,1:2],Method="ABC MCMC")
testTotal2<-rbind(testGibbs,testABCMCMC)

pl<-ggplot(data=as.data.frame(testTotal2),aes(x=lambda,y=beta,col=Method))+geom_point()+theme_bw()+geom_density2d()+xlab(expression(lambda)) +ylab(expression(beta))
print(pl)






