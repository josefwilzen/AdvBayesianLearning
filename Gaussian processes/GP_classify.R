rm(list=ls())
graphics.off()
setwd("/home/joswi05/Dropbox/Josef/Advanced Bayesian Learning/AdvBayesianLearning/Gaussian processes")
source("GPfunc.R")

#dir()
# get data!


# variables:
# sbp-  	systolic blood pressure
# tobacco-		cumulative tobacco (kg)
# ldl-		low densiity lipoprotein cholesterol
# adiposity-
# famhist-		family history of heart disease (Present, Absent)
# typea-		type-A behavior
# obesity-
# alcohol-		current alcohol consumption
# age-		age at onset
# chd-		response, coronary heart disease


#----------------------------------------------------
# Part a)
#----------------------------------------------------
library(ggplot2)
heart<-read.csv(file="SAheart.csv")
heartGLM<-heart
heart$chd<-ifelse(heart$chd==0,-1,1)  # change coding of the binary respose to -1/1


head(heart)
heart2<-heart
heart2$chd2<-as.factor(heart2[,11])

names(heart)
testData1<-heart[,c(8,10,11)]
qplot(x=age,y=obesity, color=chd2,data=heart2)

#----------------glm------------------------------------


modelA<-glm(formula=chd~age+obesity,data=heartGLM,family=binomial(link = "logit"))
summary(modelA)
modelA$coefficients
names(heart)
heart
as.matrix(modelA$coefficients)

hist(heart2$age,30)
hist(heart2$obesity,30)
quantile(heart2$obesity,probs=seq(0,1,0.1))

set.seed(123144)
a<-sample(heart2$obesity,size=100,prob=density(heart2$obesity,n=462)$y)
b<-sample(heart2$age,size=100,prob=density(heart2$age,n=462)$y)
xyz<-expand.grid(a,b)
xyz<-cbind(constant=1,xyz,z=0)
names(xyz)[1:2]<-c("x","y")
for(i in 1:dim(xyz)[1]){
  xyz[i,4]<-as.matrix(xyz[i,1:3])%*%as.matrix(modelA$coefficients)
}
image(x=xyz[,2],y=xyz[,3],xyz[,4])
modelA$coefficients
binomial()
heart2<-heart[,c()]
library(reshape2) # for melt
volcano3d <- melt(volcano)
names(volcano3d) <- c("x", "y", "z")

# Basic plot
v <- ggplot(volcano3d, aes(x, y, z = z))
v + stat_contour()
melt
heart[,1:2]

g<-ggplot(heart2,mapping=aes(age,obesity))
p<-g+geom_point(aes(color=chd))+theme_bw()
p


#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# Part b)
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------


#install.packages("optimx")
#library(optimx)
#optimx



system.time(test1Theta<-covMatrixVect(theta=c(0,1,1,1),myData=testData1[,1:2]))
#system.time(test1<-covMatrixVect(theta=c(1,1,1),myData=testData1[,1:2]))
#system.time(test2<-covMatrixVect2(theta=c(1,1,1,1),myData=testData1[,1:2]))
system.time(test2Theta<-covMatrixVect2(theta=c(0,1,1,1),myData=testData1[,1:2]))
all(test1Theta==test2Theta)


#-----------------------------------------------------------------------------------
# laplace approx:
#-----------------------------------------------------------------------------------
install.packages("VGAM")
library(VGAM)
library(boot)
logit
inv.logit(x=)
str(binomial(link="logit"))
testX<-t(as.matrix(testData1[5,1:2]))
theta=c(1,0,1,1)
SqrExpCovVect(theta=c(1,0,1,1),x1=testX,x2=testX)
theta
mu<-0

pnorm(0,mean=mu,sd=sum(theta[1:2]))


round(inv.logit(as.matrix(testData1[,1:2])%*%as.matrix(c(-3,3))),4)
dim(t(testData1[,1:2]))
modelA$coefficients
inv.logit(as.matrix(cbind(1,testData1[,1:2]))%*%modelA$coefficients)
fitted(modelA)
str(modelA)
glm


# log(likelihoodGPclassify(1,y=1,LOG=FALSE))
# log(likelihoodGPclassify(1,y=-1,LOG=FALSE))
# likelihoodGPclassify(1,y=1,LOG=TRUE)
# likelihoodGPclassify(1,y=-1,LOG=TRUE)

#-----------------------------------------------------------------------------------
# posterior

#install.packages("Matrix")
library(Matrix)

a<-logisticHessian(f=fitted(object=modelA)[1:10],y=heart$chd[1:10])
logisticGradient(f=fitted(object=modelA)[1:10],y=heart$chd[1:10])

class(a)
a@x
attributes(a)
b<-sqrt(a)%*%sqrt(a)
all(b==a)
sqrt(diag(a))
a<-diag(9,nrow=3)
sqrt(a)%*%sqrt(a)
d<-sqrt(a)*sqrt(a)
a<-as.matrix(testData1[,1:2])
a<-Matrix(as.matrix(testData1[,1:2]))
class(a)
testNew<-covMatrixVect2(theta=c(1,1,1,1),myData=testData1[,1:2])

a.eig <- eigen(a)
a.sqrt <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
Matrix(data=)
all(abs(a@x-sqrt(a@x)*sqrt(a@x))<1e-16)

a<-covMatrixVect2(theta=c(1,1,1,1),myData=testData1[,1:2])

h<-matrix(data=sample(100,size=9),3,3)
h<-chol(a)
h
str(h)


# y<-testData1$chd
# K<-covMatrixVect2(theta=c(1,1,1,1),myData=testData1[,1:2])

# noObs<-length(y)
# f<-rep(0,noObs)
# W<-logisticHessian(f=f,y=y)
# sqrtW<-Diagonal(x=sqrt(W@x))
# L<-chol(diag(noObs)+sqrtW%*%K%*%sqrtW)
# 
# fGrad<-logisticGradient(f=f,y=y)
# b<-W%*%f+fGrad
# sqrtWKb<-sqrtW%*%K%*%b
# LsolveSqrtWKb<-solve(a=L,b=sqrtWKb)
# WsqrtLt<-sqrtW%*%t(L)
# a<-b-solve(a=WsqrtLt,b=LsolveSqrtWKb)
# f<-K%*%a
# 
# str(diag(L))
# sum(log(diag(L)))

logPosterior<-function(f,y,K){
  logPost<-sum(likelihoodGPclassify(f=f,y=y))-0.5*t(f)%*%solve(K)%*%f-0.5*det(x=K,log=TRUE)-0.5*length(y)*log(2*pi)  
  return(as.vector(logPost))
}

posteriorMode<-function(K,y,nIter=100,tol=c(1e-20,1e-20),f=NULL,adder=0,step=1,...){
  noObs<-length(y)
  objectiveFunc<-rep(x=0,length.out=nIter)
  if(is.null(f)) f<-rep(0,noObs)
  objectiveFunc[1]<-logPosterior(f=f,y=y,K=K)
  stopCondition="iter"
  K<-K+adder*diag(noObs)
  
  for(i in 2:nIter){
    print(i)
    W<-logisticHessian(f=f,y=y)
    sqrtW<-Diagonal(x=sqrt(W@x))
    
    # try if the Choleski factorization shall be inverted with respect to t()s
    L<-chol(Diagonal(noObs)+as.matrix(sqrtW%*%K%*%sqrtW))
    L<-t(L)
    fGrad<-logisticGradient(f=f,y=y)
    b<-W%*%f+fGrad*step
    sqrtWKb<-sqrtW%*%K%*%b
    print(paste("a: Class of L:   ",as.character(class(L))," rcond:  ",rcond(x=as.matrix(L))))
    print(paste("b: Class of sqrtWKb:   ",class(sqrtWKb)," rcond:  ",rcond(sqrtWKb)))
    LsolveSqrtWKb<-solve(a=L,b=sqrtWKb)
    #LsolveSqrtWKb<-solve(a=L,b=as.vector(sqrtWKb),system="L")
    WsqrtLt<-sqrtW%*%t(L)
    print(paste("a: Class of WsqrtLt:   ",class(WsqrtLt)," rcond:  ",rcond(WsqrtLt)))
    print(paste("b: Class of LsolveSqrtWKb:   ",class(LsolveSqrtWKb)," rcond:  ",rcond(LsolveSqrtWKb)))
    a<-b-solve(a=WsqrtLt,b=LsolveSqrtWKb)
    #a<-b-solve(a=WsqrtLt,b=as.vector(LsolveSqrtWKb),system="L")
    fnew<-K%*%a
    objectiveFunc[i]<-logPosterior(f=fnew,y=y,K=K)
    #if(i==4) browser()
    #browser()
    if(i>10&abs(objectiveFunc[i]-objectiveFunc[(i-1)])<tol[1]){
      stopCondition<-"objectiveFunc"
      break
    }else if(sum(abs(fnew-f))<tol[2]){
      stopCondition<-"f"
      break
    }
    f<- fnew
  }
  logMargLike<- -0.5*t(a)%*%f+sum(likelihoodGPclassify(f=f,y=y))-sum(log(diag(L)))
  res<-list(f=as.vector(f),logMargLike=logMargLike,postMode=objectiveFunc[1:i],stopCondition=stopCondition)
  return(res)
  #return(f)
}


y<-testData1$chd
K<-covMatrixVect2(theta=c(1,1,1,1),myData=testData1[,1:2])
#hej<-posteriorMode(K=K,y=y,nIter=5,f=fitted(modelA))
hej<-posteriorMode(K=K,y=y,nIter=100,tol=c(1e-100,1e-10),adder=0,f=NULL,step=2)
hej
hej$postMode[11]-hej$postMode[10]
round(hej,3)
"%*%"
chol

data(iris)
iris2<-iris[,c(1,2,5)]
iris2[,3]<-ifelse(iris2[,3]=="setosa",1,-1)

Kiris<-covMatrixVect2(theta=c(1,2,1,1,1,1),myData=iris[,1:4])
irisTest<-posteriorMode(K=Kiris,y=iris2[,3],nIter=100,tol=c(1e-100,1e-100),adder=0,f=NULL)

#-----------------------------------------------------------------------------------
# marginal likelihood

logMarginalLikelihood<-function(theta,y,X,fMode,...){
  W<-logisticHessian(f=f,y=y)
  sqrtW<-Diagonal(x=sqrt(W@x))
  det(Diagonal(noObs)+as.matrix(sqrtW%*%K%*%sqrtW),log=TRUE)
  return()
}




















