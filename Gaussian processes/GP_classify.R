rm(list=ls())
graphics.off()
setwd("/home/joswi05/Dropbox/Josef/Advanced Bayesian Learning/AdvBayesianLearning/Gaussian processes")
source("GPfunc.R")

dir()
# get data!
heart<-read.csv(file="SAheart.csv")

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
head(heart)
heart2<-heart
heart2$chd2<-as.factor(heart2[,11])

names(heart)
testData1<-heart[,c(8,10,11)]
qplot(x=age,y=obesity, color=chd2,data=heart2)

#----------------glm------------------------------------
modelA<-glm(formula=chd~age+obesity,data=heart,family=binomial(link = "logit"))
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



#----------------------------------------------------
# Part b)
#----------------------------------------------------


#install.packages("optimx")
#library(optimx)
#optimx

SqrExpCovVect<-function(theta,x1,x2){
  sigma<-theta[1]
  l<-theta[2:length(theta)]
  x1<-as.matrix(x1)
  x2<-as.matrix(x2)
  M<-diag(l^(-2))
  y<-sigma^2*exp(-0.5*t(x1-x2)%*%M%*%(x1-x2))
  # possiblity to add the sigma model term above
  return(y)
}


SqrExpCovVect(theta=c(4,2,3),x1=testData1[1,1:2],x2=testData1[2,1:2])
SqrExpCovVect(theta=c(4,2,3),x1=c(1,2),x2=c(2,8))

SqrExpCov<-function(x1,x2,sigma=1,lengthScale=1){
  if(lengthScale<=0) stop("Values of lengthScale not allowed")
  if(sigma<=0) stop("Values of sigma not allowed")
  covVal<-(sigma^2)*exp(((x1-x2)^2)/(-2*lengthScale^2))
  return(covVal)
}



createCovMatrix<-function(xVal,yVal=NULL,covFunc,...){

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





