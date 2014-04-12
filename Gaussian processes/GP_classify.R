rm(list=ls())
graphics.off()
setwd("/home/joswi05/Dropbox/Josef/Advanced Bayesian Learning/AdvBayesianLearning/Gaussian processes")
source("GPfunc.R")

dir()
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
heart$chd<-ifelse(heart$chd==0,-1,1)  # change coding of the binary respose to -1/1


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




system.time(test1<-covMatrixVect(theta=c(1,1,1,1),myData=testData1[,1:2]))
system.time(test2<-covMatrixVect2(theta=c(1,1,1,1),myData=testData1[,1:2]))
all(test1==test2)





