rm(list=ls())
graphics.off()
setwd("/home/josef/Dropbox/Josef/Advanced Bayesian Learning/AdvBayesianLearning/Gaussian processes/")
dir()
source("GPfunc.R")
# data:
wages<-read.delim(file="Gaussian processes/CanadianWages2.csv",sep=" ")
# standardize data:
wages$age<-scale(wages$age)
wages$logWage<- wages$logWage-mean(wages$logWage)
plot(wages$age,wages$logWage,pch=20,col="blue",xlab="age",ylab="log(wages)")


library(soma)
#install.packages("Rcpp")












