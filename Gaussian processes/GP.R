

#draw<-sampleGP(noDraw=20,xVal=seq(-6,6,length.out=100),covFunc=SqrExpCov,lengthScale=2,sigma=1,func=sin)
#draw<-sampleGP(noDraw=20,xVal=seq(-6,6,length.out=100),covFunc=SqrExpCov,lengthScale=2,sigma=1,func=sin,onlyPara=T)
#round(draw$covMat,4)
#image(draw$covMat)
#plotGP(drawMatrix=draw$draws)

# data:
wages<-read.delim(file="Gaussian processes/CanadianWages2.csv",sep=" ")
# standardize data:
wages$age<-scale(wages$age)
plot(wages$age,wages$logWage,pch=20,ylim=c(10,17),col="blue")
# calculate posterior

# prior:
# let l=0.5 and sigma=0.4 (change those later)
# assume sigmaError to be sd of data
sigmaError<-sd(wages$logWage)
xGrid<-seq(-1.5,2.2,length=30)
draws<-sampleGP(noDraw=30,xVal=xGrid,covFunc=SqrExpCov,func=meanWages,sigma=0.4,lengthScale=0.5)
plotGP(draws$draws)
curve(expr=meanWages,from=-1.7,to=2.3,ylim=c(10,15),add=TRUE,lwd=4,col="red")

# prior mean function:
plot(wages$age,wages$logWage,pch=20,ylim=c(10,17),,xlim=c(-1.7,2.3),col="blue",xlab="age",ylab="log(wages)")
curve(expr=meanWages,from=-1.7,to=2.3,ylim=c(10,15),add=TRUE,lwd=3)

library(lubridate)
# posterior mean
# grid of xValues
xGrid<-seq(-1.5,2.2,length=30)

meanXgrid<-as.matrix(meanFunc(xValues=xGrid,func=meanWages))
xKxGrid<-createCovMatrix(xVal=wages$age,yVal=xGrid,covFunc=SqrExpCov,sigma=0.4,lengthScale=0.5)
sqrtKy<-sqrt(createCovMatrix(xVal=wages$age,covFunc=SqrExpCov,sigma=0.4,lengthScale=0.5)
             +((sigmaError)^2)*diag(1,nrow=dim(wages)[1]))
yDiff<-wages$logWage-meanFunc(xValues=wages$age,func=meanWages)
#fBar<-xKxGrid%*%sqrtKy%*%wages$logWage
fBar<-meanXgrid+xKxGrid%*%sqrtKy%*%wages$logWage
#plot(x=wages$age,y=meanFunc(xValues=wages$age,func=meanWages),type="l")
#plot(x=xGrid,y=meanXgrid,type="l")
graphics.off()

plot(x=xGrid,y=fBar,type="o",main=as.character(now()))


temp3<-createCovMatrix(xVal=xGrid,yVal=wages$age,covFunc=SqrExpCov)
temp2<-createCovMatrix(xVal=wages$age,yVal=xGrid,covFunc=SqrExpCov)
temp1<-createCovMatrix(xVal=wages$age,covFunc=SqrExpCov)
temp4<-createCovMatrix(xVal=xGrid,covFunc=SqrExpCov)
dim(temp2)
dim(temp3)
dim(temp1)
dim(temp4)







