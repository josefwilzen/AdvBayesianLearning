

#install.packages("reshape")
library(reshape)


# covFunc- covariance function
# para - list with all parameters for covFunc on the form:
#   xVal - x values to compute covariance for
#   yVal - if cross correlations shall be computed otherwise NULL
#   addional hyperparamters for covFunc

createCovMatrix2<-function(covFunc,para){
  # index<-1:length(xVal)
  #   indexMat<-as.matrix(expand.grid(index,index))
  #   covValues<-covFunc(x1=indexMat[,1],x2=indexMat[,2],...)
  #   covMat<-matrix(0,length(xVal),length(xVal))
  #   covMat[indexMat]<-covValues
  xVal<-para$xVal
  yVal<-para$yVal
  if(is.null(yVal)){
    #browser()
    rows<-length(xVal)
    y_index <- rev(abs(sequence(seq.int(rows)) - rows) + 1)
    x_index <- rep.int(seq.int(rows), rev(seq.int(rows)))
    xy_index<-cbind(x_index,y_index)
    new_x1<-xVal[xy_index[,1]]
    new_x2<-xVal[xy_index[,2]]
    para_func<-c(list(x1=new_x1,x2=new_x2),para[-(1:2)])
    cov_values<-do.call(covFunc,para_func)
    covMat<-matrix(0,rows,rows)
    covMat[xy_index]<-cov_values
    covMat2<-covMat
    covMat2<-covMat2+t(covMat)
    diag(covMat2)<-diag(covMat2)-diag(covMat)
    return(covMat2)
  }else {
    yl<-length(yVal)  # no of rows
    xl<-length(xVal)  # no of cols
    covMat<-matrix(0,nrow=yl,ncol=xl)
    if(yl>=xl){  # loop over cols
      for(i in 1:xl){
        x_index <-rep(i,xl) 
        y_index <-1:yl
        yx_index<-cbind(y_index,x_index)
        new_y<-xVal[yx_index[,1]]
        new_x<-yVal[yx_index[,2]]
        para_func<-c(list(x1=new_x,x2=new_y),para[-(1:2)])
        covMat[,i]<-do.call(covFunc,para_func)
      }
    }else{  # loop over rows
      for(i in 1:yl){
        x_index <-rep(i,xl) 
        y_index <-1:yl
        xy_index<-cbind(x_index,y_index)
        new_x<-xVal[xy_index[,1]]
        new_y<-yVal[xy_index[,2]]
        para_func<-c(list(x1=new_x,x2=new_y),para[-(1:2)])
        covMat[i,]<-do.call(covFunc,para_func)
      }
    }
#     for(i in 1:length(yVal)){  # rows 
#       for(j in 1:length(xVal)){ # cols
#         covMat[i,j]<-covFunc(x1=yVal[i],x2=xVal[j],...)  # do.call(covFunc,para)
#       } 
#     }
#     colnames(covMat)<-paste("xVal",round(xVal,3))
#     rownames(covMat)<-paste("yVal",round(yVal,3))
#     return(covMat)
  }
}


SqrExpCov2<-function(x1,x2,sigma=1,lengthScale=1){
  if(lengthScale<=0) stop("Values of lengthScale not allowed")
  if(sigma<=0) stop("Values of sigma not allowed")
  covVal<-(sigma^2)*exp(((x1-x2)^2)/(-2*lengthScale^2))
  return(covVal)
}

rows<-2000
para<-list(xVal=1:rows,yVal=NULL,lengthScale = 1,sigma = 1)
library(microbenchmark)
para1<-list(covFunc = SqrExpCov,para = para)
para2<-list(xVal = 1:rows,yVal = NULL,covFunc = SqrExpCov,lengthScale = 1,sigma = 1)
compare<-microbenchmark(do.call(createCovMatrix2,para1),do.call(createCovMatrix,para2),times = 10,unit = "s")
compare<-microbenchmark(do.call(createCovMatrix2,para1),do.call(createCovMatrix,para2),times = 10)
compare


system.time(a<-createCovMatrix2(covFunc = SqrExpCov,para = para))
system.time(b<-createCovMatrix(xVal = 1:rows,yVal = NULL,covFunc = SqrExpCov,lengthScale = 1,sigma = 1))
all(a==b)
Matrix::image(Matrix(a),colorkey=TRUE)
det(a)
isSymmetric(a)
library(MASS)
plot(mvrnorm(n=1,mu=10*cos(c(1:rows)^2),Sigma=1*a),type="l")

SqrExpCov2(x1 =1:10 ,x2 = rep(1,10),lengthScale = 1,sigma = 1)

rows <- 10 # number of rows and columns (20)

y <- rev(abs(sequence(seq.int(rows - 1)) - rows) + 1)
x <- rep.int(seq.int(rows - 1), rev(seq.int(rows - 1)))

y2 <- rev(abs(sequence(seq.int(rows)) - rows) + 1)
x2 <- rep.int(seq.int(rows), rev(seq.int(rows)))

#idx <- cbind(x, y)
idx<-cbind(x2,y2)
#idx <- cbind(y, x)

a<-1:rows
b<-SqrExpCov2(x1 = a[idx[,1]],x2 = a[idx[,2]],sigma = 1,lengthScale = 1)

d<-matrix(0,rows,rows)
d
idx
d[idx]<-b
d
Matrix::image(Matrix(d),colorkey=TRUE)
#Matrix::image(Matrix(d)[1:20,1:20],colorkey=TRUE)
Matrix(d)
#d[upper.tri(d)]
t(upper.tri(d))
e<-d
e<-e+t(d)
diag(e)<-diag(e)-diag(d)
round(e,3)
d[t(upper.tri(d))]
d[t(upper.tri(d))]<-d[upper.tri(d)]
round(d,3)
Matrix::image(Matrix(e),colorkey=TRUE)
isSymmetric(e)
det(e)
round(d[1:10,1:10],3)






