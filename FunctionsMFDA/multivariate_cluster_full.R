# Import data
setwd("C:/Users/Admin/Documents/R/FDA/FunctionsMFDA")
source("univPCA_function.R")

dfX <- read.csv("./data/xdata_norm.csv",header = FALSE)
dfY <- read.csv("./data/ydata_norm.csv",header = FALSE)
#Z data is in the wrong format

library(fda)
?gait #show gait data

temp <- gait

dataX <- data.matrix(dfX)
dataY <- data.matrix(dfY)

#3D array with dX and dY
dataXY <- array(c(dataX,dataY), dim=c(200,20,2))

#FDA parameters
breakpoint_spacing <- 20
npoints <- 200
time <- 1:npoints
knots <- c(seq(0,npoints,breakpoint_spacing))
n_knots <- length(knots)
n_order <- 4 #BSpline order
n_basis <- length(knots) + n_order - 2

dev.new()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(time,dataXY[,,1],type='l',cex.lab=1.5,cex.axis=1.5,ylab='Xdata')

dev.new()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(time,dataXY[,,2],type='l',cex.lab=1.5,cex.axis=1.5,ylab='Ydata')

#Set up general fd object
basis = create.bspline.basis(rangeval = c(0,npoints), n_basis)
Lfd <- int2Lfd(2)
fdParam <- fdPar(basis,Lfd,1e4)
df_fda = smooth.basis(1:npoints,dataXY,fdParam)

fd <- df_fda$fd
fd$fdnames[[3]] <- c("Xdata","Ydata")
fdvals <- eval.fd(time,fd)

#Plot X values from smoothed data (compare to unsmoothed) (why not in 1D)
dev.new()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(fdvals[,,1],fdvals[,,2],type='l',xlab='time',cex.lab=1.5,cex.axis=1.5)

#Plot just the first data
dev.new()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(array(time,dim = c(200,20)),fdvals[,,1],type='l',xlab='time',cex.lab=1.5,cex.axis=1.5)

#pca multidim
XY.pca = pca.fd(fd,nharm = 3)

#plot mean for the two dimensions
dev.new()
par(mfrow=c(2,1),mar = c(8, 8, 4, 2))
plot(XY.pca$meanfd,cex.lab=1.5,cex.axis=1.5)

#plot variance proportion
dev.new()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(XY.pca$values,cex.lab=1.5,cex.axis=1.5,xlab='PC',ylab='eigenvalue',col=4,cex=2)
XY.pca$varprop
sum(XY.pca$varprop)

#plot variance proportion
dev.new()
par(mfrow=c(2,1),mar = c(8, 8, 4, 2))
plot(XY.pca$harmonics,lwd=2,cex.lab=1.5,cex.axis=1.5)

#multidimensional clustering - use list of functional objects
library(funHDDC)
xfd <- univPCA(dfX)$fda
yfd <- univPCA(dfY)$fda
res.uni <- funHDDC(list(xfd,yfd),K=4,model="AkBkQkDk",init="kmeans",threshold=0.2)
print(res.uni$class)