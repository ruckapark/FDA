# Import data
setwd("C:/Users/Admin/Documents/R/FDA/FunctionsMFDA")
source("univPCA_function.R")

dfX <- read.csv("./data/E_X_m_data.csv",header = FALSE)
dfY <- read.csv("./data/G_Y_m_data.csv",header = FALSE)
dfZ <- read.csv("./data/R_Z_m_data.csv",header = FALSE)

breakpoint_spacing <- 30
npoints = breakpoint_spacing * (dim(dataXYZ)[1] %/% breakpoint_spacing)

#remove unecessary rows
dfX <- dfX[1:npoints,]
dfY <- dfY[1:npoints,]
dfZ <- dfZ[1:npoints,]

library(fda)

dataX <- data.matrix(dfX)
dataY <- data.matrix(dfY)
dataZ <- data.matrix(dfY)

#3D array with dX and dY
dataXYZ <- array(c(dataX,dataY,dataZ), dim=c(dim(dataX)[1],dim(dataX)[2],3))

#FDA parameters
breakpoint_spacing <- 30
npoints <- dim(dataXYZ)[1]
time <- 1:npoints
knots <- c(seq(0,npoints,breakpoint_spacing))
n_knots <- length(knots)
n_order <- 4 #BSpline order
n_basis <- length(knots) + n_order - 2

dev.new()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(time,dataXYZ[,,1],type='l',cex.lab=1.5,cex.axis=1.5,ylab='Xdata')

dev.new()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(time,dataXYZ[,,2],type='l',cex.lab=1.5,cex.axis=1.5,ylab='Ydata')

#Set up general fd object
basis = create.bspline.basis(rangeval = c(0,npoints), n_basis)
Lfd <- int2Lfd(2)
fdParam <- fdPar(basis,Lfd,1e4)
df_fda = smooth.basis(1:npoints,dataXYZ,fdParam)

fd <- df_fda$fd
fd$fdnames[[3]] <- c("Xdata","Ydata","Zdata")
fdvals <- eval.fd(time,fd)

#Plot X values from smoothed data (compare to unsmoothed) (why not in 1D)
dev.new()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(fdvals[,,1],fdvals[,,2],type='l',xlab='time',cex.lab=1.5,cex.axis=1.5)

#Plot just the first data
dev.new()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(array(time,dim = c(870,12)),fdvals[,,1],type='l',xlab='time',cex.lab=1.5,cex.axis=1.5)

#pca multidim
XYZ.pca = pca.fd(fd,nharm = 3)

#plot mean for the three dimensions
dev.new()
par(mfrow=c(3,1),mar = c(8, 8, 4, 2))
plot(XYZ.pca$meanfd,cex.lab=1.5,cex.axis=1.5)

#plot variance proportion
dev.new()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(XYZ.pca$values,cex.lab=1.5,cex.axis=1.5,xlab='PC',ylab='eigenvalue',col=4,cex=2)
XYZ.pca$varprop
sum(XY.pca$varprop)

#plot variance proportion
dev.new()
par(mfrow=c(3,1),mar = c(8, 8, 4, 2))
plot(XYZ.pca$harmonics,lwd=2,cex.lab=1.5,cex.axis=1.5)

#multidimensional clustering - use list of functional objects
library(funHDDC)
xfd <- univPCA(dfX)$fda
yfd <- univPCA(dfY)$fda
zfd <- univPCA(dfZ)$fda
res.uni <- funHDDC(list(xfd,yfd,zfd),K=3,model="AkBkQkDk",init="kmeans",threshold=0.2)
print(res.uni$class)