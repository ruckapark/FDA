# Import data
setwd("C:/Users/Admin/Documents/R/FDA/FunctionsMFDA")
dfX <- read.csv("./data/xdata_norm.csv",header = FALSE)
dfY <- read.csv("./data/ydata_norm.csv",header = FALSE)
dfZ <- read.csv("./data/zdata_norm.csv",header = FALSE)

# Make functional data objects of x data
library(fda)

#generate new df for x with each time series in vector form
ncurves <- ncol(dfX)
npoints <- nrow(dfX)
ID <- 1:ncurves
t <- vector(mode = "list", length = ncurves)
x <- vector(mode = "list", length = ncurves)

for (i in 1:ncurves){
  col <- colnames(dfX)[i]
  t[i] <- list(1:npoints)
  x[i] <- list(dfX[,col])
}

columns <- c("id","time","xvalue")
df <- data.frame(matrix(nrow = ncurves,ncol = length(columns)))
colnames(df) <- columns
df$id <- ID
df$time <- t
df$xvalue <- x

# FDA parameters
knots    = c(seq(0,npoints,20)) #Location of knots
n_knots   = length(knots) #Number of knots
n_order   = 4 # order of basis functions: for cubic b-splines: order = 3 + 1
n_basis   = length(knots) + n_order - 2;
basis = create.bspline.basis(rangeval = c(0,npoints), n_basis)

#FDA
Lfd <- int2Lfd(2)
fdParam <- fdPar(basis,Lfd,1e4)
df_fd <- smooth.basis(1:200,temp,fdParam)

#visualise
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(df_fd$fd,xlab='time',ylab='Y',cex.lab=1.5,cex.axis=1.5)

#variance and covariance of the functional data + correlation coefficient
library(fields)
yvariance <- var.fd(df_fd$fd)
var_values <- eval.bifd(1:npoints,1:npoints,yvariance)
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
contour(1:npoints,1:npoints,var_values,xlab='time',ylab='Y',cex.lab=1.5,cex.axis=1.5)
#image.plot(1:365,1:365,tvvals,xlab='day',ylab='day',cex.lab=1.5,cex.axis=1.5)

corr = cor.fd(1:npoints,df_fd$fd)
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
contour(1:npoints,1:npoints,corr,xlab='time',ylab='Y',cex.lab=1.5,cex.axis=1.5)

#perform PCA
X_pca <- pca.fd(df_fd$fd,nharm = 2)
names(X_pca)
X_pca$varprop #shows only two are necessary here

#Plot the eigenvalues (proportional to explained variance of each component)
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(X_pca$values,xlab='component',ylab='variance',col="red",
     cex.lab=1.5,cex.axis=1.5,cex=2)
#Plot cumulative percentage explained with added components
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(cumsum(X_pca$values[1:10])/sum(X_pca$values),xlab='Number of Components',
     ylab='cumulative variance explained',col=2,cex.lab=2,
     cex.axis=2,cex=2)
abline(h=0.99)

#Mean curves - meanfd (not sure the lines function is working)
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(df_fd$fd,xlab='time',ylab='Y',cex.lab=1.5,cex.axis=1.5,col=4)
lines(df_fd$meanfd,lwd=2.5,col=2)

#Functional principal components
harmfd = X_pca$harmonics
harmvals = eval.fd(1:npoints,harmfd)
dim(harmvals) # The top 2 FPCs of time dimension

#plot the first FPC
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(1:npoints,harmvals[,1],xlab='time',ylab='PC1',
     lwd=4,lty=1,cex.lab=1,cex.axis=1,type='l',ylim=c(-0.09,0.16),col='blue')
lines(1:npoints,harmvals[,2],xlab='time',ylab='PC2',
     lwd=4,lty=1,cex.lab=1,cex.axis=1,type='l',col='red')

#plot scores plot
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(X_pca$scores[,1:2],xlab='PC Score 1',ylab='PC Score 2',col=4,
     cex.lab=1.5,cex.axis=1.5,cex=1)
text(temppca$scores[,1],temppca$scores[,2],labels=daily$place,cex=1)

# Use funHDDC for MFDA with single with x values
res.uni <- funHDDC(df_fd$fd,K=3,model="AkBkQkDk",init="kmeans",threshold=0.2)
print(res.uni$class)
