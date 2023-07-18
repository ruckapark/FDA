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

# Use funHDDC for MFDA