# Import data
setwd("C:/Users/Admin/Documents/R/FDA/FunctionsMFDA")
source("univPCA_function.R")

dfX <- read.csv("./data/xdata_norm.csv",header = FALSE)
dfY <- read.csv("./data/ydata_norm.csv",header = FALSE)
dfZ <- read.csv("./data/zdata_norm.csv",header = FALSE)

#Perform PCA for each of the data frames - could put this in a function
pca.X <- univPCA(dfX)
plot_smoothed(pca.X$fda)
plotPCA_comps(pca.X$PCA)
plotPCA_scores(pca.X$PCA)

#Depending on PCA plot, decide on number of clusters
res.X <- univCluster(pca.X$fda,3)
print(res.X$class)

#repeat for Y
pca.Y <- univPCA(dfY)
plot_smoothed(pca.Y$fda)
plotPCA_comps(pca.Y$PCA)
plotPCA_scores(pca.Y$PCA)
res.Y <- univCluster(pca.Y$fda,3)

#repeat for z - should be two distinct pulse clusters
#the cluster works, however it does not succeed in capture the shortest pulse signals

#there should always be a visual check on the goodness of fit of the splines
pca.Z <- univPCA(dfZ,breakpoint_spacing = 5,base = 'f')
plot_smoothed(pca.Z$fda)
plotPCA_comps(pca.Z$PCA)
plotPCA_scores(pca.Z$PCA)
res.Z <- univCluster(pca.Z$fda,2)

#Print class results for three dimensions
classes <- data.frame(Xclasses = res.X$class, Yclasses = res.Y$class, Zclasses = res.Z$class)
print(classes)


#Perform multivariate PCA