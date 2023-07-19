# Import data
setwd("C:/Users/Admin/Documents/R/FDA/FunctionsMFDA")
dfX <- read.csv("./data/xdata_norm.csv",header = FALSE)
dfY <- read.csv("./data/ydata_norm.csv",header = FALSE)
dfZ <- read.csv("./data/zdata_norm.csv",header = FALSE)

library(fda)
