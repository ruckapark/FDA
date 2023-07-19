library(fda)
library(funHDDC)

univPCA <-function(df,n_breakpoints=20,spline_order=4,Lfd_order = 2,fPC_no = 2) {
  
  #define no. of points and timeseires
  ncurves <- ncol(df)
  npoints <- nrow(df)
  
  # FDA parameters
  knots    = c(seq(0,npoints,n_breakpoints)) #Location of knots
  n_knots   = length(knots) #Number of knots
  n_order   = spline_order # order of basis functions: for cubic b-splines: order = 3 + 1
  n_basis   = length(knots) + n_order - 2;
  basis = create.bspline.basis(rangeval = c(0,npoints), n_basis)
  
  #FDA
  Lfd <- int2Lfd(Lfd_order)
  fdParam <- fdPar(basis,Lfd,1e4)
  df_fda <- smooth.basis(1:npoints,data.matrix(df),fdParam)
  fd <- df_fda$fd
  
  #PCA
  df_PCA <- pca.fd(df_fda$fd,nharm = fPC_no)
  
  return(list(fda = fd,PCA = df_PCA))
  
}

plot_smoothed <-function(fd,ylab = 'Y'){
  par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
  plot(fd,xlab='time',ylab=ylab,cex.lab=1.5,cex.axis=1.5)
}

plotPCA_comps <-function(df_PCA) {
  
  harmfd = df_PCA$harmonics
  harmvals = eval.fd(1:npoints,harmfd)
  dim(harmvals)
  
  #par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
  plot(1:npoints,harmvals[,1],xlab='time',ylab='PC1',
       lwd=4,lty=1,cex.lab=1,cex.axis=1,type='l',ylim=c(-0.09,0.16),col='blue')
  lines(1:npoints,harmvals[,2],xlab='time',ylab='PC2',
        lwd=4,lty=1,cex.lab=1,cex.axis=1,type='l',col='red')
  
}

plotPCA_scores <-function(df_PCA){
  #plot scores plot
  par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
  plot(df_PCA$scores[,1:2],xlab='PC Score 1',ylab='PC Score 2',col=4,
       cex.lab=1.5,cex.axis=1.5,cex=1)
}

#threshold is what?
univCluster <- function(fd_obj,k_clus,mod="AkBkQkDk"){
  
  res.uni <- funHDDC(fd_obj,K=k_clus,model=mod,init="kmeans",threshold=0.2)
  
  return(res.uni)
  
}