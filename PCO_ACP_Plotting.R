## This script contains the functions for performing 2-D PCO given a dissimilarity matrix as well as the ACP 

#var.cols <- c(1:20)
library(colorspace)

# Predefined dark colors
var.cols <- c("darkred", "darkblue", "darkgreen", "darkcyan", "darkmagenta", "darkorange", "darkviolet", "darkgoldenrod", "darkslategray")


pco.fun <- function(dissim.mat, ACP=FALSE){

  if(ACP){
    tmp.pco.out <- pco.fun(dissim.mat=dissim.mat, ACP = FALSE)
    c0 <- -2*rev(tmp.pco.out$eig.values)[1]
    A.mat <- -0.5*(dissim.mat^2+ c0)
  } else {
    A.mat <- -0.5*dissim.mat^2
  }

  col.centers <- colMeans(A.mat)
  row.centers <- colMeans(t(A.mat))
  center <-  mean(A.mat)
  centred_mat <- apply(apply(A.mat,1,function(x)x-col.centers),
                       1,function(x)x-row.centers)+center
  
  
  
  eig.dat <- eigen(centred_mat)
  # The coordinates for the factual observations
  coords <- eig.dat$vectors[,1:2]%*%diag(eig.dat$values[1:2]^0.5)#%*%(rot.mat)
  
  if(ACP){
    mds.obj=list(dissim.mat=dissim.mat, col.centers=col.centers, center=center, coords=coords, eig.values=eig.dat$values, n=nrow(coords),eig.vectors=eig.dat$vectors,ACP=c0)
  } else{
    mds.obj=list(dissim.mat=dissim.mat, col.centers=col.centers, center=center, coords=coords, eig.values=eig.dat$values, n=nrow(coords),eig.vectors=eig.dat$vectors)
  }
  
  return(mds.obj)
}

rot.plot.fun <- function(mds.obj,rot.angle=0){
  rot.mat <- matrix(c(cos(rot.angle/180*pi),sin(rot.angle/180*pi),-sin(rot.angle/180*pi),cos(rot.angle/180*pi)), ncol=2)
  mds.obj$coords <- mds.obj$coords%*%rot.mat
  return(mds.obj)
}


refl.plot.fun <- function(mds.obj,x.refl=FALSE,y.refl=FALSE){
  if(!(x.refl)& !(y.refl)){
    return(mds.obj)    
  }
    
  if(x.refl & y.refl){
    refl.angle <- 45
  } else if(x.refl){
    refl.angle <- 90
  } else if(y.refl){
    refl.angle <- 180
  }
  refl.mat <- matrix(c(cos(2*refl.angle/180*pi),
                       sin(2*refl.angle/180*pi),
                       sin(2*refl.angle/180*pi),
                       -cos(2*refl.angle/180*pi)), ncol=2)
  mds.obj$coords <- mds.obj$coords%*%refl.mat
  
  return(mds.obj)
}

add.labels <- function(mds.obj,labels.vec=1:mds.obj$n, pos.vec=rep(1,mds.obj$n)){
  if(is.null(mds.obj$labels.list)){
    ind <- 1+length(mds.obj)
    mds.obj[[ind]] <- list(labels.vec=labels.vec,pos.vec=pos.vec)
    names(mds.obj)[ind] <- "labels.list"
  } else {
    ind <- which(names(mds.obj)=="labels.list")
    mds.obj[[ind]] <- list(labels.vec=labels.vec,pos.vec=pos.vec)
    names(mds.obj)[ind] <- "labels.list"
  }
  return(mds.obj)
}

mds.plot <- function(mds.obj, col.vec=NULL, exp.factor=NULL){

  coords <- mds.obj$coords
  
  if(is.null(col.vec)){
    col.vec <- rep(1,nrow(coords))
  }
  if(is.null(exp.factor)){
    graphics::plot(coords,asp=1,
         xlab="",#paste("Dim 1 (",round(eig.dat$values[1]/sum(eig.dat$values)*100,2),"%)",sep=""),
         ylab="",#paste("Dim 2 (",round(eig.dat$values[2]/sum(eig.dat$values)*100,2),"%)",sep=""),
         pch=16, 
         axes =FALSE,
         frame.plot=TRUE, 
         col=col.vec
    )
  } else {
    graphics::plot(coords,asp=1,
                   xlab="",#paste("Dim 1 (",round(eig.dat$values[1]/sum(eig.dat$values)*100,2),"%)",sep=""),
                   ylab="",#paste("Dim 2 (",round(eig.dat$values[2]/sum(eig.dat$values)*100,2),"%)",sep=""),
                   pch=16, 
                   axes =FALSE,
                   frame.plot=TRUE,ylim=range(coords)+exp.factor*c(-1,1)*diff(range(coords)),xlim=range(coords)+exp.factor*c(-1,1)*diff(range(coords)),
                   col=col.vec
    ) 
    }
  if(!is.null(mds.obj$labels.list)){
    labels.vec <- mds.obj$labels.list$labels.vec
    pos.vec <- mds.obj$labels.list$pos.vec
    text(coords,labels = labels.vec, pos=pos.vec,cex=1.2)
  }
  invisible(mds.obj)
}

#pco.fun(dissim.mat, ACP=TRUE)|> 
#  add.labels(labels.vec = rownames(dissim.mat) )|>
#  mds.plot()
#dev.new()
#pco.fun(dissim.mat, ACP=FALSE)|> 
#  add.labels(labels.vec = rownames(dissim.mat) )|>
#  refl.plot.fun(y.refl=TRUE,x.refl=FALSE)|>
#  rot.plot.fun(rot.angle = -35)|>
#  mds.plot()



###############################################################################
##  Predictive Axes
###############################################################################


# Adding linear axes to MDS config with linear regression.
#     Takes a centered MDS configuration
#     Takes a uncentered dataset with column names
mds.linreg <- function(Z,dat, Class=1, Adj=c(0.5,-0.5), samp.cols=NULL, Var.cols=rep("grey",ncol(dat)), axes.labels=TRUE){
  Y <- as.matrix(dat)
  centering.mat <- matrix(apply(Y,2,mean),byrow = TRUE, nrow=nrow(Y), ncol=ncol(Y))
  Y.cen <- Y - centering.mat
  A <- solve(t(Z)%*%Z)%*%t(Z)%*%Y.cen
  pred.data <- Z%*%A+centering.mat
  
  if(is.null(samp.cols)){samp.cols <- Class}
  colnames_dat <- colnames(dat)
  n <- nrow(dat)
  G <- Z
  G.lim <- max(abs(min(G)),abs(max(G)))+0.01
  H <- t(A)
  plot(G, asp=1, col=samp.cols, pch=16, ann=FALSE, 
       axes=FALSE,ylim=c(-G.lim,G.lim),xlim=c(-G.lim,G.lim),
       frame.plot=TRUE, xpd=TRUE)
  
  
  #apply(H,1, function(x)abline(a=0,b=x[2]/x[1],col='gray'))
  sapply(1:nrow(H), function(x)abline(a=0,b=H[x,2]/H[x,1],col=Var.cols[x]))
  
  for (i in 1:ncol(dat)){
    marker.val <- pretty(dat[,i])
    c.val <- (marker.val-mean(dat[,i]))/sum(H[i,]^2)
    c.lim <- (G.lim-mean(dat[,i]))/sum(H[i,]^2)
    tick.coords<- t(sapply(c.val,function(x)c(H[i,1]*x,H[i,2]*x)))
    
    #   points(tick.coords, pch=16, col='gray')
    
    # Draw tick marks perpendicular to the axes
    # Steepness = S, calculated as minus change in x over change in y
    S <- -H[i,1]/H[i,2]
    # Tick length = tl
    tl <- 0.01*(par('usr')[2]-par('usr')[1])#0.05
    # These were calculated off of the pythogorean theorem with 
    #     tl=sqrt(x.len^2+y.len^2) where x.len/y.len=S
    x.len <- tl/sqrt(S^2+1)
    y.len <- S*x.len
    apply(tick.coords,1,function(x)lines(x=c(-x.len,x.len)+x[1],
                                         y=x[2]+c(-y.len,y.len), col=Var.cols[i]))
    
    text(tick.coords,adj=Adj,labels=marker.val,col=Var.cols[i])
    #    par(xpd=TRUE)
    
    if(abs(tick.coords[1,1]-tick.coords[2,1])-0.03
       < abs(tick.coords[1,2]-tick.coords[2,2])){
      # side 1 or 3
      index <- 2
      if(tick.coords[1,2]>tick.coords[2,2]){
        sidez <- 1
        ADJS <- (G.lim-H[i,1]/H[i,2]*G.lim)/(2*G.lim)
      } else {
        sidez <- 3
        ADJS <- (G.lim+H[i,1]/H[i,2]*G.lim)/(2*G.lim)
      }
    }else {
      index<- 2
      if(tick.coords[1,1]>tick.coords[2,1]){
        sidez <- 2
        ADJS <- (G.lim-H[i,2]/H[i,1]*G.lim)/(2*G.lim)
      } else {
        sidez <- 4
        ADJS <- (G.lim+H[i,2]/H[i,1]*G.lim)/(2*G.lim)
      }
    }
    if(axes.labels){
      mtext(colnames_dat[i],side=sidez
            #,adj=(G.lim+H[i,index]*c.lim)/(2*G.lim)
            ,adj = ADJS
      )
    }
  }
  SSE.linreg <- sapply(1:ncol(Y),function(x)sum((pred.data[,x]-Y[,x])^2))
  R2 <- 1-(SSE.linreg/apply(Y,2,function(x)var(x)*(length(x)-1)))
  return(R2)
}
#mds.linreg(Z=coords[,1:2], dat=tmp.dat)



ortho.procrustes <- function(Target, Testee){
  p <- ncol(Target)
  q <- ncol(Testee)
  if(!p==q){
    n.vec <- c("Target", "Testee")
    n <- nrow(Target)
    num.pad <- abs(p-q)
    pad.mat <- matrix(0,ncol=num.pad, nrow=n)
    index <- as.numeric(p>q)+1
    assign(n.vec[index],cbind(get(n.vec[index]),pad.mat))
    print(list(Target, Testee))
  }
  
  # Takes care of translation
  Y <- as.matrix(apply(Target,2,function(x)x-mean(x)))
  X <- as.matrix(apply(Testee,2,function(x)x-mean(x)))
  
  # Find the rho that takes care of Dilation
  svd.y.x <- svd(t(Y)%*%X)
  A <- svd.y.x$v%*%t(svd.y.x$u)
  rho <- sum(diag(X%*%A%*%t(Y)))/sum(diag(t(X)%*%X))
  
  # Find fitted Testee
  fit.X <- rho*X%*%A
  
  # Find the Procrustes Statistic
  proc.stat <- 1- sum(svd.y.x$d)^2/
    (sum(diag(t(X)%*%X))%*%sum(diag(t(Y)%*%Y)))
  
  return(list(fit.X, proc.stat) )
}

ortho.procrustes2 <- function(Target, Testee){
  p <- ncol(Target)
  q <- ncol(Testee)
  if(!p==q){
    n.vec <- c("Target", "Testee")
    n <- nrow(Target)
    num.pad <- abs(p-q)
    pad.mat <- matrix(0,ncol=num.pad, nrow=n)
    index <- as.numeric(p>q)+1
    assign(n.vec[index],cbind(get(n.vec[index]),pad.mat))
    print(list(Target, Testee))
  }
  
  # Takes care of translation
  Y <- as.matrix(apply(Target,2,function(x)x-mean(x)))
  X <- as.matrix(apply(Testee,2,function(x)x-mean(x)))
  
  # Find the rho that takes care of Dilation
  svd.y.x <- svd(t(Y)%*%X)
  #A <- svd.y.x$v%*%t(svd.y.x$u)
  
  #A <- svd.y.x$v%*%diag(svd.y.x$d^0.5)%*%t(svd.y.x$v)
  A <- svd.y.x$v%*%diag(svd.y.x$d)%*%t(svd.y.x$v)%*%solve(t(Y)%*%X)
  rho <- sum(diag(X%*%A%*%t(Y)))/sum(diag(t(X)%*%X))
  
  # Find fitted Testee
  fit.X <- rho*X%*%A
  
  # Find the Procrustes Statistic
  proc.stat <- 1- sum(svd.y.x$d)^2/
    (sum(diag(t(X)%*%X))%*%sum(diag(t(Y)%*%Y)))
  
  return(list(fit.X, proc.stat,A,rho) )
}


#Adapted from https://chitchatr.wordpress.com/2011/12/30/convex-hull-around-scatter-plot-in-r/

# INPUTS:
# xcoords: x-coordinates of point data
# ycoords: y-coordinates of point data
# lcolor: line color

# OUTPUTS:
# convex hull around data points in a particular color (specified by lcolor)

# FUNCTION:
Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
}  
# END OF FUNCTION


#############################################################################################
#####   Interpolation biplots
#############################################################################################

## Distance measures

sqrt.manhattan.dist <- function(x,y){as.matrix(proxy::dist(x,y,method = "Minkowski",p=1))^0.5}

eucl.dist <- function(x,y){as.matrix(proxy::dist(x,y,method = "Minkowski",p=2)^(2))^0.5}

pal.rf.dist <- function(x, model){(1-extract_proximity_full(model,x))^0.5}
pal.rf.dist.new <- function(x, y, model,obs_id = NULL){(1-extract_proximity_new_full(model,x,newdata = y,obs_id = NULL))^0.5}


add.nonlinear.axis.dist <- function(mds.obj, dist.fun="sqrt.manhattan.dist",orig.dat,scaled.dat=scale(orig.dat),which=c(1:ncol(orig.dat)),exp.factor=0, num.ticks=5){
  pretty.tick.vals <- list()
  tick.coords <- list()
  axes.coords <- list()
  m <- ncol(orig.dat)
  n <- nrow(orig.dat)
  for(i in 1:length(which)){
    j <- which[i]
    pretty.tick.vals[[i]] <- pretty(range(orig.dat[,j])+0.3*c(-1,1)*diff(range(orig.dat[,j])),n=num.ticks)
    tick.vals <- (pretty.tick.vals[[i]]-mean(orig.dat[,j]))/sd(orig.dat[,j])
    tick.mat <- matrix(0,ncol=m,nrow=length(tick.vals),byrow = TRUE)
    tick.mat[,j] <- tick.vals
    colnames(tick.mat) <- colnames(orig.dat)
    #tick.dists <- as.matrix(proxy::dist(tick.mat,scaled.dat,method = dist.fun,p=p))^0.5
    tick.dists <- eval(call(dist.fun,quote(tick.mat),quote(scaled.dat)))
    
    col.centers <- mds.obj$col.centers
    b.col.means <- matrix(col.centers,nrow = nrow(tick.dists),ncol=n,byrow = TRUE)
    tick.ddists <- -0.5*tick.dists^2
    b.row.means <- matrix(apply(tick.ddists,1,mean),nrow = nrow(tick.dists),ncol=n,byrow = FALSE)
    tick.ddists <- tick.ddists - b.col.means #- b.row.means - mds.obj$center
    tick.coords.tmp <- tick.ddists%*%mds.obj$eig.vectors[,1:2]%*%diag(mds.obj$eig.values^(-0.5))[1:2,1:2]
    tick.coords[[i]] <- tick.coords.tmp
  }
  for(i in 1:length(which)){
    j <- which[i]
    axes.range <- (range(pretty.tick.vals[[i]])-mean(orig.dat[,j]))/sd(orig.dat[,j])
    axes.vals <- seq(axes.range[1],axes.range[2], length.out=50)
    axes.mat <- matrix(0,ncol=m,nrow=length(axes.vals),byrow = TRUE)
    axes.mat[,j] <- axes.vals
    colnames(axes.mat) <- colnames(orig.dat)
    axes.dists <- eval(call(dist.fun,quote(axes.mat),quote(scaled.dat)))
    
    col.centers <- mds.obj$col.centers
    b.col.means <- matrix(col.centers,nrow = nrow(axes.dists),ncol=n,byrow = TRUE)
    axes.ddists <- -0.5*axes.dists^2
    b.row.means <- matrix(apply(axes.ddists,1,mean),nrow = nrow(axes.dists),ncol=n,byrow = FALSE)
    axes.ddists <- axes.ddists - b.col.means #- b.row.means - mds.obj$center
    axes.coords.tmp <- axes.ddists%*%mds.obj$eig.vectors[,1:2]%*%diag(mds.obj$eig.values^(-0.5))[1:2,1:2]
    axes.coords[[i]] <- axes.coords.tmp
  }
  if(is.null(mds.obj$axes.coords)){
    ind <- 1+length(mds.obj)
    mds.obj[[ind]] <- list(pretty.tick.vals=pretty.tick.vals,
                           tick.coords=tick.coords,
                           axes.coords = axes.coords)
    names(mds.obj)[ind] <- "axes.coords"
  } else {
    ind <- which(names(mds.obj)=="axes.coords")
    mds.obj[[ind]] <- list(pretty.tick.vals=pretty.tick.vals,
                           tick.coords=tick.coords,
                           axes.coords = axes.coords)
    names(mds.obj)[ind] <- "axes.coords"
  }
  return(mds.obj)
}


add.interpolated.points.dist <- function(mds.obj, dist.fun="sqrt.manhattan.dist",orig.dat,scaled.dat=scale(orig.dat),new.data){
  tick.coords <- list()
  orig.dat <- as.matrix(orig.dat)
  m <- ncol(orig.dat)
  n <- nrow(orig.dat)
  new.vals.mat <- sapply(1:m,function(x)(new.data[,x]-mean(orig.dat[,x]))/sd(orig.dat[,x]))
  colnames(new.vals.mat) <- colnames(orig.dat)
  new.dists <- eval(call(dist.fun,quote(new.vals.mat),quote(scaled.dat)))
  #    print(new.vals.mat)
  col.centers <- mds.obj$col.centers
  b.col.means <- matrix(col.centers,nrow = nrow(new.dists),ncol=n,byrow = TRUE)
  if(is.null(mds.obj$ACP)){
    new.ddists <- -0.5*new.dists^2
  } else{
    new.ddists <- -0.5*(new.dists^2+mds.obj$ACP)
  }
  b.row.means <- matrix(apply(new.ddists,1,mean),nrow = nrow(new.dists),ncol=n,byrow = FALSE)
  new.ddists <- new.ddists - b.col.means# - b.row.means - mds.obj$center
  new.coords.tmp <- new.ddists%*%mds.obj$eig.vectors[,1:2]%*%diag(mds.obj$eig.values^(-0.5))[1:2,1:2]
  new.coords <- new.coords.tmp
  if(is.null(mds.obj$new.coords)){
    ind <- 1+length(mds.obj)
    mds.obj[[ind]] <- list(new.data=new.data,
                           new.coords=new.coords)
    names(mds.obj)[ind] <- "new.coords"
  } else {
    ind <- which(names(mds.obj)=="new.coords")
    mds.obj[[ind]] <- list(new.data=new.data,
                           new.coords=new.coords)
    names(mds.obj)[ind] <- "new.coords"
  }
  return(mds.obj)
}


# For transformation on variables prior to scaling eg. log transform

add.interpolated.transformed.points.dist <- function(mds.obj,
                                                     dist.fun="eucl.dist",
                                                     orig.dat,trans.fun="log",
                                                     trans.dat=eval(call(trans.fun,orig.dat)),
                                                     scaled.dat = scale(trans.dat),
                                                     new.data){
  tick.coords <- list()
  orig.dat <- as.matrix(orig.dat)
  m <- ncol(orig.dat)
  n <- nrow(orig.dat)
  new.vals.mat <- (sapply(1:m,function(x){
    (eval(call(trans.fun,new.data[,x,drop=FALSE]))-mean(trans.dat[,x]))/sd(trans.dat[,x])
  }))
  colnames(new.vals.mat) <- colnames(orig.dat)
  new.dists <- eval(call(dist.fun,quote(new.vals.mat),quote(scaled.dat)))
  #    print(new.vals.mat)
  col.centers <- mds.obj$col.centers
  b.col.means <- matrix(col.centers,nrow = nrow(new.dists),ncol=n,byrow = TRUE)
  
  new.ddists <- -0.5*new.dists^2
  
  b.row.means <- matrix(apply(new.ddists,1,mean),nrow = nrow(new.dists),ncol=n,byrow = FALSE)
  new.ddists <- new.ddists - b.col.means# - b.row.means - mds.obj$center
  new.coords.tmp <- new.ddists%*%mds.obj$eig.vectors[,1:2]%*%diag(mds.obj$eig.values^(-0.5))[1:2,1:2]
  new.coords <- new.coords.tmp
  if(is.null(mds.obj$new.coords)){
    ind <- 1+length(mds.obj)
    mds.obj[[ind]] <- list(new.data=new.data,
                           new.coords=new.coords)
    names(mds.obj)[ind] <- "new.coords"
  } else {
    ind <- which(names(mds.obj)=="new.coords")
    mds.obj[[ind]] <- list(new.data=new.data,
                           new.coords=new.coords)
    names(mds.obj)[ind] <- "new.coords"
  }
  return(mds.obj)
}

add.transformed.axis <- function(mds.obj,orig.dat,trans.fun="log",trans.dat=eval(call(trans.fun,orig.dat)), scaled.dat = scale(trans.dat),which=c(1:ncol(orig.dat)),exp.factor=0.3){
  pretty.tick.vals <- list()
  tick.coords <- list()
  m <- ncol(orig.dat)
  n <- nrow(orig.dat)
  for(i in 1:length(which)){
    j <- which[i]
    pretty.tick.vals[[i]] <- pretty(range(orig.dat[,j])+exp.factor*c(-1,1)*diff(range(orig.dat[,j])))
    tick.vals <- (eval(call(trans.fun,pretty.tick.vals[[i]]))-mean(trans.dat[,j]))/sd(trans.dat[,j])
    tick.mat <- matrix(0,ncol=m,nrow=length(tick.vals))
    tick.mat[,j] <- tick.vals
    tick.dists <- as.matrix(proxy::dist(tick.mat,scaled.dat))
    col.centers <- mds.obj$col.centers
    b.col.means <- matrix(col.centers,nrow = nrow(tick.dists),ncol=n,byrow = TRUE)
    tick.ddists <- -0.5*tick.dists^2
    b.row.means <- matrix(colMeans(t(tick.ddists)),nrow = nrow(tick.dists),ncol=n,byrow = FALSE)
    tick.ddists <- tick.ddists  - b.col.means #- b.row.means + mds.obj$center
    tick.coords.tmp <- tick.ddists%*%mds.obj$coords[,1:2]%*%diag(mds.obj$eig.values^(-1))[1:2,1:2]
    tick.coords[[i]] <- tick.coords.tmp
  }
  if(is.null(mds.obj$axes.coords)){
    ind <- 1+length(mds.obj)
    mds.obj[[ind]] <- list(pretty.tick.vals=pretty.tick.vals,
                           tick.coords=tick.coords)
    names(mds.obj)[ind] <- "axes.coords"
  } else {
    ind <- which(names(mds.obj)=="axes.coords")
    mds.obj[[ind]] <- list(pretty.tick.vals=pretty.tick.vals,
                           tick.coords=tick.coords)
    names(mds.obj)[ind] <- "axes.coords"
  }
  return(mds.obj)
}


add.nonlinear.axis <- function(mds.obj,model=penguin.mod,orig.dat,which=c(1:ncol(orig.dat)),exp.factor=0.3, num.ticks=5, obs.ind = 11,num.axes.ticks = 50){
  pretty.tick.vals <- list()
  tick.coords <- list()
  axes.coords <- list()
  m <- ncol(orig.dat)
  n <- nrow(orig.dat)
  #  obs.ind <- 185
  mean.vec <- orig.dat[obs.ind,]
  for(i in 1:length(which)){
    j <- which[i]
    pretty.tick.vals[[i]] <- pretty(range(orig.dat[,j])+exp.factor*c(-1,1)*diff(range(orig.dat[,j])),n=num.ticks)
    tick.vals <- pretty.tick.vals[[i]]
    tick.mat <- matrix(mean.vec,ncol=m,nrow=length(tick.vals),byrow = TRUE)#matrix(0,ncol=m,nrow=length(tick.vals),byrow = TRUE)#matrix(mean.vec,ncol=m,nrow=length(tick.vals),byrow = TRUE)
    tick.mat[,j] <- tick.vals
    colnames(tick.mat) <- colnames(orig.dat)
    tmp.dissim.array <- pal.rf.dist.new(x = orig.dat,y=tick.mat,model=model)
    #print(dim(tmp.dissim.array))
    tick.dists <- tmp.dissim.array#apply(tmp.dissim.array,1,function(x)apply(x,1,sum))
    col.centers <- mds.obj$col.centers
    b.col.means <- matrix(col.centers,nrow = nrow(tick.dists),ncol=n,byrow = TRUE)
    tick.ddists <- -0.5*tick.dists^2
    b.row.means <- matrix(apply(tick.ddists,1,mean),nrow = nrow(tick.dists),ncol=n,byrow = FALSE)
    tick.ddists <- tick.ddists - b.col.means - b.row.means + mds.obj$center
    tick.coords.tmp <- tick.ddists%*%mds.obj$eig.vectors[,1:2]%*%diag(mds.obj$eig.values^(-0.5))[1:2,1:2]
    tick.coords[[i]] <- tick.coords.tmp
  }
  for(i in 1:length(which)){
    j <- which[i]
    axes.range <- range(pretty.tick.vals[[i]])
    axes.vals <- seq(axes.range[1],axes.range[2], length.out=num.axes.ticks)
    axes.mat <- matrix(mean.vec,ncol=m,nrow=length(axes.vals),byrow = TRUE)
    axes.mat[,j] <- axes.vals
    colnames(axes.mat) <- colnames(orig.dat)
    axes.dists <- pal.rf.dist.new(x=orig.dat,y=axes.mat,model=model)
    
    col.centers <- mds.obj$col.centers
    b.col.means <- matrix(col.centers,nrow = nrow(axes.dists),ncol=n,byrow = TRUE)
    axes.ddists <- -0.5*axes.dists^2
    b.row.means <- matrix(apply(axes.ddists,1,mean),nrow = nrow(axes.dists),ncol=n,byrow = FALSE)
    axes.ddists <- axes.ddists - b.col.means #- b.row.means - mds.obj$center
    axes.coords.tmp <- axes.ddists%*%mds.obj$eig.vectors[,1:2]%*%diag(mds.obj$eig.values^(-0.5))[1:2,1:2]
    axes.coords[[i]] <- axes.coords.tmp
  }
  if(is.null(mds.obj$axes.coords)){
    ind <- 1+length(mds.obj)
    mds.obj[[ind]] <- list(pretty.tick.vals=pretty.tick.vals,
                           tick.coords=tick.coords,
                           axes.coords = axes.coords)
    names(mds.obj)[ind] <- "axes.coords"
  } else {
    ind <- which(names(mds.obj)=="axes.coords")
    mds.obj[[ind]] <- list(pretty.tick.vals=pretty.tick.vals,
                           tick.coords=tick.coords,
                           axes.coords = axes.coords)
    names(mds.obj)[ind] <- "axes.coords"
  }
  # if(is.null(mds.obj$axes.coords)){
  #   ind <- 1+length(mds.obj)
  #   mds.obj[[ind]] <- list(pretty.tick.vals=pretty.tick.vals,
  #                          tick.coords=tick.coords)
  #   names(mds.obj)[ind] <- "axes.coords"
  # } else {
  #   ind <- which(names(mds.obj)=="axes.coords")
  #   mds.obj[[ind]] <- list(pretty.tick.vals=pretty.tick.vals,
  #                          tick.coords=tick.coords)
  #   names(mds.obj)[ind] <- "axes.coords"
  # }
  return(mds.obj)
}

# Isomap interp biplot

eig.isomap.fun.interp <- function(dist.mat,eps=Opt.eps.fun(dist.mat,p=2)$opt.eps,
                                  shortest=TRUE,
                                  Trace=FALSE,
                                  points.to.plot=1:nrow(dist.mat)){
  if(shortest){
    s.paths <- vegan::stepacross(dist.mat,
                                 toolong=eps, 
                                 trace = Trace)
    path.meth <- "shortest"
  } else {
    s.paths <- vegan::stepacross(dist.mat,
                                 path='extended',
                                 toolong=eps, 
                                 trace = Trace)
    path.meth <- "extended"
  }
  dissim.mat <- as.matrix(s.paths)[points.to.plot,points.to.plot]
  
  A.mat <- -0.5*dissim.mat^2
  col.centers <- colMeans(A.mat)
  row.centers <- colMeans(t(A.mat))
  center <-  mean(A.mat)
  centred_mat <- apply(apply(A.mat,1,function(x)x-col.centers),
                       1,function(x)x-row.centers)+center
  
  
  eig.dat <- eigen(centred_mat)
  # The coordinates for the factual observations
  coords <- eig.dat$vectors[,1:2]%*%diag(eig.dat$values[1:2]^0.5)#%*%(rot.mat)
  
  mds.obj=list(dissim.mat=dissim.mat, col.centers=col.centers, center=center, coords=coords, eig.values=eig.dat$values, n=nrow(coords),eig.vectors=eig.dat$vectors,path.meth=path.meth, eps=eps )
  
  return(mds.obj)
}

add.isomap.axis <- function(mds.obj,model=penguin.mod,orig.dat,which=c(1:ncol(orig.dat)),exp.factor=0.3, num.ticks=5, obs.ind = 11, num.axes.points = 50){
  pretty.tick.vals <- list()
  tick.coords <- list()
  axes.coords <- list()
  m <- ncol(orig.dat)
  n <- nrow(orig.dat)
  eps <- mds.obj$eps
  #  obs.ind <- 185
  mean.vec <- orig.dat[obs.ind,]
  isomap.dissim <- mds.obj$dissim.mat
  for(i in 1:length(which)){
    j <- which[i]
    pretty.tick.vals[[i]] <- pretty(range(orig.dat[,j])+exp.factor*c(-1,1)*diff(range(orig.dat[,j])),n=num.ticks)
    tick.vals <- pretty.tick.vals[[i]]
    tick.mat <- matrix(mean.vec,ncol=m,nrow=length(tick.vals),byrow = TRUE)#matrix(0,ncol=m,nrow=length(tick.vals),byrow = TRUE)#matrix(mean.vec,ncol=m,nrow=length(tick.vals),byrow = TRUE)
    tick.mat[,j] <- tick.vals
    colnames(tick.mat) <- colnames(orig.dat)
    tmp.dissim.array <- pal.rf.dist.new(x = orig.dat,y=tick.mat,model=model)
    tmp.dissim.array <- base::ifelse(tmp.dissim.array>=eps,NA,tmp.dissim.array)
    tick.dists <- tmp.dissim.array
    
    for(w in 1:nrow(tmp.dissim.array)){
      tmp.appended.array <- cbind(rbind(isomap.dissim,tmp.dissim.array[w,]),c(tmp.dissim.array[w,],0))
      tmp.isomap.dissim <- as.matrix(vegan::stepacross(tmp.appended.array,path = mds.obj$path.meth,toolong = max(tmp.appended.array)+2,trace = FALSE))
      tick.dists[w,] <-  tmp.isomap.dissim[n+1,1:n]
    }
    #print(dim(tmp.dissim.array))
    #tick.dists <- tmp.dissim.array
    
    col.centers <- mds.obj$col.centers
    b.col.means <- matrix(col.centers,nrow = nrow(tick.dists),ncol=n,byrow = TRUE)
    tick.ddists <- -0.5*tick.dists^2
    b.row.means <- matrix(apply(tick.ddists,1,mean),nrow = nrow(tick.dists),ncol=n,byrow = FALSE)
    tick.ddists <- tick.ddists - b.col.means - b.row.means + mds.obj$center
    tick.coords.tmp <- tick.ddists%*%mds.obj$eig.vectors[,1:2]%*%diag(mds.obj$eig.values^(-0.5))[1:2,1:2]
    tick.coords[[i]] <- tick.coords.tmp
  }
  for(i in 1:length(which)){
    j <- which[i]
    axes.range <- range(pretty.tick.vals[[i]])
    axes.vals <- seq(axes.range[1],axes.range[2], length.out=num.axes.points)
    axes.mat <- matrix(mean.vec,ncol=m,nrow=length(axes.vals),byrow = TRUE)
    axes.mat[,j] <- axes.vals
    colnames(axes.mat) <- colnames(orig.dat)
    axes.dists <- pal.rf.dist.new(x=orig.dat,y=axes.mat,model=model)
    tmp.dissim.array <- ifelse(axes.dists>=eps,NA,axes.dists)
    axes.dists <- tmp.dissim.array
    
    for(w in 1:nrow(tmp.dissim.array)){
      tmp.appended.array <- cbind(rbind(isomap.dissim,tmp.dissim.array[w,]),c(tmp.dissim.array[w,],0))
      tmp.isomap.dissim <- as.matrix(vegan::stepacross(tmp.appended.array,path = mds.obj$path.meth,toolong = max(tmp.appended.array)+2,trace = FALSE))
      axes.dists[w,] <-  tmp.isomap.dissim[n+1,1:n]
    }
    
    col.centers <- mds.obj$col.centers
    b.col.means <- matrix(col.centers,nrow = nrow(axes.dists),ncol=n,byrow = TRUE)
    axes.ddists <- -0.5*axes.dists^2
    b.row.means <- matrix(apply(axes.ddists,1,mean),nrow = nrow(axes.dists),ncol=n,byrow = FALSE)
    axes.ddists <- axes.ddists - b.col.means #- b.row.means - mds.obj$center
    axes.coords.tmp <- axes.ddists%*%mds.obj$eig.vectors[,1:2]%*%diag(mds.obj$eig.values^(-0.5))[1:2,1:2]
    axes.coords[[i]] <- axes.coords.tmp
  }
  if(is.null(mds.obj$axes.coords)){
    ind <- 1+length(mds.obj)
    mds.obj[[ind]] <- list(pretty.tick.vals=pretty.tick.vals,
                           tick.coords=tick.coords,
                           axes.coords = axes.coords)
    names(mds.obj)[ind] <- "axes.coords"
  } else {
    ind <- which(names(mds.obj)=="axes.coords")
    mds.obj[[ind]] <- list(pretty.tick.vals=pretty.tick.vals,
                           tick.coords=tick.coords,
                           axes.coords = axes.coords)
    names(mds.obj)[ind] <- "axes.coords"
  }
  return(mds.obj)
}



## Spline based axes functions 
# The original code is attributable to Prof. Niel le Roux and appeared in the 
#   code for the unpublished Understanding Biplots package. An implementation 
#   of spline biplots have since become available in the R package biplotEZ 
#   available on CRAN.

biplot.spline.axis <- function(j, X, Ytilde, means, sd, n.int, spline.control, dmeth=0, ... )
{
  n<-nrow(X)
  p<-ncol(X)
  if (n > 103)
  {  my.sample <- sample (1:n, size=103, replace=F)
  X <- X[my.sample,]
  Ytilde <- Ytilde[my.sample,]
  n <- nrow(X)
  }
  
  tau <- spline.control$tau
  nmu <- spline.control$nmu
  u <- spline.control$u
  v <- spline.control$v
  lambda <- spline.control$lambda
  smallsigma <- spline.control$smallsigma
  bigsigma <- spline.control$bigsigma
  gamma <- spline.control$gamma
  bigsigmaactivate <- spline.control$bigsigmaactivate
  eps <- spline.control$eps
  tiny <- spline.control$tiny
  itmax <- spline.control$itmax
  ftol <- spline.control$ftol
  
  cat ("Calculation spline axis for variable", j, "\n")
  require(splines)
  if(dmeth==1)stop("dmeth should be equal to zero or integer greater than 1 \n")  
  Y<-scale(Ytilde,center=means,scale=sd)
  
  ytilde<-Ytilde[,j]
  mutilde<-seq(from=min(ytilde),to=max(ytilde),length.out=nmu)
  y<-Y[,j]
  rangey<-max(y)-min(y)
  mu<-seq(from=min(y)-.3*rangey,to=max(y)+.3*rangey,length.out=nmu)
  markers <- (pretty(ytilde)-means[j])/sd[j]
  mu <- sort(c(mu,markers))
  mu <- unique(mu)
  nmu <- length(mu)
  
  if (v>0)
  {
    knots<-seq.int(from=0,to=1,length.out=v+2)[-c(1,v+2)]
    knots<-quantile(y,knots)
    M<-bs(mu,knots=knots,degree=u,intercept=FALSE)
  }
  else M<-bs(mu,df=u+v,degree=u,intercept=FALSE)
  M<-scale(M,scale=FALSE,center=M[which.min(abs(mu)),]) # To ensure that the spline passes through the origin at the calibration which represents the mean of the variable
  Breg<-t(solve(t(X)%*%X)%*%t(X)%*%y)
  Zreg<-mu%*%Breg/sum(Breg^2)
  Bvec<-as.vector(solve(t(M)%*%M)%*%t(M)%*%Zreg)  # Closest to regression biplot
  
  const1<-sum(y^2)
  const2<-sum(X^2)/(n*p)
  TotalNumberOfLossFunctionCalls<-0
  
  optimtouse<-function(Bvec)
  {
    timetemp<-proc.time()[3]
    LOSS<-1.0
    LOSS1<-1.0
    Ind<-rep(1,n)
    pred<-rep(0,nmu)
    deltmp<-0
    tau<-tau
    #.5 # the choice of tau seems to affect perfomance quite substantially.
    # tau is used to specify the points on the inital simplex.
    Ay<-rep(0,(u+v)*p+1)
    TEMPVK<-rep(0,(u+v)*p)
    iter1<-0
    iter<-0
    ERRO <-0
    
    # Prepare for Fortran subroutine
    storage.mode(X)<- "double"
    storage.mode(Ind)<- "integer"
    storage.mode(mu)<- "double"
    storage.mode(pred)<- "double"
    storage.mode(y)<- "double"
    storage.mode(M)<- "double"
    storage.mode(Bvec)<- "double"
    storage.mode(Ay)<- "double"
    storage.mode(TEMPVK)<- "double"
    
    returned_data<-.Fortran('L',LOSS=as.double(LOSS),X=X,n=as.integer(n),p=as.integer(p),nmu=as.integer(nmu),Ind=Ind,
                            mu=mu,pred=pred,lambda=as.double(lambda),y=y,const1=as.double(const1),const2=as.double(const2),u=as.integer(u),
                            v=as.integer(v),M=M,Bvec=Bvec,tau=as.double(tau),Ay=Ay,TEMPVEK=TEMPVK,iter=as.integer(iter),
                            ftol=as.double(ftol),LOSS1=as.double(LOSS1),iter1=as.integer(iter1),fout = as.integer(ERRO),
                            const3=as.double(tiny), itmax=as.integer(itmax))
    if(returned_data$fout > 0)
    {
      cat("Fout is: ", returned_data$fout, "\n")
      warning("Increase itmax for Fortran \n")
    }
    
    B<-matrix(returned_data$Bvec,ncol=p)
    Z<-M%*%B 
    
    aa<-list(BestValue=returned_data$LOSS,BestSolution=returned_data$Bvec,ConvergenceCode=returned_data$fout, iter1=returned_data$iter1,
             iter=returned_data$iter,TimeTaken=proc.time()[3]-timetemp)
    aa
  }
  EuclidDist2 <- function (X, Y) 
  {
    n <- nrow(X)
    m <- nrow(Y)
    bx <- rowSums(X^2)
    by <- rowSums(Y^2)
    outer(bx, by, FUN = "+") - 2 * X %*% t(Y)
  }
  
  ### Variable initialisation
  outBestValues<-rep(NA,gamma+1)
  outBestSolutions<-matrix(nrow=2*(u+v),ncol=gamma+1)
  outTimeTaken<-rep(NA,gamma+1) # Is made one element longer at each iteration.
  BestSolutionsFrequency<-rep(NA,gamma+1)
  BestSolutionsIndices<-rep(NA,gamma+1) # Is made one element longer at each iteration.
  SquaredDistancesBetweenBestSolutions<-matrix(nrow=gamma+1,ncol=gamma+1)
  
  ### Initial coefficients closest to regression biplot
  temp<-optimtouse(Bvec)
  outBestValues[1]<-temp$BestValue
  outBestSolutions[,1]<-temp$BestSolution
  outTimeTaken[1]<-temp$TimeTaken
  BestSolutionsFrequency[1]<-1
  BestSolutionsIndices[1]<-1
  DistinctSolutions<-1
  PreviousBestSolution<-NA
  nSameSolutionConsecutively<-0
  BigSigmaActivations<-NULL
  
  test.iter <- temp$iter
  test.iter1 <- temp$iter1
  
  ### Last best coefficients perturbed
  for (gammacounter in 2:(gamma+1))
  {
    if (nSameSolutionConsecutively>=bigsigmaactivate)
    {
      temp<-optimtouse(outBestSolutions[,which.min(outBestValues)]+rnorm((u+v)*2,mean=0,sd=bigsigma))
      BigSigmaActivations<-c(BigSigmaActivations,gammacounter)
    }
    else temp<-optimtouse(outBestSolutions[,which.min(outBestValues)]+rnorm((u+v)*2,mean=0,sd=smallsigma))
    outTimeTaken[gammacounter]<-temp$TimeTaken
    tempSquaredDistances<-EuclidDist2(matrix(temp$BestSolution,nrow=1),t(outBestSolutions[,1:DistinctSolutions]))
    if (any(tempSquaredDistances<eps))
    {
      BestSolutionsIndices[gammacounter]<-tempAA<-which.min(tempSquaredDistances)
      BestSolutionsFrequency[tempAA]<-BestSolutionsFrequency[tempAA]+1
      if (!is.na(PreviousBestSolution) && tempAA==PreviousBestSolution) nSameSolutionConsecutively<-nSameSolutionConsecutively+1
      else
      {
        PreviousBestSolution<-tempAA
        nSameSolutionConsecutively<-0
      }
    }
    else
    {
      DistinctSolutions<-DistinctSolutions+1
      outBestValues[DistinctSolutions]<-temp$BestValue
      outBestSolutions[,DistinctSolutions]<-temp$BestSolution
      BestSolutionsFrequency[DistinctSolutions]<-1
      BestSolutionsIndices[gammacounter]<-DistinctSolutions
      SquaredDistancesBetweenBestSolutions[1:(DistinctSolutions-1),DistinctSolutions]<-tempSquaredDistances
      nSameSolutionConsecutively<-0
    }
  }
  axis.points <- cbind(M%*%matrix(outBestSolutions[,which.min(outBestValues)],ncol=2), mu, 0)
  
  for (i in 1:nrow(axis.points)) if (any(zapsmall(axis.points[i,3] - markers) == 0)) axis.points[i, 4] <- 1
  axis.points[,3] <- axis.points[,3]*sd[j] + means[j]
  TTT[[length(TTT)+1]] <<- axis.points
  axis.points
}

