# Isomap Functions


# Find minimum eps for connected graph
## Takes: Dissimilarity matrix and a step size for the search for the minimum 
##        epsilon value
## Outputs: Minimum Epsilon and Maximum Epsilon as numeric scalars and numeric vector of epsilons between 
##          the eps.min and eps.max
min.eps.fun <- function(dist.mat,step.size=(max(dist.mat)-min(dist.mat))/200){
  library('igraph')
  eps.vec <- c()
  groups.vec <- c()
  num.groups = 0
  eps <- max.eps <-  max(dist.mat)
  i <- 0
  g.list <- list()
  while(num.groups <2){
    if(i!=0){
      g.list[[length(g.list)+1]] <- g
    }
    i <- i+1
    eps.vec[i] <- eps
    star.mat <- ifelse( dist.mat < eps, dist.mat, 0)
    g <- graph.adjacency(star.mat, weighted = TRUE, mode="undirected")
    num.groups <- components(g)$no
    groups.vec[i] <- num.groups
    eps <- eps - step.size
  }
  #plot(eps.vec, groups.vec)
  min.eps <- min(eps.vec[groups.vec==1])
  output <- list(min.eps=min.eps, max.eps=max.eps ,eps.vec=eps.vec[groups.vec==1])
  return(output)
}

# Normalized Stress Function
Stress.fun <- function(Dis1,Dis2, Normalize=TRUE){
  H <- diag(nrow(Dis1))-1/nrow(Dis1)  # Centering matrix
  if(Normalize){
    A <- ((-H%*%(Dis1^2)%*%H/2) - (-H%*%(Dis2^2)%*%H/2))^2/sum((H%*%(Dis1^2)%*%H/2)^2)
  } else {
    A <- ((-H%*%(Dis1^2)%*%H/2) - (-H%*%(Dis2^2)%*%H/2))^2
  }  
  return(sqrt(sum(A[upper.tri(A)])))
}


# Find the epsilon value that corresponds to the lowest Normalized Stress
## Takes: Dissimilarity matrix, a vector of allowable epsilon values and
##        the number of dimensions of the solution.
## Outputs: The optimal epsilon value, vector of epsilons evaluated,
##          vector of loss values, 
##          vector of 1-cor(geodesic distances,solution distances)
Opt.eps.fun <- function(dist.mat,
                        eps.vec=min.eps.fun(dist.mat=dist.mat)$eps.vec,
                        p=2,
                        shortest=TRUE,
                        Trace=FALSE,
                        Normalize=TRUE){
  library('vegan')
  loss.vec <- c() # Store the loss values
  dist.correl.vec <- c() # Store the correlations between the distances
  for(i in 1:length(eps.vec)){
    if(shortest){
      s.paths <- stepacross(dist.mat,
                            toolong=eps.vec[i], 
                            trace = Trace)
    } else {
      s.paths <- stepacross(dist.mat,
                            path='extended',
                            toolong=eps.vec[i], 
                            trace = Trace)
    }

    
    dissim.mat <- as.matrix(s.paths)
    
    A.mat <- -0.5*dissim.mat^2
    col.centers <- colMeans(A.mat)
    row.centers <- colMeans(t(A.mat))
    center <-  mean(A.mat)
    centred_mat <- apply(apply(A.mat,1,function(x)x-col.centers),
                         1,function(x)x-row.centers)+center
    
    
    eig.dat <- eigen(centred_mat)
    tmp.coords <- eig.dat$vectors[,1:p,drop=FALSE]%*%diag(eig.dat$values[1:p]^(0.5),nrow=p)
    tmp.dist <- as.matrix(dist(tmp.coords))
    E.val <- Stress.fun(dissim.mat,tmp.dist, Normalize=Normalize)
    loss.vec[i] <- E.val
    dist.correl <- cor(dissim.mat[upper.tri(dissim.mat)],
                       tmp.dist[upper.tri(tmp.dist)])
    dist.correl.vec[i] <- dist.correl
  }
#  opt.eps <- eps.vec[which.min(loss.vec)]
#  print(eps.vec)
  if(length(eps.vec)>3){
    window.size <- ceiling(length(loss.vec)*0.01)
    loss.vec.expanded <- c(rep(r[1]+0.01,window.size),loss.vec,rep(r[length(loss.vec)]+0.01,window.size))
    local.min.vec <- sapply((1+window.size):(length(loss.vec)),function(x)all(loss.vec.expanded[x]<c(loss.vec.expanded[c((x-1):(x-window.size),(x+1):(x+window.size))])))
    opt.eps <- eps.vec[which(local.min.vec)][which.max(dist.correl.vec[which(local.min.vec)])]
  } else { opt.eps <- max(eps.vec)}
  
  output <- list(opt.eps=opt.eps,
                 eps.vec=eps.vec,
                 loss.vec=loss.vec, 
                 dist.correl.vec=dist.correl.vec)
}


# Find the Eigenvalues and Vectors for the Isomap algorithm at a given 
#   epsilon value
## Takes: Dissimilarity matrix, epsilon value
## Outputs: Eigenvectors and Eigenvalues
eig.isomap.fun <- function(dist.mat,eps=Opt.eps.fun(dist.mat,p=2)$opt.eps,
                           shortest=TRUE,
                           Trace=FALSE,
                           points.to.plot=1:nrow(dist.mat)){
  if(shortest){
    s.paths <- stepacross(dist.mat,
                          toolong=eps, 
                          trace = Trace)
  } else {
    s.paths <- stepacross(dist.mat,
                          path='extended',
                          toolong=eps, 
                          trace = Trace)
  }
  dissim.mat <- as.matrix(s.paths)[points.to.plot,points.to.plot]

  A.mat <- -0.5*dissim.mat^2
  col.centers <- colMeans(A.mat)
  row.centers <- colMeans(t(A.mat))
  center <-  mean(A.mat)
  centred_mat <- apply(apply(A.mat,1,function(x)x-col.centers),
                       1,function(x)x-row.centers)+center
  
  
  eig.dat <- eigen(centred_mat)
#  return(list(eig.obj=eig.dat,eps=eps))
  return(list(eig.obj=eig.dat,eps=eps,dissim.mat=dissim.mat))
}

eigs.to.coords.fun <- function(eig.obj,p=2){
  coords <- eig.obj$vectors[,1:p,drop=FALSE]%*%diag(eig.obj$values[1:p],nrow=p)^0.5
  output <- list(coords)
  return(output)
}







