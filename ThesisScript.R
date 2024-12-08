# Extracting the full proximities between every pair of observations
extract_proximity_full = function(fit, olddata) {
  pred = predict(fit, olddata, type = "terminalNodes")$predictions
  prox = matrix(NA, nrow(pred), nrow(pred))
  ntree = ncol(pred)
  n = nrow(prox)
  
  if (is.null(fit$inbag.counts)) {
    stop("call ranger with keep.inbag = TRUE")
  }
  
  # Get inbag counts
  inbag = simplify2array(fit$inbag.counts)
  
  for (i in 1:n) {
    for (j in 1:n) {
      prox[i, j] = sum(pred[i, ] == pred[j, ]) / ncol(pred)
    }
  }
  
  prox
}

# Extracting the GAP proximities between every pair of observations
extract_proximity_gap_old = function(fit, olddata) {
  pred = predict(fit, olddata, type = "terminalNodes")$predictions
  prox = matrix(NA, nrow(pred), nrow(pred))
  ntree = ncol(pred)
  n = nrow(prox)
  
  if (is.null(fit$inbag.counts)) {
    stop("call ranger with keep.inbag = TRUE")
  }
  
  # Get inbag counts
  inbag = simplify2array(fit$inbag.counts)
  
  for (i in 1:n) {
    for (j in 1:n) {
      # Assymetric so the proximities must be averaged
      tree_idx1 = inbag[i, ] == 0 & inbag[j, ] > 0
      prox1 <- sum(inbag[j,tree_idx1][pred[i, tree_idx1] == pred[j, tree_idx1]]) / n#sum(tree_idx1)
      tree_idx2 = inbag[i, ] == 0 & inbag[j, ] > 0
      prox2 <- sum(inbag[i,tree_idx2][pred[i, tree_idx2] == pred[j, tree_idx2]]) / n#sum(tree_idx2)
      prox[i, j] <-  (prox1 + prox2)/2
    }
  }
  
  prox
}


# Do not use # Must take the full dataset, cannot yet do any interpolation
extract_proximity_gap = function(fit, olddata) {
  pred = predict(fit, olddata, type = "terminalNodes")$predictions
  #prox = matrix(NA, nrow(pred), nrow(pred))
  prox = matrix(0, nrow(pred), nrow(pred))
  ntree = ncol(pred)
  n = nrow(prox)
  
  if (is.null(fit$inbag.counts)) {
    stop("call ranger with keep.inbag = TRUE")
  }
  
  # Get inbag counts
  inbag = simplify2array(fit$inbag.counts)
  
  for (i in 1:n) {
    for (j in 1:i) {
      # Assymetric so the proximities must be averaged
      tree_idx1 = inbag[i, ] == 0 #& inbag[j, ] > 0 # Trees for which i is oob and j ib
      share.trees <- pred[i,] == pred[j, ]# Trees where i and j share a terminal node      
      size.terminal.node1 <- sapply(1:ntree,function(x)sum(c(pred[i,x]==pred[,x])*c(inbag[,x]))) # The number of training observations that share a terminal node with i
      prox1 <- sum(inbag[j,tree_idx1&share.trees]/ size.terminal.node1[tree_idx1&share.trees]) #sum(tree_idx1)
      tree_idx2 = inbag[j, ] == 0 #& inbag[i, ] > 0
      size.terminal.node2 <- sapply(1:ntree,function(x)sum(c(pred[j,x]==pred[,x])*c(inbag[,x]))) # The number of training observations that share a terminal node with j
      prox2 <- sum(inbag[i,tree_idx2&share.trees]/ size.terminal.node2[tree_idx2&share.trees]) #sum(tree_idx2)
      prox[i, j] <-  (prox1 + prox2)/2
    }
  }
  
  prox <- prox + t(prox)
  prox
}


# Must take the full dataset, cannot yet do any interpolation. Self-similarities are zero
extract_gap.prox <- function(fit, olddata){
  prox <- matrix(0, nrow=nrow(olddata), ncol=nrow(olddata))
  # Get Terminal Node Predictions
  pred = predict(fit, olddata, type = "terminalNodes")$predictions
  # Get inbag counts
  inbag = simplify2array(fit$inbag.counts)
#  Si <- 
#  Si.size <- 
#  cj(t) <- 
#  share.tree <- 
#  in.bag <- 
#  Ji <- share.tree*in.bag
#  node.size <-
  for (i in 1:n) {
    for (j in 1:n) {    
      oob.trees <- inbag[i, ] == 0 # Trees for which i is oob
      num.oob.trees <- sum(oob.trees)
      node.size <- sapply(which(oob.trees),function(x)sum(c(pred[i,x]==pred[,x])*c(inbag[,x])))
      cj <- inbag[j,oob.trees]
      share.pred <- pred[i,oob.trees]==pred[j,oob.trees]
      prox_ij <- sum(cj*share.pred/node.size)/num.oob.trees
      prox[i,j] <- prox_ij
    }  
  }
  prox <- (prox + t(prox))/2
}
#aa <- extract_gap.prox(fit=mod,olddata = tmp.dat)
#aa2 <- extract_proximity_gap(fit=mod,olddata = tmp.dat)
# Extracting the OOB proximities between every pair of observations
extract_proximity_oob = function(fit, olddata) {
  pred = predict(fit, olddata, type = "terminalNodes")$predictions
  prox = matrix(NA, nrow(pred), nrow(pred))
  ntree = ncol(pred)
  n = nrow(prox)
  
  if (is.null(fit$inbag.counts)) {
    stop("call ranger with keep.inbag = TRUE")
  }
  
  # Get inbag counts
  inbag = simplify2array(fit$inbag.counts)
  
  for (i in 1:n) {
    for (j in 1:n) {
      # Use only trees where both obs are OOB
      tree_idx = inbag[i, ] == 0 & inbag[j, ] == 0
      prox[i, j] = sum(pred[i, tree_idx] == pred[j, tree_idx]) / sum(tree_idx)
    }
  }
  
  prox
}


# Extracting the proximities between a psuedo point and every existing observation
extract_proximity_new = function(fit, olddata, newdata, obs_id=NULL) {
  pred.old = predict(fit, olddata, type = "terminalNodes")$predictions
  pred.new = predict(fit, newdata, type = "terminalNodes")$predictions
  prox = matrix(NA, nrow=nrow(pred.new), ncol=nrow(pred.old))
  ntree = ncol(pred.old)
  n1 = nrow(pred.new)
  n2 = nrow(pred.old)
  
  if (is.null(fit$inbag.counts)) {
    stop("call ranger with keep.inbag = TRUE")
  }
  
  # Get inbag counts
  inbag = simplify2array(fit$inbag.counts)
  
  for (i in 1:n1) {
    for (j in 1:n2) {
      # Use only trees where both obs are OOB
      if(is.null(obs_id)){
        tree_idx = inbag[j, ] == 0
      } else {
        tree_idx = inbag[j, ] == 0 & inbag[obs_id, ] == 0
      }
      prox[i, j] = sum(pred.new[i,tree_idx]==pred.old[j,tree_idx])/sum(tree_idx)
    }
  }
  
  prox
}


# Extracting the proximities between a psuedo point and every existing observation with the full proximity measure
extract_proximity_new_full = function(fit, olddata, newdata,...) {
  pred.old = predict(fit, olddata, type = "terminalNodes")$predictions
  pred.new = predict(fit, newdata, type = "terminalNodes")$predictions
  prox = matrix(NA, nrow=nrow(pred.new), ncol=nrow(pred.old))
  ntree = fit$num.trees
  n1 = nrow(pred.new)
  n2 = nrow(pred.old)
  

  
  for (i in 1:n1) {
    for (j in 1:n2) {
      prox[i, j] = sum(pred.new[i,]==pred.old[j,])/ntree
    }
  }
  
  prox
}

# Extracting the Monte Carlo realizations of the full proximities between every pair of observations
extract_proximity_full_monte_carlo = function(fit, olddata, K=10) {
  pred = predict(fit, olddata, type = "terminalNodes")$predictions
  prox  <- prox.tmp <-  matrix(NA, nrow(pred), nrow(pred))
  ntree = ncol(pred)
  n = nrow(prox)
  B = ncol(pred)
  
  if (is.null(fit$inbag.counts)) {
    stop("call ranger with keep.inbag = TRUE")
  }
  
  # Get inbag counts
  #  inbag = simplify2array(fit$inbag.counts)
  
  for (i in 1:n) {
    for (j in 1:n) {
      prox[i, j] = sum(pred[i, ] == pred[j, ]) / B
    }
  }
  
  prox.monte.array <- array(NA, c(n,n,K))
  for (k in 1:K){
    pred.tmp <- pred[,sample(1:B,replace=TRUE)]
    for (i in 1:n) {
      for (j in 1:n) {
        prox.tmp[i, j] = sum(pred.tmp[i, ] == pred.tmp[j, ]) / B
      }
    }
    prox.monte.array[,,k] <- prox.tmp
  }
  
  
  output <- list(original=prox, monte.carlo=prox.monte.array)
  return(output)
  
}




