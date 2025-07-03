perf.PLSda <- function(object, K=nrow(object$X), ncomp = object$ncomp, method = "max.dist"){
  
  X <- object$X
  Y <- map(object$Y)
  n <- nrow(X)
  p <- ncol(X)
  
  if(method == "max.dist"){dist=1}else if(method == "centroids.dist"){dist=2}else{dist=3}
  
  # prediction analysis
  err <- matrix(NA, nrow = K, ncol = ncomp)

  b <- floor(n/K) # block size
  ind <- sample(n)
  
  ## Indices for beginning and end of each block
  bloc.ind <- matrix(NA, nrow = K, ncol = 2)
  for (k in seq_len(K)) {
    bloc.ind[k,1] <- (k-1) * b + 1 # block k beginning
    bloc.ind[k,2] <- k * b # block k end
  }
  
  for(h in 1:ncomp){

    for(k in seq_len(K)){
      
      ind.test <- ind[bloc.ind[k,1]:bloc.ind[k,2]]
      nk <- length(ind.test)
      X.train <- X[-ind.test,]
      Y.train <- Y[-ind.test]
      X.test <- X[ind.test,]
      Y.test <- Y[ind.test]
      
      # model created
      modele <- PLSda(X = X.train,Y = Y.train, ncomp = h)
      pred <- predict.PLSda(modele, newdata = X.test, methode = methode)$class[[dist]]
      equal <- Y.test == pred[,h]
      err[k,h] <- sum(1-equal)

      
    }
  }
  
  err.moy <- colSums(err)/b/K

  h.best <- min(which.min(err.moy))
  plot(err.moy, col="blue", pch = 16, type = "b", main = "Error rate of the model", xlab = "number of components", ylab = "Error")
  
  return(setNames(list(err.moy,h.best),c("error","h.best")))
}
