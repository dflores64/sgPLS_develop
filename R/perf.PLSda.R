perf.PLSda <- function(object, ncomp = object$ncomp, method = "max.dist"){
  
  X <- object$X
  Y <- map(object$Y)
  n <- nrow(X)
  p <- ncol(X)
  
  if(method == "max.dist"){k=1}else if(method == "centroids.dist"){k=2}else{k=3}
  
  # prediction analysis
  err <- matrix(NA, nrow = n, ncol = ncomp)
  
  matconf <- list()
  
  for(h in 1:ncomp){

    for(i in 1:n){
      
      X.train <- X[-i,]
      Y.train <- Y[-i]
      X.test <- X[i,]
      Y.test <- Y[i]
      
      # model created
      modele <- PLSda(X = X.train,Y = Y.train, ncomp = ncomp)
      pred <- predict.PLSda(modele, newdata = X.test, methode = methode)$class[[k]]
      equal <- Y.test == pred[,h]
      err[i,h] <- 1-equal
      
    }
  }
  
  err.moy <- colMeans(err)
  h.best <- min(which.min(err.moy))
  plot(err.moy, col="blue", pch = 16, type = "b", main = "Error rate of the model", xlab = "number of components", ylab = "Error")
  
  return(setNames(list(err.moy,h.best),c("error","h.best")))
}

perf2.PLSda <- function(object, ncomp = object$ncomp, method = "max.dist"){
  
  X <- object$X
  Y <- map(object$Y)
  n <- nrow(X)
  p <- ncol(X)
  
  if(method == "max.dist"){k=1}else if(method == "centroids.dist"){k=2}else{k=3}
  
  # prediction analysis
  err.moy <- numeric(ncomp)
  matconf <- list()
  
  for(h in 1:ncomp){
    err.moy[h] <- CV(X,Y,mode = "class",ncomp = h,K=n)
  }
  
  h.best <- min(which.min(err.moy))
  plot(err.moy, col="blue", pch = 16, type = "b", main = "Error rate of the model", xlab = "number of components", ylab = "Error")
  
  return(setNames(list(err.moy,h.best),c("error","h.best")))
}
