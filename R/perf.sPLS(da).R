# perf.sPLS function ----

perf.sPLS <- function(object, K=nrow(object$X), ncomp = object$ncomp){
  
  X <- object$X
  Y <- object$Y
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  keepX = object$keepX
  keepY = object$keepY
  
  # PART 1 : number of components selection
  
  ## prediction analysis
  err <- matrix(NA, nrow = K, ncol = ncomp)
  
  b <- floor(n/K) # block size
  ind <- sample(n)
  
  for(k in seq_len(K)){
    
    ## bloc definition
    ind.beg <- (k-1) * b + 1 # block k beginning
    ind.end <- k * b # block k end
    ind.test <- ind[ind.beg:ind.end]
    nk <- length(ind.test)
    X.train <- X[-ind.test,]
    Y.train <- Y[-ind.test,]
    X.test <- X[ind.test,]
    Y.test <- Y[ind.test,]
    modele <- sPLS(X = X.train, Y = Y.train, ncomp = ncomp, mode = "regression", keepX = keepX, keepY = keepY)  
    
    for(h in 1:ncomp){  
      
      ## predictions
      pred <- predict.sPLS(modele, newdata = X.test)$predict[,,h]
      err[k,h] <- sum(colSums(as.matrix((Y.test - pred)^2)))
      
    }
  }
  
  err.moy <- colSums(err)/b/K
  
  h.best <- min(which.min(err.moy))
  
  par(mfrow = c(1,2))
  plot(err.moy, col="blue", pch = 16, type = "b", main = "MSEP of the model", xlab = "number of components", ylab = "MSEP")
  
  
  # PART 2 : number of variables selection
  
  ## prediction analysis
  err <- matrix(0, nrow = q, ncol = p)
  
  for(k in seq_len(K)){
    
    ## bloc definition
    ind.beg <- (k-1) * b + 1 # block k beginning
    ind.end <- k * b # block k end
    ind.test <- ind[ind.beg:ind.end]
    X.train <- X[-ind.test,]
    Y.train <- Y[-ind.test,]
    X.test <- X[ind.test,]
    Y.test <- Y[ind.test,]
    
    
    for(i in 1:q){  
      
      for(j in 1:p){
        
        modele <- sPLS(X = X.train, Y = Y.train, ncomp = h.best, mode = "regression", keepX = rep(j,h.best), keepY = rep(i,h.best))  
        
        ## predictions
        pred <- predict.sPLS(modele, newdata = X.test)$predict[,,h.best]
        err[i,j] <- err[i,j] + sum(colSums(as.matrix((Y.test - pred)^2)))
        
      }
    }
  }
  
  err.moy2 <- err/b/K
  err.min <- min(err.moy2)
  
  # min error coordinates identification
  ind.best <- numeric(2)
  for(i in 1:q){
    for(j in 1:p){
      if(err.moy2[i,j] == err.min){ind.best <- c(j,i)}
    }
  }
  
  nb.keepX <- ind.best[1]
  nb.keepY <- ind.best[2]
  keepX.best <- rep(nb.keepX,h.best)
  keepY.best <- rep(nb.keepY,h.best)
  
  round.err <- round(err.moy2,3)
  
  xmin <- max(nb.keepX-2.5,0.5)
  xmax <- min(nb.keepX+2.5,p+0.5)
  ymin <- max(nb.keepY-4.5,0.5)
  ymax <- min(nb.keepY+4.5,q+0.5)
  
  if(q==1){
    plot(as.vector(err.moy2), col="darkorange", pch = 16, type = "b", main = "MSEP of the model", ylab = "MSEP", xlab = "number of variables selected in X (keepX)")
  }else{
    plot(0, main = "MSEP of the model", xlim = c(xmin,xmax), ylim = c(ymin,ymax), ylab = "number of variables selected in Y (keepY)", xlab = "number of variables selected in X (keepX)")
    text(gl(p,q), rep(1:q,p), labels = as.vector(round.err), col = (round.err==min(round.err))+1)
  }
  
  return(setNames(list(err.moy,h.best,err.moy2,keepX.best,keepY.best),c("MSEP.h","h.best","MSEP.q.p","keepX.best","keepY.best")))
  
}


# perf.sPLS function ------

perf2.sPLS <- function(object, K=nrow(object$X), ncomp = object$ncomp){
  
  X <- object$X
  Y <- object$Y
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  keepX = object$keepX
  keepY = object$keepY
  
  # PART 1 : number of components selection
  
  ## prediction analysis
  err <- matrix(NA, nrow = K, ncol = ncomp)
  
  b <- floor(n/K) # block size
  ind <- sample(n)
  
  for(k in seq_len(K)){
    
    ## bloc definition
    ind.beg <- (k-1) * b + 1 # block k beginning
    ind.end <- k * b # block k end
    ind.test <- ind[ind.beg:ind.end]
    nk <- length(ind.test)
    X.train <- X[-ind.test,]
    Y.train <- Y[-ind.test,]
    X.test <- X[ind.test,]
    Y.test <- Y[ind.test,]
    modele <- sPLS(X = X.train, Y = Y.train, ncomp = ncomp, mode = "regression", keepX = keepX, keepY = keepY)  
    
    for(h in 1:ncomp){  
      
      ## predictions
      pred <- predict.sPLS(modele, newdata = X.test)$predict[,,h]
      err[k,h] <- sum(colSums(as.matrix((Y.test - pred)^2)))
      
    }
  }
  
  err.moy <- colSums(err)/b/K
  
  h.best <- min(which.min(err.moy))
  
  par(mfrow = c(1,2))
  plot(err.moy, col="blue", pch = 16, type = "b", main = "MSEP of the model", xlab = "number of components", ylab = "MSEP")
  
  
  # PART 2 : number of variables selection
  
  ## prediction analysis
  err <- matrix(0, nrow = p, ncol = q)
  
  for(k in seq_len(K)){
    
    ## bloc definition
    ind.beg <- (k-1) * b + 1 # block k beginning
    ind.end <- k * b # block k end
    ind.test <- ind[ind.beg:ind.end]
    X.train <- X[-ind.test,]
    Y.train <- Y[-ind.test,]
    X.test <- X[ind.test,]
    Y.test <- Y[ind.test,]
    
    
    for(i in 1:p){  
      
      for(j in 1:q){
        
        modele <- sPLS(X = X.train, Y = Y.train, ncomp = h.best, mode = "regression", keepX = rep(i,h.best), keepY = rep(j,h.best))  
        
        ## predictions
        pred <- predict.sPLS(modele, newdata = X.test)$predict[,,h.best]
        err[i,j] <- err[i,j] + sum(colSums(as.matrix((Y.test - pred)^2)))
        
      }
    }
  }
  
  err.moy2 <- err/b/K
  err.min <- min(err.moy2)
  
  # min error coordinates identification
  ind.best <- numeric(2)
  for(i in 1:p){
    for(j in 1:q){
      if(err.moy2[i,j] == err.min){ind.best <- c(i,j)}
    }
  }
  
  
  keepX.best <- rep(ind.best[1],h.best)
  keepY.best <- rep(ind.best[2],h.best)
  round.err <- round(err.moy2,round(12/q))
  
  if(q==1){
    plot(err.moy2, col="orange", pch = 16, type = "b", main = "MSEP of the model", ylab = "MSEP", xlab = "number of variables selected in X (keepX)")
  }else{
    plot(0, main = "MSEP of the model", xlim = c(0.5,q+0.5), ylim = c(0.5,p+0.5), xlab = "number of variables selected in Y (keepY)", ylab = "number of variables selected in X (keepX)")
    text(gl(q,p), rep(1:p,q), labels = as.vector(round.err), col = (round.err==min(round.err))+1)
  }
  
  return(setNames(list(err.moy,h.best,err.moy2,keepX.best,keepY.best),c("MSEP.h","h.best","MSEP.p.q","keepX.best","keepY.best")))
}


# perf.sPLSda function -------

perf.sPLSda <- function(object, K=nrow(object$X), ncomp = object$ncomp, method = "max.dist", keepX = object$keepX){
  
  X <- object$X
  Y <- map(object$Y)
  n <- nrow(X)
  p <- ncol(X)
  
  if(method == "max.dist"){dist=1}else if(method == "centroids.dist"){dist=2}else{dist=3}
  
  # prediction analysis
  err <- matrix(NA, nrow = K, ncol = ncomp)
  
  b <- floor(n/K) # block size
  ind <- 1:n
  
  for(k in seq_len(K)){
    
    # bloc definition
    ind.beg <- (k-1) * b + 1 # block k beginning
    ind.end <- k * b # block k end
    ind.test <- ind[ind.beg:ind.end]
    X.train <- X[-ind.test,]
    Y.train <- Y[-ind.test]
    X.test <- X[ind.test,]
    Y.test <- Y[ind.test]
    
    for(h in 1:ncomp){
      # model created
      modele <- sPLSda(X = X.train,Y = Y.train, ncomp = h, keepX = keepX)
      pred <- predict.sPLSda(modele, newdata = X.test, methode = methode)$class[[dist]]
      equal <- Y.test == pred[,h]
      err[k,h] <- sum(1-equal)
      
    }
  }
  
  err.moy <- colSums(err)/b/K
  
  h.best <- min(which.min(err.moy))
  
  par(mfrow = c(1,2))
  plot(err.moy, col="blue", pch = 16, type = "b", main = "Error rate of the model", xlab = "number of components", ylab = "Error")
  
  # PART 2 : number of variables selection
  
  ## prediction analysis
  err <- numeric(p)
  
  for(k in seq_len(K)){
    
    ## bloc definition
    ind.beg <- (k-1) * b + 1 # block k beginning
    ind.end <- k * b # block k end
    ind.test <- ind[ind.beg:ind.end]
    X.train <- X[-ind.test,]
    Y.train <- Y[-ind.test]
    X.test <- X[ind.test,]
    Y.test <- Y[ind.test]
    
    
    for(i in 1:p){  
        
      modele <- sPLS(X = X.train, Y = Y.train, ncomp = h.best, mode = "regression", keepX = rep(i,h.best))  
        
      ## predictions
      pred <- predict.sPLS(modele, newdata = X.test)$predict[,,h.best]
      err[i] <- err[i] + sum(colSums(as.matrix((Y.test - pred)^2)))

    }
  }
  
  err.moy2 <- err/b/K
  err.min <- min(err.moy2)
  
  # min error coordinates identification
  
  nb.keepX <- which.min(err)
  keepX.best <- rep(nb.keepX,h.best)
  
  plot(err.moy2, col="orange", pch = 16, type = "b", main = "MSEP of the model", ylab = "MSEP", xlab = "number of variables selected in X (keepX)")
  
  return(setNames(list(err.moy,h.best,err.moy2,keepX.best),c("error.h","h.best","error.p","keepX.best")))
}
