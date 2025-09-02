tuning.sPLS.XY <- function(X,Y,folds=10,validation=c("Mfold","loo"), ncomp = ncol(X)){
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  # conditions check-up
  
  if(ncomp > p || ncomp <= 0){ stop(paste("ncomp.max must be a value between 0 and",p,"which is the total number of variables in the X matrix."))}
  
  if(validation[1] == "Mfold"){
    if(folds < 2 || folds > n){ stop(paste("folds must be a value between 2 and",n))}
    K = folds
  }else{K = n}
  
  # PART 1 : number of components selection
  
  model <- sPLS(X,Y,ncomp = ncomp)
  perf <- perf2.sPLS(model, criterion = "MSEP")
  err.moy <- colSums(perf$MSEP)
  h.best <- which.min(err.moy)
  
  par(mfrow = c(1,2))
  plot(err.moy, col="blue", pch = 16, type = "b", main = "MSEP of the model", xlab = "number of components", ylab = "MSEP")
  
  
  # PART 2 : number of variables selection

  ## prediction analysis
  err <- matrix(0, nrow = q, ncol = p)
  
  b <- floor(n/K) # block size
  ind <- 1:n
  
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
    plot(as.vector(err.moy2), col="darkorange", pch = 16, type = "b", main = "MSEP of the model", ylab = "MSEP", xlab = "Number of variables selected in X (keepX)")
  }else{
    plot(0, main = "MSEP of the model", xlim = c(xmin,xmax), ylim = c(ymin,ymax), ylab = "Number of variables selected in Y (keepY)", xlab = "Number of variables selected in X (keepX)")
    text(gl(p,q), rep(1:q,p), labels = as.vector(round.err), col = (round.err==min(round.err))+1)
    abline(v = xmin:xmax, h = ymin:ymax, lty = 3, col = "yellow")
  }
  
  return(setNames(list(err.moy,h.best,err.moy2,keepX.best,keepY.best),c("MSEP.h","h.best","MSEP.q.p","keepX.best","keepY.best")))
 
}

tuning0.sPLS.XY <- function(object, K=nrow(object$X), ncomp = object$ncomp){
  
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
    plot(as.vector(err.moy2), col="darkorange", pch = 16, type = "b", main = "MSEP of the model", ylab = "MSEP", xlab = "Number of variables selected in X (keepX)")
  }else{
    plot(0, main = "MSEP of the model", xlim = c(xmin,xmax), ylim = c(ymin,ymax), ylab = "Number of variables selected in Y (keepY)", xlab = "Number of variables selected in X (keepX)")
    text(gl(p,q), rep(1:q,p), labels = as.vector(round.err), col = (round.err==min(round.err))+1)
    abline(v = xmin:xmax, h = ymin:ymax, lty = 3, col = "yellow")
  }
  
  return(setNames(list(err.moy,h.best,err.moy2,keepX.best,keepY.best),c("MSEP.h","h.best","MSEP.q.p","keepX.best","keepY.best")))
  
  
}
