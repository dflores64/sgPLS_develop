q2.pls8 <- function(object, mode = "regression", ncomp.max = object$ncomp){
  
  # object attributes recovery
  X <- scale(object$X)
  Y <- scale(object$Y)
  c <- object$mat.c
  d <- object$mat.d
  p <- ncol(X)
  q <- ncol(Y)
  
  if(nrow(X)!=nrow(Y)){ stop("X and Y have not the same number of rows.")}
  n <- nrow(X)
  
  #  variables definition 
  q2 <- numeric(ncomp.max)
  PRESS <- numeric(ncomp.max)
  RSS <- numeric(ncomp.max)
  Y.test <- Y
  
  if(mode == "regression"){ # CASE OF REGRESSION MODE 
    
    # RSS0 computing
    Y.mean <- t(matrix(colMeans(Y), nrow = q, ncol = n)) #mean for each column
    RSS0 <- sum(colSums((Y-Y.mean)^2))
    
    for(h in 1:ncomp.max){ # 
      
      Y.test <- Y  
      
      for(i in 1:n){
        
        model <- PLS(X = X[-i,], Y = Y[-i,], ncomp = h, mode = "regression")
        a <- model$loadings$X[, 1, drop = FALSE]
        d <- model$mat.d[, 1, drop = FALSE]
        Y.test[i, ] <- Y.test[i, ] - X[i, , drop = FALSE] %*% a %*% t(d)
        
      }
      
      # deflation matrices for RSS
      u <- object$loadings$X[,h]
      v <- object$loadings$Y[,h]
      res.deflat <- step2.spls(X=X,Y=Y,u,v,mode="regression")
      X = res.deflat$X.h
      Y = res.deflat$Y.h
      
      
      # RSS computing 
      RSS[h] <- mean(colSums((Y)^2)) # why not the sum ?
      
      # PRESSh computing
      PRESS[h] <- mean(colSums((Y.test)^2)) # why not the sum ?
      
      # Q2 
      q2[h] <- 1-PRESS[h]/RSS[max(h-1,1)]
      
    }# end h loop
    
    # first value correction 
    q2[1] <- 1-PRESS[1]/RSS0
    
  }else if(mode == "canonical"){ # CASE OF CANONICAL MODE 
    
    # RSS0 computing
    X.mean <- t(matrix(colMeans(X), nrow = p, ncol = n)) #mean for each column
    
    RSS0 <- sum(colSums((X-X.mean)^2))
    
    for(h in 1:ncomp.max){ # (3)
      
      X.test <- X  
      
      for(i in 1:n){
        
        model <- PLS(X = X[-i,], Y = Y[-i,], ncomp = h, mode = "regression")
        a <- model$loadings$X[, 1, drop = FALSE]
        c <- model$mat.c[, 1, drop = FALSE]
        X.test[i, ] <- X.test[i, ] - X[i, , drop = FALSE] %*% a %*% t(c)
        
      }
      
      # deflation matrices for RSS
      u <- object$loadings$X[,h]
      v <- object$loadings$Y[,h]
      res.deflat <- step2.spls(X=X,Y=Y,u,v,mode="regression")
      X = res.deflat$X.h
      Y = res.deflat$Y.h
      
      
      # RSS computing 
      RSS[h] <- mean(colSums((X)^2)) # why not the sum ?
      
      # PRESSh computing
      PRESS[h] <- mean(colSums((X.test)^2)) # why not the sum ?
      
      # Q2 
      q2[h] <- 1-PRESS[h]/RSS[max(h-1,1)]
      
    }# end h loop
    
    # first value correction 
    q2[1] <- 1-PRESS[1]/RSS0
    
  }else{
    stop("The mode must be either regression or canonical.")
  }# end if loop
  
  lim <- 0.0975
  h.best <- min(which(q2 > lim))
  
  # Plot
  plot(q2, type = "b", col = "blue", pch = 16,
       main = "Q² performance according to the number of components",
       xlab = "Number of components", ylab = "Q²")
  abline(h = 0.0975, col = "red", lty = 2)
  
  suggestion <- paste("best number of components : H =",h.best)
  
  q2.pls.results <- list(q2 = q2, suggestion = suggestion, h.best = h.best, PRESS = PRESS, RSS = RSS)
  return(q2.pls.results)
  
}
