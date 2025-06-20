## 9e version (cf 8e version) : scaling ----------------

q2.pls9 <- function(object, mode = "regression", ncomp.max = object$ncomp){
  
  # object attributes recovery (1)
  X <- scale(object$X0)
  Y <- scale(object$Y0)
  c <- object$mat.c
  d <- object$mat.d
  p <- ncol(X)
  q <- ncol(Y)
  
  if(nrow(X)!=nrow(Y)){ stop("X and Y have not the same number of rows.")}
  n <- nrow(X)
  
  #  variables definition (2)
  q2 <- numeric(ncomp.max)
  MSEP <- numeric(ncomp.max) 
  PRESS <- numeric(ncomp.max)
  RSS <- numeric(ncomp.max)
  Y.test <- Y
  
  if(mode == "regression"){ 
    
    # RSS0 computing
    Y.mean <- t(matrix(colMeans(Y), nrow = q, ncol = n)) #mean for each column
    RSS0 <- sum(colSums((Y-Y.mean)^2))
    
    for(h in 1:ncomp.max){ # (3)
      
      for(i in 1:n){
        
        # Training data
        X_train <- X[-i, , drop = FALSE]
        Y_train <- Y[-i, , drop = FALSE]
        
        # Test sample
        X_test <- X[i, , drop = FALSE]
        Y_test <- Y[i, , drop = FALSE]
        
        # Center and scale training data
        X_train_scaled <- scale(X_train, center = TRUE, scale = TRUE)
        Y_train_scaled <- scale(Y_train, center = TRUE, scale = TRUE)
        
        # Save means and sds
        X_means <- attr(X_train_scaled, "scaled:center")
        X_sds   <- attr(X_train_scaled, "scaled:scale")
        Y_means <- attr(Y_train_scaled, "scaled:center")
        Y_sds   <- attr(Y_train_scaled, "scaled:scale")
        
        # Scale test point using training stats
        X_test_scaled <- scale(X_test, center = X_means, scale = X_sds)
        
        model <- PLS(X = X_train_scaled, Y = Y_train_scaled, ncomp = h, mode = "regression")
        
        a <- model$loadings$X[, 1, drop = FALSE]
        d <- model$mat.d[, 1, drop = FALSE]
        
        Y_test_scaled_pred <- X_test_scaled %*% a %*% t(d)
        Y_pred <- Y_test_scaled_pred * Y_sds + Y_means
        Y.test[i, ] <- Y.test[i, ] - Y_pred
        
      }
      
      # deflation matrices for RSS(9)
      u <- object$loadings$X[,h]
      v <- object$loadings$Y[,h]
      res.deflat <- step2.spls(X=X,Y=Y,u,v,mode="regression")
      X = res.deflat$X.h
      Y = res.deflat$Y.h
      
      
      # RSS computing (10)
      RSS[h] <- sum(colSums((Y)^2)) # why not the sum ?
      
      # PRESSh computing (11)
      PRESS[h] <- sum(colSums((Y.test)^2)) # why not the sum ?
      
      # Q2 (12)
      q2[h] <- 1-PRESS[h]/RSS[max(h-1,1)]
      
      # MSEP 
      MSEP[h] <- RSS[h]/(n*q)
      
      
    }# end h loop
    
    # first value correction 
    q2[1] <- 1-PRESS[1]/RSS0
    
  }else if(mode == "canonical"){ # CASE OF CANONICAL MODE 
    
    X.test <- X
    
    # RSS0 computing
    X.mean <- t(matrix(colMeans(X), nrow = p, ncol = n)) #mean for each column
    RSS0 <- sum(colSums((X-X.mean)^2))
    
    for(h in 1:ncomp.max){ # (3)
      
      for(i in 1:n){
        
        # Training data
        X_train <- X[-i, , drop = FALSE]
        Y_train <- Y[-i, , drop = FALSE]
        
        # Test sample
        X_test <- X[i, , drop = FALSE]
        Y_test <- Y[i, , drop = FALSE]
        
        # Center and scale training data
        X_train_scaled <- scale(X_train, center = TRUE, scale = TRUE)
        Y_train_scaled <- scale(Y_train, center = TRUE, scale = TRUE)
        
        # Save means and sds
        X_means <- attr(X_train_scaled, "scaled:center")
        X_sds   <- attr(X_train_scaled, "scaled:scale")
        Y_means <- attr(Y_train_scaled, "scaled:center")
        Y_sds   <- attr(Y_train_scaled, "scaled:scale")
        
        # Scale test point using training stats
        X_test_scaled <- scale(X_test, center = X_means, scale = X_sds)
        
        model <- PLS(X = X_train_scaled, Y = Y_train_scaled, ncomp = h, mode = "regression")
        
        a <- model$loadings$X[, 1, drop = FALSE]
        c <- model$mat.c[, 1, drop = FALSE]
        
        X_test_scaled_pred <- X_test_scaled %*% a %*% t(c)
        X_pred <- X_test_scaled_pred * X_sds + X_means
        X.test[i, ] <- X.test[i, ] - X_pred
        
      }
      
      # deflation matrices for RSS(9)
      u <- object$loadings$X[,h]
      v <- object$loadings$Y[,h]
      res.deflat <- step2.spls(X=X,Y=Y,u,v,mode="regression")
      X = res.deflat$X.h
      Y = res.deflat$Y.h
      
      # RSS computing (10)
      RSS[h] <- sum(colSums((X)^2)) # why not the sum ?
      
      # PRESSh computing (11)
      PRESS[h] <- sum(colSums((X.test)^2)) # why not the sum ?
      
      # Q2 (12)
      q2[h] <- 1-PRESS[h]/RSS[max(h-1,1)]
      
      # MSEP 
      MSEP[h] <- RSS[h]/(n*q)
      
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
  
  q2.pls.results <- list(q2 = q2, suggestion = suggestion, h.best = h.best, PRESS = PRESS, RSS = RSS, MSEP = MSEP)
  return(q2.pls.results)
  
}

q2pls9 <- q2.pls9(result.pls)
q2pls9
