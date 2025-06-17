
# q2.pls5a ----------------

q2.pls5a <- function(object, mode = "regression", ncomp.max = object$ncomp){
  
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
  PRESS <- numeric(ncomp.max)
  RSS <- numeric(ncomp.max)
  Y.test <- Y
  
  if(mode == "regression"){ 
    
    # RSS0 computing
    Y.mean <- t(matrix(colMeans(Y), nrow = q, ncol = n)) #mean for each column
    RSS0 <- sum(colSums((Y-Y.mean)^2))
    
    # prediction computing for press1
    Y.hat.press1 <- matrix(nrow = n, ncol = q) # prediction matrix for PRESS1
    
    
    
    for(h in 1:ncomp.max){ # (3)
      
      # h-dependant variables definition (4)
      Y.hat.rss <- matrix(nrow = n, ncol = q) # prediction matrix for RSS
      
      for(i in 1:n){ # (5)
        
        # training on the dataset without the ith individual (6)
        model <- PLS(X = X[-i,], Y = Y[-i,], ncomp = 1, mode = "regression")
        
        # predictions on the ith individual (7)
        Y.hat.press1[i,] <- predict.PLS(model, newdata = X[i,])$predict[,,1]
        
      } # (8)
      
      
      # deflation matrices (9)
      u <- object$loadings$X[,h]
      v <- object$loadings$Y[,h]
      res.deflat <- step2.spls(X=X,Y=Y,u,v,mode="regression")
      X = res.deflat$X.h
      Y = res.deflat$Y.h
      Y.test <- Y.test - Y.hat.press1
      
      # RSS computing (10)
      RSS[h] <- mean(colSums((Y)^2)) # why not the sum ?
      
      # PRESSh computing (11)
      PRESS[h] <- mean(colSums((Y.test)^2)) # why not the sum ?
      
      # Q2 (12)
      q2[h] <- 1-PRESS[h]/RSS[max(h-1,1)]

    }# end h loop
    
    # first value correction 
    q2[1] <- 1-PRESS[1]/RSS0
    
  }else if(mode == "canonical"){ # CASE OF CANONICAL MODE 
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

q2pls <- q2.pls5a(result.pls)


# q2.pls5b -------------

q2.pls5b <- function(object, mode = "regression", ncomp.max = object$ncomp){
  
  # object attributes recovery (1)
  X <- scale(object$X0)
  Y <- scale(object$Y0)
  c <- object$mat.c
  d <- object$mat.d
  p <- ncol(X)
  q <- ncol(Y)
  
  if(nrow(X)!=nrow(Y)){ stop("X and Y have not the same number of rows.")}
  n <- nrow(X)
  
  # variables definition (2)
  q2 <- numeric(ncomp.max)
  PRESS <- numeric(ncomp.max)
  RSS <- numeric(ncomp.max)
  Y.test <- Y
  
  if(mode == "regression"){
    
    # RSS0 computing
    Y.mean <- t(matrix(colMeans(Y), nrow = q, ncol = n)) #mean for each column
    RSS0 <- sum(colSums((Y-Y.mean)^2))
    
    # prediction computing for press1
    Y.hat.press1 <- matrix(nrow = n, ncol = q) # prediction matrix for PRESS1
    
    
    
    for(h in 1:ncomp.max){ # (3)
      
      # h-dependant variables definition (4)
      Y.hat.rss <- matrix(nrow = n, ncol = q) # prediction matrix for RSS
      
      
      for(i in 1:n){
        
        # training on the dataset without the ith individual
        model <- PLS(X = X[-i,], Y = Y[-i,], ncomp = 1, mode = "regression")
        
        # predictions on the ith individual
        Y.hat.press1[i,] <- predict.PLS(model, newdata = X[i,])$predict[,,1]
        
      }
      
      
      # deflation matrices
      u <- object$loadings$X[,h]
      v <- object$loadings$Y[,h]
      res.deflat <- step2.spls(X=X,Y=Y,u,v,mode="regression")
      X = res.deflat$X.h
      Y = res.deflat$Y.h
      Y.test <- Y.test - Y.hat.press1
      
      # RSS computing
      RSS[h] <- sum(colSums((Y)^2)) # why not the sum ? So we do the sum.
      
      # PRESSh computing
      PRESS[h] <- sum(colSums((Y.test)^2)) # why not the sum ? So we do the sum.
      
      # Q2
      q2[h] <- 1-PRESS[h]/RSS[max(h-1,1)]
      
    }# end h loop
    
    # first value correction 
    q2[1] <- 1-PRESS[1]/RSS0
    
  }else if(mode == "canonical"){ # CASE OF CANONICAL MODE 
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

q2pls <- q2.pls5b(result.pls)