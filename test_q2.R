
# function creating dataset ---------------------------------

create.data <- function(n = 40, p = 10, q = 1){
  X <- matrix(data = runif(n*p),n,p)
  U <- matrix(data = runif(q*p,-10,10), nrow = p, ncol = q)
  E <- matrix(data = rnorm(n*q),n,q)
  Y <- X%*%U + E
  D <- data.frame(X,Y)
  return(list(D = D,X = X,Y = Y))
}

jeu <- create.data(n=50, q=3)

X <- jeu$X
Y <- jeu$Y

result.pls <- PLS(X = X, Y = Y, ncomp = 10, mode = "regression")
result2.pls <- pls(X = X, Y = Y, ncomp = 10, mode = "regression")

q2pls8 <- q2.pls8(result.pls)
q2pls8$q2

# MIX OMICS ----------------

perf <- perf(result2.pls, validation = "loo")
q2 <- perf$measures$Q2.total$values
press <- perf$measures$PRESS$values
rss <- perf$measures$RSS$values

d <- data.frame(q2pls8$q2, q2$value)
matplot(d, ylab = "q2", ylim = c(-10,2), pch = 16, type ="b", col = c("blue", "green"))























q2.pls8 <- function(object, mode = "regression", ncomp.max = object$ncomp){
  
  # object attributes recovery (1)
  X <- scale(object$X)
  Y <- scale(object$Y)
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
    
    # h-dependant variables definition (4)
    Y.hat.rss <- matrix(nrow = n, ncol = q) # prediction matrix for RSS
    
    for(h in 1:ncomp.max){ # (3)
      
      Y.test <- Y  
      
      for(i in 1:n){
        
        model <- PLS(X = X[-i,], Y = Y[-i,], ncomp = h, mode = "regression")
        a <- model$loadings$X[, 1, drop = FALSE]
        d <- model$mat.d[, 1, drop = FALSE]
        Y.test[i, ] <- Y.test[i, ] - X[i, , drop = FALSE] %*% a %*% t(d)
        
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

q2.pls8(result.pls)$q2


