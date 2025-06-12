# Q2 function 

q2.pls <- function(X,Y, mode = "regression", ncomp.max = 10){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  p <- ncol(X)
  q <- ncol(Y)
  if(nrow(X)!=nrow(Y)){ stop("X and Y have not the same number of rows.")}
  n <- nrow(X)
  X0 <- X
  Y0 <- Y
  q2 <- numeric(ncomp.max)
  
  if(mode == "regression"){
    
    # RSS0 computing
    Y.mean <- t(matrix(colMeans(Y), nrow = q, ncol = n)) #mean for each column
    RSS0 <- sum(colSums((Y0-Y.mean)^2))
    model.all <- PLS(X = X, Y = Y, ncomp = ncomp.max, mode = mode)
    
  
    for(h in 1:ncomp.max){
      
      Y.hat.rss <- matrix(nrow = n, ncol = q) # prediction matrix for RSS
      Y.hat.press <- matrix(nrow = n, ncol = q) # prediction matrix for PRESS
      
      for(i in 1:n){
        
        
        # training on the dataset without the ith individual
        model <- PLS(X = X[-i,], Y = Y[-i,], ncomp = h, mode = "regression")

        # predictions on the ith individual
        Y.hat.rss[i,] <- predict.PLS(model.all, newdata = X[i,])$predict[,,max(h-1,1)]
        Y.hat.press[i,] <- predict.PLS(model, newdata = X[i,])$predict[,,h]
        
        
        
      }
      # RSS computing
      RSS <- sum(colSums((Y-Y.hat.rss)^2))
      
      # PRESS computing
      PRESS <- sum(colSums((Y-Y.hat.press)^2))
      
      if(h == 1){
          PRESS1 <- PRESS
      }
      
      # new Y matrice
      Y <- Y - Y.hat.rss
      
      # Q2
      q2[h] <- 1-PRESS/RSS
      
    }# end h loop
    
  }else if(mode == "canonical"){ # CASE OF CANONICAL MODE
    
    # RSS0 computing
    X.mean <- t(matrix(colMeans(X), nrow = p, ncol = n)) #mean for each column
    RSS0 <- sum(colSums((X0-X.mean)^2))
    model.all <- PLS(X = X, Y = Y, ncomp = ncomp.max, mode = "canonical")
    mat.t.all <- as.matrix(model.all$variates$X)
    mat.c.all <- as.matrix(model.all$mat.c)
    
    for(h in 1:ncomp.max){
      
      X.hat.rss <- matrix(nrow = n, ncol = p) # prediction matrix for RSS
      X.hat.press <- matrix(nrow = n, ncol = p) # prediction matrix for PRESS
      
      for(i in 1:n){
        
        # training on the dataset without the ith individual
        model <- PLS(X = X[-i,], Y = Y[-i,], ncomp = h, mode = "canonical")
        mat.t <- as.matrix(model$variates$X)
        mat.c <- as.matrix(model$mat.c)
        load.u <- as.matrix(model$loadings$X)

        # predictions on the ith individual
        X.hat.rss <- t(t(mat.t.all[,max(h-1,1)]))%*%t(mat.c.all[,max(h-1,1)])
        X.hat.press[i,] <- X[i,]%*%load.u[,h]%*%t(mat.c[,h])
        
      } # end i loop
      
      # RSS computing
      RSS <- sum(colSums((X-X.hat.rss)^2))
      
      # PRESS computing
      PRESS <- sum(colSums((X-X.hat.press)^2))
      
      if(h == 1){
        RSS <- RSS0
      }
      
      # new Y matrice
      X <- X - X.hat.rss
      
      # Q2
      q2[h] <- 1-PRESS/RSS

      
    }# end h loop
    
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
  
  q2.pls.results <- list(q2 = q2, suggestion = suggestion, h.best = h.best)
  return(q2.pls.results)
  
}





