# function creating dataset ---------------------------------

library(sgPLSdevelop)
library(mixOmics)
library(pls)

data.create <- function(n = 40, p = 10, q = 1){
  X <- matrix(data = runif(n*p),n,p)
  U <- matrix(data = runif(q*p,-10,10), nrow = p, ncol = q)
  E <- matrix(data = rnorm(n*q),n,q)
  Y <- X%*%U + E
  D <- data.frame(X,Y)
  return(list(D = D,X = X,Y = Y))
}

jeu <- data.create(n=50, q=1)


train <- 1:40
test <- 41:50

X <- jeu$X[train,]
Y <- jeu$Y[train,]
X.test <- jeu$X[test,]
Y.test <- jeu$Y[test,]


result.pls <- PLS(X = X, Y = Y, ncomp = 10, mode = "regression")


result2.pls <- pls(X = X, Y = Y, ncomp = 10, mode = "regression")
colnames(X.test) <- colnames(result2.pls$X)


result3.pls <- plsr(Y ~ X, ncomp=10, validation = "LOO")
result3.pls$validation$PRESS

# Q2 function ---------

q2pls13 <- q2.pls13(result.pls)
q2pls13$q2

# MIX OMICS ----------------

perf <- perf(result2.pls, validation = "loo")
q2 <- perf$measures$Q2.total$values
press <- perf$measures$PRESS$values
rss <- perf$measures$RSS$values

perf2 <- perf.pls(result2.pls, validation = "loo")
q2 <- perf$measures$Q2.total$values
press <- perf$measures$PRESS$values
rss <- perf$measures$RSS$values

d <- data.frame(q2pls13$q2, q2$value)
matplot(d, ylab = "q2", ylim = c(-10,2), pch = 16, type ="b", col = c("blue", "green"))


# Ancienne version de mixOmics ---------

n <- nrow(X)
q2.pls(X,Y,ncomp = 10,mode = "regression",M = n,fold=1:n,max.iter = 500,tol = 10^(-6))


## 13e version ----------------

q2.pls13 <- function(object, mode = "regression", ncomp.max = object$ncomp){
  
  X <- object$X
  Y <- object$Y
  c <- object$mat.c
  d <- object$mat.d
  p <- ncol(X)
  q <- ncol(Y)
  
  if(nrow(X)!=nrow(Y)){ stop("X and Y have not the same number of rows.")}
  n <- nrow(X)
  
  q2 <- numeric(ncomp.max)
  PRESS <- numeric(ncomp.max)
  RSS <- numeric(ncomp.max)
  
  
  if(mode == "regression"){
    
    # RSS0 computing
    Y.mean <- t(matrix(colMeans(Y), nrow = q, ncol = n)) #mean for each column
    RSS0 <- sum(colSums((Y-Y.mean)^2))
    
    # prediction computing for press1
    Y.hat.press1 <- matrix(nrow = n, ncol = q) # prediction matrix for PRESS1
    
    Y_test <- list()
    Y_test[[1]] <- Y
    
    for(h in 1:ncomp.max){
      Y_test[[h+1]] <- matrix(nrow = n, ncol = q)
    }
    
    for(i in 1:n){
      
      # Training data
      X_train <- X[-i, , drop = FALSE]
      Y_train <- Y[-i, , drop = FALSE]
      Y_train_scaled <- scale(Y_train, center = TRUE, scale = TRUE)
      
      X_means <- apply(X_train,2,mean)
      X_sds <- apply(X_train,2,sd)
      Y_means <- apply(Y_train,2,mean)
      Y_sds <- apply(Y_train,2,sd)
      
      
      X_test <- X[i,]
      X_test_scaled <- (X_test - X_means) / X_sds
      
      for(h in 1:ncomp.max){
        
        # Center and scale training data
        X_train_scaled <- scale(X_train, center = TRUE, scale = TRUE)
        
        # training on the dataset without the ith individual
        model <- PLS(X = X_train_scaled, Y = Y_train_scaled, ncomp = 1, mode = "regression")
        
        # predictions on the ith individual
        Y_pred <- predict.PLS(model, newdata = X_test_scaled)$predict[,,1]
        Y_pred <- Y_pred * Y_sds + Y_means
        
        # deflation matrices
        ui <- model$loadings$X[,1]
        vi <- model$loadings$Y[,1]
        res.deflat <- step2.spls(X=X_train,Y=Y_train,ui,vi,mode="regression")
        X_train = res.deflat$X.h
        Y_train = res.deflat$Y.h
        Y_test[[h+1]][i,] <- Y_test[[h]][i,] - Y_pred
        
      }
      
    }
    
    for(h in 1:ncomp.max){
      
      # deflation matrices
      u <- object$loadings$X[,h]
      v <- object$loadings$Y[,h]
      res.deflat <- step2.spls(X=X,Y=Y,u,v,mode="regression")
      X = res.deflat$X.h
      Y = res.deflat$Y.h
      
      
      # RSS computing
      RSS[h] <- mean(colSums((Y)^2)) # why not the sum ?
      
      # PRESSh computing
      PRESS[h] <- mean(colSums((Y_test[[h+1]])^2)) # why not the sum ?
      
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


q2pls13 <- q2.pls13(result.pls)
q2pls13








