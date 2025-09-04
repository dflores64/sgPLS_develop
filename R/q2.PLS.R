q2.PLS <- function(object, ncomp.max = object$ncomp, mode = "regression", plot = TRUE){
  
  X <- object$X
  Y <- object$Y
  c <- object$mat.c
  d <- object$mat.d
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  ncomp.object <- object$ncomp
  
  # conditions check-up
  if(!("pls" %in% class(object)) && class(object) != "mixo_pls"){ stop("object class must either contain pls class or be mixo_pls class."); print(class(object))}
  
  if(ncomp.max > object$ncomp || ncomp.max <= 0){ stop(paste("ncomp.max must be set up between 0 and",object$ncomp,"which is the total number of components computed in the object model."))}
  
  if(mode != "regression" && mode != "canonical"){ stop("mode must be either << regression >> or << canonical >>")}
  
  q2 <- numeric(ncomp.max)
  PRESS <- numeric(ncomp.max)
  RSS <- numeric(ncomp.max)
  
  
  if(mode == "regression"){
    
    RSSj <- matrix(nrow = ncomp.max, ncol = q)
    PRESSj <- matrix(nrow = ncomp.max, ncol = q)
    
    # RSS0 computing
    Y.mean <- t(matrix(colMeans(Y), nrow = q, ncol = n)) #mean for each column
    RSS0 <- sum(colSums((Y-Y.mean)^2))
    
    
    Y_test <- list()
    Y_test[[1]] <- Y
    
    for(h in 1:ncomp.max){
      Y_test[[h+1]] <- matrix(nrow = n, ncol = q)
    }
    
    for(i in 1:n){
      
      # Training data
      X_train <- X[-i, , drop = FALSE]
      Y_train <- Y[-i, , drop = FALSE]
      
      # Center and scale training data
      X_train_scaled <- scale(X_train, center = TRUE, scale = TRUE)
      Y_train_scaled <- scale(Y_train, center = TRUE, scale = TRUE)
      
      X_means <- apply(X_train,2,mean)
      X_sds <- apply(X_train,2,sd)
      Y_means <- apply(Y_train,2,mean)
      Y_sds <- apply(Y_train,2,sd)
      
      
      X_test <- X[i,]
      X_test_scaled <- (X_test - X_means) / X_sds
      
      for(h in 1:ncomp.max){
        
        # training on the dataset without the ith individual
        model <- PLS(X = X_train_scaled, Y = Y_train_scaled, ncomp = 1, mode = "regression")
        
        # predictions on the ith individual
        Y_pred <- predict.PLS(model, newdata = X_test_scaled)$predict[,,1]
        Y_pred <- Y_pred * Y_sds + Y_means
        
        # deflation matrices
        ui <- model$loadings$X[,1]
        vi <- model$loadings$Y[,1]
        res.deflat <- step2.spls(X=X_train_scaled,Y=Y_train_scaled,ui,vi,mode="regression")
        X_train_scaled = res.deflat$X.h
        Y_train_scaled = res.deflat$Y.h
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
      
      for(j in 1:q){
        RSSj[h,j] <- sum((Y[,j])^2)
        PRESSj[h,j] <- sum((Y_test[[h+1]][,j])^2)
      }
      colnames(RSSj) <- paste0("Y",1:q)
      colnames(PRESSj) <- paste0("Y",1:q)
      
      # RSS computing
      RSS[h] <- sum(RSSj[h,]) 
      
      # PRESSh computing
      PRESS[h] <- sum(PRESSj[h,]) 
      
      # Q2
      q2[h] <- 1-PRESS[h]/RSS[max(h-1,1)]
      
      
    }# end h loop
    
    # first value correction 
    q2[1] <- 1-PRESS[1]/RSS0
    
  }else if(mode == "canonical"){ # CASE OF CANONICAL MODE 
    
    RSSj <- matrix(nrow = ncomp.max, ncol = p)
    PRESSj <- matrix(nrow = ncomp.max, ncol = p)
    
    # RSS0 computing
    X.mean <- t(matrix(colMeans(X), nrow = p, ncol = n)) #mean for each column
    RSS0 <- sum(colSums((X-X.mean)^2))

    X_press <- list()
    X_press[[1]] <- X
    
    for(h in 1:ncomp.max){
      X_press[[h+1]] <- matrix(nrow = n, ncol = p)
    }
    
    for(i in 1:n){
      
      # Training data
      X_train <- X[-i, , drop = FALSE]
      Y_train <- Y[-i, , drop = FALSE]
      
      # Center and scale training data
      X_train_scaled <- scale(X_train, center = TRUE, scale = TRUE)
      Y_train_scaled <- scale(Y_train, center = TRUE, scale = TRUE)
      
      X_means <- apply(X_train,2,mean)
      X_sds <- apply(X_train,2,sd)
      Y_means <- apply(Y_train,2,mean)
      Y_sds <- apply(Y_train,2,sd)
      
      
      X_test <- X[i,]
      X_test_scaled <- (X_test - X_means) / X_sds
      
      for(h in 1:ncomp.max){
        
        model <- PLS(X = X_train_scaled, Y = Y_train_scaled, ncomp = h, mode = "regression")
        
        a <- model$loadings$X[, 1, drop = FALSE]
        c <- model$mat.c[, 1, drop = FALSE]
        
        X_test_scaled_pred <- X_test_scaled %*% a %*% t(c)
        X_pred <- X_test_scaled_pred * X_sds + X_means
        X_press[[h+1]][i, ] <- X_press[[h]][i, ] - X_pred
        
        # deflation matrices
        ui <- model$loadings$X[,1]
        vi <- model$loadings$Y[,1]
        res.deflat <- step2.spls(X=X_train_scaled,Y=Y_train_scaled,ui,vi,mode="regression")
        X_train_scaled = res.deflat$X.h
        Y_train_scaled = res.deflat$Y.h

      }
      
    }
    
    for(h in 1:ncomp.max){
      
      # deflation matrices
      u <- object$loadings$X[,h]
      v <- object$loadings$Y[,h]
      res.deflat <- step2.spls(X=X,Y=Y,u,v,mode="regression")
      X = res.deflat$X.h
      Y = res.deflat$Y.h
      
      for(j in 1:p){
        RSSj[h,j] <- sum((X[,j])^2)
        PRESSj[h,j] <- sum((X_press[[h+1]][,j])^2)
      }
      colnames(RSSj) <- paste0("X",1:p)
      colnames(PRESSj) <- paste0("X",1:p)
      
      # RSS computing
      RSS[h] <- sum(RSSj[h,]) 
      
      # PRESSh computing
      PRESS[h] <- sum(PRESSj[h,]) 
      
      # Q2
      q2[h] <- 1-PRESS[h]/RSS[max(h-1,1)]
      
      
    }# end h loop
    
    # first value correction 
    q2[1] <- 1-PRESS[1]/RSS0
    
  }else{
    stop("The mode must be either regression or canonical.")
  }# end if loop
  
  lim <- 0.0975

  h <- 1
  while(!is.na(q2[h]) && q2[h] > lim){
    h <- h + 1
  }
  h.best <- max(h-1,1)
  
  # Plot
  if(plot){
    plot(q2, type = "b", col = "blue", pch = 16,
       main = "Model Q² performance",
       xlab = "Number of components", ylab = "Q²")
    abline(h = lim, col = "red", lty = 2)
  }
  
  suggestion <- paste("best number of components : H =",h.best)
  
  q2.pls.results <- list(suggestion = suggestion, h.best = h.best, q2 = q2, PRESS = PRESS, RSS = RSS, PRESSj = PRESSj, RSSj = RSSj)
  return(q2.pls.results)
  
}


