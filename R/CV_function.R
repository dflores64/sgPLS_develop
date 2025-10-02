# CV function

CV <- function(X,Y,type = "loo",mode = "regression",ncomp){
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  #if(nrow(X) != nrow(Y)) stop("X and Y must have the same number of rows.")
  
  err <- 0

  if(mode == "regression"){
      
    for(i in 1:n){
        
      X.train <- X[-i,]
      Y.train <- Y[-i,]
      X.test <- X[i,]
      Y.test <- Y[i,]
        
      modele <- PLS(X = X.train,Y = Y.train, ncomp = ncomp, mode = "regression")
      pred <- predict.PLS(modele, newdata = X.test)$predict
      err <- err + sum((Y.test - pred)^2)
        
    }
      
  }else if (mode == "class"){
      
    for(i in 1:n){
        
      X.train <- X[-i,]
      Y.train <- Y[-i]
      X.test <- X[i,]
      Y.test <- Y[i]
        
      #Y.train <- as.factor(Y.train)
      modele <- PLSda(X = X.train,Y = Y.train, ncomp = ncomp)
      pred <- predict.PLSda(modele, newdata = X.test)$class[[1]][, ncomp]
      equal <- Y.test == pred
      err <- err + 1-equal
        
    }
  }else{stop("mode must be either 'regression' or 'classif'")}
  return(err/n)
}


# seconde version -----------

CV2 <- function(X,Y,K,mode = "regression",ncomp){
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  if(K > n) stop("K mustn't exceed n = nrow(X)")
  if(ncomp > p) stop("ncomp mustn't exceed p = ncol(X)")
  
  err <- 0
  b <- floor(n/K) # block size
  
  ind <- sample(n)
  
  ## Indices for beginning and end of each block
  bloc.ind <- matrix(NA, nrow = K, ncol = 2)
  for (k in seq_len(K)) {
    bloc.ind[k,1] <- (k-1) * b + 1 # block k beginning
    bloc.ind[k,2] <- k * b # block k end
  }
  
  if(mode == "regression"){
    
    for (k in seq_len(K)){
      
      ind.test <- ind[bloc.ind[k,1]:bloc.ind[k,2]]
      nk <- length(ind.test)
      X.train <- X[-ind.test,]
      Y.train <- Y[-ind.test,]
      X.test <- X[ind.test,]
      Y.test <- Y[ind.test,]
      
      modele <- PLS(X = X.train, Y = Y.train, ncomp = ncomp, mode = "regression")
      pred <- predict.PLS(modele, newdata = X.test)$predict[,,ncomp]
      
      err <- err + sum(colSums(as.matrix((Y.test - pred)^2)))
      
    }
    
  }else if (mode == "class"){
    
    for (k in seq_len(K)){
      
      ind.test <- ind[bloc.ind[k,1]:bloc.ind[k,2]]
      nk <- length(ind.test)
      X.train <- X[-ind.test,]
      Y.train <- Y[-ind.test]
      X.test <- X[ind.test,]
      Y.test <- Y[ind.test]
      
      modele <- PLSda(X = X.train,Y = Y.train, ncomp = ncomp)
      pred <- predict.PLSda(modele, newdata = X.test)$class[[1]][, ncomp]
      equal <- Y.test == pred
      err <- err + sum(1-equal)
      
    }
  }else{stop("mode must be either 'regression' or 'classif'")}
  return(err/K/b)
}

# 3e version (object) -----------

CV3 <- function(object,K,mode = "regression",ncomp){
  
  X <- object$X
  Y <- object$Y
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  if(K > n) stop("K mustn't exceed n = nrow(X)")
  if(ncomp > p) stop("ncomp mustn't exceed p = ncol(X)")
  
  err <- 0
  b <- floor(n/K) # block size
  
  ind <- sample(n)
  
  ## Indices for beginning and end of each block
  bloc.ind <- matrix(NA, nrow = K, ncol = 2)
  for (k in seq_len(K)) {
    bloc.ind[k,1] <- (k-1) * b + 1 # block k beginning
    bloc.ind[k,2] <- k * b # block k end
  }
  
  if(mode == "regression"){
    
    for (k in seq_len(K)){
      
      ind.test <- ind[bloc.ind[k,1]:bloc.ind[k,2]]
      nk <- length(ind.test)
      X.train <- X[-ind.test,]
      Y.train <- Y[-ind.test,]
      X.test <- X[ind.test,]
      Y.test <- Y[ind.test,]
      
      modele <- PLS(X = X.train, Y = Y.train, ncomp = ncomp, mode = "regression")
      pred <- predict.PLS(modele, newdata = X.test)$predict[,,ncomp]
      
      err <- err + sum(colSums(as.matrix((Y.test - pred)^2)))
      
    }
    
  }else if (mode == "class"){
    
    Y <- map(Y)
    
    for (k in seq_len(K)){
      
      ind.test <- ind[bloc.ind[k,1]:bloc.ind[k,2]]
      nk <- length(ind.test)
      X.train <- X[-ind.test,]
      Y.train <- Y[-ind.test]
      X.test <- X[ind.test,]
      Y.test <- Y[ind.test]
      
      modele <- PLSda(X = X.train,Y = Y.train, ncomp = ncomp)
      pred <- predict.PLSda(modele, newdata = X.test)$class[[1]][, ncomp]
      equal <- Y.test == pred
      err <- err + sum(1-equal)
      
    }
    
  }else{stop("mode must be either 'regression' or 'classif'")}
  return(err/K/b)
}

