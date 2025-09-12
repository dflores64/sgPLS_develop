# msep.PLS function ------

msep.PLS <- function(object, ncomp = object$ncomp, K=nrow(object$X), plot = TRUE){
  
  X <- object$X
  Y <- object$Y
  n <- nrow(X)
  p <- ncol(X)
  
  # conditions check-up
  if(!("pls" %in% class(object)) && class(object) != "mixo_pls"){ stop("object class must either contain pls class or be mixo_pls class."); print(class(object))}
  
  if(ncomp > object$ncomp || ncomp <= 0){ stop(paste("ncomp.max must be a value between 0 and",object$ncomp,"which is the total number of components computed in the object model."))}
  
  if(K < 2 || K > n){ stop(paste("K must be a value between 2 and",n))}
  
  # MSEP TRAIN -----
  
  Y.pred <- predict(object, newdata = X, ncomp = ncomp)$predict
  msep.train <- numeric(ncomp)
  for(h in 1:ncomp){msep.train[h] <- sum(colSums((Y - Y.pred[,,h])**2))}
  
  # MSEP CV ------
  
  ## prediction analysis
  err <- matrix(NA, nrow = K, ncol = ncomp)
  
  b <- floor(n/K) # block size
  ind <- sample(n)

  for(k in seq_len(K)){
    
    # bloc definition
    ind.beg <- (k-1) * b + 1 # block k beginning
    ind.end <- k * b # block k end
    ind.test <- ind[ind.beg:ind.end]
    nk <- length(ind.test)
    X.train <- X[-ind.test,]
    Y.train <- Y[-ind.test,]
    X.test <- X[ind.test,]
    Y.test <- Y[ind.test,]
    modele <- PLS(X = X.train, Y = Y.train, ncomp = ncomp, mode = "regression")  
    
    for(h in 1:ncomp){  
      
      # predictions
      pred <- predict.PLS(modele, newdata = X.test)$predict[,,h]
      err[k,h] <- sum(colSums(as.matrix((Y.test - pred)^2)))
      
    }
  }
  msep.cv <- colSums(err)/b/K
  
  h.best <- min(which.min(msep.cv))
  if(plot){
    plot(msep.cv, col="blue", pch = 16, type = "b", main = "Model MSEP", xlab = "number of components", ylab = "MSEP", axes = FALSE)
    axis(1, at = 1:ncomp)
    axis(2, labels = TRUE)
  }
  
  return(setNames(list(msep.train, msep.cv, h.best),c("MSEP.train","MSEP.cv","h.best")))
}



