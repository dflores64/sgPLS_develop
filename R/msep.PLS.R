# msep.PLS function ------




perf.PLSda <- function(object,
                        method.predict = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
                        validation = c("Mfold", "loo"), 
                        folds = 10, progressBar = TRUE){
  
  ncomp <- object$ncomp
  X <- object$X
  Y <- map(object$Y)
  n <- nrow(X)
  p <- ncol(X)
  method <- method.predict
  
  # conditions check-up
  if(!("plsda" %in% class(object)) && !("mixo_plsda" %in% class(object))){ stop("object class must either contain plsda class or contain mixo_plsda class."); print(class(object))}
  
  #if(ncomp > object$ncomp || ncomp <= 0){ stop(paste("ncomp.max must be a value between 0 and",object$ncomp,"which is the total number of components computed in the object model."))}
  
  if(method[1] == "max.dist"||method[2] == "max.dist"){dist=1}else if(method[1] == "centroids.dist"){dist=2}else{dist=3}
  
  if(validation[1] == "Mfold"){
    if(folds < 2 || folds > n){ stop(paste("folds must be a value between 2 and",n))}
    K = folds
  }else{K = n}
  
  if (progressBar == TRUE) pb <- txtProgressBar(style = 3)
  setTxtProgressBar(pb,1)
  cat('\n')
  
  b <- floor(n/K) # block size
  ind <- 1:n
  
  # prediction analysis
  err <- matrix(NA, nrow = K, ncol = ncomp)
  
  for(k in seq_len(K)){
    
    # bloc definition
    ind.beg <- (k-1) * b + 1 # block k beginning
    ind.end <- k * b # block k end
    ind.test <- ind[ind.beg:ind.end]
    X.train <- X[-ind.test,]
    Y.train <- Y[-ind.test]
    X.test <- X[ind.test,]
    Y.test <- Y[ind.test]
    modele <- PLSda(X = X.train,Y = Y.train, ncomp = ncomp)
    
    for(h in 1:ncomp){
      # model created
      pred <- predict.PLSda(modele, newdata = X.test, methode = methode)$class[[dist]][,h]
      equal <- Y.test == pred
      err[k,h] <- sum(1-equal)
      
    }
  }
  
  err.moy <- colSums(err)/b/K
  
  h.best <- min(which.min(err.moy))
  plot(err.moy, col="blue", pch = 16, type = "b", main = "Error rate of the model", xlab = "number of components", ylab = "Error")
  abline(h = (1:9)/10, lty = 3, col = "grey")
  
  return(setNames(list(err.moy,h.best),c("error.rate","h.best")))
}




