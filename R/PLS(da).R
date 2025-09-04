# PLS function ----
PLS <- function(X,Y,ncomp,mode = "regression"){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  q <- ncol(Y)
  p <- ncol(X)
  n <- nrow(X)
  
  X.names = dimnames(X)[[2]]
  if (is.null(X.names)) 
    X.names = paste("X", 1:p, sep = "")
  if (dim(Y)[2] == 1) 
    Y.names = "Y"
  else {
    Y.names = dimnames(Y)[[2]]
    if (is.null(Y.names)) 
      Y.names = paste("Y", 1:q, sep = "")
  }
  ind.names = dimnames(X)[[1]]
  if (is.null(ind.names)) {
    ind.names = dimnames(Y)[[1]]
    rownames(X) = ind.names
  }
  if (is.null(ind.names)) {
    ind.names = 1:n
    rownames(X) = rownames(Y) = ind.names
  }
  
  
  X.s <- scale(X,scale=TRUE)
  Y.s <- scale(Y,scale=TRUE)
  
  
  mat.c <- matrix(nrow = p, ncol = ncomp)
  mat.d <- matrix(nrow = q, ncol = ncomp)
  mat.e <- matrix(nrow = q, ncol = ncomp)
  mat.t <- matrix(nrow = n, ncol = ncomp)
  mat.u <- matrix(nrow = n, ncol = ncomp)
  
  svd <- svd(t(X.s)%*%Y.s)
  load.u <- svd$u[,1]
  load.v <- svd$v[,1]
  res.deflat <- step2.spls(X=X.s,Y=Y.s,load.u,load.v,mode=mode)
  mat.c[,1] <- res.deflat$c
  
  if (mode=="regression") mat.d[,1] <- res.deflat$d else mat.e[,1] <- res.deflat$e
  
  mat.t[, 1] <- X.s%*%load.u
  mat.u[, 1] <- Y.s%*%load.v
  
  if(ncomp>1) {
    
    for (h in 2:ncomp) {
      X.h <- res.deflat$X.h
      Y.h <- res.deflat$Y.h
      svd <- svd(t(X.h)%*%Y.h)
      load.u.new <- svd$u[,1]
      load.v.new <- svd$v[,1]
      load.u <- cbind(load.u,load.u.new)
      load.v <- cbind(load.v,load.v.new)
      mat.t[, h] <- X.h%*%load.u.new
      mat.u[, h] <- Y.h%*%load.v.new   
      res.deflat <- step2.spls(X=res.deflat$X.h,Y=res.deflat$Y.h,load.u.new,load.v.new,mode=mode)
      mat.c[,h] <- res.deflat$c
      if (mode=="regression") mat.d[,h] <- res.deflat$d else mat.e[,h] <- res.deflat$e
    }
  }else{
    load.u <- matrix(load.u,ncol=1)
    load.v <- matrix(load.v,ncol=1)
  }
  colnames(load.u) <- NULL
  colnames(load.v) <- NULL  
  cl = match.call()
  result <- list(X=X.s,Y=Y.s,ncomp=ncomp,loadings=list(X = load.u, Y = load.v),variates=list(X = mat.t, Y = mat.u),mat.c=mat.c,mat.d=mat.d,mat.e=mat.e,mode=mode)
  class(result) = c("sPLS", "spls","pls")
  return(invisible(result))
  
}

# PLSda function ------------------------------------------
PLSda <- function(X,Y,ncomp = 2){
  
  # Testing the input Y
  if (is.null(dim(Y)))
  {
    Y = as.factor(Y)	
    ind.mat = unmap(as.numeric(Y))					
  }else {
    stop("'Y' should be a factor or a class vector.")						
  }		
  
  result = PLS(X, ind.mat, ncomp = ncomp, mode = "regression")
  
  cl = match.call()
  cl[[1]] = as.name('PLSda')
  result$call = cl
  
  result$ind.mat = ind.mat
  result$names$Y = levels(Y)
  
  class(result) = c("sPLSda","splsda","plsda")
  return(invisible(result))	
}

