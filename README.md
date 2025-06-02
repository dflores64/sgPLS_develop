# q2_function

# prediction function ----

predict.PLS <- predict.sPLS <- predict.gPLS <- predict.sgPLS <- function(object, newdata,  ...)
{
  
  #-- validation des arguments --#
  if (missing(newdata))
    stop("No new data available.")
  
  X = object$X
  Y = object$Y
  q = ncol(Y)
  p = ncol(X)
  
  if (length(dim(newdata)) == 2) {
    if (ncol(newdata) != p)
      stop("'newdata' must be a numeric matrix with ncol = ", p,
           " or a vector of length = ", p, ".")
  }
  
  if (length(dim(newdata)) == 0) {
    if (length(newdata) != p)
      stop("'newdata' must be a numeric matrix with ncol = ", p,
           " or a vector of length = ", p, ".")
    dim(newdata) = c(1, p)
  }
  
  #-- initialisation des matrices --#
  ncomp = object$ncomp
  a = object$loadings$X
  b = object$loadings$Y
  c = object$mat.c
  
  means.X = attr(X, "scaled:center")
  means.Y = attr(Y, "scaled:center")
  sigma.X = attr(X, "scaled:scale")
  sigma.Y = attr(Y, "scaled:scale")
  
  newdata = as.matrix(newdata)
  ones = matrix(rep(1, nrow(newdata)), ncol = 1)
  ##- coeff de regression
  B.hat = array(0, dim = c(p, q, ncomp))
  ##- prediction
  Y.hat = array(0, dim = c(nrow(newdata), q, ncomp))
  Y.hat2 = array(0, dim = c(nrow(newdata), q, ncomp))
  ##- variates
  t.pred = array(0, dim = c(nrow(newdata), ncomp))
  
  variates.X = object$variates$X
  betay = list()
  
  #-- prediction --#
  for(h in 1:ncomp){
    
    dd= coefficients(lm(Y~variates.X[,1:h,drop=FALSE])) #regression of Y on variates.global.X => =loadings.global.Y at a scale factor
    if(q==1){betay[[h]]=(dd[-1])}
    if(q>=2){betay[[h]]=(dd[-1,])}
    
    W = a[, 1:h,drop=FALSE] %*% solve(t(c[, 1:h,drop=FALSE]) %*% a[, 1:h,drop=FALSE])
    B = W %*% drop(betay[[h]])
    
    Y.temp=scale(newdata,center=means.X,scale=sigma.X) %*% as.matrix(B) #so far: gives a prediction of Y centered and scaled
    Y.temp2=scale(Y.temp,center=FALSE,scale=1/sigma.Y) #so far: gives a prediction of Y centered, with the right scaling
    Y.temp3=scale(Y.temp2,center=-means.Y,scale=FALSE) #so far: gives a prediction of Y with the right centering and scaling
    
    Y.hat[, , h] = Y.temp3 # we add the variance and the mean of Y used in object to predict
    t.pred[, h] = scale(newdata, center = means.X, scale = sigma.X) %*% W[, h]
    B.hat[, , h] = B
  }  #end h
  
  #-- valeurs sortantes --#
  rownames(t.pred) = rownames(newdata)
  colnames(t.pred) = paste("dim", c(1:ncomp), sep = " ")
  rownames(Y.hat) = rownames(newdata)
  colnames(Y.hat) = colnames(Y)
  
  return(invisible(list(predict = Y.hat, variates = t.pred, B.hat = B.hat,betay=betay)))
}



# function creating dataset ------------------------------------

D1 <- function(n, min = -3, max = 2) {
  x <- runif(n, min = min, max = max)
  f <- function(x) -0.2*x^5 + x^3 - x^2 + 3*x + 1
  eps <- rnorm(n, sd = 1)
  y <- f(x) + eps
  data.frame(x = x, y = y)
}


D2 <- function(q, n = 50,p = 45){
  x <- matrix(data = rnorm(n*p),n,p)
  y <- 1 + rowSums(x[,1:q]) + rnorm(n)
  d <- data.frame(x,y)
  list(x,y)
}

D3 <- function(q, n = 50,p = 45){
  x <- matrix(data = rnorm(n*p),n,p)
  y1 <- 1 + rowSums(x[,1:q]) + rnorm(n)
  y2 <- 1 + 5*x[,1] + -10*x[,2] + 2*x[,q] + rnorm(n)
  y3 <- 1 + -2*x[,1] + 8*x[,2] + 7*x[,q] + rnorm(n)
  y4 <- 1 + 4*x[,2] - 12*x[,q] + rnorm(n)
  y <- cbind(y1,y2,y3,y4)
  d <- data.frame(x,y)
  list(x,y)
}


jeu <- D3(q=3,n=30,p=10)
X <- jeu[[1]]
Y <- jeu[[2]]

jeu.test <- D3(q=3,n=10,p=10)
X.test <- jeu[[1]]
Y.test <- jeu[[2]]

# necessaries functions for the PLS function -----------------------------

# Fonction norme / norm function
normv <- function(x) sqrt(sum(x**2))

# Fonction de déflation / deflation function 
step2.spls <- function(X,Y,u.tild.new,v.tild.new,mode){
  ### Step d
  xi.h <- X%*% matrix(u.tild.new,ncol=1)/((normv(u.tild.new))**2)
  w.h  <- Y%*% matrix(v.tild.new,ncol=1)/((normv(v.tild.new))**2)
  
  ### Step e
  c.h <- t(X)%*%matrix(xi.h,ncol=1)/((normv(xi.h))**2)
  
  d.rh <- t(Y)%*%matrix(xi.h,ncol=1)/(sum(xi.h*xi.h))
  
  d.h <- t(Y)%*%matrix(w.h,ncol=1)/(sum(w.h*w.h))
  
  ###Step f and g
  X.h <- X - xi.h%*%t(c.h)
  if (mode=="regression") Y.h <- Y - xi.h%*%t(d.rh) else Y.h <- Y - w.h%*%t(d.h)
  
  res <- list(X.h=X.h,Y.h=Y.h,c=c.h,d=d.rh,e=d.h)
  return(res)
}


# PLS function ---------------------------------------------------

PLS <- function(X,Y,ncomp,mode){
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
  
  result <- list(X = X.s, Y = Y.s, ncomp = ncomp, loadings = list(X = load.u, Y = load.v),variates = list(X = mat.t, Y = mat.u), mat.c = mat.c, mat.d = mat.d, mat.e = mat.e)
  return(result)

}

# TEST ---------------------------------------------------------


result.pls <- PLS(X = X, Y = Y, ncomp = 1, mode = "regression")
Y.chap <- predict.PLS(result.pls, newdata = X.test)

result.pls2 <- PLS(X = X, Y = Y, ncomp = 2, mode = "regression")
Y.chap2 <- predict.PLS(result.pls2, newdata = X.test)

result.pls3 <- PLS(X = X, Y = Y, ncomp = 3, mode = "regression")
Y.chap3 <- predict.PLS(result.pls3, newdata = X.test)



head(Y)
head(Y.chap$predict)
head(Y.chap2$predict)
head(Y.chap3$predict)

# Q2 function -------------------------------------------

# This function is aimed to select the best number of components

q2.pls <- function(X,Y, mode = "regression"){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  p <- ncol(X)
  q <- ncol(Y)
  ncomp.max <- p
  if(nrow(X)!=nrow(Y)){ stop("X and Y have not the same number of rows.")}
  n <- nrow(X)
  X0 <- X
  Y0 <- Y
  q2 <- numeric(ncomp.max)
  
  # RSS0 computing
  Y.mean <- t(matrix(colMeans(Y), nrow = q, ncol = n)) #mean for each column
  RSS0 <- sum(colSums((Y0-Y.mean)^2))
  model.all <- PLS(X = X, Y = Y, ncomp = ncomp.max, mode = "regression")
  

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
    #Y <- Y - Y.hat.rss
    
    # Q2
    q2[h] <- 1-PRESS/RSS
    
    print(head(Y.hat.rss))
    
  }# end h loop
  
  
  q2[1] <- 1-PRESS1/RSS0 # correction of the first value in Q2 vector
  lim <- 0.0975
  h.best <- max(which(q2 > lim))
  
  # Plot
  plot(q2, type = "b", col = "blue", pch = 16,
       main = "Q² performance according to the number of components",
       xlab = "Number of components", ylab = "Q²")
  abline(h = 0.0975, col = "red", lty = 2)
  
  suggestion <- paste("best number of components : H =",h.best)
  
  q2.pls.results <- list(q2 = q2, suggestion = suggestion, h.best = h.best)
  return(q2.pls.results)
  
}

# test q2
q2.pls(X,Y)
