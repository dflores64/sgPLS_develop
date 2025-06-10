
# function creating dataset ---------------------------------

create.data <- function(n = 40, p = 10, q = 1){
  X <- matrix(data = runif(n*p),n,p)
  U <- matrix(data = runif(q*p,-10,10), nrow = p, ncol = q)
  E <- matrix(data = rnorm(n*p),n,q)
  Y <- X%*%U + E
  D <- data.frame(X,Y)
  return(list(X = X,Y = Y))
}

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

# 


# PLS function -----------------------------------------
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
  
  
  mat.c <-matrix(nrow = p, ncol = ncomp)
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
  
  result <- list(X = X.s, Y = Y.s, ncomp = ncomp, loadings = list(X = load.u, Y = load.v),variates = list(X = mat.t, Y = mat.u), mat.c = mat.c)
  return(result)

}

# prediction function -------------------------------------------

predict.PLS <- function(object, newdata,  ...)
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


# TEST ----

# data simulation

train <- 1:40
test <- 41:50

jeu <- create.data(n=50)

X <- jeu$X[train,]
Y <- jeu$Y[train,]
X.test <- jeu$X[test,]
Y.test <- jeu$Y[test,]

result.pls <- PLS(X = X, Y = Y, ncomp = 10, mode = "regression")
Y.chap <- predict.PLS(result.pls, newdata = X.test)

head(Y.test)
head(Y.chap$predict)
sum((Y.test - Y.chap$predict[,,10])^2)

## data yarn

library(pls)
data(yarn)

train <- sample(x = 1:28, size = 20)
X <- yarn$NIR[train,]
Y <- yarn$density[train]
X.test <- yarn$NIR[-train,]
Y.test <- yarn$density[-train]

result.pls3 <- PLS(X = X, Y = Y, ncomp = 19, mode = "regression")
Y.chap3 <- predict.PLS(result.pls3, newdata = X.test)
head(Y.chap3$predict)
head(Y.test)

# PLSda function ------------------------------------------

PLSda <- function(X,Y,ncomp = 2, 
                  keepX = rep(ncol(X), ncomp),
                  max.iter = 500,		 
                  tol = 1e-06){
  return(sPLSda(X,Y,ncomp = ncomp, 
                keepX = keepX,
                max.iter = max.iter,		 
                tol = tol))
}

predict.PLSda <- 
  function(object, newdata, 
           method = c("all", "max.dist", "centroids.dist", "mahalanobis.dist")){
    #-- validation des arguments --#
    if (missing(newdata))
      stop("No new data available.")
    
    X = object$X
    Y = object$Y 
    Yprim = object$ind.mat   
    q = ncol(Yprim)          
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
    
    G = matrix(0, nrow = q, ncol = ncomp)
    cls = list()
    
    for (i in 1:q) {
      if(ncomp > 1) {
        
        G[i, ] = apply(object$variates$X[Yprim[, i] == 1, , drop = FALSE], 2, mean)
      }
      else {
        G[i, ] = mean(object$variates$X[Yprim[, i] == 1, ])
      }
    }	
    
    # ----    max distance -----------------
    
    if (any(method == "all") || any(method == "max.dist")) {
      
      function.pred = function(x){
        nr = nrow(x)
        tmp = vector("numeric", nr)
        for(j in 1:nr){
          tmp[j] = (which(x[j, ] == max(x[j, ]))[1])
        }
        return(tmp)
      }
      cls$max.dist = matrix(apply(Y.hat, 3, function.pred), ncol = ncomp)
      colnames(cls$max.dist) = paste(rep("comp", ncomp), 1:ncomp, sep = " ")
    }
    
    # ----    centroids distance -----------------
    
    if (any(method == "all") || any(method == "centroids.dist")) {
      
      cl = matrix(nrow = nrow(newdata), ncol = ncomp)
      
      centroids.fun = function(x, G, h) {
        q = nrow(G)
        x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
        if (h > 1) {
          d = apply((x - G[, 1:h])^2, 1, sum)
        }
        else {
          d = (x - G[, 1])^2
        }
        cl.id = which.min(d)
      }
      
      for (h in 1:ncomp) {
        cl.id = apply(matrix(t.pred[, 1:h], ncol = h), 1, centroids.fun, G = G, h = h)
        cl[, h] = cl.id		
      }
      colnames(cl) = paste(rep("comp", ncomp), 1:ncomp, sep = " ")
      cls$centroids.dist = cl
    }	
    
    # ----    mahalanobis distance -----------------
    
    if (any(method == "all") || any(method == "mahalanobis.dist")) {
      
      cl = matrix(nrow = nrow(newdata), ncol = ncomp)
      
      Sr.fun = function(x, G, Yprim, h) {
        q = nrow(G)
        Xe = Yprim %*% G[, 1:h]
        Xr = object$variates$X[, 1:h] - Xe
        Sr = t(Xr) %*% Xr / nrow(Y)
        Sr.inv = solve(Sr)
        x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
        if (h > 1) {
          mat = (x - G[, 1:h]) %*% Sr.inv %*% t(x - G[, 1:h])
          d = apply(mat^2, 1, sum)
        }
        else {
          d = drop(Sr.inv) * (x - G[, 1])^2
        }
        cl.id = which.min(d)
      }
      
      for (h in 1:ncomp) {
        cl.id = apply(matrix(t.pred[, 1:h], ncol = h), 1, Sr.fun, G = G, Yprim = Yprim, h = h)
        cl[, h] = cl.id		
      }
      colnames(cl) = paste(rep("comp", ncomp), 1:ncomp, sep = " ")
      cls$mahalanobis.dist = cl
    }
    
    #-- valeurs sortantes --#
    if (any(method == "all")) method = "all"
    rownames(t.pred) = rownames(newdata)
    colnames(t.pred) = paste("dim", c(1:ncomp), sep = " ")
    rownames(Y.hat) = rownames(newdata)
    colnames(Y.hat) = colnames(Y)
    colnames(G) = paste("dim", c(1:ncomp), sep = " ")
    
    return(invisible(list(predict = Y.hat, variates = t.pred, B.hat = B.hat, 
                          centroids = G, method = method, class = cls)))
  }


# TEST PLSda function --------------------------------------

## 2classes --------
library(sgPLS)
library(mvtnorm)
rAD2 <- function(n, prob, mu1, mu2, Sigma1, Sigma2) {
  y <- numeric(n)
  x <- matrix(NA, nrow = n, ncol = length(mu1))
  for (i in seq_len(n)) {
    y[i] <- rbinom(1, size = 1, prob = 1 - prob) + 1 
    x[i,] <- if (y[i] == 1) rmvnorm(1, mean = mu1, sigma = Sigma1) else
      rmvnorm(1, mean = mu2, sigma = Sigma2)
  }
  data.frame(x, y)
}

prob <- 0.6
mu1 <- c(1, 1)
mu2 <- c(-1, -1)
Sigma1 <- matrix(c(2, 0, 0, 2), nrow = 2, ncol = 2)
Sigma2 <- matrix(c(0.5, 0, 0, 0.5), nrow = 2, ncol = 2)

train_12 <- rAD2(n = 80, prob = prob, mu1 = mu1, mu2 = mu2,
                 Sigma1 = Sigma1, Sigma2 = Sigma2)
test_12 <- rAD2(n = 20, prob = prob, mu1 = mu1, mu2 = mu2,
                Sigma1 = Sigma1, Sigma2 = Sigma2)

modele <- sPLSda(X = train_12[,1:2],Y = train_12$y, ncomp = 2)
pred <- predict.PLSda(modele, newdata = test_12[,1:2])$class$max.dist
table(pred[,2],test_12$y)

## 3 classes ------------

rAD3 <- function(n,mu1, mu2, mu3, Sigma1, Sigma2, Sigma3) {
  y <- numeric(n)
  x <- matrix(NA, nrow = n, ncol = length(mu1))
  for (i in seq_len(n)) {
    y[i] <- sample(seq(3),1) # tirage aléatoire uniforme de la classe
    x[i,] <- if (y[i] == 1) rmvnorm(1, mean = mu1, sigma = Sigma1) else if (y[i] == 2)
      rmvnorm(1, mean = mu2, sigma = Sigma2) else
        rmvnorm(1, mean = mu3, sigma = Sigma3)
  }
  data.frame(x, y)
}

mu1 <- c(1, 1)
mu2 <- c(-1, -1)
mu3 <- c(0,0)
Sigma1 <- matrix(c(2, 0, 0, 2), nrow = 2, ncol = 2)
Sigma2 <- matrix(c(0.5, 0, 0, 0.5), nrow = 2, ncol = 2)
Sigma3 <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)

train_123 <- rAD3(n = 80, mu1 = mu1, mu2 = mu2, mu3 = mu3,
                  Sigma1 = Sigma1, Sigma2 = Sigma2, Sigma3 = Sigma3)
test_123 <- rAD3(n = 20, mu1 = mu1, mu2 = mu2, mu3 = mu3,
                 Sigma1 = Sigma1, Sigma2 = Sigma2, Sigma3 = Sigma3)

modele <- PLSda(X = train_123[,1:2],Y = train_123$y, ncomp = 6)
pred <- predict.PLSda(modele, newdata = test_123[,1:2])$class$max.dist


# PLSda performance ------------------------------------

## first function ---------------
perf.PLSda <- function(X,Y,X.test,Y.test){
  modele <- PLSda(X = X,Y = Y)
  pred <- predict.PLSda(modele, newdata = X.test)$class$max.dist
  
  erreur <- numeric(2)
  matconf <- list()
  matconf[[1]] <- table(pred[,1],Y.test)
  matconf[[2]] <- table(pred[,2],Y.test)
  erreur[1] <- 1 - sum(diag(matconf[[1]])) / sum(matconf[[1]])
  erreur[2] <- 1 - sum(diag(matconf[[2]])) / sum(matconf[[2]])
  h.best <- min(which.min(erreur))
  
  plot(erreur, col="blue", pch = 16, type = "b", main = "Error rate of the model", xlab = "number of components", ylab = "Error")
  return(setNames(list(pred, erreur, matconf,h.best),c("predictions","erreur","confusion","h.best")))
}

perf.PLSda(X = train_12[,1:2],Y = train_12$y, X.test = test_12[,1:2],Y.test = test_12$y)
perf.PLSda(X = train_123[,1:2],Y = train_123$y, X.test = test_123[,1:2],Y.test = test_123$y)

## second function ---------------------
perf.PLSda.bis <- function(X,Y,X.test,Y.test){
  modele <- PLSda(X = X,Y = Y)
  pred <- predict.PLSda(modele, newdata = X.test)$class$max.dist
  
  erreur <- numeric(2)
  matconf <- list()
  matconf[[1]] <- table(pred[,1],Y.test)
  matconf[[2]] <- table(pred[,2],Y.test)
  erreur[1] <- 1 - sum(diag(matconf[[1]])) / sum(matconf[[1]])
  erreur[2] <- 1 - sum(diag(matconf[[2]])) / sum(matconf[[2]])
  barplot(erreur, col  ="blue", pch = 16, width = c(0.2,0.2), main = "Error rate of the model", xlab = "number of components", ylab = "Error")
  return(setNames(list(pred, erreur, matconf),c("predictions","erreur","confusion")))
}

perf.PLSda.bis(X = train_12[,1:2],Y = train_12$y, X.test = test_12[,1:2],Y.test = test_12$y)
perf.PLSda.bis(X = train_123[,1:2],Y = train_123$y, X.test = test_123[,1:2],Y.test = test_123$y)

# Avec MixOmics -----------------------------------------------------------

library(mixOmics) # import the mixOmics library
data(srbct) # extract the small round bull cell tumour data
X <- srbct$gene # use the gene expression data as the X matrix
Y <- srbct$class # use the class data as the Y matrix

result.plsda.srbct <- plsda(X, Y) # run the method
plotIndiv(result.plsda.srbct) # plot the samples
plotVar(result.plsda.srbct) # plot the variables

splsda.result <- splsda(X, Y, keepX = c(50,30)) # run the method
plotIndiv(splsda.result) # plot the samples
plotVar(splsda.result) # plot the variables

# extract the variables used to construct the first latent component
selectVar(splsda.result, comp = 1)$name 
# depict weight assigned to each of these variables
plotLoadings(splsda.result, method = 'mean', contrib = 'max')  
