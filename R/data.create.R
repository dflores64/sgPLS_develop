data.create <- function(n = 40, p = 10, q = 1, list = TRUE){
  
  X <- matrix(data = rnorm(n*p),n,p)
  B <- matrix(data = runif(q*p,-1,1), nrow = p, ncol = q)
  E <- matrix(data = rnorm(n*q, sd = 0.1),n,q)
  Y <- X%*%B + E
  
  if(p > 1){colnames(X) <- paste0(rep("X",p),1:p)}
  if(q > 1){colnames(Y) <- paste0(rep("Y",q),1:q)}
  D <- data.frame(X,Y)
  
  if(list){return(list(B = B,D = D,X = X,Y = Y))}else{return(D)}
}

data.cl.create <- function(n = 40, p = 10, classes = 2, list = TRUE){
  
  X <- matrix(data = rnorm(n*p),n,p)
  B <- matrix(data = runif(p,-1,1), nrow = p, ncol = 1)
  E <- matrix(data = rnorm(n, sd = 0.1),n,1)
  Y0 <- X%*%B + E
  Y <- numeric(n)
  
  for(i in 1:n){
    for(cl in 1:classes){
      if(Y0[i] >= quantile(Y0, prob = (cl-1)/classes) && Y0[i] <= quantile(Y0, prob = cl/classes)){Y[i] <- cl}
    }
  }
  
  if(p > 1){colnames(X) <- paste0(rep("X",p),1:p)}
  D <- data.frame(X,Y)
  
  if(list){return(list(B = B,D = D,X = X,Y = Y))}else{return(D)}
  
}

data.spls.create <- function(n = 100, p100 = 4, q100 = 5, list = TRUE){
  
  # noise standard deviation
  p <- 100*p100
  q <- 100*q100
  sigma.gamma <- 1
  sigma.e <- 1.5
  
  theta.x1 <- c(rep(1, 4*p100), rep(0, p100), rep(-1, 4*p100), rep(0, p100), rep(1.5,2*p100), 
                rep(0, p100), rep(-1.5, 2*p100), rep(0, 85*p100))
  theta.x2 <- c(rep(0, 84*p100), rep(1, 4*p100), rep(0, p100), rep(-1, 4*p100), rep(0, p100),
                rep(1.5, 2*p100), rep(0, p100), rep(-1.5, 2*p100), rep(0, p100))
  
  theta.y1 <- c(rep(1, 4*q100), rep(0, q100), rep(-1, 4*q100), rep(0, q100), rep(1.5, 2*q100),
                rep(0, q100), rep(-1.5, 2*q100), rep(0, 85*q100))
  theta.y2 <- c(rep(0, 84*q100), rep(1, 4*q100), rep(0, q100), rep(-1, 4*q100), rep(0, q100),
                rep(1.5, 2*q100), rep(0, q100), rep(-1.5, 2*q100), rep(0, q100))
  
  # covariance matrices
  Sigmax <- matrix(0, nrow = p, ncol = p)
  diag(Sigmax) <- sigma.e ^ 2
  Sigmay <- matrix(0,nrow = q, ncol = q)
  diag(Sigmay) <- sigma.e ^ 2
  
  set.seed(125)
  
  gam1 <- rnorm(n)
  gam2 <- rnorm(n)
  
  GAM <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE)
  Thetax <- matrix(c(theta.x1, theta.x2),
                   nrow = 2, byrow = TRUE)
  Thetay <- matrix(c(theta.y1, theta.y2),
                   nrow = 2, byrow = TRUE)
  E1 <- rmvnorm(n, mean = rep(0, p), sigma = Sigmax, method = "svd")
  E2 <- rmvnorm(n, mean = rep(0, q), sigma = Sigmay, method = "svd")
  
  X <- GAM %*% Thetax + E1                                                
  Y <- GAM %*% Thetay + E2
  
  ind.block.x <- seq(20, p-20, 20)
  ind.block.y <- seq(20, q-20, 20)
  
  colnames(X) <- paste0(rep("X",p),1:p)
  colnames(Y) <- paste0(rep("Y",q),1:q)
  D <- data.frame(X,Y)
  
  if(list){return(list(D = D,X = X,Y = Y, ind.block.x = ind.block.x, ind.block.y =ind.block.y))}else{return(D)}
               
}


