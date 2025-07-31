data.create <- function(n = 40, p = 3, q = 1, list = FALSE){  
  X <- matrix(data = rnorm(n*p),n,p)
  B <- matrix(data = runif(q*p,-1,1), nrow = p, ncol = q)
  E <- matrix(data = rnorm(n*q, sd = 0.05),n,q)
  Y <- X%*%B + E
  
  if(p > 1){colnames(X) <- paste0(rep("X",p),1:p)}
  if(q > 1){colnames(Y) <- paste0(rep("Y",q),1:q)}
  D <- data.frame(X,Y)
  
  if(list){return(list(D = D,X = X,Y = Y))}else{return(D)}
}

data.cl.create <- function(n = 40, p = 3, classes = 2, list = FALSE){
  X <- matrix(nrow = n, ncol = p)
  Y <- numeric(n)
  for(i in seq_len(n)){
    Y[i] <- sample(x = seq_len(classes), size = 1)
    X[i,] <- rnorm(p) + Y[i] + rnorm(p, sd = 0.05)
  }
  
  if(p > 1){colnames(X) <- paste0(rep("X",p),1:p)}
  D <- data.frame(X,Y)
  if(list){return(list(D = D,X = X,Y = Y))}else{return(D)}
}
