data.create <- function(n = 40, p = 10, q = 1){
  X <- matrix(data = runif(n*p),n,p)
  U <- matrix(data = runif(q*p,-10,10), nrow = p, ncol = q)
  E <- matrix(data = rnorm(n*q),n,q)
  Y <- X%*%U + E
  D <- data.frame(X,Y)
  colnames(X) <- paste0(rep("X",p),1:p)
  return(list(D = D,X = X,Y = Y))
}

data.cl.create <- function(n = 40, p = 10, classes = 2){
  X <- matrix(nrow = n, ncol = p)
  Y <- numeric(n)
  for(i in seq_len(n)){
    Y[i] <- sample(x = seq_len(classes), size = 1)
    X[i,] <- runif(p) + Y[i]
  }
  D <- data.frame(X,Y)
  colnames(X) <- paste0(rep("X",p),1:p)
  return(list(D = D,X = X,Y = Y))
}
