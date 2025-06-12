create.data <- function(n = 40, p = 10, q = 1){
  X <- matrix(data = runif(n*p),n,p)
  U <- matrix(data = runif(q*p,-10,10), nrow = p, ncol = q)
  E <- matrix(data = rnorm(n*p),n,q)
  Y <- X%*%U + E
  D <- data.frame(X,Y)
  return(list(X = X,Y = Y))
}
