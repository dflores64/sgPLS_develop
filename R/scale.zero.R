scale.zero <- function(X, center = TRUE, scale = TRUE){
  
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  X.s <- scale(X, center = TRUE, scale = TRUE)
  
  for(j in 1:p){
    if(unique(X.s[,j])[1] == "NaN"){X.s[,j] <- X[,j]}
  }
  
  return(X.s)
}
