plot.indiv <- function(object, compX = c(1,2), compY = NULL){
  
  class.obj <- class(object)
  X <- object$X
  Y <- object$Y
  n <- nrow(X)
  ncomp <- object$ncomp
  mat.t <- as.data.frame(object$variates$X)

  compX1 <- compX[1]
  compX2 <- compX[2]
  
  # check arguments
  if(is.na(compX2)){stop("'comp' must be a numeric vector of length 2.")}
  if(compX1 < 1 || compX2 < 1){stop("compX and compY must be strictly positive.")}
  if(compX1 > ncomp || compX2 > ncomp){stop("compX elements mustn't exceed ncomp.")}
  
  t1 <- mat.t[,compX1]
  t2 <- mat.t[,compX2]
  data <- data.frame(t1,t2)
  t1.name <- paste("X-variate",compX1)
  t2.name <- paste("X-variate",compX2)
  colnames(data) <- c(t1.name,t2.name)
  
  # PLS/PLSda condition 
  if(class.obj[3] == "plsda"){
    Y <- map(Y)
    Classes <- data$class <- as.factor(Y)
    graphX <- ggplot(mat.t, aes(t1,t2, colour = Classes)) + geom_point() + 
      stat_ellipse(level = 0.95, type = "norm") +
      labs(x = paste("X-variate",compX1), y = paste("X-variate ",compX2)) + 
      geom_text(aes(label = 1:n), vjust = -1) +
      ggtitle("Individuals projected on X latent variables") +
      theme(plot.title = element_text(hjust = 0.5))
    
  }else{
    graphX <- ggplot(mat.t, aes(t1,t2)) + geom_point(colour = "blue") + 
      labs(x = paste("X-variate",compX1), y = paste("X-variate ",compX2)) + 
      geom_text(aes(label = 1:n), vjust = -1) +
      ggtitle("Individuals projected on X latent variables") +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  # case of compY input and PLS
  if(!is.null(compY) && class.obj[3] == "pls"){
    
    compY1 <- compY[1]
    compY2 <- compY[2]
    
    #check arguments
    if(is.na(compY2)){stop("'comp' must be a numeric vector of length 2.")}
    if(compY1 < 1 || compY2 < 1){stop("compX and compY must be strictly positive.")}
    if(compY1 > ncomp || compY2 > ncomp){stop("compY elements mustn't exceed ncomp.")}

    mat.s <- as.data.frame(object$variates$Y)
    s1 <- mat.t[,compY1]
    s2 <- mat.t[,compY2]
    
    data <- data.frame(t1,t2,s1,s2)
    s1.name <- paste("Y-variate",compY1)
    s2.name <- paste("Y-variate",compY2)
    colnames(data) <- c(t1.name,t2.name,s1.name,s2.name)
    
    graphY <- ggplot(mat.t, aes(s1,s2)) + geom_point(colour = "blue") + 
      labs(x = paste("Y-variate",compY1), y = paste("Y-variate ",compY2)) + 
      geom_text(aes(label = 1:n), vjust = -1) +
      ggtitle("Individuals projected on Y latent variables") +
      theme(plot.title = element_text(hjust = 0.5))
    
    return(list(data = data, graphX = graphX, graphY = graphY))
  }else{return(list(data = data, graphX = graphX))}
  
}
