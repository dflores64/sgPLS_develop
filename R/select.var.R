select.var <- function(object, comp = 1){
  
  if(!("spls" %in% class(object)) && class(object) != "mixo_spls"){ stop("object class must either contain spls class or be mixo_spls class."); print(class(object))}
  
  if(comp > object$ncomp || comp <= 0){ stop(paste("comp must be a value between 0 and",object$ncomp,"which is the total number of components computed in the object model."))}
  
  ncomp <- object$ncomp
  keepX <- object$keepX
  keepY <- object$keepY
  loads.X <- object$loadings$X
  loads.Y <- object$loadings$Y
  colnames(loads.X) <- paste0("comp",1:ncomp)
  colnames(loads.Y) <- paste0("comp",1:ncomp)
  out <- list()
  
  # X loads --------------
  nb.var.select <- keepX[comp]
  loads.comp <- loads.X[,comp]
  loads.comp.abs <- abs(loads.X[,comp])
  loads.comp.abs.sort <- sort(loads.comp.abs,decreasing = TRUE)[1:nb.var.select]
  loads.comp.sort <- sort(loads.comp,decreasing = TRUE)
  loads.name <- rownames(as.matrix(loads.comp.abs.sort))
  
  loads.value <- numeric(nb.var.select)
  for(i in 1:nb.var.select){
    if(-loads.comp.abs.sort[i] %in% loads.comp.sort){ 
      # means that this value is an absolute value and its true value is the opposite value
      loads.value[i] <- -loads.comp.abs.sort[i]
    }else{loads.value[i] <- loads.comp.abs.sort[i]}
  }
  
  loads.value <- data.frame(loads.value)
  rownames(loads.value) <- loads.name
  
  out$X$name <- loads.name
  out$X$value <- data.frame(loads.value)
  
  # Y loads --------------
  nb.var.select <- keepY[comp]
  loads.comp <- loads.Y[,comp]
  loads.comp.abs <- abs(loads.Y[,comp])
  loads.comp.abs.sort <- sort(loads.comp.abs,decreasing = TRUE)[1:nb.var.select]
  loads.comp.sort <- sort(loads.comp,decreasing = TRUE)
  loads.name <- rownames(as.matrix(loads.comp.abs.sort))
  
  loads.value <- numeric(nb.var.select)
  for(i in 1:nb.var.select){
    if(-loads.comp.abs.sort[i] %in% loads.comp.sort){ 
      # means that this value is an absolute value and its true value is the opposite value
      loads.value[i] <- -loads.comp.abs.sort[i]
    }else{loads.value[i] <- loads.comp.abs.sort[i]}
  }
  
  loads.value <- data.frame(loads.value)
  rownames(loads.value) <- loads.name
  
  out$Y$name <- loads.name
  out$Y$value <- data.frame(loads.value)
  
  out$comp <- comp
  
  return(out)
  
}
