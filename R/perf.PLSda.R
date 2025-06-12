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
