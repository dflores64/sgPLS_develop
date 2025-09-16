# ---------------------------------------------------
# perf for PLS object ----
# ---------------------------------------------------

# perf.PLS function -----

perf.PLS <- function(object, criterion = c("all","MSEP","Q2"), validation = c("Mfold","loo"),
                     folds = 10, progressBar = TRUE, setseed = 1, plot = FALSE){
  
  X <- object$X
  Y <- object$Y
  c <- object$mat.c
  d <- object$mat.d
  ncomp <- object$ncomp
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  # conditions check-up
  if(!("pls" %in% class(object)) && class(object) != "mixo_pls"){ stop("object class must either contain pls class or be mixo_pls class."); print(class(object))}

  if(validation[1] == "Mfold"){
    if(folds < 2 || folds > n){ stop(paste("folds must be a value between 2 and",n))}
    K = folds
  }else{K = n}
  
  if (progressBar == TRUE) pb <- txtProgressBar(style = 3)
  setTxtProgressBar(pb,1)
  cat('\n')
  
  set.seed(setseed)
  ind <- sample(1:n,n,replace = FALSE)
  b <- floor(n/K) # block size
  res = list()
  
  # MSEP CRITERION
  if(any(criterion %in% c("all","MSEP"))){
    
    # prediction analysis
    err <- matrix(NA, nrow = K, ncol = ncomp)
    MSEPj <- matrix(0, nrow = ncomp, ncol = q)
    colnames(MSEPj) <- paste0("Y",1:q)
    rownames(MSEPj) <- paste0("Comp",1:ncomp)
    
    for(k in seq_len(K)){
      
      # bloc definition
      ind.beg <- (k-1) * b + 1 # block k beginning
      ind.end <- k * b # block k end
      ind.test <- ind[ind.beg:ind.end]
      nk <- length(ind.test)
      X.train <- X[-ind.test,]
      Y.train <- Y[-ind.test,]
      X.test <- X[ind.test,]
      Y.test <- Y[ind.test,]
      modele <- PLS(X = X.train, Y = Y.train, ncomp = ncomp, mode = "regression")  
      
      for(h in 1:ncomp){  
        
        # predictions
        pred <- predict.PLS(modele, newdata = X.test)$predict[,,h]
        MSEPj[h,] <- MSEPj[h,] + colMeans(as.matrix((Y.test - pred)^2))/K
        
      }
    }
    err.moy <- rowMeans(MSEPj)
    
    h.best.msep <- min(which.min(err.moy))
    if(plot){
      plot(err.moy, col="blue", pch = 16, type = "b", main = "MSEP of the model", xlab = "number of components", ylab = "MSEP", axes = FALSE)
      axis(1, at = 1:ncomp)
      axis(2, labels = TRUE)
    }
    
    res$MSEP <- err.moy
    res$MSEPj <- MSEPj
    res$h.best.msep <- h.best.msep

  }
  
  if(any(criterion %in% c("all","Q2"))){
    
    q2 <- numeric(ncomp)
    PRESS <- numeric(ncomp)
    RSS <- numeric(ncomp)
    
    RSSj <- matrix(nrow = ncomp, ncol = q)
    PRESSj <- matrix(nrow = ncomp, ncol = q)
    
    # RSS0 computing
    Y.mean <- t(matrix(colMeans(Y), nrow = q, ncol = n)) # mean for each column
    RSS0 <- sum(colSums((Y-Y.mean)^2))
    
    
    Y_test <- list()
    Y_test[[1]] <- Y
    
    for(h in 1:ncomp){
      Y_test[[h+1]] <- matrix(nrow = n, ncol = q)
    }
    
    for(k in seq_len(K)){
      
      # bloc definition
      ind.beg <- (k-1) * b + 1 # block k beginning
      ind.end <- k * b # block k end
      ind.test <- ind[ind.beg:ind.end]
      nk <- length(ind.test)
      X_train <- X[-ind.test,, drop = FALSE]
      Y_train <- Y[-ind.test,, drop = FALSE]
      
      #Y.test <- Y[ind.test,]
      
      
      # Center and scale training data
      X_train_scaled <- scale(X_train, center = TRUE, scale = TRUE)
      Y_train_scaled <- scale(Y_train, center = TRUE, scale = TRUE)
      
      X_means <- apply(X_train,2,mean)
      X_sds <- apply(X_train,2,sd)
      Y_means <- apply(Y_train,2,mean)
      Y_sds <- apply(Y_train,2,sd)
      
      
      X_test <- X[ind.test,]
      X_test_scaled <- (X_test - X_means) / X_sds
      
      for(h in 1:ncomp){
        
        # training on the dataset without the ith individual
        model <- PLS(X = X_train_scaled, Y = Y_train_scaled, ncomp = 1, mode = "regression")
        
        # predictions on the ith individual
        Y_pred <- predict.PLS(model, newdata = X_test_scaled)$predict[,,1]
        Y_pred <- Y_pred * Y_sds + Y_means
        
        # deflation matrices
        ui <- model$loadings$X[,1]
        vi <- model$loadings$Y[,1]
        res.deflat <- step2.spls(X=X_train_scaled,Y=Y_train_scaled,ui,vi,mode="regression")
        X_train_scaled = res.deflat$X.h
        Y_train_scaled = res.deflat$Y.h
        Y_test[[h+1]][ind.test,] <- Y_test[[h]][ind.test,] - Y_pred
        
      }
      
    }
    
    for(h in 1:ncomp){
      
      # deflation matrices
      u <- object$loadings$X[,h]
      v <- object$loadings$Y[,h]
      res.deflat <- step2.spls(X=X,Y=Y,u,v,mode="regression")
      X = res.deflat$X.h
      Y = res.deflat$Y.h
      
      for(j in 1:q){
        RSSj[h,j] <- sum((Y[,j])^2)
        PRESSj[h,j] <- sum((Y_test[[h+1]][,j])^2)
      }
      colnames(RSSj) <- paste0("Y",1:q)
      colnames(PRESSj) <- paste0("Y",1:q)
      
      # RSS computing
      RSS[h] <- sum(RSSj[h,]) 
      
      # PRESSh computing
      PRESS[h] <- sum(PRESSj[h,]) 
      
      # Q2
      q2[h] <- 1-PRESS[h]/RSS[max(h-1,1)]
      
      
    }# end h loop
    
    # first value correction 
    q2[1] <- 1-PRESS[1]/RSS0
    
    lim <- 0.0975
    
    h <- 1
    while(!is.na(q2[h]) && q2[h] > lim){
      h <- h + 1
    }
    h.best.q2 <- max(h-1,1)
    
    # Plot
    if(plot){
      plot(q2, type = "b", col = "blue", pch = 16,
           main = "Model Q² performance",
           xlab = "Number of components", ylab = "Q²", axes = FALSE)
      abline(h = lim, col = "red", lty = 2)
      axis(1, at = 1:ncomp)
      axis(2, labels = TRUE)
    }
    
    res$q2 <- q2
    res$PRESS = PRESS
    res$RSS = RSS
    res$PRESSj = PRESSj
    res$RSSj = RSSj
    res$h.best.q2 = h.best.q2
    
    #return(res)
    
  } # end criterion loop
  
  return(res)
  
}

# ---------------------------------------------------
# perf for sPLS object (0) ----
# ---------------------------------------------------

perf0.sPLS <-
  function(object,
           criterion = c("all","MSEP", "R2", "Q2"), 
           validation = c("Mfold", "loo"),
           folds = 10,
           progressBar = TRUE,setseed=1,
           ...)
  {
    set.seed(setseed)
    #-- validation des arguments --#
    # these are the centered and scaled matrices output from spls 
    X = object$X
    Y = object$Y
    
    tol = object$tol
    max.iter = object$max.iter
    
    # tells which variables are selected in X and in Y:
    keepX = object$keepX   
    keepY = object$keepY   
    
    mode = object$mode
    ncomp = object$ncomp
    n = nrow(X)
    p = ncol(X)
    q = ncol(Y)
    res = list()
    
    validation = match.arg(validation)
    
    # initialize new objects:= to record feature stability
    featuresX  = featuresY =  list()
    for(k in 1:ncomp){
      featuresX[[k]] = featuresY[[k]] = NA
    }
    
    
    if (length(dim(X)) != 2) 
      stop("'X' must be a numeric matrix for validation.")
    
    if(object$mode == 'canonical') stop('sPLS mode should be set to regression, invariant or classic')
    
    if (any(criterion == "Q2") & ncomp == 1)
      stop("'ncomp' must be > 1 for Q2 criterion.")
    
    if (any(is.na(X)) || any(is.na(Y))) 
      stop("Missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.")
    
    
    #-- M fold or loo cross validation --#
    #- define the folds
    if (validation == "Mfold") {
      if (is.list(folds)) {
        if (length(folds) < 2 | length(folds) > n)
          stop("Invalid number of folds.")
        if (length(unique(unlist(folds))) != n)
          stop("Invalid folds.")
        
        M = length(folds)
      }
      else {
        if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
          stop("Invalid number of folds.")
        else {
          M = round(folds)
          folds = split(sample(1:n), rep(1:M, length = n)) 
        }
      }
    } 
    else { 
      folds = split(1:n, rep(1:n, length = n)) 
      M = n
    }
    
    
    #-- compute criteria ------------ --#
    RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
    RSS.indiv = array(NA, c(n, q, ncomp+1))
    PRESS.inside = Q2.inside = matrix(nrow = ncomp, ncol = q)
    
    #KA: all criteria were included.
    if (any(criterion %in% c("all", "MSEP", "R2", "Q2"))) {
      press.mat = Ypred = array(NA, c(n, q, ncomp))
      MSEP = R2 = matrix(NA, nrow = q, ncol = ncomp)
      
      # set up dimnames
      rownames(MSEP) = rownames(R2) = colnames(Q2.inside) = colnames(Y)
      dimnames(press.mat)[[2]] = colnames(Y)
      
      # in case the test set only includes one sample, it is better to advise the user to
      # perform loocv
      stop.user = FALSE
      if (progressBar == TRUE) pb <- txtProgressBar(style = 3)
      
      for (i in 1:M) {
        if (progressBar == TRUE) setTxtProgressBar(pb, i/M)
        
        omit = folds[[i]]
        # see below, we stop the user if there is only one sample drawn on the test set using MFold
        if(length(omit) == 1) stop.user = TRUE
        
        # the training set is NOT scaled
        X.train = X[-omit, ]
        Y.train = Y[-omit, ]
        X.test = matrix(X[omit, ], nrow = length(omit))
        Y.test = matrix(Y[omit, ], nrow = length(omit))
        
        
        #-- spls --#
        spls.res = sPLS(X.train, Y.train, ncomp, mode, max.iter, tol, keepX=keepX, keepY=keepY)     ## change
        
        # added: record selected features in each set
        for(k in 1:ncomp){
          featuresX[[k]] = c(unlist(featuresX[[k]]), select.var(spls.res, comp = k)$X$name)
          featuresY[[k]] = c(unlist(featuresY[[k]]), select.var(spls.res, comp = k)$Y$name)
        }
        
        
        #if (!is.null(spls.res$nzv$Position)) X.test = X.test[, -spls.res$nzv$Position]
        Y.hat = predict(spls.res, X.test)$predict
        
        #compute prediction of Y
        for (h in 1:ncomp) {
          Ypred[omit, , h] = Y.hat[, , h]
          
          # KA: this bunch was added:
          # compute the press and the RSS
          # by definition (tenenhaus), RSS[h+1,] = (y_i - y.hat_(h-1)_i)^2
          press.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
          RSS.indiv[omit, ,h+1] = (Y.test - Y.hat[, , h])^2
        } # end h
      } #end i (cross validation)
      
      # KA added
      # warn the user that at least test set had a length of 1
      if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
      
      # these criteria are computed across all folds.
      for (h in 1:ncomp) { 
        MSEP[, h] = apply(as.matrix(press.mat[, , h]), 2, mean, na.rm = TRUE)
        R2[, h] = (diag(cor(Y, Ypred[,,h], use = "pairwise")))^2
        #KA:  PRESS is also computed as well as Q2 inside this procedure 
        if(q>1){
          RSS[h+1,] = t(apply(RSS.indiv[,,h+1], 2, sum))
          PRESS.inside[h, ] = colSums(press.mat[, , h], na.rm = TRUE)
        }else{
          RSS[h+1,q] = sum(RSS.indiv[,q,h+1])
          PRESS.inside[h, q] = sum(press.mat[,q, h], na.rm = TRUE)
        }
        
        Q2.inside[h, ] = 1 - PRESS.inside[h, ]/RSS[h, ]
        
      }
      
      colnames(MSEP) = colnames(R2) = rownames(Q2.inside) = paste('ncomp', c(1:ncomp), sep = " ")
      rownames(MSEP) = rownames(R2) = colnames(Q2.inside)  = colnames(Y)
      
      if (q == 1) rownames(MSEP) = rownames(R2) = ""
      
      if(ncomp>1){
        # compute Q2 total
        if(q>1){
          Q2.total = 1 - rowSums(PRESS.inside, na.rm = TRUE)/rowSums(RSS[-(ncomp+1), ], na.rm = TRUE)
        }else{ # for q == 1
          Q2.total = t(1 - PRESS.inside/RSS[-(ncomp+1), ])
        }
      }else{
        Q2.total = NA
      }
      
      names(Q2.total) = paste('comp', 1:ncomp, sep = " ")
      
    } # end all, MSEP, R2
    
    if (progressBar == TRUE) cat('\n')
    
    
    # ---- extract stability of features ----- # NEW
    list.featuresX = list.featuresY =list()
    for(k in 1:ncomp){
      #remove the NA value that was added for initialisation
      remove.naX = which(is.na(featuresX[[k]]))
      remove.naY = which(is.na(featuresY[[k]]))
      # then summarise as a factor and output the percentage of appearance
      list.featuresX[[k]] = sort(summary(as.factor(featuresX[[k]][-remove.naX]))/M, decreasing = TRUE)
      list.featuresY[[k]] = sort(summary(as.factor(featuresY[[k]][-remove.naY]))/M, decreasing = TRUE)
    }
    
    # extract features selected from the full model ---------
    features.finalX = features.finalY =list()
    for(k in 1:ncomp){
      features.finalX[[k]] = select.var(object, comp = k)$X$value
      features.finalY[[k]] = select.var(object, comp = k)$Y$value
    }
    
    names(features.finalX)  = names(features.finalY) = names(list.featuresX) = names(list.featuresY) = paste('comp', 1:ncomp)
    
    
    # ----------  final outputs:
    if (any(criterion %in% c("all", "MSEP"))) res$MSEP = MSEP
    if (any(criterion %in% c("all", "R2"))) res$R2 = R2
    if (any(criterion %in% c("all", "Q2"))) res$Q2 = t(Q2.inside)
    if (any(criterion %in% c("all", "Q2"))) res$Q2.total = Q2.total
    
    
    #features
    res$features$stable.X = list.featuresX
    res$features$stable.Y = list.featuresY
    res$features$final.X = features.finalX
    res$features$final.Y = features.finalY
    
    res$press.mat=press.mat
    res$RSS.indiv=RSS.indiv
    res$PRESS.inside=PRESS.inside
    res$RSS=RSS
    
    method = "pls.mthd"
    class(res) = c("perf", method)
    return(invisible(res))
  }

# ---------------------------------------------------
# perf for sPLS object ----
# ---------------------------------------------------

perf.sPLS <- function(object, criterion = c("all","MSEP","Q2"), validation = c("Mfold","loo"),
                      folds = 10, progressBar = TRUE, setseed = 1, plot = FALSE){
  
  X <- object$X
  Y <- object$Y
  c <- object$mat.c
  d <- object$mat.d
  ncomp <- object$ncomp
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  # sPLS model attributes
  tol <- object$tol
  max.iter <- object$max.iter
  keepX <- object$keepX   
  keepY <- object$keepY   
  
  # conditions about sPLS model attributes
  if(is.null(tol)){tol = 10^(-6)}
  if(is.null(max.iter)){max.iter = 500}
  if(is.null(keepX)){keepX = rep(ncol(X),ncomp)}
  if(is.null(keepY)){keepY = rep(ncol(Y),ncomp)}
  
  # others conditions check-up
  if(!("pls" %in% class(object)) && class(object) != "mixo_pls"){ stop("object class must either contain pls class or be mixo_pls class."); print(class(object))}
  
  if(validation[1] == "Mfold"){
    if(folds < 2 || folds > n){ stop(paste("folds must be a value between 2 and",n))}
    K = folds
  }else{K = n}
  
  if (progressBar == TRUE){
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb,1)
    cat('\n')
  }
  
  set.seed(setseed)
  ind <- sample(1:n,n,replace = FALSE)
  b <- floor(n/K) # block size
  res = list()
  
  # MSEP CRITERION
  if(any(criterion %in% c("all","MSEP"))){
    
    # prediction analysis
    err <- matrix(NA, nrow = K, ncol = ncomp)
    MSEPj <- matrix(0, nrow = ncomp, ncol = q)
    colnames(MSEPj) <- paste0("Y",1:q)
    rownames(MSEPj) <- paste0("Comp",1:ncomp)
    
    for(k in seq_len(K)){
      
      # bloc definition
      ind.beg <- (k-1) * b + 1 # block k beginning
      ind.end <- k * b # block k end
      ind.test <- ind[ind.beg:ind.end]
      nk <- length(ind.test)
      X.train <- X[-ind.test,]
      Y.train <- Y[-ind.test,]
      X.test <- X[ind.test,]
      Y.test <- Y[ind.test,]
      modele <- sPLS(X = X.train, Y = Y.train, ncomp = ncomp, mode = "regression", keepX = keepX, keepY = keepY, tol = tol, max.iter = max.iter)  
      
      for(h in 1:ncomp){  
        
        # predictions
        pred <- predict.PLS(modele, newdata = X.test)$predict[,,h]
        MSEPj[h,] <- MSEPj[h,] + colMeans(as.matrix((Y.test - pred)^2))/K
        
      }
    }
    err.moy <- rowMeans(MSEPj)
    
    h.best.msep <- min(which.min(err.moy))
    if(plot){
      plot(err.moy, col="blue", pch = 16, type = "b", main = "MSEP of the model", xlab = "number of components", ylab = "MSEP", axes = FALSE)
      axis(1, at = 1:ncomp)
      axis(2, labels = TRUE)
    }
    
    res$MSEP <- err.moy
    res$MSEPj <- MSEPj
    res$h.best.msep <- h.best.msep
    
  }
  
  if(any(criterion %in% c("all","Q2"))){
    
    q2 <- numeric(ncomp)
    PRESS <- numeric(ncomp)
    RSS <- numeric(ncomp)
    
    RSSj <- matrix(nrow = ncomp, ncol = q)
    PRESSj <- matrix(nrow = ncomp, ncol = q)
    
    # RSS0 computing
    Y.mean <- t(matrix(colMeans(Y), nrow = q, ncol = n)) # mean for each column
    RSS0 <- sum(colSums((Y-Y.mean)^2))
    
    
    Y_test <- list()
    Y_test[[1]] <- Y
    
    for(h in 1:ncomp){
      Y_test[[h+1]] <- matrix(nrow = n, ncol = q)
    }
    
    for(k in seq_len(K)){
      
      # bloc definition
      ind.beg <- (k-1) * b + 1 # block k beginning
      ind.end <- k * b # block k end
      ind.test <- ind[ind.beg:ind.end]
      nk <- length(ind.test)
      X_train <- X[-ind.test,, drop = FALSE]
      Y_train <- Y[-ind.test,, drop = FALSE]
      
      #Y.test <- Y[ind.test,]
      
      
      # Center and scale training data
      X_train_scaled <- scale(X_train, center = TRUE, scale = TRUE)
      Y_train_scaled <- scale(Y_train, center = TRUE, scale = TRUE)
      
      X_means <- apply(X_train,2,mean)
      X_sds <- apply(X_train,2,sd)
      Y_means <- apply(Y_train,2,mean)
      Y_sds <- apply(Y_train,2,sd)
      
      
      X_test <- X[ind.test,]
      X_test_scaled <- (X_test - X_means) / X_sds
      
      for(h in 1:ncomp){
        
        # training on the dataset without the ith individual
        model <- sPLS(X = X_train_scaled, Y = Y_train_scaled, ncomp = 1, mode = "regression", keepX = keepX, keepY = keepY, tol = tol, max.iter = max.iter)
        
        # predictions on the ith individual
        Y_pred <- predict.PLS(model, newdata = X_test_scaled)$predict[,,1]
        Y_pred <- Y_pred * Y_sds + Y_means
        
        # deflation matrices
        ui <- model$loadings$X[,1]
        vi <- model$loadings$Y[,1]
        res.deflat <- step2.spls(X=X_train_scaled,Y=Y_train_scaled,ui,vi,mode="regression")
        X_train_scaled = res.deflat$X.h
        Y_train_scaled = res.deflat$Y.h
        Y_test[[h+1]][ind.test,] <- Y_test[[h]][ind.test,] - Y_pred
        
      }
      
    }
    
    for(h in 1:ncomp){
      
      # deflation matrices
      u <- object$loadings$X[,h]
      v <- object$loadings$Y[,h]
      res.deflat <- step2.spls(X=X,Y=Y,u,v,mode="regression")
      X = res.deflat$X.h
      Y = res.deflat$Y.h
      
      for(j in 1:q){
        RSSj[h,j] <- sum((Y[,j])^2)
        PRESSj[h,j] <- sum((Y_test[[h+1]][,j])^2)
      }
      colnames(RSSj) <- paste0("Y",1:q)
      colnames(PRESSj) <- paste0("Y",1:q)
      
      # RSS computing
      RSS[h] <- sum(RSSj[h,]) 
      
      # PRESSh computing
      PRESS[h] <- sum(PRESSj[h,]) 
      
      # Q2
      q2[h] <- 1-PRESS[h]/RSS[max(h-1,1)]
      
      
    }# end h loop
    
    # first value correction 
    q2[1] <- 1-PRESS[1]/RSS0
    
    lim <- 0.0975
    
    h <- 1
    while(!is.na(q2[h]) && q2[h] > lim){
      h <- h + 1
    }
    h.best.q2 <- max(h-1,1)
    
    # Plot
    if(plot){
      plot(q2, type = "b", col = "blue", pch = 16,
           main = "Model Q² performance",
           xlab = "Number of components", ylab = "Q²", axes = FALSE)
      abline(h = lim, col = "red", lty = 2)
      axis(1, at = 1:ncomp)
      axis(2, labels = TRUE)
    }
    
    res$q2 <- q2
    res$PRESS = PRESS
    res$RSS = RSS
    res$PRESSj = PRESSj
    res$RSSj = RSSj
    res$h.best.q2 = h.best.q2
    
  } # end criterion loop
  
  method = "pls.mthd"
  class(res) = c("perf", method)
  return(invisible(res))
  
}


# ---------------------------------------------------
# perf for gPLS object (0) -----
# ---------------------------------------------------

perf0.gPLS <-
  function(object,
           criterion = c("all","MSEP", "R2", "Q2"), 
           validation = c("Mfold", "loo"),
           folds = 10,
           progressBar = TRUE,setseed=1,
           ...)
  {
    set.seed(setseed)
    #-- validation des arguments --#
    # these are the centered and scaled matrices output from spls 
    X = object$X
    Y = object$Y
    
    tol = object$tol
    max.iter = object$max.iter
    
    # tells which variables are selected in X and in Y:
    keepX = object$keepX   
    keepY = object$keepY   
    
    mode = object$mode
    ncomp = object$ncomp
    ind.block.x = object$ind.block.x
    ind.block.y = object$ind.block.y

    
    n = nrow(X)
    p = ncol(X)
    q = ncol(Y)
    res = list()
    
    validation = match.arg(validation)
    
    # initialize new objects:= to record feature stability
    featuresX  = featuresY =  list()
    for(k in 1:ncomp){
      featuresX[[k]] = featuresY[[k]] = NA
    }
    
    
    if (length(dim(X)) != 2) 
      stop("'X' must be a numeric matrix for validation.")
    
    if(object$mode == 'canonical') stop('sPLS mode should be set to regression, invariant or classic')
    
    if (any(criterion == "Q2") & ncomp == 1)
      stop("'ncomp' must be > 1 for Q2 criterion.")
    
    if (any(is.na(X)) || any(is.na(Y))) 
      stop("Missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.")
    
    
    #-- M fold or loo cross validation --#
    #- define the folds
    if (validation == "Mfold") {
      if (is.list(folds)) {
        if (length(folds) < 2 | length(folds) > n)
          stop("Invalid number of folds.")
        if (length(unique(unlist(folds))) != n)
          stop("Invalid folds.")
        
        M = length(folds)
      }
      else {
        if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
          stop("Invalid number of folds.")
        else {
          M = round(folds)
          folds = split(sample(1:n), rep(1:M, length = n)) 
        }
      }
    } 
    else { 
      folds = split(1:n, rep(1:n, length = n)) 
      M = n
    }
    
    
    #-- compute criteria ------------ --#
    RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
    RSS.indiv = array(NA, c(n, q, ncomp+1))
    PRESS.inside = Q2.inside = matrix(nrow = ncomp, ncol = q)
    
    #KA: all criteria were included.
    if (any(criterion %in% c("all", "MSEP", "R2", "Q2"))) {
      press.mat = Ypred = array(NA, c(n, q, ncomp))
      MSEP = R2 = matrix(NA, nrow = q, ncol = ncomp)
      
      # set up dimnames
      rownames(MSEP) = rownames(R2) = colnames(Q2.inside) = colnames(Y)
      dimnames(press.mat)[[2]] = colnames(Y)
      
      # in case the test set only includes one sample, it is better to advise the user to
      # perform loocv
      stop.user = FALSE
      if (progressBar == TRUE) pb <- txtProgressBar(style = 3)
      
      for (i in 1:M) {
        if (progressBar == TRUE) setTxtProgressBar(pb, i/M)
        
        omit = folds[[i]]
        # see below, we stop the user if there is only one sample drawn on the test set using MFold
        if(length(omit) == 1) stop.user = TRUE
        
        # the training set is NOT scaled
        X.train = X[-omit, ]
        Y.train = Y[-omit, ]
        X.test = matrix(X[omit, ], nrow = length(omit))
        Y.test = matrix(Y[omit, ], nrow = length(omit))
        
        #-- gpls --#
        spls.res = gPLS(X.train, Y.train, ncomp, mode, max.iter, tol, keepX=keepX, keepY=keepY,ind.block.x=ind.block.x,ind.block.y=ind.block.y)     ## change
        
        # added: record selected features in each set
        for(k in 1:ncomp){
          featuresX[[k]] = c(unlist(featuresX[[k]]), select.var(spls.res, comp = k)$X$name)
          featuresY[[k]] = c(unlist(featuresY[[k]]), select.var(spls.res, comp = k)$Y$name)
        }
        
        
        #if (!is.null(spls.res$nzv$Position)) X.test = X.test[, -spls.res$nzv$Position]
        Y.hat = predict.gPLS(spls.res, X.test)$predict
        
        #compute prediction of Y
        for (h in 1:ncomp) {
          Ypred[omit, , h] = Y.hat[, , h]
          
          # KA: this bunch was added:
          # compute the press and the RSS
          # by definition (tenenhaus), RSS[h+1,] = (y_i - y.hat_(h-1)_i)^2
          press.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
          RSS.indiv[omit, ,h+1] = (Y.test - Y.hat[, , h])^2
        } # end h
      } #end i (cross validation)
      
      # KA added
      # warn the user that at least test set had a length of 1
      if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
      
      # these criteria are computed across all folds.
      for (h in 1:ncomp) { 
        MSEP[, h] = apply(as.matrix(press.mat[, , h]), 2, mean, na.rm = TRUE)
        R2[, h] = (diag(cor(Y, Ypred[,,h], use = "pairwise")))^2
        #KA:  PRESS is also computed as well as Q2 inside this procedure 
        if(q>1){
          RSS[h+1,] = t(apply(RSS.indiv[,,h+1], 2, sum))
          PRESS.inside[h, ] = colSums(press.mat[, , h], na.rm = TRUE)
        }else{
          RSS[h+1,q] = sum(RSS.indiv[,q,h+1])
          PRESS.inside[h, q] = sum(press.mat[,q, h], na.rm = TRUE)
        }
        
        Q2.inside[h, ] = 1 - PRESS.inside[h, ]/RSS[h, ]
        
      }
      
      colnames(MSEP) = colnames(R2) = rownames(Q2.inside) = paste('ncomp', c(1:ncomp), sep = " ")
      rownames(MSEP) = rownames(R2) = colnames(Q2.inside)  = colnames(Y)
      
      if (q == 1) rownames(MSEP) = rownames(R2) = ""
      
      if(ncomp>1){
        # compute Q2 total
        if(q>1){
          Q2.total = 1 - rowSums(PRESS.inside, na.rm = TRUE)/rowSums(RSS[-(ncomp+1), ], na.rm = TRUE)
        }else{ # for q == 1
          Q2.total = t(1 - PRESS.inside/RSS[-(ncomp+1), ])
        }
      }else{
        Q2.total = NA
      }
      
      names(Q2.total) = paste('comp', 1:ncomp, sep = " ")
      
    } # end all, MSEP, R2
    
    if (progressBar == TRUE) cat('\n')
    
    
    # ---- extract stability of features ----- # NEW
    list.featuresX = list.featuresY =list()
    for(k in 1:ncomp){
      #remove the NA value that was added for initialisation
      remove.naX = which(is.na(featuresX[[k]]))
      remove.naY = which(is.na(featuresY[[k]]))
      # then summarise as a factor and output the percentage of appearance
      list.featuresX[[k]] = sort(summary(as.factor(featuresX[[k]][-remove.naX]))/M, decreasing = TRUE)
      list.featuresY[[k]] = sort(summary(as.factor(featuresY[[k]][-remove.naY]))/M, decreasing = TRUE)
    }
    
    # extract features selected from the full model ---------
    features.finalX = features.finalY =list()
    for(k in 1:ncomp){
      features.finalX[[k]] = select.var(object, comp = k)$X$value
      features.finalY[[k]] = select.var(object, comp = k)$Y$value
    }
    
    names(features.finalX)  = names(features.finalY) = names(list.featuresX) = names(list.featuresX) = paste('comp', 1:ncomp)
    
    
    # ----------  final outputs:
    if (any(criterion %in% c("all", "MSEP"))) res$MSEP = MSEP
    if (any(criterion %in% c("all", "R2"))) res$R2 = R2
    if (any(criterion %in% c("all", "Q2"))) res$Q2 = t(Q2.inside)
    if (any(criterion %in% c("all", "Q2"))) res$Q2.total = Q2.total
    
    
    #features
    res$features$stable.X = list.featuresX
    res$features$stable.Y = list.featuresY
    res$features$final.X = features.finalX
    res$features$final.Y = features.finalY
    
    res$press.mat=press.mat
    res$RSS.indiv=RSS.indiv
    res$PRESS.inside=PRESS.inside
    res$RSS=RSS
    
    method = "pls.mthd"
    class(res) = c("perf", method)
    return(invisible(res))
  }

# ---------------------------------------------------
# perf for gPLS object ----
# ---------------------------------------------------

perf.gPLS <- function(object, criterion = c("all","MSEP","Q2"), validation = c("Mfold","loo"),
                      folds = 10, progressBar = TRUE, setseed = 1, plot = FALSE){
  
  X <- object$X
  Y <- object$Y
  c <- object$mat.c
  d <- object$mat.d
  ncomp <- object$ncomp
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  # sPLS model attributes
  tol <- object$tol
  max.iter <- object$max.iter
  keepX <- object$keepX   
  keepY <- object$keepY 
  ind.block.x = object$ind.block.x
  ind.block.y = object$ind.block.y
  
  # conditions about sPLS model attributes
  if(is.null(tol)){tol = 10^(-6)}
  if(is.null(max.iter)){max.iter = 500}
  if(is.null(keepX)){keepX = rep(ncol(X),ncomp)}
  if(is.null(keepY)){keepY = rep(ncol(Y),ncomp)}
  if(is.null(ind.block.x)){stop("Model object doesn't have ind.block.x attribute.")}
  
  # others conditions check-up
  if(!("pls" %in% class(object)) && class(object) != "mixo_pls"){ stop("object class must either contain pls class or be mixo_pls class."); print(class(object))}
  
  if(validation[1] == "Mfold"){
    if(folds < 2 || folds > n){ stop(paste("folds must be a value between 2 and",n))}
    K = folds
  }else{K = n}
  
  if (progressBar == TRUE){
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb,1)
    cat('\n')
  }
  
  set.seed(setseed)
  ind <- sample(1:n,n,replace = FALSE)
  b <- floor(n/K) # block size
  res = list()
  
  # MSEP CRITERION
  if(any(criterion %in% c("all","MSEP"))){
    
    # prediction analysis
    err <- matrix(NA, nrow = K, ncol = ncomp)
    MSEPj <- matrix(0, nrow = ncomp, ncol = q)
    colnames(MSEPj) <- paste0("Y",1:q)
    rownames(MSEPj) <- paste0("Comp",1:ncomp)
    
    for(k in seq_len(K)){
      
      # bloc definition
      ind.beg <- (k-1) * b + 1 # block k beginning
      ind.end <- k * b # block k end
      ind.test <- ind[ind.beg:ind.end]
      nk <- length(ind.test)
      X.train <- X[-ind.test,]
      Y.train <- Y[-ind.test,]
      X.test <- X[ind.test,]
      Y.test <- Y[ind.test,]
      modele <- gPLS(X = X.train, Y = Y.train, ncomp = ncomp, mode = "regression", keepX = keepX, keepY = keepY, tol = tol, max.iter = max.iter, ind.block.x = ind.block.x, ind.block.y = ind.block.y)  
      
      for(h in 1:ncomp){  
        
        # predictions
        pred <- predict.PLS(modele, newdata = X.test)$predict[,,h]
        MSEPj[h,] <- MSEPj[h,] + colMeans(as.matrix((Y.test - pred)^2))/K
        
      }
    }
    err.moy <- rowMeans(MSEPj)
    
    h.best.msep <- min(which.min(err.moy))
    if(plot){
      plot(err.moy, col="blue", pch = 16, type = "b", main = "MSEP of the model", xlab = "number of components", ylab = "MSEP", axes = FALSE)
      axis(1, at = 1:ncomp)
      axis(2, labels = TRUE)
    }
    
    res$MSEP <- err.moy
    res$MSEPj <- MSEPj
    res$h.best.msep <- h.best.msep
    
  }
  
  if(any(criterion %in% c("all","Q2"))){
    
    q2 <- numeric(ncomp)
    PRESS <- numeric(ncomp)
    RSS <- numeric(ncomp)
    
    RSSj <- matrix(nrow = ncomp, ncol = q)
    PRESSj <- matrix(nrow = ncomp, ncol = q)
    
    # RSS0 computing
    Y.mean <- t(matrix(colMeans(Y), nrow = q, ncol = n)) # mean for each column
    RSS0 <- sum(colSums((Y-Y.mean)^2))
    
    
    Y_test <- list()
    Y_test[[1]] <- Y
    
    for(h in 1:ncomp){
      Y_test[[h+1]] <- matrix(nrow = n, ncol = q)
    }
    
    for(k in seq_len(K)){
      
      # bloc definition
      ind.beg <- (k-1) * b + 1 # block k beginning
      ind.end <- k * b # block k end
      ind.test <- ind[ind.beg:ind.end]
      nk <- length(ind.test)
      X_train <- X[-ind.test,, drop = FALSE]
      Y_train <- Y[-ind.test,, drop = FALSE]
      
      #Y.test <- Y[ind.test,]
      
      
      # Center and scale training data
      X_train_scaled <- scale(X_train, center = TRUE, scale = TRUE)
      Y_train_scaled <- scale(Y_train, center = TRUE, scale = TRUE)
      
      X_means <- apply(X_train,2,mean)
      X_sds <- apply(X_train,2,sd)
      Y_means <- apply(Y_train,2,mean)
      Y_sds <- apply(Y_train,2,sd)
      
      
      X_test <- X[ind.test,]
      X_test_scaled <- (X_test - X_means) / X_sds
      
      for(h in 1:ncomp){
        
        # training on the dataset without the ith individual
        model <- gPLS(X = X_train_scaled, Y = Y_train_scaled, ncomp = 1, mode = "regression", keepX = keepX, keepY = keepY, tol = tol, max.iter = max.iter, ind.block.x = ind.block.x, ind.block.y = ind.block.y)
        
        # predictions on the ith individual
        Y_pred <- predict.PLS(model, newdata = X_test_scaled)$predict[,,1]
        Y_pred <- Y_pred * Y_sds + Y_means
        
        # deflation matrices
        ui <- model$loadings$X[,1]
        vi <- model$loadings$Y[,1]
        res.deflat <- step2.spls(X=X_train_scaled,Y=Y_train_scaled,ui,vi,mode="regression")
        X_train_scaled = res.deflat$X.h
        Y_train_scaled = res.deflat$Y.h
        Y_test[[h+1]][ind.test,] <- Y_test[[h]][ind.test,] - Y_pred
        
      }
      
    }
    
    for(h in 1:ncomp){
      
      # deflation matrices
      u <- object$loadings$X[,h]
      v <- object$loadings$Y[,h]
      res.deflat <- step2.spls(X=X,Y=Y,u,v,mode="regression")
      X = res.deflat$X.h
      Y = res.deflat$Y.h
      
      for(j in 1:q){
        RSSj[h,j] <- sum((Y[,j])^2)
        PRESSj[h,j] <- sum((Y_test[[h+1]][,j])^2)
      }
      colnames(RSSj) <- paste0("Y",1:q)
      colnames(PRESSj) <- paste0("Y",1:q)
      
      # RSS computing
      RSS[h] <- sum(RSSj[h,]) 
      
      # PRESSh computing
      PRESS[h] <- sum(PRESSj[h,]) 
      
      # Q2
      q2[h] <- 1-PRESS[h]/RSS[max(h-1,1)]
      
      
    }# end h loop
    
    # first value correction 
    q2[1] <- 1-PRESS[1]/RSS0
    
    lim <- 0.0975
    
    h <- 1
    while(!is.na(q2[h]) && q2[h] > lim){
      h <- h + 1
    }
    h.best.q2 <- max(h-1,1)
    
    # Plot
    if(plot){
      plot(q2, type = "b", col = "blue", pch = 16,
           main = "Model Q² performance",
           xlab = "Number of components", ylab = "Q²", axes = FALSE)
      abline(h = lim, col = "red", lty = 2)
      axis(1, at = 1:ncomp)
      axis(2, labels = TRUE)
    }
    
    res$q2 <- q2
    res$PRESS = PRESS
    res$RSS = RSS
    res$PRESSj = PRESSj
    res$RSSj = RSSj
    res$h.best.q2 = h.best.q2
    
  } # end criterion loop
  
  method = "pls.mthd"
  class(res) = c("perf", method)
  return(invisible(res))
  
}

# ---------------------------------------------------
# perf for sgPLS object ----
# ---------------------------------------------------

perf.sgPLS <-
  function(object,
           criterion = c("all","MSEP", "R2", "Q2"), 
           validation = c("Mfold", "loo"),
           folds = 10,
           progressBar = TRUE,setseed=1,
           ...)
  {
    set.seed(setseed)
    #-- validation des arguments --#
    # these are the centered and scaled matrices output from spls 
    X = object$X
    Y = object$Y
    
    tol = object$tol
    max.iter = object$max.iter
    
    # tells which variables are selected in X and in Y:
    keepX = object$keepX   
    keepY = object$keepY   
    
    mode = object$mode
    ncomp = object$ncomp
    ind.block.x = object$ind.block.x
    ind.block.y = object$ind.block.y
    alpha.x = object$alpha.x
    alpha.y = object$alpha.y
    upper.lambda = object$upper.lambda
    
    n = nrow(X)
    p = ncol(X)
    q = ncol(Y)
    res = list()
    
    validation = match.arg(validation)
    
    # initialize new objects:= to record feature stability
    featuresX  = featuresY =  list()
    for(k in 1:ncomp){
      featuresX[[k]] = featuresY[[k]] = NA
    }
    
    
    if (length(dim(X)) != 2) 
      stop("'X' must be a numeric matrix for validation.")
    
    if(object$mode == 'canonical') stop('sPLS mode should be set to regression, invariant or classic')
    
    if (any(criterion == "Q2") & ncomp == 1)
      stop("'ncomp' must be > 1 for Q2 criterion.")
    
    if (any(is.na(X)) || any(is.na(Y))) 
      stop("Missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.")
    
    
    #-- M fold or loo cross validation --#
    #- define the folds
    if (validation == "Mfold") {
      if (is.list(folds)) {
        if (length(folds) < 2 | length(folds) > n)
          stop("Invalid number of folds.")
        if (length(unique(unlist(folds))) != n)
          stop("Invalid folds.")
        
        M = length(folds)
      }
      else {
        if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
          stop("Invalid number of folds.")
        else {
          M = round(folds)
          folds = split(sample(1:n), rep(1:M, length = n)) 
        }
      }
    } 
    else { 
      folds = split(1:n, rep(1:n, length = n)) 
      M = n
    }
    
    
    #-- compute criteria ------------ --#
    RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
    RSS.indiv = array(NA, c(n, q, ncomp+1))
    PRESS.inside = Q2.inside = matrix(nrow = ncomp, ncol = q)
    
    #KA: all criteria were included.
    if (any(criterion %in% c("all", "MSEP", "R2", "Q2"))) {
      press.mat = Ypred = array(NA, c(n, q, ncomp))
      MSEP = R2 = matrix(NA, nrow = q, ncol = ncomp)
      
      # set up dimnames
      rownames(MSEP) = rownames(R2) = colnames(Q2.inside) = colnames(Y)
      dimnames(press.mat)[[2]] = colnames(Y)
      
      # in case the test set only includes one sample, it is better to advise the user to
      # perform loocv
      stop.user = FALSE
      if (progressBar == TRUE) pb <- txtProgressBar(style = 3)
      
      for (i in 1:M) {
        if (progressBar == TRUE) setTxtProgressBar(pb, i/M)
        
        omit = folds[[i]]
        # see below, we stop the user if there is only one sample drawn on the test set using MFold
        if(length(omit) == 1) stop.user = TRUE
        
        # the training set is NOT scaled
        X.train = X[-omit, ]
        Y.train = Y[-omit, ]
        X.test = matrix(X[omit, ], nrow = length(omit))
        Y.test = matrix(Y[omit, ], nrow = length(omit))
        
        #-- sgpls --#
        spls.res = sgPLS(X.train, Y.train, ncomp, mode, max.iter, tol, keepX=keepX, keepY=keepY,ind.block.x=ind.block.x,ind.block.y=ind.block.y,alpha.x=alpha.x,alpha.y=alpha.y,upper.lambda=upper.lambda)     ## change
        #  Sparse.Group.spls.BP(X,Y,ncomp=1,mode="regression",max.iter=500,tol=1e-06,keepX=c(4),keepY=c(4),ind.block.y=ind.block.y,ind.block.x=ind.block.x,alpha.x=0.05,alpha.y=0.95,upper.lambda=1000000000)
        #res.sparse <- s
        # added: record selected features in each set
        for(k in 1:ncomp){
          featuresX[[k]] = c(unlist(featuresX[[k]]), select.var(spls.res, comp = k)$X$name)
          featuresY[[k]] = c(unlist(featuresY[[k]]), select.var(spls.res, comp = k)$Y$name)
        }
        
        
        #if (!is.null(spls.res$nzv$Position)) X.test = X.test[, -spls.res$nzv$Position]
        Y.hat = predict.sgPLS(spls.res, X.test)$predict
        
        #compute prediction of Y
        for (h in 1:ncomp) {
          Ypred[omit, , h] = Y.hat[, , h]
          
          # KA: this bunch was added:
          # compute the press and the RSS
          # by definition (tenenhaus), RSS[h+1,] = (y_i - y.hat_(h-1)_i)^2
          press.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
          RSS.indiv[omit, ,h+1] = (Y.test - Y.hat[, , h])^2
        } # end h
      } #end i (cross validation)
      
      # KA added
      # warn the user that at least test set had a length of 1
      if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
      
      # these criteria are computed across all folds.
      for (h in 1:ncomp) { 
        MSEP[, h] = apply(as.matrix(press.mat[, , h]), 2, mean, na.rm = TRUE)
        R2[, h] = (diag(cor(Y, Ypred[,,h], use = "pairwise")))^2
        #KA:  PRESS is also computed as well as Q2 inside this procedure 
        if(q>1){
          RSS[h+1,] = t(apply(RSS.indiv[,,h+1], 2, sum))
          PRESS.inside[h, ] = colSums(press.mat[, , h], na.rm = TRUE)
        }else{
          RSS[h+1,q] = sum(RSS.indiv[,q,h+1])
          PRESS.inside[h, q] = sum(press.mat[,q, h], na.rm = TRUE)
        }
        
        Q2.inside[h, ] = 1 - PRESS.inside[h, ]/RSS[h, ]
        
      }
      
      colnames(MSEP) = colnames(R2) = rownames(Q2.inside) = paste('ncomp', c(1:ncomp), sep = " ")
      rownames(MSEP) = rownames(R2) = colnames(Q2.inside)  = colnames(Y)
      
      if (q == 1) rownames(MSEP) = rownames(R2) = ""
      
      if(ncomp>1){
        # compute Q2 total
        if(q>1){
          Q2.total = 1 - rowSums(PRESS.inside, na.rm = TRUE)/rowSums(RSS[-(ncomp+1), ], na.rm = TRUE)
        }else{ # for q == 1
          Q2.total = t(1 - PRESS.inside/RSS[-(ncomp+1), ])
        }
      }else{
        Q2.total = NA
      }
      
      names(Q2.total) = paste('comp', 1:ncomp, sep = " ")
      
    } # end all, MSEP, R2
    
    if (progressBar == TRUE) cat('\n')
    
    
    # ---- extract stability of features ----- # NEW
    list.featuresX = list.featuresY =list()
    for(k in 1:ncomp){
      #remove the NA value that was added for initialisation
      remove.naX = which(is.na(featuresX[[k]]))
      remove.naY = which(is.na(featuresY[[k]]))
      # then summarise as a factor and output the percentage of appearance
      list.featuresX[[k]] = sort(summary(as.factor(featuresX[[k]][-remove.naX]))/M, decreasing = TRUE)
      list.featuresY[[k]] = sort(summary(as.factor(featuresY[[k]][-remove.naY]))/M, decreasing = TRUE)
    }
    
    # extract features selected from the full model ---------
    features.finalX = features.finalY =list()
    for(k in 1:ncomp){
      features.finalX[[k]] = select.var(object, comp = k)$X$value
      features.finalY[[k]] = select.var(object, comp = k)$Y$value
    }
    
    names(features.finalX)  = names(features.finalY) = names(list.featuresX) = names(list.featuresX) = paste('comp', 1:ncomp)
    
    
    # ----------  final outputs:
    if (any(criterion %in% c("all", "MSEP"))) res$MSEP = MSEP
    if (any(criterion %in% c("all", "R2"))) res$R2 = R2
    if (any(criterion %in% c("all", "Q2"))) res$Q2 = t(Q2.inside)
    if (any(criterion %in% c("all", "Q2"))) res$Q2.total = Q2.total
    
    
    #features
    res$features$stable.X = list.featuresX
    res$features$stable.Y = list.featuresY
    res$features$final.X = features.finalX
    res$features$final.Y = features.finalY
    
    res$press.mat=press.mat
    res$RSS.indiv=RSS.indiv
    res$PRESS.inside=PRESS.inside
    res$RSS=RSS
    
    method = "pls.mthd"
    class(res) = c("perf", method)
    return(invisible(res))
  }

#------------------------------------------------------------------------
  # perf for PLSda object ----
#------------------------------------------------------------------------

perf.PLSda <- function(object,
                        method.predict = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
                        validation = c("Mfold", "loo"), 
                        folds = 10, progressBar = TRUE){
  
  ncomp <- object$ncomp
  X <- object$X
  Y <- map(object$Y)
  n <- nrow(X)
  p <- ncol(X)
  method <- method.predict
  
  # conditions check-up
  if(!("plsda" %in% class(object)) && !("mixo_plsda" %in% class(object))){ stop("object class must either contain plsda class or contain mixo_plsda class."); print(class(object))}
  
  #if(ncomp > object$ncomp || ncomp <= 0){ stop(paste("ncomp.max must be a value between 0 and",object$ncomp,"which is the total number of components computed in the object model."))}
  
  if(method[1] == "max.dist"||method[2] == "max.dist"){dist=1}else if(method[1] == "centroids.dist"){dist=2}else{dist=3}
  
  if(validation[1] == "Mfold"){
    if(folds < 2 || folds > n){ stop(paste("folds must be a value between 2 and",n))}
    K = folds
  }else{K = n}
  
  if (progressBar == TRUE){
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb,1)
  }
  cat('\n')
  
  b <- floor(n/K) # block size
  ind <- 1:n
  
  # prediction analysis
  err <- matrix(NA, nrow = K, ncol = ncomp)
  
  for(k in seq_len(K)){
    
    # bloc definition
    ind.beg <- (k-1) * b + 1 # block k beginning
    ind.end <- k * b # block k end
    ind.test <- ind[ind.beg:ind.end]
    X.train <- X[-ind.test,]
    Y.train <- Y[-ind.test]
    X.test <- X[ind.test,]
    Y.test <- Y[ind.test]
    modele <- PLSda(X = X.train,Y = Y.train, ncomp = ncomp)
    
    for(h in 1:ncomp){
      # model created
      pred <- predict.PLSda(modele, newdata = X.test, methode = methode)$class[[dist]][,h]
      equal <- Y.test == pred
      err[k,h] <- sum(1-equal)
      
    }
  }
  
  err.moy <- colSums(err)/b/K
  
  h.best <- min(which.min(err.moy))
  plot(err.moy, col="blue", pch = 16, type = "b", main = "Error rate of the model", xlab = "number of components", ylab = "Error")
  abline(h = (1:9)/10, lty = 3, col = "grey")
  
  return(setNames(list(err.moy,h.best),c("error.rate","h.best")))
}

# ---------------------------------------------------
# perf for splsda object ----
# ---------------------------------------------------
perf.sPLSda <- function(object,                                               
                        method.predict = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),     
                        validation = c("Mfold", "loo"),                                          
                        folds = 10,                                                              
                        progressBar = TRUE,...)                                                        
{
  
  #-- initialising arguments --#
  # these data are the centered and scaled X output or the unmapped(Y) scaled and centered
  X = object$X
  level.Y = object$names$Y  #to make sure the levels are ordered
  Y = object$ind.mat
  Y = map(Y)
  Y = factor(Y,labels = level.Y)
  ncomp = object$ncomp
  n = nrow(X)
  keepX = object$keepX  
  
  tol = object$tol
  max.iter = object$max.iter
  
  # initialize new objects:
  features <- list()
  for(k in 1:ncomp){
    features[[k]] = NA
  }
  
  method.predict = match.arg(method.predict, choices = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"), several.ok = TRUE)
  if (any(method.predict == "all")) nmthdd = 3 
  else nmthdd = length(method.predict)  
  
  
 
  
  
  error.fun = function(x, y) {
    error.vec = sweep(x, 1, y, FUN = "-")
    error.vec = (error.vec != 0)
    error.vec = apply(error.vec, 2, sum) / length(y)
    return(error.vec)
  }
  
  #-- define the folds --#
  if (validation == "Mfold") {
    if (is.list(folds)) {
      if (length(folds) < 2 | length(folds) > n)
        stop("Invalid number of folds.")
      if (length(unique(unlist(folds))) != n)
        stop("Invalid folds.")
      
      M = length(folds)
    }
    else {
      if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
        stop("Invalid number of folds.")
      else {
        M = round(folds)
        folds = split(sample(1:n), rep(1:M, length = n)) 
      }
    }
  } 
  else { 
    folds = split(1:n, rep(1:n, length = n)) 
    M = n
  }
  
  
  error.mat = array(0, dim = c(ncomp, nmthdd, M))
  
  # in case the test set only includes one sample, it is better to advise the user to perform loocv
  stop.user = FALSE
  # set up a progress bar
  if (progressBar == TRUE) pb <- txtProgressBar(style = 3)
  
  for (i in 1:M) {
    if (progressBar == TRUE) setTxtProgressBar(pb, i/M)
    
    #set up leave out samples.
    omit = folds[[i]]
    
    # see below, we stop the user if there is only one sample drawn on the test set using MFold
    if(length(omit) == 1) stop.user = TRUE
    
    # the training set is NOT scaled
    X.train = X[-omit, ]
    Y.train = Y[-omit]
    X.test = matrix(X[omit, ], nrow = length(omit))
    
    spls.res = sPLSda(X.train, Y.train, ncomp, max.iter = max.iter , tol = tol, keepX=keepX)  
    # added: record selected features
    for(k in 1:ncomp){
      features[[k]] = c(unlist(features[[k]]), select.var(spls.res, comp = k)$name)
    }
    
    if (!is.null(spls.res$nzv$Position)) X.test = X.test[, -spls.res$nzv$Position]
    Y.predict = predict(spls.res, X.test, method = method.predict)$class
    error.mat[, , i] = sapply(Y.predict, error.fun, y = as.numeric(Y[omit]))
    
    
  } # end loop on i
  
  # warn the user that at least test set had a length of 1
  if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
  
  if (progressBar == TRUE) cat('\n')
  
  #-- compute the error --#
  res = apply(error.mat, 1:2, mean)
  
  rownames(res) = paste('ncomp', 1:ncomp, sep = " ")
  colnames(res) = names(Y.predict)
  
  # ---- extract stability of features ----- # NEW
  list.features = list()
  for(k in 1:ncomp){
    #remove the NA value that was added for initialisation
    remove.na = which(is.na(features[[k]]))
    # then summarise as a factor and output the percentage of appearance
    list.features[[k]] = sort(table(as.factor(features[[k]][-remove.na]))/M, decreasing = TRUE)
  }
  

  
  names(list.features) = paste('comp', 1:ncomp)
  
  result = list()
  result$error.rate = res
  result$features$stable = list.features
  
  # added
 # result$nzvX = nzv$Position
  
  
  method = "plsda.mthd"
  result$meth = "splsda.mthd"
  class(result) = c("perf", method)
  #updated outputs
  return(invisible(result))
}


# ---------------------------------------------------
# perf for gplsda object ----
# ---------------------------------------------------

perf.gPLSda <- function(object,                                               
                        method.predict = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),     
                        validation = c("Mfold", "loo"),                                          
                        folds = 10,                                                              
                        progressBar = TRUE,...)                                                        
{
  
  #-- initialising arguments --#
  # these data are the centered and scaled X output or the unmapped(Y) scaled and centered
  X = object$X
  level.Y = object$names$Y  #to make sure the levels are ordered
  Y = object$ind.mat
  Y = map(Y)
  Y = factor(Y,labels = level.Y)
  ncomp = object$ncomp
  n = nrow(X)
  keepX = object$keepX  
  
  # tells which variables are selected in X and in Y:
  
  
  ind.block.x = object$ind.block.x
  ind.block.y = object$ind.block.y
  
  
  
  
  tol = object$tol
  max.iter = object$max.iter
  
  # initialize new objects:
  features <- list()
  for(k in 1:ncomp){
    features[[k]] = NA
  }
  
  method.predict = match.arg(method.predict, choices = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"), several.ok = TRUE)
  if (any(method.predict == "all")) nmthdd = 3 
  else nmthdd = length(method.predict)  
  
  
  
  
  error.fun = function(x, y) {
    error.vec = sweep(x, 1, y, FUN = "-")
    error.vec = (error.vec != 0)
    error.vec = apply(error.vec, 2, sum) / length(y)
    return(error.vec)
  }
  
  #-- define the folds --#
  if (validation == "Mfold") {
    if (is.list(folds)) {
      if (length(folds) < 2 | length(folds) > n)
        stop("Invalid number of folds.")
      if (length(unique(unlist(folds))) != n)
        stop("Invalid folds.")
      
      M = length(folds)
    }
    else {
      if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
        stop("Invalid number of folds.")
      else {
        M = round(folds)
        folds = split(sample(1:n), rep(1:M, length = n)) 
      }
    }
  } 
  else { 
    folds = split(1:n, rep(1:n, length = n)) 
    M = n
  }
  
  
  error.mat = array(0, dim = c(ncomp, nmthdd, M))
  
  # in case the test set only includes one sample, it is better to advise the user to perform loocv
  stop.user = FALSE
  # set up a progress bar
  if (progressBar == TRUE) pb <- txtProgressBar(style = 3)
  
  for (i in 1:M) {
    if (progressBar == TRUE) setTxtProgressBar(pb, i/M)
    
    #set up leave out samples.
    omit = folds[[i]]
    
    # see below, we stop the user if there is only one sample drawn on the test set using MFold
    if(length(omit) == 1) stop.user = TRUE
    
    # the training set is NOT scaled
    X.train = X[-omit, ]
    Y.train = Y[-omit]
    X.test = matrix(X[omit, ], nrow = length(omit))
    
    spls.res = gPLSda(X = X.train,Y = Y.train, ncomp = ncomp , max.iter = max.iter, tol = tol, keepX=keepX, ind.block.x=ind.block.x)     ## change
    
    
    
    
    # added: record selected features
    for(k in 1:ncomp){
      features[[k]] = c(unlist(features[[k]]), select.var(spls.res, comp = k)$name)
    }
    
    if (!is.null(spls.res$nzv$Position)) X.test = X.test[, -spls.res$nzv$Position]
    Y.predict = predict(spls.res, X.test, method = method.predict)$class
    error.mat[, , i] = sapply(Y.predict, error.fun, y = as.numeric(Y[omit]))
    
    
  } # end loop on i
  
  # warn the user that at least test set had a length of 1
  if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
  
  if (progressBar == TRUE) cat('\n')
  
  #-- compute the error --#
  res = apply(error.mat, 1:2, mean)
  
  rownames(res) = paste('ncomp', 1:ncomp, sep = " ")
  colnames(res) = names(Y.predict)
  
  # ---- extract stability of features ----- # NEW
  list.features = list()
  for(k in 1:ncomp){
    #remove the NA value that was added for initialisation
    remove.na = which(is.na(features[[k]]))
    # then summarise as a factor and output the percentage of appearance
    list.features[[k]] = sort(table(as.factor(features[[k]][-remove.na]))/M, decreasing = TRUE)
  }
  
  
  
  names(list.features) = paste('comp', 1:ncomp)
  
  result = list()
  result$error.rate = res
  result$features$stable = list.features
  
  # added
  #result$nzvX = nzv$Position
  
  
  method = "plsda.mthd"
  result$meth = "splsda.mthd"
  class(result) = c("perf", method)
  #updated outputs
  return(invisible(result))
}

# ---------------------------------------------------
# perf for sgplsda object ----
# ---------------------------------------------------

perf.sgPLSda <- function(object,                                               
                         method.predict = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),     
                         validation = c("Mfold", "loo"),                                          
                         folds = 10,                                                              
                         progressBar = TRUE,...)                                                        
{
  
  #-- initialising arguments --#
  # these data are the centered and scaled X output or the unmapped(Y) scaled and centered
  X = object$X
  level.Y = object$names$Y  #to make sure the levels are ordered
  Y = object$ind.mat
  Y = map(Y)
  Y = factor(Y,labels = level.Y)
  ncomp = object$ncomp
  n = nrow(X)
  keepX = object$keepX  
  
  # tells which variables are selected in X and in Y:
  
  
  ind.block.x = object$ind.block.x
  ind.block.y = object$ind.block.y
  alpha.x = object$alpha.x
  alpha.y = object$alpha.y
  upper.lambda = object$upper.lambda
  
  
  
  tol = object$tol
  max.iter = object$max.iter
  
  # initialize new objects:
  features <- list()
  for(k in 1:ncomp){
    features[[k]] = NA
  }
  
  method.predict = match.arg(method.predict, choices = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"), several.ok = TRUE)
  if (any(method.predict == "all")) nmthdd = 3 
  else nmthdd = length(method.predict)  
  
  
  
  
  error.fun = function(x, y) {
    error.vec = sweep(x, 1, y, FUN = "-")
    error.vec = (error.vec != 0)
    error.vec = apply(error.vec, 2, sum) / length(y)
    return(error.vec)
  }
  
  #-- define the folds --#
  if (validation == "Mfold") {
    if (is.list(folds)) {
      if (length(folds) < 2 | length(folds) > n)
        stop("Invalid number of folds.")
      if (length(unique(unlist(folds))) != n)
        stop("Invalid folds.")
      
      M = length(folds)
    }
    else {
      if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
        stop("Invalid number of folds.")
      else {
        M = round(folds)
        folds = split(sample(1:n), rep(1:M, length = n)) 
      }
    }
  } 
  else { 
    folds = split(1:n, rep(1:n, length = n)) 
    M = n
  }
  
  
  error.mat = array(0, dim = c(ncomp, nmthdd, M))
  
  # in case the test set only includes one sample, it is better to advise the user to perform loocv
  stop.user = FALSE
  # set up a progress bar
  if (progressBar == TRUE) pb <- txtProgressBar(style = 3)
  
  for (i in 1:M) {
    if (progressBar == TRUE) setTxtProgressBar(pb, i/M)
    
    #set up leave out samples.
    omit = folds[[i]]
    
    # see below, we stop the user if there is only one sample drawn on the test set using MFold
    if(length(omit) == 1) stop.user = TRUE
    
    # the training set is NOT scaled
    X.train = X[-omit, ]
    Y.train = Y[-omit]
    X.test = matrix(X[omit, ], nrow = length(omit))
    
    spls.res = sgPLSda(X = X.train,Y = Y.train, ncomp = ncomp , max.iter = max.iter, tol = tol, keepX=keepX, ind.block.x=ind.block.x, alpha.x = alpha.x, upper.lambda = upper.lambda)     ## change
    
    
    
    
    # added: record selected features
    for(k in 1:ncomp){
      features[[k]] = c(unlist(features[[k]]), select.var(spls.res, comp = k)$name)
    }
    
    if (!is.null(spls.res$nzv$Position)) X.test = X.test[, -spls.res$nzv$Position]
    Y.predict = predict(spls.res, X.test, method = method.predict)$class
    error.mat[, , i] = sapply(Y.predict, error.fun, y = as.numeric(Y[omit]))
    
    
  } # end loop on i
  
  # warn the user that at least test set had a length of 1
  if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
  
  if (progressBar == TRUE) cat('\n')
  
  #-- compute the error --#
  res = apply(error.mat, 1:2, mean)
  
  rownames(res) = paste('ncomp', 1:ncomp, sep = " ")
  colnames(res) = names(Y.predict)
  
  # ---- extract stability of features ----- # NEW
  list.features = list()
  for(k in 1:ncomp){
    #remove the NA value that was added for initialisation
    remove.na = which(is.na(features[[k]]))
    # then summarise as a factor and output the percentage of appearance
    list.features[[k]] = sort(table(as.factor(features[[k]][-remove.na]))/M, decreasing = TRUE)
  }
  
  
  
  names(list.features) = paste('comp', 1:ncomp)
  
  result = list()
  result$error.rate = res
  result$features$stable = list.features
  
  # added
  # X = nzv$Position
  
  
  method = "plsda.mthd"
  result$meth = "splsda.mthd"
  class(result) = c("perf", method)
  #updated outputs
  return(invisible(result))
}













