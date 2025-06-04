---
title: "Q2 indicator for PLS model assessment"
output: html_document
date: "2025-06-03"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(sgPLS)
```

## What is Q² indicator ?

The $Q^2$ is a assessment indicator for PLS models; for each new component $h$, a new matrix $Y^{(h)}$ is obtained by deflation and compared to the corresponding prediction matrix $\hat{Y}^{(h)}$. The Q2 therefore takes this comparison into account. A $Q^2$ value close to $1$ indicates a good performance. To compute this figure, we must compute two more indicators : the $RSS$ and the $PRESS$.


$RSS_h = \sum_{i=1}^{n} \sum_{j=1}^{q} (Y^{(h)} - \hat{Y})^2 =  \sum_{i=1}^{n} \sum_{j=1}^{q} (Y^{(h+1)}_{i,j})^2$

$PRESS_h = \sum_{i\in test}^{n} \sum_{j=1}^{q} (Y_{i,j}^{(h)} - \hat{Y_{i,j}})^2 =  \sum_{i\in test}^{n} \sum_{j=1}^{q} (Y^{(h+1)}_{i,j})^2$

Then, $Q^2$ is defined by this formula :

$Q^2_h = 1 - \frac{PRESS_h}{RSS_{h-1}}$

## How to use Q² ?

We compare the value of this criterion to a certain limit $l$ ; this limit is conventionally equal to $1-0.95^2 = 0.0975$. As long as we have this inequality : $Q^2_h \geq l$, we keep on following iteration ; therefore we stop when we have this inequality : $Q^2_h < l$.

## Using Q² with R

The $Q^2$ function, available below and named as `q2.pls()`, takes four parameters : 

- the $X$ matrix which is the predictors values

- the $Y$ matrix which is the response values

- the mode : we must choose between "regression" or "canonical"

- the number of maximal components

Note that $X$ and $Y$ can be dataframes : they are automatically converted into matrix.

```{r}
q2.pls <- function(X,Y, mode = "regression", ncomp.max = 10){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  p <- ncol(X)
  q <- ncol(Y)
  if(nrow(X)!=nrow(Y)){ stop("X and Y have not the same number of rows.")}
  n <- nrow(X)
  X0 <- X
  Y0 <- Y
  q2 <- numeric(ncomp.max)
  
  if(mode == "regression"){
    
    # RSS0 computing
    Y.mean <- t(matrix(colMeans(Y), nrow = q, ncol = n)) #mean for each column
    RSS0 <- sum(colSums((Y0-Y.mean)^2))
    model.all <- sPLS(X = X, Y = Y, ncomp = ncomp.max, mode = mode)
    
  
    for(h in 1:ncomp.max){
      
      Y.hat.rss <- matrix(nrow = n, ncol = q) # prediction matrix for RSS
      Y.hat.press <- matrix(nrow = n, ncol = q) # prediction matrix for PRESS
      
      for(i in 1:n){
        
        
        # training on the dataset without the ith individual
        model <- sPLS(X = X[-i,], Y = Y[-i,], ncomp = h, mode = "regression")
        
        # predictions on the ith individual
        Y.hat.rss[i,] <- predict.sPLS(model.all, newdata = X[i,])$predict[,,max(h-1,1)]
        Y.hat.press[i,] <- predict.sPLS(model, newdata = X[i,])$predict[,,h]
        
        
        
      }
      # RSS computing
      RSS <- sum(colSums((Y-Y.hat.rss)^2))
      
      # PRESS computing
      PRESS <- sum(colSums((Y-Y.hat.press)^2))
      
      if(h == 1){
          PRESS1 <- PRESS
      }
      
      # new Y matrice
      Y <- Y - Y.hat.rss
      
      # Q2
      q2[h] <- 1-PRESS/RSS
      
    }# end h loop
    
  }else if(mode == "canonical"){ # CASE OF CANONICAL MODE
    
    # RSS0 computing
    X.mean <- t(matrix(colMeans(X), nrow = p, ncol = n)) #mean for each column
    RSS0 <- sum(colSums((X0-X.mean)^2))
    model.all <- sPLS(X = X, Y = Y, ncomp = ncomp.max, mode = "canonical")
    
    for(h in 1:ncomp.max){
      
      X.hat.rss <- matrix(nrow = n, ncol = p) # prediction matrix for RSS
      X.hat.press <- matrix(nrow = n, ncol = p) # prediction matrix for PRESS
      
      for(i in 1:n){
        
        
        # training on the dataset without the ith individual
        model <- sPLS(X = X[-i,], Y = Y[-i,], ncomp = h, mode = "canonical")
        
        # predictions on the ith individual
        mat.t.all <- as.matrix(model.all$variates$X)
        mat.c.all <- as.matrix(model.all$mat.c)
        mat.t <- as.matrix(model$variates$X)
        mat.c <- as.matrix(model$mat.c)
        mat.u <- as.matrix(model$loadings$X)
        #mat.td <- t(t(mat.t[,h]))
        
        X.hat.rss <- t(t(mat.t.all[,max(h-1,1)]))%*%t(mat.c.all[,max(h-1,1)])
        X.hat.press[i,] <- X[i,]%*%mat.u[,h]%*%t(mat.c[,h])
        
      } # end i loop
      
      # RSS computing
      RSS <- sum(colSums((X-X.hat.rss)^2))
      
      # PRESS computing
      PRESS <- sum(colSums((X-X.hat.press)^2))
      
      if(h == 1){
        PRESS1 <- PRESS
      }
      
      # new Y matrice
      X <- X - X.hat.rss
      
      # Q2
      q2[h] <- 1-PRESS/RSS

      
    }# end h loop
    
  }else{
    stop("The mode must be either regression or canonical.")
  }# end if loop

  
  q2[1] <- 1-PRESS1/RSS0 # correction of the first value in Q2 vector
  lim <- 0.0975
  h.best <- min(which(q2 > lim))
  
  # Plot
  plot(q2, type = "b", col = "blue", pch = 16,
       main = "Q² performance according to the number of components",
       xlab = "Number of components", ylab = "Q²")
  abline(h = 0.0975, col = "red", lty = 2)
  
  suggestion <- paste("best number of components : H =",h.best)
  
  q2.pls.results <- list(q2 = q2, suggestion = suggestion, h.best = h.best)
  return(q2.pls.results)
  
}
```


Let's give some simple examples : we will create and use two datasets:

- one is a dataset with only one response variable $Y$.

- the other is a dataset with five response variables $Y = (Y1,Y2,Y3,Y4,Y5)$.

- the last dataset contains real data about NIR spectra.

The function below allow to create the first two datasets.

```{r pressure, echo=FALSE}
create.data <- function(n = 40, p = 10, q = 1){
  X <- matrix(data = runif(n*p),n,p)
  U <- matrix(data = runif(q*p,-10,10), nrow = p, ncol = q)
  E <- matrix(data = rnorm(n*p),n,p)
  Y <- X%*%U + E
  D <- data.frame(X,Y)
  return(list(X = X,Y = Y))
}
```

By default, the population is set to $n = 40$ which is close to actual conditions. In this case, we have $p < n$ ; a not very large value of $p$ avoids a long time of execution.
Let's also notice that, on average, the response $Y$ is a linear combination from the predictors $X$. Indeed, the function include a matrix product $Y = XU + E$ with $U$ the weight matrix and $E$ matrix the gaussian noise. This linearity condition is important in order to have a good performance of the model, PLS method using linearity combinaison. 

### First data set

```{r}
data <- create.data()
X <- data$X
Y <- data$Y

print("X matrix")
print(head(X))
print("Y matrix")
print(head(Y))
```

This first dataset contains only one $Y$ variable. Now, let's compute q2 values.

```{r}
q2.pls(X,Y)
```

According to the graph, the $Q^2$ is rapidly close to $1$ from the second component. And even with the first component, it is greater than the limit. So here, one component may be enough.

### Second data set

```{r}
data <- create.data(q = 5)
X <- data$X
Y <- data$Y

print("X matrix")
print(head(X))
print("Y matrix")
print(head(Y))
```

This second dataset contains the five $Y$ variables as announced previously. Now, let's compute q2 values. 

```{r}
q2.pls(X,Y)
```

According to the graph, the $Q^2$ is rapidly close to $1$ from the second component. And even with the first component, it is greater than the limit. So here, one component may be enough.



### Second data set

```{r}
data <- create.data(q = 5)
X <- data$X
Y <- data$Y

print("X matrix")
print(head(X))
print("Y matrix")
print(head(Y))
```

This second dataset contains the five $Y$ variables announced previously. Now, let's compute q2 values. 

```{r}
q2.pls(X,Y)
```

According to the graph, the $Q^2$ is rapidly close to $1$ from the second component. And even with the first component, it is greater than the limit. So here, one component may be enough.


### Third data set

This real dataset deals with NIR spectra and density measurements of PET yarns. It contains $28$ rows ans $268$ $X$ variables for only $1$ $Y$ variable.

```{r}
library(pls)
data(yarn)
X <- yarn$NIR
Y <- yarn$density

print("Y matrix")
print(head(Y))
```

This second dataset contains the five $Y$ variables announced previously. Now, let's compute q2 values. 

```{r}
q2.pls(X,Y)
```

