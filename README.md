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

The $Q^2$ is an assessment indicator for PLS models; for each new component $h$, a new matrix $Y^{(h)}$ is obtained by deflation and compared to the corresponding prediction matrix $\hat{Y}^{(h)}$. The Q2 therefore takes this comparison into account. A $Q^2$ value close to $1$ indicates a good performance. To compute this figure, we must compute two more indicators : the $RSS$ and the $PRESS$.


$RSS_h = \sum_{i=1}^{n} \sum_{j=1}^{q} (Y^{(h)} - \hat{Y})^2 =  \sum_{i=1}^{n} \sum_{j=1}^{q} (Y^{(h+1)}_{i,j})^2$

$PRESS_h = \sum_{i\in test}^{n} \sum_{j=1}^{q} (Y_{i,j}^{(h)} - \hat{Y_{i,j}})^2 =  \sum_{i\in test}^{n} \sum_{j=1}^{q} (Y^{(h+1)}_{i,j})^2$

Then, $Q^2$ is defined by this formula :

$Q^2_h = 1 - \frac{PRESS_h}{RSS_{h-1}}$

## How to use Q² ?

We compare the value of this criterion to a certain limit $l$ ; this limit is conventionally equal to $1-0.95^2 = 0.0975$. As long as we have the inequality $Q^2_h \geq l$, we keep on following iteration ; therefore we stop when we have $Q^2_h < l$.

## Using Q² with R

The $Q^2$ function, available below and named as `q2.pls()`, takes four parameters : 

- the $X$ matrix which is the predictors values

- the $Y$ matrix which is the response values

- the mode : we must choose between "regression" or "canonical"

- the number of maximal components

Note that $X$ and $Y$ can be dataframes : they are automatically converted into matrix.



Let's give some simple examples : we will create and use three datasets:

- one is a dataset with only one response variable $Y$.

- the other is a dataset with five response variables $Y = (Y1,Y2,Y3,Y4,Y5)$.

- the last dataset contains real data about NIR spectra.

To access to predefined functions from sgPLSdevelop package, run this line :

```{r pressure, echo=FALSE}
library(sgPLSdevelop)
```

By default, the population is set to $n = 40$ which is close to actual conditions. In this case, we have $p < n$ ; a not very large value of $p$ avoids a long time of execution.
Let's also notice that, on average, the response $Y$ is a linear combination from the predictors $X$. Indeed, the function includes a matrix product $Y = XU + E$ with $U$ the weight matrix and $E$ matrix the gaussian noise. This linearity condition is important in order to have a good performance of the model, the PLS method using linearity combinaison. 

### First data set

```{r}
data <- data.create()
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

This second dataset contains the five $Y$ variables announced previously. Now, let's compute q2 values. 

```{r}
q2.pls(X,Y)
```

According to the graph, the $Q^2$ is rapidly close to $1$ from the second component. And even with the first component, it is greater than the limit. So here, one component may be enough.


### Third data set

This real dataset deals with NIR spectra and density measurements of PET yarns. It contains $n = 28$ rows and $p = 268$ $X$ variables for only $1$ $Y$ variable.

```{r}
library(pls)
data(yarn)
X <- yarn$NIR
Y <- yarn$density

print("Y matrix")
print(head(Y))
```

This second dataset contains the five $Y$ variables as announced previously. Now, let's compute the q2 values. 

```{r}
q2.pls(X,Y)
```

