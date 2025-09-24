## How to install and load sgPLSdevelop package ?

This package can be installed by runing the following lines

```{r install}
library(devtools) # if 'devtools' is already installed
install_github("dflores64/sgPLS_develop")
```

Once the package installed, it is necessary to load it by runing :

```{r loading}
library(sgPLSdevelop)
```

You can now try to execute some specific functions...

```{r}
# data building
data <- data.create(p = 10) 
X <- data$X
Y <- data$Y
print(head(X))

# model building
ncomp.max <- 8
model <- PLS(X,Y,mode = "regression", ncomp = ncomp.max)
print(head(model$loadings$X))
```


