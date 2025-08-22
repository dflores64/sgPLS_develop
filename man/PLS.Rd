\name{PLS}
\encoding{latin1}
\alias{PLS}

\title{Partial Least Squares (PLS)}

\description{A function for performing a partial least squares (OLS) method on two data sets with more columns than rows.
The PLS approach is therefore useful for data sets where OLS regression is not possible.
}

\usage{
PLS(X, Y, ncomp, mode = "regression",scale=TRUE)
}	

\arguments{
  \item{X}{numeric matrix of predictors.}
  \item{Y}{numeric vector or matrix of responses (for multi-response models).}
  \item{ncomp}{the number of components to include in the model (see Details).}
  \item{mode}{character string. What type of algorithm to use, (partially) matching 
    one of \code{"regression"} or \code{"canonical"}. See Details.}
  \item{scale}{a logical indicating if the orignal data set need to be scaled. By default \code{scale}=TRUE}
}

\details{
\code{PLS} function fits PLS models with \eqn{1, \ldots ,}\code{ncomp} components. 
Multi-response models are fully supported. 

The type of algorithm to use is specified with the \code{mode} argument. Two PLS 
algorithms are available: PLS regression \code{("regression")} and PLS canonical analysis 
\code{("canonical")} (see References). 

}

\value{
\code{PLS} returns an object of class \code{"PLS"}, a list 
that contains the following components:

  \item{X}{the centered and standardized original predictor matrix.}
  \item{Y}{the centered and standardized original response vector or matrix.}
  \item{ncomp}{the number of components included in the model.}
  \item{mode}{the algorithm used to fit the model.}
  \item{mat.c}{matrix of coefficients to be used internally by \code{predict}.}
  \item{variates}{list containing the variates.}
  \item{loadings}{list containing the estimated loadings for the \eqn{X} and 
	\eqn{Y} variates.}
  \item{names}{list containing the names to be used for individuals and variables.}
}

\references{
Liquet Benoit, Lafaye de Micheaux Pierre, Hejblum Boris, Thiebaut Rodolphe. Group and sparse group partial least square
approaches applied in genomics context. \emph{bioinformatics}
}

\author{Daniel FLORES.}

\seealso{\code{\link{sPLS}}, \code{\link{gPLS}}, \code{\link{sgPLS}}, \code{\link{predict}}, \code{\link{perf}}, \code{\link{cim}} and functions from \code{mixOmics} package: \code{summary}, \code{plotIndiv}, \code{plotVar}, \code{plot3dIndiv}, \code{plot3dVar}.}


\examples{
# Simulation of datasets X and Y with group variables

library(sgPLSdevelop)

## First example 

### paramaters
n <- 100
sigma.gamma <- 1
sigma.e <- 1.5
p <- 400
q <- 500

theta.x1 <- c(rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5), rep(1.5,15), 
              rep(0, 5), rep(-1.5, 15), rep(0, 325))
theta.x2 <- c(rep(0, 320), rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5),
              rep(1.5, 15), rep(0, 5), rep(-1.5, 15), rep(0, 5))

theta.y1 <- c(rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5), rep(1.5, 15),
              rep(0, 5), rep(-1.5, 15), rep(0, 425))
theta.y2 <- c(rep(0, 420), rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5),
              rep(1.5, 15), rep(0, 5), rep(-1.5, 15), rep(0, 5))                            

### covariance matrices
Sigmax <- matrix(0, nrow = p, ncol = p)
diag(Sigmax) <- sigma.e ^ 2
Sigmay <- matrix(0,nrow = q, ncol = q)
diag(Sigmay) <- sigma.e ^ 2

set.seed(125)

gam1 <- rnorm(n)
gam2 <- rnorm(n)

GAM <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE)
Thetax <- matrix(c(theta.x1, theta.x2), nrow = 2, byrow = TRUE)
Thetay <- matrix(c(theta.y1, theta.y2), nrow = 2, byrow = TRUE)
E1 <- rmvnorm(n, mean = rep(0, p), sigma = Sigmax, method = "svd")
E2 <- rmvnorm(n, mean = rep(0, q), sigma = Sigmay, method = "svd")
  
X <- GAM \%*\% Thetax + E1                                                
Y <- GAM \%*\% Thetay + E2

### PLS model
model.PLS <- PLS(X, Y, ncomp = 2, mode = "regression")


## Second example

train <- 1:40
test <- 41:50
n.test <- length(test)

d <- data.create(n = 50, p = 10, q = 2, list = TRUE)

X <- d$X[train,]
Y <- d$Y[train,]
X.test <- d$X[test,]
Y.test <- d$Y[test,]

ncompmax <- 10
model.pls <- PLS(X = X, Y = Y, ncomp = ncompmax, mode = "regression")
pred <- predict.PLS(model.pls, newdata = X.test)$predict

}


\keyword{regression}
\keyword{multivariate}
