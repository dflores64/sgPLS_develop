\name{perf.PLS}
\alias{perf.PLS}

\title{
  PLS function performance assessment using \eqn{MSEP} indicator
}
\description{
  The \code{perf.PLS} function allows to assess PLS models using\eqn{MSEP} criterion. 
  This function returns a list including \eqn{MSEP} values for each number of components.
  \code{perf.PLS} gives also a suggestion about the number of components selection.
  A plot allows to visualize the model MSEP according to the number of components.
}
\usage{
  perf.PLS(object, ncomp = object$ncomp, K=nrow(object$X))
}

\arguments{
  \item{object}{Object of class inheriting from "pls".}
  \item{ncomp}{The number of components desired for MSEP computing (no more than the number of components computed in the model).}
  \item{K}{The number of blocks for cross-validation (leave-one-out by default).}
}

\details{
  \eqn{K} must be a value between 2 and the the number of rows of the dataset used for training model.
  
  MSEP is defined by mean of squared error between the true \eqn{Y} values and their associated predictions.
}



%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{PLS}, \link{q2.PLS}}
}
\examples{

library(sgPLSdevelop)

## data and model creation 
d <- data.create(p = 10, list = TRUE)
X <- d$X
Y <- d$Y
n <- nrow(X)
ncomp.max <- 10
modele <- PLS(X,Y,ncomp = ncomp.max, mode = "regression")

## using perf.PLS function
perf.res <- perf.PLS(modele)
perf.res$MSEP 
perf.res$h.best #optimal components number suggestion 

}
