\name{msep.PLS}
\alias{msep.PLS}

\title{
  PLS function performance assessment using \eqn{MSEP} indicator
}

\description{
  The \code{msep.PLS} function allows to assess PLS models using\eqn{MSEP} criterion. 
  This function returns a list including \eqn{MSEP} values for each number of components.
  \code{msep.PLS} gives also a suggestion about the number of components selection.
  A plot allows to visualize the model performance according to the number of components.
}

\usage{
  msep.PLS(object, ncomp = object$ncomp, K=nrow(object$X), plot = TRUE)
}

\arguments{
  \item{object}{Object of class inheriting from "pls".}
  \item{ncomp}{The number of components desired for MSEP computing (no more than the number of components computed in the model).}
  \item{K}{The number of blocks for cross-validation (leave-one-out by default).}
  \item{plot}{Boolean indicating if a plot must be output or not}
}

\details{
  A leave-one-out cross-validation is used by default.

  \eqn{K} must be a value between 2 and the the number of rows of the dataset used for training model.
  
  MSEP is defined by mean of squared error between the true \eqn{Y} values and their associated predictions.

  If \code{plot = TRUE}, the plot displays a curve of the $Q^2$ values.
}

\value{
  \item{\code{MSEP.train}}{\eqn{MSEP} of the training set WITHOUT cross-validation}
  \item{\code{MSEP.cv}}{\eqn{MSEP} of the training set WITH cross-validation}
  \item{\code{h.best}}{Indice of the minimum \code{MSEP.cv}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{PLS}, \link{perf.PLS}, \link{q2.PLS}}
}
\examples{
library(sgPLSdevelop)

## data and model creation 
d <- data.create(p = 10, list = TRUE)
n <- nrow(d$X)
ncomp.max <- 10
X <- d$X
Y <- d$Y
model <- PLS(X,Y,ncomp = ncomp.max, mode = "regression")

## using msep.PLS function
msep.res <- msep.PLS(model, ncomp = ncomp.max, K = n)
msep.res$MSEP 
msep.res$h.best #number of components suggestion 


}


