\name{tuning.sPLS.XY}
\alias{tuning.sPLS.XY}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tuning of the parameters for sPLS models
}
\description{
This function calculate the best number of components and the best number of variables to select for sPLS models.
}
\usage{
tuning.sPLS.XY(X,Y,folds=10,validation=c("Mfold","loo"), ncomp = ncol(X))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{The centered and standardized original predictor matrix.}
  \item{Y}{The centered and standardized original response vector or matrix.}
  \item{validation}{Character. What kind of (internal) validation to use, matching one of "Mfold" or "loo". Default is "Mfold".}
  \item{folds}{The folds in the Mfold cross-validation.}
  \item{ncomp}{The number of components}
}

\details{
This function outputs the best tuning parameters by using MSEP criterion.
There is also two MSEP plots : the one shows the MSEP according to the number of components and the other shows the MSEP according to the number of selected variables in X and Y. 
Particularly, in multivariate sPLS (q>1), the "plot" is a table of MSEP values ; the best (minimum) value of MSEP is then colored in red.
}
\value{
  \item{MSEP.h}{Vector of MSEP values according to the number of components}
  \item{h.best}{Index of the minimum value of \code{MSEP.h}}
  \item{MSEP.q.p}{Matrix of MSEP values according the number of X selected variables (by column) and the number of Y selected variables (by row)}
  \item{keepX.best}{A vector of the best number of X selected variables repeated \code{h.best} times}
  \item{keepY.best}{A vector of the best number of Y selected variables repeated \code{h.best} times}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Daniel FLORES
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{tuning.sPLS.X}}, \code{\link{tuning.gPLS.X}}, \code{\link{tuning.sgPLS.X}}
}
\examples{
library(sgPLSdevelop)

# number of components
ncomp.max <- 8

# dataset
data <- data.create(p = 10, list = TRUE)

# model
X <- data$X
Y <- data$Y

# tuning
perf.res <- tuning.sPLS.XY(X,Y,ncomp = ncomp.max)
}