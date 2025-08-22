\name{PLSda}
\alias{PLSda}
\title{Sparse Partial Least Squares Discriminant Analysis (sPLS-DA)}

\description{Function to perform Partial Least Squares to classify samples (supervised analysis).
}

\usage{
PLSda(X, Y, ncomp = 2)
}	

\arguments{
  \item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
  \item{Y}{a factor or a class vector for the discrete outcome.}
  \item{ncomp}{the number of components to include in the model (see Details).}
  }

\details{
\code{PLSda} function fit PLS models with \eqn{1, \ldots ,}\code{ncomp} components
to the factor or class vector \code{Y}. The appropriate indicator (dummy)
matrix is created.
}

\value{
\code{PLSda} returns an object of class \code{"PLSda"}, a list 
that contains the following components:

  \item{X}{the centered and standardized original predictor matrix.}
  \item{Y}{the centered and standardized indicator response vector or matrix.}
  \item{ind.mat}{the indicator matrix.}
  \item{ncomp}{the number of components included in the model.}
  \item{mat.c}{matrix of coefficients to be used internally by \code{predict}.}
  \item{variates}{list containing the variates.}
  \item{loadings}{list containing the estimated loadings for the \code{X} and 
	\code{Y} variates.}
  \item{names}{list containing the names to be used for individuals and variables.}
}


\author{Benoit Liquet and Pierre Lafaye de Micheaux.}

\seealso{\code{\link{sPLS}}, \code{\link{summary}}, 
\code{\link{plotIndiv}}, \code{\link{plotVar}}, 
\code{\link{cim}}, \code{\link{network}}, \code{\link{predict}}, \code{\link{perf}} and http://www.mixOmics.org for more details.}

\examples{
library(sgPLSdevelop)

ncompmax <- 5
d <- data.cl.create(classes = 4, n = 30, list = TRUE)
X <- d$X
Y <- d$Y

modele <- PLSda(X = X, Y = Y, ncomp = ncompmax)
}

\keyword{regression}
\keyword{multivariate}