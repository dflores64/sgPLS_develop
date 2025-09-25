\name{data.create}
\alias{data.create}
\alias{data.cl.create}

\title{
Dataset simulation
}

\description{
These functions allow to generate a dataset with linear dependance between \eqn{Y} and \eqn{X}.
\code{data.create} is used for quantitative response while \code{data.cl.create} is used for qualitative response.
}

\usage{
data.create(n = 40, p = 10, q = 1, seed = NULL, list = TRUE)
data.cl.create(n = 40, p = 10, classes = 2, seed = NULL, list = TRUE)
}

\arguments{
  \item{n}{Number of dataset rows.}
  \item{p}{Number of \eqn{X} variables.}
  \item{q}{Number of \eqn{Y} variables with \code{data.create} function.} 
  \item{classes}{Number of classes to generate with \code{data.cl.create} function}
  \item{seed}{Uses \code{set.seed} function to generate data.}
  \item{list}{By default, returns a list including the dataframe. If the argument \code{list} is set to \code{FALSE}, a dataframe is returned.}
}

\details{
By default, the population is set to \eqn{n=40} which is close to actual conditions. In this case, we have 
\eqn{p<n}. 

\eqn{Y} is a linear combinaison from each gaussian variable \eqn{X_j} of \eqn{X}.
Indeed, the function includes a matrix product to compute the response : \eqn{Y = XB+E} with \eqn{B} the weight (coefficients) matrix and \eqn{E} matrix the gaussian noise.
\eqn{B} matrix can be found in the list returned by the function (if \code{list = TRUE}).

With \code{data.cl.create} function, \code{Y} object returned is a factor class.
}


\examples{
library(sgPLSdevelop, warn.conflicts = FALSE, verbose = FALSE, quietly = TRUE)

# data.create
data <- data.create(n = 20, p = 5, q = 2, list = TRUE)
X <- data$X
Y <- data$Y

# data.cl.create
data.cl <- data.cl.create(n = 20, p = 5, classes = 3, list = TRUE)
X <- data.cl$X
Y <- data.cl$Y
}

