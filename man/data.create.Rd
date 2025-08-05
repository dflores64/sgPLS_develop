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
data.create(n = 40, p = 10, q = 1, list = FALSE)
data.cl.create(n = 40, p = 10, classes = 2, list = FALSE)
}

\arguments{
  \item{n}{number of dataset rows}
  \item{p}{number of \eqn{X} variables}
  \item{q}{number of \eqn{Y} variables with \code{data.create} function} 
  \item{classes}{number of classes to generate with \code{data.cl.create} function}
  \item{list}{By default, return a dataframe. If the argument \code{list} is set to \code{TRUE}, a list including the dataframe is returned.}
}

\details{
By default, the population is set to \eqn{n=40} which is close to actual conditions. In this case, we have 
\eqn{p<n}. 

With \code{data.create} function, \eqn{Y} is a linear combinaison from each gaussian variable \eqn{X_j} of \eqn{X}.
Indeed, the function includes a matrix product to compute the response : \eqn{Y = XB+E} with \eqn{B} the weight (coefficients) matrix and \eqn{E} matrix the gaussian noise.
\eqn{B} matrix can be found in the list returned by the function (if \code{list = TRUE}).

With \code{data.cl.create} function, there is a link between \eqn{X} and the classes of \eqn{Y}.
The list returns also \code{Y.f} an other version of \code{Y} but from factor class.
}


\examples{

library(sgPLS)
library(sgPLSdevelop)

# data.create
data <- data.create(n = 20, p = 5, q = 2, list = TRUE)
X <- data$X
Y <- data$Y

# data.cl.create
data.cl <- data.cl.create(n = 20, p = 5, classes = 3, list = TRUE)
X <- data.cl$X
Y <- data.cl$Y
