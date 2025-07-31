\name{q2.PLS_doc}
\alias{q2.PLS}

\title{
PLS function performance assessment using Q2 indicator.
}
\description{
The \code{q2.PLS} function allows to assess PLS models using \eqn{Q2} criterion. 
This function returns a list including \eqn{Q2} values for each number of components.
A plot allows to visualize the model performance according to the number of components.
}

\usage{
q2.PLS(object, mode = "regression", ncomp.max = object$ncomp)
}

\arguments{
  \item{object}{Object of class inheriting from \code{"pls"}.}
  \item{mode}{Character string. What type of algorithm to use, (partially) matching one of "regression" or "canonical".}
  \item{ncomp.max}{The number of components desired for q2 computing (no more than the number of components computed in the model).}
}

\details{
A suggestion about the number of components selection is also given by Q2 criterion.
Endly, other indicators are given such as PRESS and RSS for each number of components 
and such as \eqn{PRESSj} and \eqn{RSSj} (matrices) for a given number of components and for a given column of Y dataset.
}

\value{
\eqn{PRESSj} and \eqn{RSSj} are matrices of size \eqn{ncomp.max} x \eqn{q}.
The row sums of \eqn{PRESSj} and \eqn{RSSj} give respectively PRESS and RSS.
}

\references{
K.-A. Le Cao, Zoe Welham, \emph{Multivariate data integration using R} (pages 172 to 174), MixOmics
}


\seealso{
\code{\link{PLS}, \link{perf.PLS}}
}

\examples{

library(sgPLSdevelop)

d <- data.create(p = 10, q = 3, list = TRUE)
X <- d$X
Y <- d$Y
n <- nrow(X)
ncomp.max <- 10
model.pls <- PLS(X = X, Y = Y, ncomp = ncomp.max, mode = "regression")

# q2 function
q2.res <- q2.PLS(model.pls, ncomp.max = ncomp.max, mode = "regression")

# q2 values
q2.res$q2

# PRESS values
q2.res$PRESS

# RSS values
q2.res$RSS

# PRESS values by column
q2.res$PRESSj

# RSS values by column
q2.res$RSSj

}
