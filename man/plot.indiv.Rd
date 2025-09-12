\name{plot.indiv}
\alias{plot.indiv}

\title{
Population projection on others dimension
}
\description{
This function displays scatterplot of the population (dataset individuals) projected on dimensions selected in parameters.
}
\usage{
plot.indiv(object, compX = c(1,2), compY = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class inheriting from "sPLS", "gPLS", "sgPLS", "sPLSda", "gPLSda" or "sgPLSda". The function will retrieve some key parameters stored in that object.}
  \item{compX}{The X components on what the population is projected.}
  \item{compY}{The Y components on what the population is projected (available only for regression).}
}
\details{
In the case of PLSDA, the scatterplot is accompanied by ellipses circling each point group.
}
\value{
  \item{data}{gives the latentes variables selected in parameters (in other words, the point coordinates).}
  \item{graphX}{scatterplot for X selected components.}
  \item{graphY}{scatterplot for Y selected components for regression case.}
}

\author{
Daniel FLORES
}

\seealso{
\code{\link{perf.PLS}},\code{\link{perf.PLSda}}
}
\examples{
library(sgPLSdevelop)

ncomp.max <- 10

# PLS 
data <- data.create(n = 20, p = 5, q = 2)
X <- data$X
Y <- data$Y
model.pls <- PLS(X,Y,ncomp = ncomp.max)
pop.model.pls <- plot.indiv(model.pls, compX = c(5,4), compY = c(1,2))

# PLSDA
data.cl <- data.cl.create(n = 20, p = 5, classes = 3)
X <- data.cl$X
Y <- data.cl$Y
model.plsda <- PLSda(X,Y,ncomp = ncomp.max)
pop.model.plsda <- plot.indiv(model.plsda, compX = c(5,4))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
