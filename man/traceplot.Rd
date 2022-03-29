\name{traceplot}
\alias{traceplot}

\title{ Trace plot of graph size }

\description{
Trace plot for graph size for the objects of \code{S3} class "\code{bdgraph}", from function \code{\link{bdgraph}}. 
It is a tool for monitoring the convergence of the sampling algorithms, BDMCMC and RJMCMC.
}

\usage{ traceplot ( bdgraph.obj, acf = FALSE, pacf = FALSE, main = NULL, ... ) }

\arguments{
  \item{bdgraph.obj}{ object of \code{S3} class "\code{bdgraph}", from function \code{\link{bdgraph}}. 
					  It also can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}} of \code{R} package \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}}.
  }
  \item{acf}{ visualize the autocorrelation functions for graph size. }
  \item{pacf}{ visualize the partial autocorrelations for graph size. }
  \item{main}{ graphical parameter (see plot). }
  \item{\dots}{ system reserved (no specific usage). }
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138, \doi{10.1214/14-BA889}

Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845, \doi{10.1214/18-AOAS1164}

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645, \doi{10.1111/rssc.12171}

Mohammadi, A. and Dobra, A. (2017). The \code{R} Package \pkg{BDgraph} for Bayesian Structure Learning in Graphical Models, \emph{ISBA Bulletin}, 24(4):11-16
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} and Ernst Wit }

\seealso{ \code{\link{plotcoda}}, \code{\link{bdgraph}}, \code{\link{bdgraph.mpl}} }

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )
   
bdgraph.obj <- bdgraph( data = data.sim, iter = 10000, burnin = 0, save = TRUE )
   
traceplot( bdgraph.obj )
   
traceplot( bdgraph.obj, acf = TRUE, pacf = TRUE )
}
}

\keyword{hplot}
