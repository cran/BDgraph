\name{traceplot}
\alias{traceplot}

\title{ Trace plot of graph size }

\description{
Trace plot for graph size for the objects of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. 
It is a tool for monitoring the convergence of the sampling algorithms, BDMCMC and RJMCMC.
}

\usage{ traceplot ( x, acf = FALSE, pacf = FALSE, main = NULL, ... ) }

\arguments{
  \item{x}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }
  \item{acf}{ Visualize the autocorrelation functions for graph size. The default is FALSE. }
  \item{pacf}{ Visualize the partial autocorrelations for graph size. The default is FALSE. }
  \item{main}{ Graphical parameter (see plot). }
  \item{\dots}{ System reserved (no specific usage). }
}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv:1501.05108} 

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Gaussian Copula Graphical Modeling for Dupuytren Disease, \emph{arXiv:1501.04849} 
}

\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
# generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )

output <- bdgraph( data = data.sim, iter = 10000, burnin = 0, save.all = TRUE )

traceplot( output )

traceplot( output, acf = TRUE, pacf = TRUE )
}
}
