\name{traceplot}
\alias{traceplot}

\title{ Trace plot of graph size }

\description{
Trace plot for graph size for the objects of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. 
It is a tool for monitoring the convergence of the sampling algorithms, BDMCMC and RJMCMC.
}

\usage{ traceplot ( bdgraph.obj, acf = FALSE, pacf = FALSE, main = NULL, ... ) }

\arguments{
  \item{bdgraph.obj}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }
  \item{acf}{ Visualize the autocorrelation functions for graph size. }
  \item{pacf}{ Visualize the partial autocorrelations for graph size. }
  \item{main}{ Graphical parameter (see plot). }
  \item{\dots}{ System reserved (no specific usage). }
}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv preprint arXiv:1501.05108} 

Dobra, A. and A. Mohammadi (2017). Loglinear Model Selection and Human Mobility, \emph{arXiv preprint arXiv:1711.02623}

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C} 

Mohammadi, A. and A. Dobra (2017). The \code{R} Package \pkg{BDgraph} for Bayesian Structure Learning in Graphical Models, \emph{ISBA Bulletin}, 24(4):11-16
}

\author{ Reza Mohammadi and Ernst Wit }

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )
   
bdgraph.obj <- bdgraph( data = data.sim, iter = 10000, burnin = 0, save.all = TRUE )
   
traceplot( bdgraph.obj )
   
traceplot( bdgraph.obj, acf = TRUE, pacf = TRUE )
}
}
