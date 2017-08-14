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

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C} 

Mohammadi, A., Massam H., and G. Letac (2017). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416} 
}

\author{ Abdolreza Mohammadi and Ernst Wit }

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
