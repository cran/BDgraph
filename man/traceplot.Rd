\name{traceplot}
\alias{traceplot}

\title{Trace plot of graph size}

\description{
Trace plot for graph size for the objects of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. 
It is a tool for monitoring the convergence of the BD-MCMC sampling algorithm.
}

\usage{ traceplot ( output, acf = FALSE, pacf = FALSE, main = NULL, ... ) }

\arguments{
  \item{output}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }
  \item{acf}{Logical: if TRUE you will see the plot of autocorrelation functions for graph size. The default is FALSE. }
  \item{pacf}{Logical: if TRUE you will see the plot of partial autocorrelations for graph size. The default is FALSE. }
  \item{main}{ Graphical parameter (see plot). }
  \item{\dots}{ System reserved (no specific usage). }
}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Arxiv preprint arXiv:1501.05108v2} 

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Modeling of Dupuytren Disease Using Gaussian Copula Graphical Models, \emph{Arxiv preprint arXiv:1501.04849v2} 
}

\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
# generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )

output   <- bdgraph( data = data.sim, iter = 10000 )

traceplot(output)

traceplot( output, acf = TRUE, pacf = TRUE )
}
}
